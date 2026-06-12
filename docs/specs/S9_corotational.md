# S9 交接規格 — Co-rotational 大位移(平面幾何非線性梁,NR + load stepping)

> 研究依據:`docs/research/WS_F_corot.md`(完整 CR 路線)+ `docs/research/WS_F2_corot3d_opensees.md`
> (3D CR 的 OpenSees `CorotCrdTransf3d` 逐行公式,供 S9b)+ `WS_R2_experiments.md §9`(elastica
> shooting 表)。**本檔描述「已實作」狀態**(平面 v1);3D 全式設計見 §⑪ 與 WS_F2。
>
> **狀態:S9 平面 v1 完成(`02cb4d3` 後一輪),經 4-agent 對抗審核(零 CRITICAL 殘留),五腿全綠** —
> standalone **F1-F50** / UE **47** / OpenSees PASS / audit **95** / CLI round-trip **9**。
> `$ExpectedUeTests` 46→47;audit 90→95;CLI 8→9。
>
> **誠實定位(本階段最重要)**:S9 v1 = **平面(global XY 平面、繞 global Z 軸彎曲、用 `Iz`、
> Euler-Bernoulli)co-rotational** 大位移 driver。**3D 通用 CR**(扭轉、雙軸彎曲、spin correction、
> SO(3) 有限旋轉)= **S9b 後續**,公式已備(WS_F2,OpenSees CorotCrdTransf3d 逐行)。理由 §⑨。

---

## ① 目標 / 不做

**做** — `runCorotational`:平面梁的幾何非線性(大位移)分析。每根梁隨其**當前弦**共旋,局部仍小應變
(small-strain large-rotation CR);局部勁度用既有 Euler-Bernoulli 梁(`EA`、`EIz`)。求解 =
**Newton-Raphson + load stepping**(載重因子 λ:0→1 分步,每步 NR 迭代到平衡)。

- **平面簡化(關鍵)**:所有旋轉繞**單一** global Z 軸 → **可交換可加** → 節點轉角直接由 `Rz` 位移
  DOF 累積,**不需** SO(3) 旋轉矩陣 / `exp`/`log` map。剛體繞 Z 大旋轉天然零內力(§⑥ 機器精度驗證)。
- **切線勁度** `K_t = Bᵀ K_l B + K_geom`(標準 2D co-rotational,Crisfield Vol.1 Ch.7):
  `K_geom = (N/Ln)·q·qᵀ + ((Mi+Mj)/Ln²)·(p·qᵀ + q·pᵀ)`,`p=[-c,-s,0,c,s,0]ᵀ`(=∂Ln/∂d)、
  `q=[s,-c,0,-s,c,0]ᵀ`(=Ln·∂β/∂d),對稱。

**不做(誠實邊界,§⑨)**:
- **3D 通用 CR**(扭轉/雙軸/spin correction):S9b,WS_F2 公式已備。
- **snap-through / 極限點後路徑**:載重控制 NR 在極限點 `det(K_T)=0` 發散 → **誠實偵測為 `diverged`、
  不靜默算錯**;弧長法(Riks/Crisfield)留後續。主 oracle(橫向端載懸臂)**單調無極限點**。
- **member UDL / prescribed 支承位移 / 殼 / grillage**:driver **明確拒絕**(回 `singular`+diagnostic,
  不靜默丟載),v1 只吃**節點力**。tension-only / 移除 / 塑鉸 / release 在 CR **不啟用**(純彈性大位移)。
- **動力大位移**、接觸、大應變(非 geometrically-exact Reissner;小應變假設下與 GE 梁收斂同解)。

---

## ② 公開 API(`Public/FrameCore/CorotationalAnalysis.h`,POD/std、零 Eigen)

```cpp
struct CorotationalOptions {
    int  loadSteps = 10;     // equal lambda increments 0->1 (elastica alpha=10 wants ~10-30)
    int  maxIter   = 50;     // NR iterations per load step (non-convergence at a NON-limit point = raise loadSteps)
    real tolR      = 1e-9;   // ||residual_free|| / ||lambda*F_ext_free||      (force residual)
    real tolU      = 1e-12;  // ||du_free|| / max(1,||u_free||)                (displacement)
    SolveOptions solve;      // pivotTol passthrough (useTimoshenko reserved; v1 is Euler-Bernoulli)
};
struct CorotationalResult {
    bool        converged = false;
    bool        diverged  = false;   // limit point: K_T not positive-definite (snap-through; needs arc-length)
    int         loadStepsCompleted = 0;
    int         totalIterations    = 0;
    real        lastResidual       = 0;
    SolveResult finalState;          // large-displacement state (u, reactions, member local end forces)
};
FRAMECORE_API CorotationalResult runCorotational(const FrameModel& model, const CorotationalOptions& opts = {});
```
**無新 model / SolveResult 欄位、不動 `modelFingerprint`**(`CorotationalOptions` 是 solve 參數,同
`PDeltaOptions` 列)。`finalState.u` 的 `Rz` 分量 = 累積平面轉角(平面繞單軸,直接可加,無歧義)。

---

## ③ 資料流 / 架構(獨立 driver,**不改 `IElement`**)

侵入性評估結論:**獨立 NR driver,`IElement` 與所有既有 element 零改動**。理由:`IElement` 幾何凍結於
初始(`prepare` 一次算 `kl_/T_/Qf_`),CR 每迭代要在變形構型重算;改介面會迫使所有 element(含殼)實作
用不到的 hook 並威脅既有逐位元行為。仿 `runPDelta` / `runProgressiveCollapse` 的「獨立 driver + const
工作副本 + 自組矩陣 + recover 對齊 `SolveResult`」模式;driver 內持輕量 `CrBeam`(**非** `IElement`),
複用 `localStiffness12` 的同款 `EA/EIz`(直接算 natural 勁度)。

**driver 狀態**(工作副本,const model 不變,安全並行):每根 `CrBeam` 存 `L0/beta0/EA/EIz/betaPrev`;
全域累積位移 `u[6N]`(平移 + `Rz` 角)。**拒絕(reject helper,皆零填 `finalState.u` 防消費端越界)**:
active 殼、active 桿上的 UDL、非零 prescribed 位移、**桿有平面外 z 延伸**(`|dz|>1e-9·L3D`,否則只用 x,y
會靜默把 3D 桿當 xy 投影)、**已形成塑鉸**(否則以全剛度算)。皆 `singular`+diagnostic(§⑨)。

**一次 NR 迭代**:對每根 active 桿用當前節點座標 `(Xi+ui, Xj+uj)` 算 `Ln/β`(atan2,相對 `betaPrev`
unwrap)→ `ū=Ln−L0`、`θ̄_{i,j}=Rz_{i,j}−(β−β0)` → natural `[N,Mi,Mj]=k_l·[ū,θ̄_i,θ̄_j]` →
`fe=Bᵀ[N,Mi,Mj]`、`Ke=BᵀK_lB+K_geom` → scatter 到 6N(平面 DOF `Ux,Uy,Rz`)→ residual
`r=λF−f_int` → reduceFF → LDLT(`vectorD` 測正定;非正定 = 極限點 → `diverged`)→ `Δu` → `u+=Δu`(含
`Rz` 累加)→ 收斂判據。load step 收斂後推進 λ;recover `finalState`。

**接點複用**:`localAxes` 概念(此處平面直接 atan2)、natural `EA/EIz`(同 `localStiffness12` 的軸/彎
係數)、`reduceFF`/LDLT/`ldltPosDef`(同 `PDeltaAnalysis.cpp`)、`gdof`/`nodeIndex`。

---

## ④ 演算法(平面 2D co-rotational,Crisfield Vol.1 Ch.7)

**local natural 變形**(3-DOF):`ū=Ln−L0`、`θ̄_i=Rz_i−α`、`θ̄_j=Rz_j−α`,`α=β−β0`(弦剛體轉角)。
**natural 勁度**(EB):`N=(EA/L0)ū`、`Mi=(EIz/L0)(4θ̄_i+2θ̄_j)`、`Mj=(EIz/L0)(2θ̄_i+4θ̄_j)`。
**B(3×6)** over `[uxi,uyi,rzi,uxj,uyj,rzj]`(`c=cosβ,s=sinβ`):
```
B0 (ū)   = [ -c,    -s,    0,  c,    s,    0 ]
B1 (θ̄_i) = [ -s/Ln, c/Ln, 1,  s/Ln, -c/Ln, 0 ]
B2 (θ̄_j) = [ -s/Ln, c/Ln, 0,  s/Ln, -c/Ln, 1 ]
```
**內力** `f_int = Bᵀ[N,Mi,Mj]`(虛功一致)。**切線** `K_t = Bᵀ k_l B + K_geom`,
`K_geom = (N/Ln)q qᵀ + ((Mi+Mj)/Ln²)(p qᵀ + q pᵀ)`,`p=[-c,-s,0,c,s,0]ᵀ`、`q=[s,-c,0,-s,c,0]ᵀ`
(`K_geom` 對稱,經 §⑥ P-Delta 退化 oracle 釘住軸力幾何項正確)。
**旋轉更新**:`Rz` 直接累加(平面繞單軸 abelian,無 2π 奇異);`β` 用 atan2 並相對 `betaPrev` ±2π
unwrap 保連續。**剛體零內力**:剛體繞 Z 轉 φ → `ū=0`、`θ̄=Rz−(β−β0)=φ−φ=0` → `f_int=0`(§⑥ 機器精度)。

---

## ⑤ 檔案

**新增**:`Public/FrameCore/CorotationalAnalysis.h`、`Private/CorotationalAnalysis.cpp`(driver + 匿名
namespace `CrBeam`/`crCompute`/`ldltPosDef`;Eigen 走 `FrameEigen.h`,零洩漏)。
**build/gate 同步**(加 1 個 `.cpp`):`build.bat` / `build_linear_audit.bat` / `build_cli.bat` /
`build_capi.bat` 源檔清單 +`CorotationalAnalysis.cpp`;`main.cpp` +F50 + `FrameTestFixtures.h`
+`cantileverPlanarTipShearN`;`Private/Tests/CorotationalTest.cpp`(UE)+ `run_gate.ps1`
`$ExpectedUeTests` 46→47;`linear_deep_audit.cpp` +`testCorotational`(5 checks,+5 → 95);`frame_cli_core.cpp`
+`COROT` 指令;`Tools/cli_roundtrip.py` +COROT check(8→9);`docs/CLI_PROTOCOL.md` `COROT`。
**不動**:`IElement.h`、所有既有 element、`FrameSolver.cpp`、`SolveResult.h`/`SolveOptions.h`(無新欄位)、
`FrameModel.cpp`(validate/fingerprint)。**`.gitignore`/`ArchSim.uproject`/`Plugins/LevelSim/` 絕不碰。**

---

## ⑥ Oracle(誠實分級,全 `[VERIFIED]`)

- **F50a elastica shooting 表(主 oracle)**:橫向端載懸臂,`α=PL²/EIz=1,2,5,10` vs WS_F F-3 獨立
  shooting 表(δv/L、δh/L)。**實測**:N16 δv/L rel **1.1e-4 ~ 6.5e-4**(優於 0.1% 目標),δh/L < 0.07%;
  N4 粗網格 0.16%~1.02%。容差 N16 δv 1.5e-3 / δh 5e-3 / N4 1.5%(皆 2~10× 餘裕)。
- **F50b 平面剛體旋轉不變性**:整模型+載重繞 Z 轉 φ → `|u_tip|` 不變。**實測 rel = 0.00e+00**
  (機器精度;= 剛體旋轉零偽內力的鐵證,共旋框架精確 frame-indifferent)。
- **F50c P-Delta 退化**:小位移軸壓柱 + 微小側載,CR sway == 線性化 `runPDelta`。**實測 rel 5.0e-3**
  (容差 1.2e-2;驗證 CR 切線是 P-Delta 幾何勁度的有限旋轉推廣)。
- **audit `testCorotational`(+5)**:elastica α=5(rel 3.5e-4)、剛體旋轉不變性(0.0)、P-Delta 退化
  (真 rel **5.0e-3**)、**殼模型拒絕守門**、**小位移平衡**(stiff 懸臂 base `Ry=-P` 精確 + member-0
  `|Mz|=P*L`,rel 3.5e-7;驗 recover 反力 + 端彎矩,補端力 oracle 缺口)。
- **CLI COROT**:`frame_cli COROT` 經文字橋跑平面 elastica α=5,dv/L=0.714138 vs 0.713792(round-trip）。
- **誠實**:elastica shooting 表 = **獨立第三方解**(與 FrameCore 完全獨立)。OpenSees `corotCrdTransf`
  對照**未加**(elastica 已是獨立第三方;3D corot 對照排 S9b,WS_F2 §F50d 注意忽略 element loads)。

---

## ⑦ Gate

F50(a/b/c);UE `FrameCore.Corotational.ElasticaCantilever`;audit `testCorotational`(5 checks,95);
`cli_roundtrip` COROT(9);`$ExpectedUeTests` **47**;四 build 腳本 +`CorotationalAnalysis.cpp`。
**commit 前五腿全綠**(`run_gate.ps1 -RequireOpenSees`)— 已驗證 PASS。

---

## ⑧ 效能

CR 是迭代(load stepping × NR),每迭代 = 一次組裝 + 一次 fresh LDLT(切線每迭代變,**無法 factorize-once
重用** — 同 `runProgressiveCollapse` 每步 fresh factor 的誠實成本)。典型 elastica α=10 N16 ≈ 20 步 ×
~6 迭代 = ~120 次 factor/solve(小模型秒級)。額外狀態 O(N)。不參與 ReSolve / 模態 / 屈曲(CR 切線非定常);
誠實標:CR 與既有分析無交互 oracle(本階段 oracle 聚焦 driver 正確性 + P-Delta 退化銜接)。

---

## ⑨ 誠實邊界 / novelty

- **CR = Crisfield / Battini 成熟文獻方法**,FrameCore 為其**實作**,非新方法。`[VERIFIED]` = 文獻方法的
  oracle 化(elastica / 剛體 patch / P-Delta 退化);`[NEW CODE]` 僅限「獨立 NR driver 整合 + 平面 CR 接線
  + `SolveResult`/CLI 對接」。對標禁裸宣稱「優於 Karamba」,附先行技術定位(CR 是業界標準,非 FrameCore 獨有)。
- **平面 v1(核心邊界)**:僅 global XY 平面、繞 Z、`Iz`、Euler-Bernoulli。caller **必須**約束平面外 DOF
  (`Uz,Rx,Ry`)否則 LDLT 報 singular;**桿不得有平面外 z 延伸**(driver **守門拒絕**,否則只用 x,y 會靜默把
  3D 桿當 xy 投影 → 審核抓到的 CRITICAL,已修)。3D 通用 CR = **S9b**(WS_F2 公式已備)。
- **無 snap-through**:極限點 `diverged`(不靜默算錯);弧長法後續。主 oracle 單調不受限。
- **只吃節點力**:member UDL / prescribed 位移 / **已形成塑鉸** 皆 **拒絕**(不靜默丟載/算錯;早返回零填
  `finalState.u` 防消費端越界 — 審核抓到的 CRITICAL segfault,已修)。tension-only/release 不啟用。
- **CR = 小應變大旋轉**(非大應變、非 GE Reissner;小應變下與 GE 收斂同解)。
- **既有 `assembleGeometric` 是 P-Delta 線性化非完整 CR**(勿混淆);S9 CR 切線是其有限旋轉推廣
  (F50c 退化即此關係之證)。**塑鉸方向耦合 → S10 必在 S9 後**(R4)。

---

## ⑩ 風險 / fallback(實作後回填)

1. **B/K_geom 推導**:三道 oracle 全過(剛體不變性機器精度 → P-Delta 退化 2.4e-4 → elastica 1e-4),
   推導正確。`K_geom` 對稱(無集中力矩問題)。
2. **平面外 free DOF → singular**:caller 須約束(fixture 已約束);LDLT 偵測,不靜默。
3. **極限點**:`K_T` 非正定 = `diverged`(誠實),測試只取單調案例。snap-through(Williams)誠實標後續。
4. **β atan2 wrap**:相對 `betaPrev` ±2π unwrap(elastica 轉角 <90° 不觸發,守門備援)。
5. **UDL/prescribed/殼/平面外 z/塑鉸誤用**:driver 拒絕 + diagnostic(不靜默;早返回零填 `finalState.u`);
   audit 有殼拒絕守門 check;CLI 實測平面外桿 + 殼皆優雅拒絕(`SINGULAR 1`,exit 0 無 crash)。
6. **build/gate 同步**:四 build 腳本 +cpp;`$ExpectedUeTests` 47;`COROT` 用 `if(ss>>tok)` 向後相容。
7. **Eigen 洩漏 / UE unity / SAL 巨集**:Eigen 走 `FrameEigen.h`;Adaptive Build 自動排除新 cpp 出 unity
   (UE build 確認);UE 測試未用 `IN/OUT`。
8. **commit 衛生**:顯式 `git add` 源碼,勿 `-A`(`.gitignore`/`uproject`/`LevelSim`/`frame_capi.*` 禁碰)。

---

## ⑪ S9b 後續 — 3D 通用 CR(公式已備,未實作)

3D co-rotational(扭轉 + 雙軸彎曲 + spin correction + SO(3) 有限旋轉)的**逐行可實作公式**已整理於
`docs/research/WS_F2_corot3d_opensees.md`(錨定 OpenSees `CorotCrdTransf3d.cpp` 開源逐行 + Crisfield 1990
/ Battini 2002)。要點:節點維護 `R_node∈SO(3)`、spatial incremental update `R←exp(skew(Δθ))R`(避 2π
奇異)、元素框架平均三元組 `E`、局部變形 `asin` 提取、內力 `f=Tᵀpb`、切線 `TᵀK_lT+K_geom`(Ksigma1/2/3)。
oracle 增:3D 任意軸剛體旋轉 patch、空間 elastica、OpenSees corot 對照。S10(N-M 塑鉸)在 S9b 後(R4)。
