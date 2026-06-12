# S9 進度 — Co-rotational 大位移(平面幾何非線性梁,NR + load stepping)

> 接續 `PROGRESS_S8.md`(S8 殼路線 `02cb4d3`)。規格 `docs/specs/S9_corotational.md`,研究依據
> `docs/research/WS_F_corot.md`(CR 路線)+ `WS_F2_corot3d_opensees.md`(3D 公式,供 S9b)+
> `WS_R2_experiments.md §9`(elastica 表)。

## 摘要

S9 v1 新增**獨立的平面 co-rotational 大位移 driver** `runCorotational`(`CorotationalAnalysis.{h,cpp}`)
— 幾何非線性、Newton-Raphson + load stepping。**不改 `IElement`、不動既有 element / dispatch /
`SolveResult` / `modelFingerprint`**(獨立 driver,仿 `runPDelta`/`runProgressiveCollapse` 工作副本模式)。

- **平面 v1**:global XY 平面、繞 global Z 軸彎曲、用 `Iz`、Euler-Bernoulli。每根梁隨當前弦共旋,局部小應變
  餵 natural `EA/EIz` 勁度。**平面繞單軸旋轉可交換可加** → 節點轉角直接由 `Rz` 位移 DOF 累積,**不需** SO(3)
  `exp`/`log`(剛體繞 Z 大旋轉天然零內力)。
- **切線** `K_t = Bᵀ k_l B + K_geom`,`K_geom=(N/Ln)q qᵀ+((Mi+Mj)/Ln²)(p qᵀ+q pᵀ)`(對稱;標準 2D
  co-rotational,Crisfield Vol.1 Ch.7)。
- **誠實守門**:active 殼 / active 桿 UDL / 非零 prescribed 位移 / **桿有平面外 z 延伸** / **已形成塑鉸** →
  拒絕(`singular`+diagnostic,不靜默算錯;早返回亦零填 `finalState.u` 防消費端越界);極限點 `K_T`
  非正定 → `diverged`(無 snap-through)。

**五腿全綠**:standalone **F1-F50** / UE **47**(+1 `Corotational.ElasticaCantilever`)/ OpenSees PASS /
audit **95**(+5 `testCorotational`)/ CLI round-trip **9**(+1 COROT)。`$ExpectedUeTests` 46→47。
**經 4-agent 對抗審核**(零 CRITICAL 殘留;2 CRITICAL + 2 MAJOR + MINOR 全修,見下「對抗式審核」)。

## 數學(平面 2D co-rotational)

natural 變形 `ū=Ln−L0`、`θ̄_{i,j}=Rz_{i,j}−α`(`α=β−β0` 弦剛體轉角)。natural 勁度 EB:
`N=(EA/L0)ū`、`Mi=(EIz/L0)(4θ̄_i+2θ̄_j)`、`Mj=(EIz/L0)(2θ̄_i+4θ̄_j)`。B(3×6)、內力 `Bᵀ[N,Mi,Mj]`、
切線 `Bᵀk_lB+K_geom`(`p=[-c,-s,0,c,s,0]ᵀ`、`q=[s,-c,0,-s,c,0]ᵀ`,`K_geom` 二次驗證對稱)。β 用 atan2
相對 `betaPrev` ±2π unwrap。

## Oracle(實測)

- **F50a elastica**(主):橫向端載懸臂 `α=PL²/EIz=1,2,5,10` vs WS_F F-3 獨立 shooting 表。N16 δv/L rel
  **1.1e-4~6.5e-4**(優於 0.1% 目標),δh/L <0.07%;N4 粗 0.16%~1.02%。
- **F50b 平面剛體旋轉不變性**:模型+載重繞 Z 轉 φ → `|u_tip|` rel **0.00e+00**(機器精度;剛體零偽內力鐵證)。
- **F50c P-Delta 退化**:小位移 CR sway == `runPDelta` rel **5.0e-3**(CR 切線 = 幾何勁度的有限旋轉推廣)。
- **audit testCorotational(+4)**:elastica(3.5e-4)、剛體不變性(0.0)、P-Delta 退化(2.4e-4)、殼模型拒絕守門。
- **CLI COROT**:文字橋平面 elastica α=5,dv/L=0.714138 vs 0.713792。

## 檔案改動

- 新增 `Public/FrameCore/CorotationalAnalysis.h`、`Private/CorotationalAnalysis.cpp`。
- `main.cpp` +F50;`FrameTestFixtures.h` +`cantileverPlanarTipShearN`;`Private/Tests/CorotationalTest.cpp`(新)。
- `linear_deep_audit.cpp` +`testCorotational`(5 checks,90→95);`frame_cli_core.cpp` +`COROT`;
  `cli_roundtrip.py` +COROT(8→9);`run_gate.ps1` `$ExpectedUeTests` 46→47。
- 四 build 腳本(build/build_linear_audit/build_cli/build_capi)源檔清單 +`CorotationalAnalysis.cpp`。
- 文檔:`docs/specs/S9_corotational.md`(重寫為平面 v1)、`docs/research/WS_F2_corot3d_opensees.md`(3D 公式,S9b)、
  本檔。`CLI_PROTOCOL.md` 待補 `COROT`。
- **不動**:`IElement.h` / 所有 element / `FrameSolver.cpp` / `SolveResult.h` / `SolveOptions.h` /
  `FrameModel.cpp`(fingerprint/validate)。**`.gitignore`/`ArchSim.uproject`/`Plugins/LevelSim/` 未碰。**

## 誠實邊界 / novelty

- **CR = Crisfield/Battini 文獻方法**,FrameCore 為其**實作**,非新方法。`[NEW CODE]` 僅限 driver 整合 + 平面接線。
- **平面 v1**:僅 XY/繞 Z/`Iz`/EB;caller 須約束平面外 DOF。**3D 通用 CR = S9b**(WS_F2 公式已備)。
- **無 snap-through**(極限點 `diverged`);**只吃節點力**(UDL/prescribed 拒絕);tension-only/塑鉸/release 不啟用。
- **CR = 小應變大旋轉**(非大應變/非 GE Reissner)。既有 `assembleGeometric` 是 P-Delta 線性化非完整 CR
  (F50c 退化即此關係之證)。**S10 N-M 塑鉸必在 S9b 後**(R4)。

## ⚠️ 踩雷(durable)

1. **平面外 free DOF → LDLT singular**:caller 必須約束 `Uz,Rx,Ry`(fixture `cantileverPlanarTipShearN` 已約束)。
2. **`K_geom` 推導對稱性**:`p,q` 向量定義 + 二次微分驗證;P-Delta 退化 oracle 釘住軸力幾何項。
3. **UE Adaptive Build 自動排除新 .cpp 出 unity**(build log `[Adaptive Build] Excluded ... CorotationalAnalysis.cpp`)
   → 匿名 namespace helper 安全;UE 測試未用 `IN/OUT`(SAL 巨集)。
4. **Bash→powershell 路徑反斜線跳脫**;UE build 用 PowerShell 前景 Build.bat;跑 Tools/git 先 `cd /e/project/ArchSim`。
5. **NR 收斂解只取決於內力 f_int(殘差),切線只影響收斂速度** — S9b 3D 可先用主項切線求位移正確,再補完整 Kgeom。

## 對抗式審核(2026-06-12,4 agent 平行)

S9 五腿全綠後做 4-agent 平行對抗審核(各錨定真實證據:讀碼行號 + 跑既有 exe + numpy/sympy 獨立重算)。
**核心數學經獨立確認正確**:B 矩陣(3×6)符號逐項一致、`K_geom` sympy 驗證對稱 + 有限差分 6.9e-10、
elastica 表獨立 shooting 驗到 5e-9(確認**非 FrameCore 自產**)、剛體零內力機器精度。findings **全修**:

- **[CRITICAL→修] 平面外幾何靜默算錯**:driver 只用節點 x,y 算弦角/長,有 z 延伸的桿被當 xy 投影
  (45° 傾斜 → 65% 誤差)無警告。→ 加**平面幾何守門**(active 桿 `|dz|>1e-9·L3D` 即 reject)。CLI 實測
  現回 `COROT 0 0 0 0`+`SINGULAR 1`+歸零 DISP(拒絕,非靜默)。
- **[CRITICAL→修] 守門/diverged 路徑 segfault**:早返回路徑未配置 `finalState.u/reactions`,frame_cli
  `printState` 越界。→ `reject`/`fillPartial` helper 在所有早返回前零填 `u/reactions/memberForces`。CLI 實測
  殼+COROT 現 exit 0 無 crash。
- **[MAJOR→修] 塑鉸靜默忽略**:`model.hinges` 未檢查,含已形成塑鉸的模型以全剛度算。→ 加塑鉸守門
  (與 UDL 拒絕政策一致;S10 在 S9b 後)。
- **[MAJOR→修] audit P-Delta relErr 誤導**:audit `relErr` 分母 `max(1,|b|)` 對 sway~0.047 退化成絕對
  誤差(印 2.4e-4),真相對誤差 5.0e-3。→ check 內顯式算 `|a-b|/|b|`(現印 5.0e-3,tol 1.5e-2)。
- **[MINOR→修] 端力/反力無 oracle**:elastica 只測位移。→ audit **+check 5 小位移平衡**(stiff 懸臂
  base `Ry=-P` 精確 + member-0 `|Mz|=P*L`,rel **3.5e-7**;驗 recover 反力 + 端彎矩)。audit 94→**95**。
- **[MINOR→修] shear 用 L0**:`V=(Mi+Mj)/L0` 改用當前弦長 `Ln`(可伸縮截面一致性)。
- **[MINOR→修] Ln=0 守門**:弦長 `max(1e-300,Ln)` 防全塌除零。
- **[OK]** FrameCore 純度(POD API、Eigen 走 `FrameEigen.h`)、build/gate 同步(四腳本 + `ExpectedUeTests=47`
  + audit 註冊)、const-correct、不碰 `IElement`/既有 element/`SolveResult`/fingerprint —— 全乾淨。
- **[誠實確認] commit 衛生**:`.gitignore`/`ArchSim.uproject`/`Plugins/LevelSim/`/`frame_capi.*` 是**禁碰
  的既有未 commit 項**(非本輪改動),S9 commit 顯式 `git add` 排除。

審核後重跑五腿 gate 全綠(standalone F50 / UE 47 / OpenSees / audit **95** / CLI 9),**零 CRITICAL 殘留**。

## 下一步

- **S9b**:3D 通用 co-rotational(扭轉/雙軸/spin/SO(3)),WS_F2 OpenSees CorotCrdTransf3d 逐行公式已備;
  oracle 增 3D 剛體旋轉 patch + 空間 elastica + OpenSees corot 對照 + 平面退化回歸。
- **S10**:N-M 互動塑鉸(必在 S9b 後,R4 方向耦合);WS_G 路線 G1 `Mp_eff=Mp·f(N/Np)`。
- **S11**:MITC9i 高階殼(殿後)。提示詞 `docs/AGENT_PROMPT_S5_S11.md`。
- 平面 CR 後續:member UDL(follower / 初始構型等效)、prescribed 大位移、弧長法(snap-through)。
