# S9b 進度 — 3D 通用 Co-rotational 大位移(扭轉 + 雙軸彎曲 + spin + SO(3) 有限旋轉)

> 接續 `PROGRESS_S9.md`(S9 平面 v1 `efd006f`)。規格 `docs/specs/S9b_corotational3d.md`,研究依據
> `docs/research/WS_F2_corot3d_opensees.md`(OpenSees CorotCrdTransf3d 逐行)+ `WS_F_corot.md`(CR 路線)。

## 摘要

S9b 把 `runCorotational` 從**平面**升級為 **3D 通用** co-rotational:任意空間取向的梁柱(扭轉 + 繞兩截面
主軸雙軸彎曲 + 有限 SO(3) 旋轉)。**獨立 driver,不改 `IElement`、不動既有 element / dispatch /
`SolveResult` / `modelFingerprint`、API 簽名零變更**(仿 S9/`runPDelta`/`runProgressiveCollapse` 工作副本)。
**不新增 `.cpp`**(在 `CorotationalAnalysis.cpp` 內把 `CrBeam`→`CrBeam3D`)→ 四 build 腳本源檔清單免動。

- **S9b vs S9 的根本差異**:3D 有限旋轉**不可加** → 每節點維護 `R_node∈SO(3)`(初始 `I`),NR 解出旋轉
  增量後**空間左乘** `R_node ← exp(skew(Δθ))·R_node`(避 total rotation vector 在 2π 奇異)。取代 S9
  「平面繞單軸 Rz 直接累加」。
- **架構**:每根 `CrBeam3D` 持 12-DOF `gmap`、`L0`、初始框架 `E0col=localAxes(pi,pj,refVec)ᵀ`、natural
  係數 `EA/EIy/EIz/GJ`。每迭代:當前弦 `e1` + 測地平均三元組 Gram-Schmidt → 共旋框架 `E`;局部變形
  `θ̄=logSO3(Eᵀ·q_node)`、`ū=Ln−L0`;natural 力 `pb=Kl·v`;內力 `f=Tᵀpb`;切線 `Kt=TᵀKlT+Ksigma1`。
- **切線(誠實)**:`TᵀKlT + Ksigma1`(嚴格軸力幾何項,reduces to S9 的 `(N/Ln)q qᵀ`,對稱 → 續用
  `SimplicialLDLT`)。**完整 spin/moment 修正 Ksigma2/3 未加**(只影響收斂速度,NR 收斂解只取決於內力
  `f_int`(完整正確);主項已收斂 elastica α=1..10)→ 留 S9c。
- **守門(沿用 S9 + 解除平面外)**:殼 / UDL / prescribed / 塑鉸 / **release** / **tension-only** 皆 reject
  (零填 `finalState.u`);**S9 的平面外 z 守門已解除**(支援任意空間取向);極限點 `K_T` 非正定 → `diverged`。

**五腿全綠**:standalone **F1-F51** / UE **48**(+1 `Corotational.SpatialElastica`)/ OpenSees PASS(+corot
leg)/ audit **98**(+3)/ CLI round-trip **10**(+1 3D)。`$ExpectedUeTests` 47→48;audit 95→98;CLI 9→10。

## 數學(3D co-rotational,逐行對齊 OpenSees CorotCrdTransf3d)

- **SO(3) helper**(手寫 `Mat3`,`.cpp` 匿名 ns,不碰 `FrameEigen.h`):`skew`、`expSO3`(Rodrigues,θ→0
  級數)、`logSO3`(θ→0 級數 / θ→π 從對稱部分 `(R+I)/2` 提軸 + `veeAsym` 定符號)。
- **框架 `E=[e1|e2|e3]`**:`e1`=當前弦;`Rbar=expSO3(½·logSO3(qJ·qIᵀ))·qI`(測地中點);`e2,e3` = `Rbar`
  列對 `e1` Gram-Schmidt(分母 `1+r1·e1`)。`q_{I,J}=Rnode·E0col`。
- **局部變形**:`θ̄_{I,J}=logSO3(Eᵀ·q_{I,J})` 取分量(繞 `e1`=扭、`e2`=My/Iy、`e3`=Mz/Iz);`ū=Ln−L0`;
  扭 `θ̄x=θ̄x_J−θ̄x_I`。natural `Kl`(6×6):`EA/L0`、`(EIz/L0)[4,2;2,4]`、`(EIy/L0)[4,2;2,4]`、`GJ/L0`。
- **T(6×12)** basic `[ū,θ̄z_I,θ̄z_J,θ̄y_I,θ̄y_J,θ̄x]`:`ū`=`[-e1,0,e1,0]`;`θ̄z_I`=`[e2/Ln,e3,-e2/Ln,0]`;
  `θ̄y_I`=`[-e3/Ln,e2,e3/Ln,0]`;`θ̄x`=`[0,-e1,0,e1]`(reduces 逐項到 S9 平面 B)。`f=Tᵀpb`。
- **Ksigma1**:`N·[[A,0,-A,0],…]`,`A=(I-e1 e1ᵀ)/Ln`(reduces 到 S9 `(N/Ln)q qᵀ`)。
- **剛體零內力**:`R_g` 下 `q_{I,J}=R_g·E0col`、`E=R_g·E0col` → `Eᵀq=E0colᵀE0col=const` → `θ̄=0,ū=0`
  → `pb=0` → `f=Tᵀpb=0`(F51a 機器精度釘住)。

## Oracle(實測)

- **F51a 3D 任意軸剛體大旋轉不變性**:X 懸臂繞 `(1,1,1)` 軸轉 φ=2.0 rad,`|u_tip|` rel **4.10e-16**
  (機器精度;驗 expSO3/logSO3/框架/T 全鏈)。
- **F51b 空間 elastica**(傾斜軸 `(1,1,1)`,對稱截面,垂直端載):`α=1,5,10` dv/L vs Mattiasson rel
  **1.1e-4~6.5e-4**;α=5 得 dv/L=0.714138(與平面 F50 完全一致)。
- **F51c 平面退化(F50 在 3D driver 下)**:elastica α=5 dv/L=**0.714138**(與 S9 記載逐字一致)、in-plane
  剛體不變性 rel **0.00e+00**、P-Delta 退化 rel **5.02e-3** — 證明 3D 推廣未破壞平面正確性。
- **F51d 純扭 + 雙軸**:純扭 `θ=TL/GJ` rel **2.60e-16**;雙軸 `dy→Iy / dz→Iz` rel **~1e-9**(雙軸分離,
  容差 1e-6)。
- **F51e 大角度純扭**(審核補強):單元素 tip 扭矩使截面扭 **0.9π(162°)**,`θ=TL/GJ` rel **1.41e-15** —
  覆蓋 `logSO3` 大角度 general 端 + `Rnode` spatial 更新在有限旋轉的正確性(vs F51c 的小角度 1e-6)。
- **OpenSees corot leg**:3D 懸臂 α=3,我方(主項切線)vs OpenSees `geomTransf Corotational`(完整切線)
  tip `|u|` **ours-vs-OS=1.22e-9**(gate tol 收緊為 `TOL_COROT_VS_OS`=1e-6)— 完全獨立第三方,證明收斂解
  一致(切線只影響速度)。
- **audit testCorotational +3**:3D 旋轉不變性 **2.22e-16**、純扭 **0.0**、雙軸 **4.56e-14**;Check 5(小位移
  平衡)改讀**橫向總彎矩** `sqrt(My²+Mz²)`(3D driver 用 localAxes 慣例報 member 端力,frame-invariant)
  rel 3.47e-7。
- **CLI COROT 3D**:X 懸臂 +Z 載重全 3D 無平面約束,文字橋 round-trip dvZ=0.714138。

## 檔案改動

- `Private/CorotationalAnalysis.cpp`:`CrBeam`→`CrBeam3D`(6→12 DOF);新增 SO(3) helper;`crCompute3D`
  3D 框架/局部變形/`Tᵀpb`/`TᵀKlT+Ksigma1`;driver 加 `Rnode` 狀態 + spatial 旋轉更新;解除平面外 z 守門 +
  加 release/tension-only 守門。`Public/.../CorotationalAnalysis.h`:SCOPE 註解平面→3D(API 簽名不變)。
- `main.cpp` +F51(a/b/c/d);`FrameTestFixtures.h` +`cantileverSpatial`(3D fixture);
  `Tests/CorotationalTest.cpp` +`Corotational.SpatialElastica`;`run_gate.ps1` `$ExpectedUeTests` 47→48;
  `linear_deep_audit.cpp` testCorotational +3 checks(+ Check 5 改總彎矩);`cli_roundtrip.py` +3D COROT;
  `opensees_compare.py` +corot leg(`corot_cantilever_model`/`run_frame_cli_corot`/`run_opensees_corot`)。
- **不動**:`IElement.h` / 所有 element / `FrameSolver.cpp` / `SolveResult.h` / `SolveOptions.h` /
  `FrameModel.cpp`。**四 build 腳本源檔清單免動**(不新增 `.cpp`)。**S9b 未修改 `.gitignore` /
  `ArchSim.uproject` / `Plugins/LevelSim/` / `frame_capi.{dll,exp,lib}`**(這些是 working tree 既有的
  **前序未 commit 項**,非本輪改動;commit 時顯式 `git add` 源碼、勿 `-A` 以排除它們)。

## 誠實邊界 / novelty

- **CR = Crisfield/Battini/OpenSees 文獻方法**,FrameCore 為其**實作**,非新方法。`[NEW CODE]` 僅限 driver
  整合 + SO(3) 節點旋轉維護 + S9→S9b 推廣。**禁裸宣稱「優於 Karamba」**(3D CR 是業界標準)。
- **切線只做主項 `Ksigma1`**(嚴格軸力幾何);完整 spin/moment `Ksigma2/3` 未加(收斂速度,非正確性)→ S9c。
- **無 snap-through**(極限點 `diverged`);**只吃節點力**(UDL/prescribed/塑鉸/release/tonly reject)。
- **小應變大旋轉**(非大應變、非 GE Reissner;小應變下與 GE 收斂同解)。**S10 N-M 塑鉸必在 S9b 後**(R4)。

## ⚠️ 踩雷(durable)

1. **`localAxes` 慣例**:`z=x×ref`、`y=z×x`。對 X 軸 + refVec=(0,0,1):local y=+Z_global、local z=-Y_global
   → global-Y 端載彎曲用 **Iy**(非直覺;F51d/audit Check 8 oracle 期望值與 member 端力分量都要照此)。
2. **平面退化 = F50 在 3D driver 下仍綠**:3D driver 跑平面 fixture,elastica/剛體不變性/P-Delta 數值與 S9
   記載逐字一致(0.714138 / 0.0 / 5.02e-3)→ 是「3D 未破壞平面」的硬保護。先讓 F50 過再開 3D oracle。
3. **主項切線足夠**:`TᵀKlT+Ksigma1`(對稱續用 LDLT)收斂 elastica α=10、OpenSees 對照 1.2e-9。完整
   Ksigma2/3 非必要(只加速)。**NR 收斂解只取決於內力 `f_int`**(S9 durable,S9b 再證)。
4. **不新增 `.cpp`** → 四 build 腳本免動(`CorotationalAnalysis.cpp` 早在清單)。
5. **audit Check 5** 在 3D driver 下原讀 `.Mz` 失效(平面 Uy 載重彎矩落 `.My`)→ 改讀總橫向彎矩
   `sqrt(My²+Mz²)`(frame-invariant)。
6. **SO(3) θ→π 分支**:`logSO3` 在 θ→π(`sinθ≈0`)用 `(R+I)/2` 對稱部分提軸 + `veeAsym` 定符號;分支
   切換閾值 **`kPi-1e-5`**(審核後從 1e-6 拓寬,拉平過渡)。**已知極窄帶**:θ∈(π-2e-6, π-1e-6) 旋轉向量
   表示精度退化至 ~1e-4(3-agent 共識),但**旋轉向量僅作輸出,所有計算用的 `R` 由 `expSO3` 產生(everywhere
   機器精度)故無害**;大角度路徑(0.9π=162°)由 **F51e 機器精度(1.41e-15)覆蓋**,θ→π 極窄帶無任何
   oracle / 工程旋轉觸及(誠實標)。
7. **UE build 用 PowerShell 前景 Build.bat**(背景增量會漏編新測試致 gate 假綠);UE 測試未用 `IN/OUT`。

## 對抗式審核(2026-06-12,3-agent 平行)

S9b 五腿全綠後做 **3-agent 平行對抗審核**(數學核心 / 誠實性與衛生 / oracle 強度),各錨定真實證據
(讀碼行號 + 跑既有 `frametest.exe` + numpy/scipy/sympy 獨立重算),不跑 build。**零 CRITICAL、零數學
MAJOR**;findings 全修後重跑五腿仍全綠。

**數學核心經三方獨立確認正確**:T 矩陣(6×12)平面退化**與 S9 B0/B1/B2 完全 exact match(誤差 0)**;
剛體零內力 3.2e-14;`expSO3` vs scipy 3e-16;`logSO3` general roundtrip 1.4e-12;測地中點等距 2e-14;
Gram-Schmidt 正交 4e-16;Kl 雙軸解耦 3e-12;spatial 左乘 + `E0col` 慣例正確。Mattiasson 表用 scipy
DOP853 shooting 獨立對拍 6e-11(確認非自產);OpenSees 1.22e-9 可信(獨立程式庫 + 唯一平衡)。

**findings(全修)**:
- **[MAJOR→修] corot leg 容差**:原用 `TOL_PDELTA`(1e-2)比實測 1.22e-9 鬆 7 個量級 → 改獨立
  `TOL_COROT_VS_OS`(strict 1e-6),spec/PROGRESS 容差宣稱對齊。
- **[MAJOR→修] spec F51c 措辭**:spec 寫「F51c vs S9 solver 逐位元<1e-10」,但 S9 solver 已被 3D 取代
  (同一 `runCorotational`)→ 改 spec:平面退化由 **F50 保留**(3D driver 跑平面 fixture,elastica/剛體
  不變性/P-Delta 與 S9 記載 + Mattiasson 表逐字一致)驗證。
- **[MAJOR→修] PROGRESS「未碰」措辭**:`.gitignore`/`uproject` 是 working tree 前序 modified(非 S9b
  改)→ 改措辭明確「S9b 未修改但 working tree 既有前序未 commit 項,commit 顯式 add 排除」。
- **[MINOR→修] logSO3 θ→π 窄帶**:閾值 1e-6→1e-5 拓寬 + 加 **F51e 大角度扭轉(0.9π,rel 1.41e-15)**
  覆蓋大角度路徑;極窄帶誠實標(僅輸出、`R` 由 expSO3 精確、無害)。
- **[MINOR→修] F51d 容差**:1e-3 vs 實測 1e-9 脫節 → 收緊 1e-6。
- **[MINOR→修] spec 措辭**:`skewMat`→`skew`(對齊 cpp);natural Kl 比對降級為「由 F51b(EI)/F51c(GJ)/
  F51d/F50(EA)間接驗各係數」(spec 原列為建議,未獨立實作)。
- **[OK]** 三方確認乾淨:FrameCore 純度(POD API、Eigen 走 `FrameEigen.h`、零洩漏、不碰 `IElement`/
  element/`SolveResult`/fingerprint)、守門完整(殼/UDL/prescribed/塑鉸/release/tonly reject + 早返回零填)、
  build/gate 同步(四腳本含 `CorotationalAnalysis.cpp`、`ExpectedUeTests=48`)、commit 衛生、無裸宣稱優於 Karamba。

審核後重跑五腿 gate 全綠(standalone F1-F51 / UE 48 / OpenSees PASS / audit 98 / CLI 10),**零 CRITICAL 殘留**。

## 下一步

- **S9c**:弧長法(Riks/Crisfield)→ snap-through;完整 Ksigma2/3 spin/moment 切線;CR member UDL +
  prescribed 大位移。
- **S10**:N-M 互動塑鉸(**必在 S9b 後**,R4;塑鉸在局部共旋座標定義,`pb` 直接餵 N-M 面)。
- **S11**:MITC9i 高階殼(殿後)。
