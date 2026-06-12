# S9b 交接規格 — 3D 通用 Co-rotational 大位移(扭轉 + 雙軸彎曲 + spin + SO(3) 有限旋轉)

> 研究依據:`docs/research/WS_F2_corot3d_opensees.md`(★ 核心:OpenSees `CorotCrdTransf3d.cpp` 逐行公式 +
> Crisfield 1990 + Battini 2002)+ `docs/research/WS_F_corot.md`(CR 全路線)+ `docs/specs/S9_corotational.md`
> (平面 v1,直接前身)+ `docs/PROGRESS_S9.md`(平面審核 + durable 踩雷)。
>
> **狀態:實作完成,五腿全綠**(standalone F1-F51 / UE 48 / OpenSees PASS +corot leg / audit 98 / CLI 10)。
> 直接前身 = S9 平面 v1(`efd006f`)。進度與實測 oracle 見 `docs/PROGRESS_S9b.md`。**實作相對本 spec 的
> 一處誠實偏離**:切線只做嚴格主項 `TᵀKlT+Ksigma1`(§④(6) 的「主項先」路徑;完整 spin/moment `Ksigma2/3`
> 留 S9c — 主項已收斂 elastica α=1..10 且 OpenSees 對照 1.2e-9)。其餘 ①-⑩ 與實作相符。
>
> **誠實定位(本階段最重要)**:CR = Crisfield / Battini / OpenSees 的**成熟文獻方法**,FrameCore 為其
> **實作**,非新方法。S9b = S9 平面 CR 的 3D 推廣,屬業界標準附先行技術定位(§⑨)。**禁裸宣稱「優於
> Karamba」**。初版誠實邊界:load-controlled NR、無 snap-through(極限點 `diverged`)、小應變大旋轉
> (非 GE Reissner)、只吃節點力(UDL/prescribed/塑鉸/release/tension-only 不啟用)。

---

## ① 目標 / 不做

**做** — 把 `runCorotational` 從平面升級為 **3D 通用** co-rotational:任意空間取向的梁柱(扭轉 + 繞兩個
截面主軸的雙軸彎曲),每根梁隨其**當前弦 + 當前截面取向**共旋,局部仍小應變(small-strain large-rotation
CR);局部勁度用既有 3D Euler-Bernoulli 梁的同款係數(`EA / EIy / EIz / GJ`,源自 `localStiffness12`)。
求解 = **Newton-Raphson + load stepping**(載重因子 λ:0→1 分步,每步 NR 迭代到平衡),與 S9 同框架。

**S9b 與 S9 的根本差異(務必理解)**:
- **S9(平面)**:所有旋轉繞單一 global Z 軸 → **可交換可加** → 節點轉角直接由 `Rz` 位移 DOF 累積,
  **不需** SO(3)。
- **S9b(3D)**:3D 有限旋轉**不可加**(`R(θ₁)R(θ₂) ≠ R(θ₁+θ₂)`,SO(3) 非交換)→ 每節點必須維護一個
  旋轉矩陣 `R_node ∈ SO(3)`(初始 `I`),NR 解出旋轉增量 `Δθ` 後**空間左乘更新**
  `R_node ← exp(skew(Δθ))·R_node`(spatial incremental,避 total rotation vector 在 `|θ|→2π` 奇異
  [Battini §2.2])。位移向量的旋轉 DOF 在每個 NR 增量代表**旋轉增量**,被吸收進 `R_node`;「累積旋轉」
  存在 `R_node` 裡,不在位移向量。

**不做(誠實邊界,§⑨)**:
- **snap-through / 極限點後路徑**:load-controlled NR 在極限點 `det(K_T)=0` 發散 → **誠實偵測為
  `diverged`、不靜默算錯**;弧長法(Riks/Crisfield)留 S9c。主 oracle(端載懸臂 elastica)單調無極限點。
- **member UDL / prescribed 支承位移 / 殼 / grillage / 已形成塑鉸 / release / tension-only**:driver
  **明確拒絕**(回 `singular`+diagnostic,零填 `finalState.u`,不靜默丟載)。v1→v2 只吃**節點力**。
- **動力大位移、接觸、大應變**(非 geometrically-exact Reissner;小應變假設下與 GE 梁收斂同解,§⑨)。
- **Timoshenko 剪變形**:`localStiffness12T` 暫不接(`CorotationalOptions::solve.useTimoshenko` 保留);
  natural 勁度走 Euler-Bernoulli(與 S9 一致)。

---

## ② 公開 API(`Public/FrameCore/CorotationalAnalysis.h`,POD/std、零 Eigen)

**API 結構不變**(`CorotationalOptions` / `CorotationalResult` / `runCorotational` 簽名與 S9 完全相同):

```cpp
struct CorotationalOptions {
    int  loadSteps = 10;     // equal lambda increments 0->1 (3D elastica / large twist may want more)
    int  maxIter   = 50;     // NR iterations per load step
    real tolR      = 1e-9;   // ||residual_free|| / ||lambda*F_ext_free||      (force residual)
    real tolU      = 1e-12;  // ||du_free|| / max(1,||u_free||)                (displacement)
    SolveOptions solve;      // pivotTol passthrough (useTimoshenko reserved; v2 is Euler-Bernoulli)
};
struct CorotationalResult {
    bool        converged = false;
    bool        diverged  = false;   // limit point: K_T not positive-definite (snap-through; needs arc-length)
    int         loadStepsCompleted = 0;
    int         totalIterations    = 0;
    real        lastResidual       = 0;
    SolveResult finalState;          // 3D large-displacement state (u, reactions, member local end forces)
};
FRAMECORE_API CorotationalResult runCorotational(const FrameModel& model, const CorotationalOptions& opts = {});
```

**無新 model / SolveResult 欄位、不動 `modelFingerprint`**(`CorotationalOptions` 是 solve 參數,同
`PDeltaOptions`/S9 列)。**API 簽名零變更** → 既有呼叫端(main.cpp F50、CLI COROT、UE 測試)源碼不動,
僅行為從「平面」擴為「3D」。`finalState.u` 的旋轉分量(`Rx,Ry,Rz`)= **節點總旋轉向量**
`logSO3(R_node)`(誠實標:3D 大旋轉下總旋轉向量非簡單可加,僅作輸出;內部狀態是 `R_node` 矩陣)。
`finalState.memberForces` 的 `endI/endJ` = 局部(共旋)座標系的 6 分量端力(N, Vy, Vz, Mx/torsion, My, Mz)。

---

## ③ 資料流 / 架構(獨立 driver,**不改 `IElement`**;沿用 S9)

侵入性結論**與 S9 相同**:獨立 NR driver,`IElement` 與所有既有 element **零改動**。理由不變:`IElement`
幾何凍結於初始(`prepare` 一次算 `kl_/T_/Qf_`),CR 每迭代要在變形構型重算框架;改介面會迫使所有 element
(含殼)實作用不到的 hook 並威脅既有逐位元行為。仿 `runPDelta`/`runProgressiveCollapse`/S9 的「獨立
driver + const 工作副本 + 自組矩陣 + recover 對齊 `SolveResult`」模式。

**driver 狀態(工作副本,const model 不變,安全並行)**:
- 每根 active 桿一個輕量 **`CrBeam3D`**(**非** `IElement`,取代 S9 的平面 `CrBeam`):
  - `e, id, ni, nj`、`gmap[12]`(12 個全域 DOF:`[uxi,uyi,uzi,rxi,ryi,rzi, uxj,uyj,uzj,rxj,ryj,rzj]`,
    對齊 `localStiffness12` 的 canonical 序);
  - `L0`(初始弦長)、`E0`(`Mat3`,**初始**元素框架 = `localAxes(pi, pj, mem.refVec)`,定義截面主軸
    `Iy/Iz` 的取向);
  - natural 勁度係數 `EA, EIy, EIz, GJ`(源自 `localStiffness12` 同款 `mat.E/mat.G + sec.A/Iy/Iz/J`);
  - recover 狀態:`pb`(natural 力 6×1)、`Ln`、當前框架 `E`。
- **每節點 `R_node ∈ SO(3)`**:`std::vector<Mat3> Rnode(nNodes, Mat3::Identity())`(driver 全域狀態,
  S9 沒有的新增)。NR 增量的旋轉分量 `Δθ_k`(node k 的 `Rx,Ry,Rz` 增量)→ `Rnode[k] ← expSO3(Δθ_k)·Rnode[k]`。
- 全域位移 `u[6N]`:**平移分量累加**(同 S9);**旋轉分量**每步存「自上次 commit 的旋轉增量」並即時清零
  到 `R_node`(即 NR 內每次 `u[rot]` 只暫存當前迭代增量,用後更新 `R_node` 即歸零;total rotation 由
  `R_node` 持有)。recover 時 `u[rot] = logSO3(R_node)`。

**拒絕(reject helper,皆零填 `finalState.u/reactions/memberForces` 防消費端越界,沿用 S9)**:active 殼、
active 桿上的 UDL、非零 prescribed 位移、**已形成塑鉸**、release(`enableReleases` 下任何 active 桿 release)、
tension-only 桿。皆 `singular`+diagnostic(§⑨)。**解除 S9 的平面外 z 守門**(改為支援任意取向);
**移除「桿須在 XY 平面」限制**。

**一次 NR 迭代**(§④ 詳述):對每根 active 桿,用當前節點座標 `(Xi+ui, Xj+uj)` 與當前 `R_node[ni]/R_node[nj]`
算當前框架 `E`、局部變形 `[ul(6); ū]` → natural 力 `pb` → 內力 `f_int=Tᵀpb`、切線 `Kt` → scatter 到 12-DOF
→ residual `r=λF−f_int` → reduceFF → LDLT(`vectorD` 測正定;非正定=極限點→`diverged`)→ `Δu` →
平移累加 + 旋轉 `R_node` 左乘 → 收斂判據。load step 收斂後推進 λ;recover `finalState`。

**接點複用**:`localAxes(pi,pj,refVec)`(初始框架 `E0`,**既有已驗證**,守垂直桿退化)、`localStiffness12`
的 `EA/EIy/EIz/GJ` 係數、`reduceFF`/LDLT/`ldltPosDef`(同 S9/`PDeltaAnalysis.cpp`)、`gdof`/`nodeIndex`/
`dofCount`。**新增**:`skewMat`/`expSO3`/`logSO3` SO(3) helper(§④、§⑤)。

---

## ④ 演算法(3D co-rotational,逐行對齊 OpenSees `CorotCrdTransf3d` + WS_F2)

> **權威來源 = `docs/research/WS_F2_corot3d_opensees.md`**(錨定 OpenSees `CorotCrdTransf3d.cpp` 行號)。
> 本節給資料流 + 關鍵公式;**實作時符號以 OpenSees 源碼逐行為準**(尤其 `getLMatrix` / `compTransfMatrixBasicGlobal`
> / `getKs2Matrix`)。下標 `I,J`=端節點;`E=[e1|e2|e3]`(共旋元素框架,列向量);`R_I,R_J`=節點旋轉。

**(0) SO(3) 基元**(手寫,`Mat3`,§⑤ 放 `.cpp` 匿名 ns,θ→0 級數避除零;實作函式名 `skew`/`expSO3`/
`logSO3` + `veeAsym`):
- `skew(v) = [[0,-v2,v1],[v2,0,-v0],[-v1,v0,0]]`。`veeAsym(R)=[R21-R12, R02-R20, R10-R01]`(= 2sinθ·axis)。
- `expSO3(θ⃗)`(Rodrigues):`θ=|θ⃗|`;`exp = I + (sinθ/θ)S + ((1−cosθ)/θ²)S²`,`S=skew(θ⃗)`;
  θ→0 級數 `sinθ/θ≈1−θ²/6`、`(1−cosθ)/θ²≈½−θ²/24`(ε<1e-7 切換)。
- `logSO3(R)`:`θ=arccos(clamp((tr R−1)/2,−1,1))`;`log = θ/(2sinθ)·veeAsym(R)`;θ→0 取
  `½·veeAsym(R)`(級數);θ→π(`θ≥kPi-1e-5`,審核後從 1e-6 拓寬以拉平過渡)用「從 `R+I` 對稱部分提軸」
  備援(最大對角列正規化,符號由 `veeAsym(R)` 定)。極窄帶 (π-2e-6,π-1e-6) ~1e-4 表示精度退化已知無害
  (旋轉向量僅輸出,計算用的 `R` 永遠由 `expSO3` 精確產生;F51e 0.9π 覆蓋大角度路徑)。

**(1) 旋轉更新(spatial incremental — 避 2π 奇異)**:NR 步末,對每節點 `Rnode[k] ← expSO3(Δθ_k)·Rnode[k]`
(`Δθ_k` = 該 NR 步 node k 的 `[Rx,Ry,Rz]` 位移增量)。**先更新 `R_node` 再進下一迭代算框架**(順序敏感)。

**(2) 元素框架 `E=[e1|e2|e3]`**(OpenSees `update()`):
- `e1 = (xJ+uJ − xI−uI)/Ln`(當前弦,`Ln`=當前弦長)。
- 平均三元組 `Rbar`(`R_I→R_J` 測地中點):`Rbar = expSO3(½·logSO3(R_J·R_Iᵀ))·R_I`(等價 WS_F2 §2 的
  四元數半角;手寫 exp/log 即可,免引入 Quaternion)。`R_I = Rnode[ni]·E0`、`R_J = Rnode[nj]·E0`
  (節點旋轉作用在初始框架)。
- `r1,r2,r3` = `Rbar` 各列;Gram-Schmidt 對 `e1`:
  `e2 = r2 − ((r2·e1)/(1+r1·e1))(e1+r1)`、`e3 = r3 − ((r3·e1)/(1+r1·e1))(e1+r1)`,
  再正規化。**分母 `1+r1·e1≈0` = 翻轉 180°**,需分支備援(§⑩;主 oracle 不觸發)。

**(3) 局部變形 `[ul(6); ū]`**(OpenSees `update()`,basic 維度 6):
`ū = Ln − L0`(natural 軸向)。6 個 natural 旋轉(扭@I/J、繞 e2/e3 彎@I/J):
`ul[k] = asin( ½(r_α^T·skew(e_β)·e_γ − r_β^T·skew(e_γ)·e_α) )`(α,β,γ 按各 DOF 循環,`r`=R_I/R_J 列)。
> basic 排序採 OpenSees 標準 6-DOF:`[ū, θ̄z_I, θ̄y_I, θ̄z_J, θ̄y_J, θ̄x]`(扭轉 `θ̄x` 為相對扭轉
> `θ̄x_J−θ̄x_I`;WS_F2 §3 的「7-DOF 分離扭轉」是等價記法)。**實作時以源碼 basic 序 + 符號為準**。
> **剛體零內力**:剛體 `R_g` 下 `R_{I,J}=R_g·E0`、`E=R_g·E0`,叉積中 `R_g` 抵消 → `ul=0, ū=0`(§⑥ F51a
> 機器精度釘住)。

**(4) natural 勁度 `Kl`(6×6,EB)+ natural 力 `pb`**:
`Kl` = blockdiag 對角塊(與 `localStiffness12` 同源係數,逐位元對齊既有梁):
- 軸:`EA/L0`(對 `ū`);扭:`GJ/L0`(對 `θ̄x`);
- 繞 e3 彎(用 `Iz`):`(EIz/L0)·[[4,2],[2,4]]`(對 `[θ̄z_I,θ̄z_J]`);
- 繞 e2 彎(用 `Iy`):`(EIy/L0)·[[4,2],[2,4]]`(對 `[θ̄y_I,θ̄y_J]`)。
> **複用 `localStiffness12` 的兩條等價路徑**(spec 容許實作擇一,audit 比對逐位元):(a) 閉式直接填上述
> 6×6(零依賴、最簡);(b) 呼叫 `localStiffness12(E,G,A,Iy,Iz,J,L0)` 取 basic 子塊。**實作採 (a) 閉式**;
> natural 係數正確性由 F51b(EI 彎)/F51c+F51e(GJ 扭)/F51d(Iy/Iz 雙軸)/F50(EA via 軸)**間接驗證**
> (未獨立加 `Kl==localStiffness12` 逐位元 check — 上述 oracle 已釘各通道係數)。
`pb = Kl·[ul; ū]`(7×1 重排為 OpenSees basic 6×1;軸力 `N=pb_axial`,端彎矩/扭矩為其餘分量)。

**(5) 內力 `f_int`(12×1)**:`f_int = Tᵀ·pb`(虛功一致;`T`=basic→global 轉換矩陣 6×12,OpenSees
`compTransfMatrixBasicGlobal()` + `getLMatrix()`)。`A=(I−e1 e1ᵀ)/Ln`;`getLMatrix(ri)`:
`L1=(ri·e1/2)A + A ri (e1+r1)ᵀ/2`、`L2=skew(ri)/2 − (ri·e1/4)skew(r1) − skew(ri)e1(e1+r1)ᵀ/4`;
`T` 各行由 `A·r{I,J}{2,3}`、`L·r2/L·r3` 與 `/(2cosθ)` 組裝(軸行 = `[−e1ᵀ,0,e1ᵀ,0]`)。**符號逐行對齊源碼**。

**(6) 一致切線 `Kt`(12×12)**(OpenSees `getGlobalStiffMatrix()`):
`Kt = Tᵀ·Kl·T + Ksigma1 + Ksigma2 + Ksigma3 + Σ_j m_j·tan(ul_j)·T_jᵀT_j`。
- `Ksigma1 = N·[[A,0,−A,0],[0,0,0,0],[−A,0,A,0],[0,0,0,0]]`(軸力幾何,**必加**;對稱)。
- `Ksigma3`(彎矩-軸耦合,`L·r2/L·r3·m`)、`Ksigma2`(spin 修正,`getKs2Matrix`,**非對稱來源**)。
- **務實兩階段(WS_F2 §6 + S9 durable 教訓:NR 收斂解只取決於內力 `f_int`,切線只影響收斂速度)**:
  - **S9b 主項階段**:`Kt = Tᵀ Kl T + Ksigma1`(對稱 → **續用 `LDLTSolver`**)。**位移即正確**,F51 oracle 全過。
  - **S9b 完整階段(可後補)**:加 `Ksigma2/3`;非對稱 → **對稱化 `(Kt+Ktᵀ)/2`**(無集中外力矩時本就近
    對稱;WS_F2 §6 認可,小步長誤差可忽略,收斂解不變)續用 LDLT。誠實標:若初版只主項 → 收斂步數較多,
    位移不受影響(§⑨⑩)。

**(7) 平面退化(回歸保護的數學基礎)**:當所有桿在 XY 平面、僅繞 Z 旋轉時,`Rnode` 退化為繞 Z 的
`expSO3([0,0,θz])`、`E` 退化為平面框架、`ul` 退化為平面 `θ̄`、`Kl` 的繞-Z 彎塊 == S9 的 `EIz·[4,2;2,4]`、
`Ksigma1` == S9 的 `(N/Ln)q qᵀ` 軸力項 → **S9b 在平面案例應收斂到 S9 同解**(§⑥ F51c 釘住)。

---

## ⑤ 檔案

**改(不新增 `.cpp` → 免動四 build 腳本源檔清單,鐵律⑥)**:
- `Private/CorotationalAnalysis.cpp`:`CrBeam`→`CrBeam3D`(平面 6-DOF → 3D 12-DOF);新增匿名 ns SO(3)
  helper `skew/veeAsym/expSO3/logSO3`;`crCompute` 重寫為 3D 框架/局部變形/`Tᵀpb`/`Kt`;driver 加 `Rnode`
  狀態 + spatial 旋轉更新;**解除平面外 z 守門**,保留殼/UDL/prescribed/塑鉸/release/tonly 守門。
  **Eigen 一律走 `FrameEigen.h`**(`Mat3/Vec3e/Mat12/Vec12/VecX/SpMat/LDLTSolver` 皆已有;**若**最終需
  `Eigen::Quaternion/AngleAxis`(預期不需,手寫 SO(3) 足夠)才加進 `FrameEigen.h`,**絕不**在 `.cpp` 直接
  `#include <Eigen/...>`)。
- `Public/FrameCore/CorotationalAnalysis.h`:**API 簽名不變**,僅更新 SCOPE 註解(平面 v1 → 3D 通用;
  移除「caller 須約束平面外 DOF」「桿須在 XY 平面」字樣;補 3D 邊界:小應變大旋轉、無 snap-through、節點力)。

**Gate 同步**:
- `main.cpp`:**保留 F50**(平面 elastica/剛體不變性/P-Delta 退化 — 現走 3D driver = 平面退化回歸);
  **新增 F51**(a 任意軸剛體旋轉不變性 / b 空間 elastica / c 純扭 `θ=TL/GJ` / d 雙軸分離 / e 大角度扭轉 0.9π)。
- `FrameTestFixtures.h`:新增 3D fixture:`cantileverSpatialTipShearN`(傾斜軸懸臂,不約束平面外旋轉)、
  `rigidBodyRotationPatch`(任意小結構 + 任意軸大旋轉 helper)、`twistCantilever`(純扭)。
- `Private/Tests/CorotationalTest.cpp`(UE):新增 `FrameCore.Corotational.SpatialElastica`(或 3D patch)
  → `run_gate.ps1` `$ExpectedUeTests` **47→48**。**UE 測試常數勿用 `IN/OUT`**(SAL 巨集)。
- `linear_deep_audit.cpp`:`testCorotational` 擴 3D checks(3D 剛體 patch=0、平面退化==S9、純扭 `θ=TL/GJ`、
  natural `Kl` 對角塊==`localStiffness12`);audit **95→ ~99**(實作後定數,run_gate 讀 `checks=N` 非寫死)。
- `frame_cli_core.cpp` / `Tools/cli_roundtrip.py`:`COROT` 指令已支援 + `MEMBER` 已帶 `rx ry rz`(refVec)
  與 3D 節點 → 加一個**傾斜取向** COROT round-trip check(`cli_roundtrip` **9→10**);`docs/CLI_PROTOCOL.md`
  COROT 註記擴 3D。
- `Tools/opensees_compare.py`:加 **corot leg**(`geomTransf Corotational $tag $vecxz`,**只用節點力、勿
  `eleLoad`** — WS_F2 §7/F-5 OpenSees corotTransf 靜默忽略 element loads;3D 需 vecxz 對齊 refVec)。

**不動**:`IElement.h`、所有既有 element、`FrameSolver.cpp`、`SolveResult.h`/`SolveOptions.h`(無新欄位)、
`FrameModel.cpp`(validate/fingerprint)。**`.gitignore`/`ArchSim.uproject`/`Plugins/LevelSim/`/
`frame_capi.{dll,exp,lib}` 絕不碰**(禁碰的既有未 commit 項;commit 顯式 `git add` 排除)。

---

## ⑥ Oracle(誠實分級,每能力先有獨立 oracle 才宣稱)

- **F51a 任意軸剛體旋轉不變性 `[VERIFIED]`**(CR 定義性質):X 懸臂 + 載重整體繞**任意軸 `n̂`(非 Z)**剛體
  旋轉 `R_g(n̂,φ)`(`n̂=(1,1,1), φ=2.0 rad`)→ `|u_tip|` 不變(剛體下 `pb=0→f=0`;**直接驗 `expSO3`/
  `logSO3`/框架 `E`/`T` 全鏈**)。實測 rel **4.10e-16**(機器精度)。
- **F51b 空間 elastica `[VERIFIED]`**(獨立第三方解):沿傾斜軸 `(1,1,1)` 大撓度懸臂(對稱截面,垂直端載),
  `α=1,5,10` dv/L vs WS_F F-3 獨立 shooting 表。容差 N16 `δv` 2e-3 / `δh` 6e-3。實測 rel 1.1e-4~6.5e-4
  (α=5 dv/L=0.714138 與平面一致)。
- **F51c 純扭 `[VERIFIED]`**:X 懸臂 tip 扭矩,tip 扭轉 `θ=T·L/GJ`(線性扭轉)。實測 rel **2.60e-16**。
- **F51d 雙軸分離 `[VERIFIED]`**:X 懸臂同時 +Y/+Z 小端載,`dy→Iy / dz→Iz`(小位移 `δ=PL³/3EI`;localAxes
  慣例 global-Y→Iy、global-Z→Iz),驗兩主軸平面獨立。容差 1e-6,實測 ~1e-9。
- **F51e 大角度純扭 `[VERIFIED]`**(審核補強):單元素 tip 扭矩使截面扭 **0.9π(162°)**,`θ=T·L/GJ` —
  覆蓋 `logSO3` 大角度 general 端 + `Rnode` spatial 更新在有限旋轉的正確性。實測 rel **1.41e-15**。
- **平面退化回歸(保護既有 F50,頭等)**:S9 平面 solver 已被 3D `runCorotational` 取代(同一函式),故由
  **保留的 F50 三子測**驗證 — 3D driver 跑平面 fixture(`cantileverPlanarTipShearN`),elastica α=5
  dv/L=**0.714138**、in-plane invariance rel **0.0**、P-Delta 退化 rel **5.02e-3**,與 S9 記載 + Mattiasson
  表逐字一致。
- **OpenSees corot 對照 `[VERIFIED]`**(`opensees_compare.py` corot leg):3D 懸臂 α=3,`geomTransf
  Corotational` + vecxz,**只節點力**(勿 eleLoad — WS_F2 §7 corotTransf 靜默忽略 element loads)。我方(主項
  切線)vs OpenSees(完整切線)tip `|u|` rel **1.22e-9**;gate 容差獨立 `TOL_COROT_VS_OS`(strict **1e-6** /
  relaxed 1e-4)。
- **CLI COROT 3D round-trip**:`frame_cli` X 懸臂 +Z 載重全 3D 無平面約束(`MEMBER … rx ry rz` refVec),
  COROT round-trip dvZ=0.714138。
- **誠實**:elastica shooting 表 = 與 FrameCore 完全獨立的第三方解(scipy DOP853 對拍 6e-11);F51a/e =
  CR 框架不變性 + 大角度旋轉鐵證;平面退化 = 對既有已驗 S9 數值的回歸;OpenSees corot 對照補「非自產第三方」第二源。

---

## ⑦ Gate

新增 **F51**(a/b/c/d/e);UE `FrameCore.Corotational.SpatialElastica`(`$ExpectedUeTests` **47→48**);
audit `testCorotational` 擴 3D(95→**98**);`cli_roundtrip` COROT 3D(9→10);`opensees_compare.py`
corot leg(`TOL_COROT_VS_OS`)。**保留 F50**(現為 3D driver 的平面退化保護)。**免動四 build 腳本源檔清單**(不新增 `.cpp`)。
**commit 前五腿全綠**(`run_gate.ps1 -RequireOpenSees`):standalone F1-F51 / UE 48 / OpenSees PASS(+corot)/
audit(擴)/ CLI round-trip 10。

---

## ⑧ 效能

同 S9:CR 是迭代(load stepping × NR),每迭代 = 一次組裝 + 一次 fresh LDLT(切線每迭代變,無法
factorize-once)。3D 每桿 12×12(vs 平面 6×6)+ 每節點 SO(3) exp/log(O(1) Mat3 運算),常數倍增。額外
狀態 O(N)(`Rnode` 每節點一個 `Mat3`)。完整 `Ksigma2/3` 較主項貴但仍 O(桿數);初版主項切線收斂步數略多
(位移不變)。不參與 ReSolve / 模態 / 屈曲(CR 切線非定常)。

---

## ⑨ 誠實邊界 / novelty

- **CR = Crisfield / Battini / OpenSees 成熟文獻方法**,FrameCore 為其**實作**,非新方法。`[VERIFIED]` =
  文獻方法的 oracle 化(elastica / 3D 剛體 patch / 平面退化 / OpenSees 對照);`[NEW CODE]` 僅限「獨立 NR
  driver 整合 + SO(3) 節點旋轉維護 + S9→S9b 平面到 3D 推廣 + `SolveResult`/CLI 對接」。**對標禁裸宣稱
  「優於 Karamba」**,附先行技術定位(3D CR 是業界標準,OpenSees/SAP2000/Karamba 皆有;FrameCore 非獨有)。
- **小應變大旋轉**:CR = small local strain + large global rotation(**非大應變、非 geometrically-exact
  Reissner**;小應變下與 GE 梁收斂同解,WS_F F-8)。
- **無 snap-through**:極限點 `det(K_T)=0` → `diverged`(不靜默算錯);弧長法(Riks/Crisfield)留 S9c。
  主 oracle(端載懸臂 elastica)單調無極限點。
- **只吃節點力**:UDL / prescribed 位移 / 已形成塑鉸 / release / tension-only 皆**拒絕**(零填 `finalState.u`)。
- **完整 Kgeom 分階段(若初版只主項)**:`Ksigma2/3` 後補時誠實標「初版主項切線 → 收斂步數較多,位移正確」
  (NR 收斂解只取決於內力,§④(6))。**非對稱切線對稱化**是 WS_F2 §6 認可的工程處理,標明。
- **既有 `assembleGeometric` 是 P-Delta 線性化非完整 CR**(勿混淆);S9b CR 切線是其有限旋轉推廣(F51c
  平面退化 + S9 F50c P-Delta 退化即此關係之證)。
- **S10 N-M 塑鉸必在 S9b 後(R4 方向耦合)**:塑鉸在局部(共旋)座標系定義,不受大轉動影響,S9b 完成才接。

---

## ⑩ 風險 / fallback(實作後回填)

1. **SO(3) `logSO3` 在 θ→π 穩健性**:θ→0 用級數、θ→π 用「從 `R` 對稱部分提軸」備援。F51a 大旋轉
   patch(`φ` 可設 >π/2 甚至接近 π)專測此路徑。**fallback**:若手寫 θ→π 不夠穩,加 `Eigen::Quaternion`
   到 `FrameEigen.h`(`<Eigen/Geometry>`,header-only/MPL2)用四元數提取——**屆時才碰 `FrameEigen.h`**
   (鐵律②:Eigen 一律走 choke point)。預期手寫足夠(oracle 旋轉量可控)。
2. **平均框架 `1+r1·e1≈0`(翻轉 180°)**:Gram-Schmidt 分母奇異 → 分支備援(WS_F2 §8);主 oracle(懸臂
   大撓度、剛體 patch φ<π)不觸發。誠實標「桿端相對旋轉接近 180° 未涵蓋」。
3. **`Kt` 非對稱(完整 `Ksigma2`)vs `LDLTSolver`**:主項切線對稱(續用 LDLT);完整則對稱化
   `(Kt+Ktᵀ)/2`(§④(6))。**不換 LU**(保持與 S9/PDelta 同 LDLT 路徑,降風險)。
4. **平面退化回歸(F51c)必須過**:這是「3D 推廣未破壞既有 S9」的硬保護。實作策略:先讓 F51c 過(平面案例
   逐位元對 S9),再開 3D oracle。若 F51c 不過 → SO(3)/框架/basic 符號有錯,優先修。
5. **`asin` 出域(`|ul|>π/2`)**:小步長保 `|ul|<π/2`;出域 → 加大 `loadSteps`(非 bug)。守門:`asin`
   參數 `clamp(±1)`。
6. **`Ln=0` / 全塌**:弦長 `max(1e-300, Ln)`(沿用 S9)。
7. **基準正確性**:`E0 = localAxes(pi,pj,refVec)` 複用既有**已驗證**函式(守垂直桿退化);natural `Kl`
   對角塊 audit 比對 `localStiffness12`(防新係數 bug);F51c 對 S9 逐位元(防平面破壞);F51a 機器精度
   (防框架/SO(3) 錯);OpenSees corot leg(第二獨立源)。
8. **build/gate 同步**:不新增 `.cpp` → 四 build 腳本免動;`$ExpectedUeTests` 47→48;audit 自報 `checks=N`;
   `COROT` CLI 已存在向後相容。**UE Adaptive Build 自動排除新 cpp 出 unity**(本檔不新增 cpp,但若擴
   `CorotationalAnalysis.cpp` 內匿名 ns helper 安全)。
9. **Eigen 洩漏 / UE / SAL**:Eigen 走 `FrameEigen.h`(POD 公開 API);UE 測試未用 `IN/OUT`;UE build 用
   **PowerShell 前景 `Build.bat`**(背景增量會漏編新測試致 gate 假綠)。
10. **commit 衛生**:顯式 `git add` 源碼,**勿 `-A`**(`.gitignore`/`uproject`/`LevelSim`/`frame_capi.*`
    禁碰);只 commit `FrameCore/docs/Tools`。

---

## ⑪ S9b 後續(未實作)

- **S9c**:弧長法(Riks/Crisfield)→ snap-through / 極限點後路徑(Williams toggle);CR member UDL
  (follower / 初始構型等效)+ prescribed 大位移。
- **S10**:N-M 互動塑鉸(**必在 S9b 後**,R4 方向耦合);塑鉸在局部共旋座標定義,`Mp_eff=Mp·f(N/Np)`
  (WS_G 路線 G1);S9b 的局部端力 `pb` 直接餵 N-M 面。
- **S11**:MITC9i 高階殼(殿後)。提示詞 `docs/AGENT_PROMPT_S5_S11.md`。
