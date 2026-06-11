# S8 交接規格 — 殼路線(8a QM6 opt-in 不協調膜 + 8b DKQ 薄板快路)

> 研究依據:`docs/research/WS_E_shell.md`(DKQ / QM6 / MITC9i / Karamba 殼查證)+ 既有膜鎖定
> P5 探索(plan `moonlit-wishing-charm.md` §P5,**已用獨立 numpy 驗 QM6**)。
> **P5 實測數據(本 spec 的 oracle 不變量來源)**:
> - 膜畸變 patch:Q4=2.8e-16、**QM6=0**(QM6 需**中心 Jacobian** 使 ∫B_inc=0 才過畸變 patch)。
> - Cook's membrane tip:Q4 @N4/8/16 = −25.8% / −9.6% / −3.2%;**QM6 = −5.9% / −2.4% / −0.9%**。
> - 細長 in-plane 懸臂:Q4 @2×1 = **−92.6%(鎖死)** → QM6 −5.8%;4×1 = −75.6% → −1.1%。
> - flat-facet 幾何誤差純 **O(1/N²)**(p=2.00);**主因判定:膜鎖定主導**(N=12/16 膜鎖定 ~16-25%
>   vs flat-facet ~0.2-0.35%,差兩數量級)。「flat-facet 才是瓶頸」被數據推翻。
>
> **架構核心(本階段最重要的設計決定)**:8a + 8b **都是 `MITC4ShellElement` 內的 `SolveOptions`
> opt-in 分支**(`useIncompatibleMembrane` / `useDKQPlate`),**共用** 24-DOF facet 框架 / 組裝 /
> 3D 旋轉 / 質量 / 壓力等效載重 / recover 骨幹。**預設兩旗標 false = 原始 MITC4(Q4 雙線性膜 +
> assumed-shear 板)逐位元不變** → **OpenSees `ShellMITC4` 平板閘門(現 ~1e-10)不破**。這是 S8
> 最大風險的根治(WS_E / P5 反覆強調:`ShellMITC4` 是原始雙線性膜,改進膜對它要 ~1e-10 會誤判退步)。
>
> 起點錨點 `ec157a0`(S7),五腿全綠 standalone F1-F47 / UE 44 / OpenSees PASS / audit 82 / CLI 8。

---

## ① 目標 / 不做

**8a 做** — QM6 不協調膜(Wilson-Taylor-Doherty-Ghaboussi 1973 Q6 + Taylor-Beresford-Wilson 1976
QM6 修正),`SolveOptions.useIncompatibleMembrane`(opt-in)。在現有 `Bmembrane`/`membraneK` 的
plane-stress Q4 膜上加 2 個不協調 bubble 模式(φ₁=1−ξ², φ₂=1−η²,各作用於 u,v → 4 個元素內部 DOF),
元素層 static condensation 消除回 12×12 膜([u,v,thz]×4),**讓 4 節點膜精確表達純彎曲、根治 in-plane
bending 膜鎖定**。**drilling(thz / Hughes-Brezzi)完全不碰**(bubble 只增強 u,v 膜應變)。

**8b 做** — DKQ 離散 Kirchhoff 薄板彎曲(Batoz & Tahar 1982),`SolveOptions.useDKQPlate`(opt-in)。
並列一條 plate 路徑取代 MITC4 assumed-shear `plateK`:8 節點 serendipity 二次旋轉場 + 沿 4 邊強制
Kirchhoff 約束(零橫向剪)消去邊中旋轉 → 12-DOF([w,θx,θy]×4)純彎曲薄板勁度,**無橫向剪自由度**。
薄板專用快路(t/L < ~1/20),**非取代 MITC4**(混合薄厚 / 中厚板回 MITC4)。

**共同**:預設(旗標 false)= 原始 MITC4,逐位元一致(同 `useTimoshenko` 的「false 保 bit-for-bit」慣例);
`SolveOptions` **不進 `modelFingerprint`**(它是 solve 參數非 model 欄位,與 `useTimoshenko`/`enableReleases`
同列)→ **不動 fingerprint**;**不新增 element 類別 / 不改 dispatch**(`FrameSolver.cpp:125` 仍只建
`MITC4ShellElement`,prepare 已收 `opts`);**不新增 `.cpp` 檔**(全寫進 `MITC4ShellElement.cpp`)→
四個 build 腳本源檔清單免動。

**不做**:
- **完整曲殼幾何**(facet 路線保留;曲面誤差隨網格 O(1/N²) 收斂,改進膜不消此幾何 floor)。
- **殼幾何勁度 / 殼挫屈**(`assembleGeometric` 殼仍空,另列後續;S8 純線性勁度)。
- **QM6 的厚殼 / 板彎曲增強**(QM6 只增強**膜** in-plane;板彎曲已由 MITC4 assumed-shear 良好處理)。
- **DKQ 作主要殼元素**(WS_E「避什麼」明列:DKQ **完全無厚板能力**,結構工程常見中厚板 t/L~1/20~1/5
  須回 MITC4;DKQ 僅薄板快路 + Kirchhoff oracle)。
- **DKQ + QM6 同時開**的交叉驗證為**次要**(兩旗標可獨立開;DKQ 板 + QM6 膜的組合本階段允許但 oracle
  以各自單獨開為主,組合留誠實標)。
- **MITC9i**(殿後 S11);**per-node 厚度 / 變厚度 / 複合層板**(單層等厚不變)。

---

## ② 公開 API

```cpp
// Public/FrameCore/SolveOptions.h — 加 2 個 opt-in bool(POD,維持 engine-agnostic)。
struct SolveOptions {
    real pivotTol       = 1e-12;
    bool enableReleases = false;
    bool useTimoshenko  = false;

    // S8-8a: opt-in QM6 incompatible-mode membrane (Wilson Q6 + Taylor 1976 correction).
    // false keeps the original bilinear Q4 membrane BIT-FOR-BIT, so the OpenSees ShellMITC4
    // plate gate (~1e-10) is preserved. true adds 2 bubble modes (1-xi^2, 1-eta^2) on the
    // in-plane u,v and condenses them out per element — defeats in-plane membrane locking.
    // The drilling (Rz / Hughes-Brezzi) block is untouched either way.
    bool useIncompatibleMembrane = false;

    // S8-8b: opt-in DKQ discrete-Kirchhoff THIN-plate bending (Batoz & Tahar 1982), replacing
    // the MITC4 assumed-shear bending block. false keeps the MITC4 Reissner-Mindlin plate
    // BIT-FOR-BIT (OpenSees ShellMITC4 gate). true = no transverse shear DOF -> thin-plate
    // only (t/L < ~1/20); mid/thick plates must keep this false (MITC4). Membrane + drilling
    // are shared (this flag ONLY swaps the bending block).
    bool useDKQPlate = false;
};
```

**無新 model 欄位、無新 result 欄位**:`ShellQuad`(matIdx/t/active)、`ShellElementForces`
(Mxx/Myy/Mxy + Qx/Qy + Nxx/Nyy/Nxy + per-corner)結構不變。DKQ 開時 `Qx=Qy=0`(Kirchhoff 無橫向剪,
誠實標於註解);膜內力 N 在 QM6 開時為**含 bubble 貢獻的一致值**(見 ④ recover)。

---

## ③ 資料流 / ④ 演算法

### 8a — QM6 不協調膜

**位移增廣**(元素層,僅膜 u,v):
```
u(ξ,η) = Σ Nᵢ uᵢ + P₁ a₁ + P₂ a₂          P₁ = 1 − ξ²,  P₂ = 1 − η²
v(ξ,η) = Σ Nᵢ vᵢ + P₁ a₃ + P₂ a₄          (a₁..a₄ = 4 個元素內部廣義 DOF,無全域對應)
```

**不協調應變** `B_inc`(3×4,對 [a₁,a₂,a₃,a₄]):
```
∂Pₖ/∂x = J₀⁻¹₀₀ ∂Pₖ/∂ξ + J₀⁻¹₀₁ ∂Pₖ/∂η       ⚠️ Jacobian 用元素中心常數 J₀ = J(0,0)
∂Pₖ/∂y = J₀⁻¹₁₀ ∂Pₖ/∂ξ + J₀⁻¹₁₁ ∂Pₖ/∂η       (∂P₁/∂ξ=−2ξ, ∂P₂/∂η=−2η, 其餘 0)
B_inc = [ ∂P₁/∂x   ∂P₂/∂x    0       0     ]     列 = [eps_x; eps_y; gamma_xy]
        [   0        0      ∂P₁/∂y  ∂P₂/∂y  ]
        [ ∂P₁/∂y   ∂P₂/∂y  ∂P₁/∂x  ∂P₂/∂x  ]
```
**中心 Jacobian = QM6 修正關鍵**(Taylor 1976):∫_Ω B_inc dΩ = J₀⁻¹·∫∫[−2ξ 或 −2η]dξdη = **0**
(∫₋₁¹ −2ξ dξ = 0)→ 常應力狀態下 bubble 不產生內力 → **畸變網格 patch test 通過**(P5 驗 QM6=0)。

**增廣勁度 + static condensation**(膜部分,Dm = plane-stress 3×3):
```
K_c  = ∫ Bm(ξ,η)ᵀ Dm Bm(ξ,η) · t dΩ  + ∫ γ Bd Bdᵀ dΩ    (12×12 協調膜+drilling = 現有 membraneK)
K_ca = ∫ Bm(ξ,η)ᵀ Dm B_inc      · t dΩ                  (12×4;僅 u,v 列非零,thz 列 0)
K_aa = ∫ B_incᵀ Dm B_inc        · t dΩ                  (4×4,正定)
K_membrane_QM6 = K_c − K_ca K_aa⁻¹ K_caᵀ                (12×12,凝聚回節點 DOF)
```
積分:2×2 Gauss(同現有);`Bm` 隨積分點(協調項不變),`B_inc` 用 J₀(常數)。drilling 項 γ=G·t **原樣保留**
(與 bubble 無耦合,K_ca 的 thz 列為 0)。prepare 末把 `K_membrane_QM6` 經 `membraneToShellMap` 注入 kl_。

**recover 一致性**(避免內力↔剛度不一致;P5 明列):
> **⚠️ 實作簡化(審核確認)**:`B_inc(0,0)=0`(∂P₁/∂ξ=−2ξ、∂P₂/∂η=−2η 在中心消失),而
> `ShellElementForces.N` 是**中心值**,故中心膜 recover 無需 bubble 貢獻——**實作不存 Hcond/Binc0
> 快取,直接用 `Bm(0,0)·dm`**(自動受益於 QM6 改善的位移解)。以下 Hcond/Binc0 路徑僅為未來
> **per-corner 膜輸出**才需(corner B_inc≠0);本階段不實作。
原始通用路徑(未來 per-corner 用):prepare 存 `Hcond = K_aa⁻¹ K_caᵀ`(4×12)與
中心 `B_inc0 = B_inc(0,0)`(3×4)。recover:
```
dm = membraneToShellMap · ul                 (12×1 局部膜位移 [u,v,thz]×4)
a  = −Hcond · dm                             (4×1 內部 DOF;僅依賴 u,v)
eps_centre = Bm(0,0)·dm + B_inc0·a           (含 bubble 貢獻)
N = t · Dm · eps_centre                      (Nxx,Nyy,Nxy 一致膜內力)
```
旗標 false 時 `a` 路徑完全略過,膜 recover 與現有逐位元相同。

### 8b — DKQ 離散 Kirchhoff 薄板(Batoz & Tahar 1982)

**旋轉場**:βx, βy 以 8 節點 serendipity 二次形函數插值(4 角 Nᵢ + 4 邊中 Pₖ);撓度 w 沿邊三次。
**Kirchhoff 約束**:沿 4 條邊在邊中點強制橫向剪應變 γ = 0 → 邊中旋轉以角點 (w,θx,θy) 唯一決定
(βₙ 沿邊 = 角點撓度與切向旋轉的 Hermite 組合的法向導數)→ 消去 4 邊中 DOF,降到 **12 DOF [w,θx,θy]×4**。
**曲率-位移** `B_DKQ`(3×12):κ = [βx,x ; βy,y ; βx,y+βy,x],由 Batoz-Tahar 標準形函數對 x,y 微分組出
(展開式見 Batoz & Tahar 1982 式(20)-(27);沿用其 `Hx,ξ / Hy,ξ / Hx,η / Hy,η` 向量結構,經 Jacobian 轉 x,y)。
**勁度**:
```
K_DKQ = ∫ B_DKQ(ξ,η)ᵀ Db B_DKQ(ξ,η) dΩ        Db = E t³/[12(1−ν²)] 板彎曲本構(同現有 plateK 的 Db)
```
積分:Batoz-Tahar 用 3 點(三角形側)或 2×2 Gauss(四邊形)——**初版 2×2 Gauss**(與現有一致,薄板足夠;
patch test 驗證);**無 Bshear 項、無 shear correction**。prepare 經 `plateToShellMap` 注入 kl_。

**recover**(DKQ 開):M = Db·B_DKQ(0,0)·dp(中心板矩),`Qx = Qy = 0`(Kirchhoff 無橫向剪);per-corner
`MxxC[k] = (Db·B_DKQ(cxi[k],ceta[k])·dp)` 同現有 MITC4 corner 機制。膜 N 共用(不受 DKQ 影響)。

### prepare 分支(`MITC4ShellElement::prepare`,讀 `opts`)
```
Mat12 Km = opts.useIncompatibleMembrane ? membraneK_QM6(...) : membraneK(...)   // 8a
Mat12 Kp = opts.useDKQPlate            ? plateK_DKQ(...)    : plateK(...)        // 8b
kl_ = Pbᵀ Kp Pb + Pmᵀ Km Pm                                                     // 組裝不變
// QM6 開時另存 qm6Hcond_/qm6Binc0_;旗標存成員 useQM6_/useDKQ_ 供 recover 分支
```

---

## ⑤ 檔案

- `Public/FrameCore/SolveOptions.h`:加 `useIncompatibleMembrane` / `useDKQPlate`(2 bool)。
- `Private/MITC4ShellElement.h`:加成員 `bool useQM6_ = false; bool useDKQ_ = false;`、QM6 recover 快取
  `Eigen::Matrix<real,4,12> qm6Hcond_;`、`Eigen::Matrix<real,3,4> qm6Binc0_;`(Private,不碰 POD 公開 API);
  標頭註解補 8a/8b opt-in 說明 + DKQ Qx=Qy=0 / 薄板限界。
- `Private/MITC4ShellElement.cpp`:匿名 namespace 加 `Binc(...)`(中心 J₀)、`membraneK_QM6(...)`、
  `plateK_DKQ(...)`(+ DKQ 形函數導數 helper);`prepare` 讀 opts 選分支 + 存快取;`recover` 依 useQM6_/useDKQ_
  分支。**Eigen 零洩漏**(本檔本就 include `FrameEigen.h`,新型別走它;**絕不**在 .cpp 直接 `#include <Eigen/...>`)。
- **不動**:`FrameSolver.cpp`(dispatch)、`Reanalysis.cpp`、`Shell.h`、`SolveResult.h`、`FrameModel.cpp`(validate /
  fingerprint)、四個 build .bat(無新 .cpp)。
- **gate 同步**:standalone `main.cpp` + `FrameTestFixtures.h` 加 F48(QM6)/F49(DKQ);UE 鏡像測試 +2;
  `linear_deep_audit.cpp` +checks;`run_gate.ps1 $ExpectedUeTests` 44→46;`Tools/opensees_compare.py` 加
  ShellDKGQ 對照(8b,若本機可用)。

---

## ⑥ Oracle(誠實分級)

### 8a QM6 — `[VERIFIED]`(QM6 是文獻方法之 FrameCore 實作)
- **F48a 膜常應變 patch test(機器精度)**:regular + parallelogram(affine,detJ 常數)網格,施加
  常應變位移場 → QM6 膜回復**精確常應力**、∫B_inc=0(中心 Jacobian)。**目標 rel ≤ 1e-12**(實測 2.2e-16)。
  **⚠️ 審核確認:∫B_inc=0 僅在 affine(平行四邊形)嚴格成立;一般非 affine 四邊形 ∫B_inc=O(h),QM6 過
  弱(Irons-Razzaque)patch、隨網格收斂(見 §⑨),非逐位元。** 仍是 QM6 中心 Jacobian 的硬獨立 oracle
  (Q6 不過此 affine patch,QM6 過)。
- **F48b Cook's membrane 收斂對照**:標準 Cook's skew membrane(梯形懸臂、tip 剪力),tip 撓度對參考
  解 ~23.96。**驗 QM6 比 Q4 顯著更接近**且符合 P5 數據量級(@N16:Q4 ≈ −3.2% / **QM6 ≈ −0.9%**;
  斷言 `|err_QM6| < |err_Q4|` 且 QM6 @N16 在 [−2%, 0] 內)。在 facet 平面內建模(膜主導)。
- **F48c 細長 in-plane 懸臂(膜鎖定釘樁)**:2×1 粗網格純 in-plane 彎曲懸臂,Q4 **鎖死 ~−90%** vs
  Euler-Bernoulli δ=PL³/3EI,QM6 大幅解鎖(P5:Q4 −92.6% → QM6 −5.8%)。斷言 `Q4 嚴重低估(<−50%)且
  QM6 接近(>−15%)`——**量化釘死膜鎖定的解除**。
- **預設 Q4 不破閘門**:`useIncompatibleMembrane=false` 時膜 K 與現有逐位元相同(F13-F16 既有 shell oracle +
  OpenSees ShellMITC4 平板 ~1e-10 不動)。

### 8b DKQ — `[VERIFIED]`(DKQ 是文獻方法之 FrameCore 實作)
- **F49a 常曲率 patch test(機器精度)**:畸變網格施加常曲率位移場(w=½(κx·x²+κy·y²+2κxy·xy),
  θ 相容)→ DKQ 回復**精確常彎矩**。**目標 rel ≤ 1e-10**(離散 Kirchhoff patch 是 DKQ 正確性必要條件)。
- **F49b 方板均布載重 vs 解析 Kirchhoff**:四邊簡支 / 固支方板,均布 q,DKQ 中心撓度對 Timoshenko-
  Woinowsky 級數解(簡支 w=0.00406 qa⁴/D;固支 0.00126 qa⁴/D),薄板 t/a=1/100。**斷言 DKQ 收斂到
  Kirchhoff 解**(N16~N24 rel < 1%);**與 MITC4 對照**:MITC4 在同網格收斂到 Mindlin(略軟,含橫剪),
  DKQ 直接 Kirchhoff → 薄板極限兩者趨近(誠實標 DKQ 無 overshoot)。
- **F49c OpenSees ShellDKGQ 對照**:`Tools/opensees_compare.py` 加 DKQ leg。**前置**:先確認本機
  openseespy 支援 `ops.element('ShellDKGQ',...)`(WS_E 標可用性兩源矛盾)。**可用** → 平板方板 DKQ vs
  ShellDKGQ 目標 rel ~1e-8~1e-10(同公式族,膜為 GQ12 非全等 → 寬容差 + 純彎曲主導案例);**不可用** →
  fallback 解析 Kirchhoff(F49b)+ vs MITC4 收斂,並在 spec/PROGRESS 誠實標「DKQ 無第三方殼 oracle,
  僅解析 + MITC4 互驗」。
- **預設 MITC4 不破閘門**:`useDKQPlate=false` 時板 K 與現有逐位元相同。

### audit(linear_deep_audit)
- **+QM6**:畸變膜 patch(a≈0、常應力)獨立重算;QM6 vs Q4 細長懸臂膜鎖定比值。
- **+DKQ**:常曲率 patch;DKQ vs MITC4 薄板方板中心撓度一致段;DKQ Qx=Qy=0 斷言(Kirchhoff)。
- **+ 預設不變**:同一模型 `useIncompatibleMembrane=false`/`useDKQPlate=false` 解 == 改旗標前 baseline
  逐位元(opt-in 安全性的執行期守門)。

---

## ⑦ Gate

F48(a/b/c)、F49(a/b/c);UE `FrameCore.Shell.IncompatibleMembrane` + `FrameCore.Shell.DKQPlate`
(**勿用 `IN`/`OUT` 常數名 = Windows SAL 巨集,用 `kIn` 等**);audit +N checks(82→~88);
`opensees_compare.py` 加 ShellDKGQ leg(若可用);`run_gate.ps1 $ExpectedUeTests` 44→**46**。
**四 build 腳本免動**(無新 .cpp)。**commit 前五腿全綠**(`run_gate.ps1 -RequireOpenSees`)。

---

## ⑧ 效能驗收

- QM6:每膜元素多一個 4×4 對稱 inverse(`K_aa`)+ 12×4 凝聚 — O(1) 小常數,組裝期一次性;recover 多一個
  4×12 mat-vec。對整體 factor/solve 可忽略(殼數量級不變)。
- DKQ:取代 assemble-shear 的等量積分;無額外全域 DOF(仍 24/facet)。
- **無新全域 DOF、無 fingerprint 變動 → ReSolve / factorize-once / 所有既有分析(modal/buckling/collapse/
  P-Delta/tension-only/BESO)沿用殼勁度,opt-in 透傳 `opts.solve`**(ReSolve `prepare(work, P.opts.solve)`
  已傳)。誠實標:QM6/DKQ 與既有殼分析的交互未逐一加 oracle(本階段 oracle 聚焦元素級正確性 + 預設不變)。

---

## ⑨ 誠實邊界 / novelty 定位

- **QM6 = Wilson-Taylor 1973 + Taylor 1976,DKQ = Batoz-Tahar 1982**:皆**成熟文獻元素**,FrameCore 為其
  **實作**,**非新方法、非新穎**。`[VERIFIED]` = 文獻方法的 oracle 化實作;`[NEW CODE]` 僅限「opt-in 旗標
  整合 + 與既有 24-DOF facet / drilling / 組裝 / ReSolve 的接線」(工程整合,非算法原創)。
- **QM6 邊界**:只幫**膜主導 / in-plane bending**問題;純板彎曲(MITC4 assumed-shear 已良好)與 flat-facet
  幾何 O(1/N²) floor **不受益**(文檔續標「flat facet, O(1/N²) geometric floor remains」);QM6 過**弱**
  patch test(modified integration,Taylor 1976),高度扭曲網格仍有殘差(WS_E F-E2);bubble 對矩形精確、
  扭曲漸降。
- **DKQ 邊界**:**完全無厚板能力**(無橫向剪;t/L > ~1/20 Kirchhoff 假設失效 → 必回 MITC4);`Qx=Qy=0` 是
  **刻意**(非遺漏);初版 2×2 Gauss(Batoz 原文部分用 3 點,薄板 patch 通過即採);DKQ 是**薄板快路 /
  Kirchhoff oracle**,**非**主殼元素(WS_E「避什麼」第一條)。**⚠️ 審核確認:DKQ per-corner 彎矩
  `MxxC[k]` 因角點 serendipity 導數大,**非該角點逐點曲率的好估計**(常曲率場角點≠中心);per-corner 仍
  線性於 dp(combine/envelope 有效),但 DKQ 設計峰值用中心值(Gauss 外插 [NOT IMPLEMENTED])。**
- **預設安全**:兩旗標 false 保 OpenSees `ShellMITC4` ~1e-10 平板閘門(此 opt-in 策略與 E2 sparse modal /
  D3 角點彎矩同模式)。**改進膜/板不得對 `ShellMITC4` 要求 ~1e-10**(會誤判退步)→ QM6/DKQ 的第三方對標
  各走 `ShellDKGQ`(DKQ)或文獻收斂表(QM6 無直接 OpenSees 對應膜)。
- **對標措辭**:FrameCore MITC4+QM6 在四邊形膜精度上對標/超越 Karamba3D TRIC(CST 膜 = 最低階),但**須附
  先行技術定位**(WS_E F-E4),禁裸宣稱「優於 Karamba」。

---

## ⑩ 風險 / fallback

1. **破 OpenSees ShellMITC4 閘門(最大風險)** → 兩旗標**預設 false 逐位元不變** + audit 加「預設 == baseline」
   守門 + F13-F16 既有 shell oracle 不動。**任一旗標 false 路徑碰到既有 K 即視為退步**,先跑 `build.bat`
   確認 F13-F16 + OpenSees 平板未動才繼續。
2. **QM6 畸變 patch 不過**(中心 Jacobian / 減均值沒做對)→ 對齊 P5 numpy 參考(已驗 QM6=0);F48a 畸變網格
   必過 1e-12,否則回查 `B_inc` 是否用 J₀ 常數 + ∫B_inc=0。
3. **K_aa 奇異 / 病態**(凝聚 inverse)→ K_aa 對正定 Dm 理論正定;加 LDLT / 條件檢查,病態則該元素退回 Q4 +
   diagnostic(誠實標,不靜默)。
4. **DKQ B 矩陣推導錯**(Batoz-Tahar 形函數冗長易錯)→ F49a 常曲率 patch(機器精度)當第一道硬閘門;旋轉等變
   性(整體旋轉 facet 解不變)當第二道;方板 vs 解析 Kirchhoff 當第三道。任一不過即 B 推導有誤。
5. **ShellDKGQ 本機不可用** → fallback 解析 Kirchhoff(F49b)+ vs MITC4 薄板收斂,PROGRESS 誠實標「無第三方
   殼 oracle」(不阻塞 8b 完成)。
6. **DKQ 用於中厚板誤用** → 文檔 + 標頭明標 t/L < ~1/20;不做自動薄/厚切換(caller 選旗標,誤用是 caller 責任,
   同 useTimoshenko)。
7. **UE SAL 巨集** → 測試常數勿命名 `IN`/`OUT`(經 CoreMinimal.h 拉入);UE build 用 PowerShell 前景 Build.bat。
8. **Bash→powershell 路徑** → 反斜線被 bash 跳脫;跑 gate 用 PowerShell 工具或正斜線路徑(本輪起點踩雷已記)。
