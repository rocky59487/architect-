> 出處:研究輪 fan-out agent 查證報告(2026-06-10,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

根據所有查證結果,整理如下:

# WS-N 三個候選差異化算法的先行技術調查(誠實定位)

## 摘要

本調查針對 FrameCore 三個候選差異化算法逐一查閱先行文獻,確認其 novelty 定位。N1(SMW 精確重分析)的核心公式由 Akgün et al. 2001 系統化建立,適合低秩修改;rank-12 梁柱元素移除雖在技術可行範圍但組合進漸進倒塌驅動器的工程整合尚無直接先例。N2(collapse-aware BESO)與 fail-safe TO(Jansen 2014)概念重疊但有具體差異:fail-safe TO 針對連續介質任意損傷域,而使用 LSP 倒塌驅動器作為 BESO 迭代評估器是尚未在文獻中明確出現的整合路線。N3(凍結分解 P-Delta)是 Wilson-Habibullah 1987 的已知算法框架,在 FrameCore 中是架構整合,非新算法。

---

## 發現

### N1 — SMW 精確重分析 / sparse Cholesky updown

**N1-a: 最接近的先行文獻**

1. **Akgün, Garcelon & Haftka (2001)** — "Fast exact linear and non-linear structural reanalysis and the Sherman–Morrison–Woodbury formulas," IJNME 50(7):1587–1606.
   - 這篇確實是做了 SMW 精確重分析的開創性論文。方法:把結構剛度改動寫為 ΔK = U C Vᵀ(低秩分解),利用 Woodbury 公式在已有分解上更新逆矩陣,不重做全分解。
   - 適用範圍:「低秩」修改(論文重點在少數幾個元素的參數變化、截面改變等),修改 DOF 數越小效益越大。論文中列舉的例子主要是截面面積/E 模量的少數元素更新,非 rank-12 梁柱元素完整加/移除的情景。
   - 限制:rank 越高(即 r 很大),需解一個 r×r 的密矩陣問題,成本隨 r 提升;對真正的「元素完整移除」(rank=12 for 3D beam),每個被移除元素帶來 12×12 塊,若一次移除多個元素代價仍可觀。[LIT:https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.87]

2. **Kirsch (2003)** — "Reanalysis and sensitivity reanalysis by combined approximations," Struct. Multidiscip. Optim.
   - Combined Approximations(CA)是另一類:以二項展開基底向量做近似重分析,low-rank 修改時可給出精確解;高 rank 修改時退化為近似。論文本身標明對「拓撲修改」(加/刪元素)的應用有精確與近似分支。[LIT:https://link.springer.com/article/10.1007/s00158-009-0369-1]

3. **Amir, Bendsoe & Sigmund (2009)** — "Approximate reanalysis in topology optimization," IJNME 78(12):1474–1491.
   - 在拓撲優化框架中引入 Kirsch CA 近似重分析,顯著降低每個優化迭代的求解成本。說明「近似誤差被靈敏度分析的不確定性掩蓋」,即拓撲優化不必精確求解。[LIT:https://onlinelibrary.wiley.com/doi/10.1002/nme.2536]

4. **Davis & Hager (2001, 2008, 2009)** — CHOLMOD sparse Cholesky updown 系列。
   - 數學後端:CHOLMOD 的 `cholmod_updown` 函式支援 A ± WW^T(rank-k update/downdate),時間正比於 L 中受影響的條目數量,多秩版本比逐步 rank-1 快(一次通過 L)。這是 N1 中 sparse LDLT downdate 方案的數值後端。
   - 文獻:Davis & Hager, "Multiple-Rank Modifications of a Sparse Cholesky Factorization," SIAM JMAA; "Algorithm 887: CHOLMOD, Supernodal Sparse Cholesky Factorization and Update/Downdate," TOMS 35(3), 2008. [LIT:https://dl.acm.org/doi/abs/10.1145/1391989.1391995][LIT:https://people.engr.tamu.edu/davis/publications_files/MultipleRank_Modifications_of_a_Sparse_Cholesky_Factorization.pdf]
   - Eigen 的 `CholmodSimplicialLDLT` 包裝了 CHOLMOD 後端,但 Eigen 公開 API **不暴露** updown 介面(Eigen 文件只有 `compute`/`factorize`/`solve`),需直接呼叫底層 `cholmod_factor` 指標或用 SuiteSparse C API。[LIT:https://eigen.tuxfamily.org/dox/classEigen_1_1CholmodSimplicialLDLT.html]

5. **Huang & Verchery(精確 SMW 元素移除)** — 基於 SMW 公式計算元素移除後的修改柔度矩陣;Ren et al. 2020 "Structural Reanalysis Based on FRFs Using Sherman-Morrison-Woodbury Formula," Shock and Vibration, 也在頻響函式語境下使用 SMW 做元素修改重分析。[LIT:https://onlinelibrary.wiley.com/doi/10.1155/2020/8212730]

**N1-b: SMW 用進漸進倒塌模擬或 BESO**

搜尋「reanalysis progressive collapse SMW」未找到直接把 SMW/sparse Cholesky downdate 嵌入 **GSA-LSP 漸進倒塌驅動器**的先行論文。Kirsch 的重分析用進拓撲優化(Amir 2009)有明確文獻,但使用 LDLT downdate 驅動 **逐步元素移除序列**以模擬倒塌動態的架構[UNKNOWN]。

**N1 結論(gap 與 novelty)**

- SMW 精確重分析本身:**已被 Akgün et al. 2001 做過**,不可宣稱算法創新。
- sparse Cholesky rank-k downdate(CHOLMOD updown):**Davis & Hager 2001 已建立**。
- 工程整合 gap:把 rank-12 梁柱元素 downdate + 倒塌驅動器序列控制(選最危元素 → downdate 舊分解 → 重解 → 再選)整合成一個架構,**尚無直接先行論文可查**。主要是「已知數學工具的新應用整合」而非數學創新。

---

### N2 — collapse-aware BESO

**N2-a: 最接近的先行文獻**

1. **Jansen, Lombaert, Schevenels et al. (2014)** — "Topology optimization of fail-safe structures using a simplified local damage model," Struct. Multidiscip. Optim. 49:657–666.
   - 這是把「構件失效」納入拓撲優化的第一篇里程碑論文。模型:對連續介質結構,在有限元網格上枚舉一組 patch(面積固定形狀)作為損傷情境,對每個損傷場景分析剛度損失,以最小最壞情況柔度(minimax formulation)為目標函式,用 SIMP 方法優化。本質是「任意一個 patch 被移除後,結構仍可承載」。[LIT:https://link.springer.com/article/10.1007/s00158-013-1001-y]

2. **Zhou & Fleury (2016)** — "Fail-safe topology optimization," Struct. Multidiscip. Optim. 54(5):1225–1243.
   - 進一步發展 fail-safe TO;延伸至大規模結構。[LIT:https://link.springer.com/article/10.1007/s00158-016-1507-1]

3. **Stolpe(2017) / Nguyen et al. (2023) — adaptive fail-safe truss TO** — 把 fail-safe 概念延伸到桁架;每個桿件可被移除,需保持承載力。[LIT:https://link.springer.com/article/10.1007/s00158-023-03585-x]

4. **Robust/robustness-based TO for progressive collapse(近年)** — "An integrated framework for optimization-based robust design to progressive collapse of RC skeleton buildings," Innov. Infrastr. Solut. 2025;及 Progressive collapse design of seismic steel frames using structural optimization, ResearchGate — 搜索結果顯示部分文獻嘗試將 GSA alternate path 作為約束放入優化框架,但論文具體使用 BESO 且搭配 LSP 倒塌驅動器的架構尚未查到直接先例。[LIT:https://link.springer.com/article/10.1007/s41062-025-02243-z]

**N2-b: fail-safe TO 與「LSP 倒塌驅動器當 BESO 評估器」的差別**

| 面向 | Jansen 2014 fail-safe TO | FrameCore 擬做的 collapse-aware BESO |
|---|---|---|
| 結構類型 | 連續介質(SIMP density) | 梁柱框架(離散桿件) |
| 損傷模型 | 有限元 patch 移除(任意形狀/位置) | 離散構件移除(GSA LSP 逐步移除序列) |
| 評估器 | 最壞情況靜態柔度 | LSP 倒塌驅動器(含 D/C 篩/碎塊連通/雙終點 Stable/Collapsed) |
| 拓撲優化方法 | SIMP | BESO(二元元素活躍/停用) |
| 目標 | 最小化最壞情況下的結構柔度 | 最大化結構在漸進倒塌下的冗餘度/延性 |

fail-safe TO = 「任一損傷後仍可承載」,與 progressive collapse resistance 概念相近但不同:fail-safe TO 是靜態單步最壞情況,而 LSP 驅動器是動態序列移除直到倒塌。把 LSP 驅動器輸出(Stable/Collapsed)當 BESO 每步評估函式,目前文獻中[UNKNOWN]。

**N2 結論(gap 與 novelty)**

- fail-safe TO 本身:Jansen 2014 **已明確做過連續介質版**。
- 把 GSA-LSP 倒塌模擬(序列線性分析/雙終點/動態載重因子)整合為 BESO 迭代評估器的框架:未查到直接先行論文,屬**方法整合 novelty**。

---

### N3 — 凍結分解 P-Delta

**N3-a: 最接近的先行文獻**

1. **Wilson & Habibullah (1987)** — "Static and Dynamic Analysis of Multi-Story Buildings, Including P-Delta Effects," Earthquake Spectra 3(2):289–298.
   - 核心貢獻:對多層建築,重力質量引起的軸力在施加側向荷載過程中**幾乎不變**,因此幾何勁度矩陣 Kg 可由重力荷載態一次計算並**線性化固定**,無需每步迭代更新。這等效於:先解 K·u₀ = f(重力),從 u₀ 算 Kg(u₀),然後對側向荷載解 (K + Kg)·u = f_lateral,一次分解搞定。這正是「凍結分解 P-Delta」的原型:Kg 凍結在重力態,K+Kg 只做一次分解。[LIT:https://pubs.geoscienceworld.org/eeri/earthquake-spectra/article-abstract/3/2/289/584240/Static-and-Dynamic-Analysis-of-Multi-Story][LIT:https://www.semanticscholar.org/paper/Static-and-Dynamic-Analysis-of-Multi-Story-P-Delta-Wilson-Habibullah/0ae7142d90889c03e1a025dd352c053ad673141c]

2. **McGuire, Gallagher & Ziemian — "Matrix Structural Analysis" (2nd ed., 2000)** — 教科書中描述 geometric stiffness 疊加 elastic stiffness 再求解的標準做法;P-Delta 的 pseudo-load 迭代法亦為教科書內容。章節:Wilson 的個人網站收藏了 Chapter 11 "Geometric Stiffness and P-Delta Effects" (Wilson 2005)。[THEORY][LIT:https://edwilson.org/bookshelf/2005/2005%20Chapter%2011.pdf]

3. **RISA-3D / SAP2000 / STAAD.Pro 技術手冊(工業標準實作)**
   - RISA「Geometric Nonlinear Stiffness Method」:assembles Kg from initial axial forces,folds into K,refactorizes once(一次)。
   - 「Standard Nodal Shear Method」(另一流派):pseudo-shear load 加回右端,不改左端矩陣,但需迭代。
   - 兩種方法均是業界標準。[LIT:https://risa.com/risahelp/adaptbuilder/Content/Analysis/2nd%20Order-PDelta-Analysis/Pdelta-Theoretical-Background.htm][LIT:https://docs.bentley.com/LiveContent/web/RAM%20Structural%20System%20Help-v4/en/GUID-98FA40E1-CC05-4897-885C-648934109607.html]

4. **Pseudo-load 迭代(另一分支)**:「Nodal shear」或「fictitious lateral load」法把 P·Δ/L 作為等效側力加回右端,用固定 K(不加 Kg)迭代直到收斂。這是 Wilson-Habibullah 方法的**對立面**,需要迭代但不修改左端矩陣。[LIT:https://risa.com/risahelp/risa3d/Content/3D_2D_Only_Topics/P-Delta%20-%20Analysis.htm]

**N3 結論**

Wilson & Habibullah 1987 的方法就是「(K + Kg)·u = f,Kg 由重力態線性化,一次分解」,這與「凍結分解 P-Delta」描述完全一致。FrameCore 若實作此方案,是在自己的 factorize-once `PreparedSystem` 架構上整合一個 **1987 年已知的算法**;差異在「與既有互動式重解架構的無縫結合」(同一 `PreparedSystem` 在不同工況重複求解,P-Delta 只需多一次 Kg assemble),屬**已知方法的架構整合**。

---

## 對 FrameCore 的含義

**N1(SMW / sparse Cholesky downdate)**
- 應採:進行「元素移除重分析」時,可用 CHOLMOD `cholmod_updown` 做 rank-12 downdate,成本正比於受影響的 L 條目數,而非全量重分解。
- Eigen 的 `SimplicialLDLT` **沒有** updown API;若要用須直接呼叫 SuiteSparse C API 或換用 CHOLMOD 後端包裝。
- 誠實邊界:對大型網格,若一步倒塌序列要移除多個元素(rank = 12×n_removed),每步仍需做 rank-12n 的 downdate,而非一次全分解;相比全分解的加速比視稀疏模式而定。對結構優化每步移除單元較多(BESO 每步可移除 1-2% 元素)時,多次 rank-12 downdate 累積成本不一定優於整批重分解。
- 避免:宣稱「SMW 重分析用於結構是新算法」——算法本身是既有 1987/2001 工作。
- 建議誠實宣稱:「利用 CHOLMOD sparse Cholesky rank-k downdate 加速 FrameCore 漸進倒塌驅動器中的元素移除重解,是既有數值工具在互動式倒塌評估場景的工程整合。」

**N2(collapse-aware BESO)**
- 應採:若要做「把漸進倒塌驅動器當 BESO 評估器」,可宣稱與 fail-safe TO(連續介質/SIMP)的概念親緣性,但強調三個具體差異:①框架結構(離散桿件);②評估器是 LSP 序列線性分析(有動態放大係數、碎塊清除等);③判定指標是 Stable/Collapsed 終點而非靜態柔度。
- 避免:把自己的方案稱為「fail-safe TO」——Jansen 2014 已確立該術語對應連續介質 SIMP 最壞情況靜態柔度,使用相同術語會造成混淆或優先權爭議。
- 建議誠實宣稱:「以 GSA 替換荷載路徑(LSP 漸進倒塌模擬)為評估器的框架結構 BESO,在概念上類似 fail-safe 拓撲優化(Jansen 2014),但針對離散框架結構並使用序列線性倒塌分析而非靜態最壞情況柔度作為評估函式。」

**N3(凍結分解 P-Delta)**
- 應採:直接在 `assembleGeometric(model, result)` 之後組成 K_eff = K_e + Kg,對 K_eff 做一次 `PreparedSystem::factorize`,後續所有工況呼叫 `solveLoad` 重用同一分解——這是 FrameCore 既有架構的自然延伸。
- 誠實邊界:Kg 依賴軸力(故依賴位移),嚴格非線性;「凍結在重力態」是線性化近似,適用於側向位移小、二階效應為修正量的情況,不適用於大側移或大位移問題。
- 避免:宣稱「凍結分解 P-Delta」是新算法。
- 建議誠實宣稱:「實作 Wilson-Habibullah (1987) 線性化 P-Delta 方法,在 FrameCore 的 factorize-once PreparedSystem 架構下一次組裝幾何勁度並整批求解多個工況。」

---

## 來源清單

- [Akgün, Garcelon & Haftka 2001 — Fast exact linear and nonlinear structural reanalysis and the Sherman–Morrison–Woodbury formulas, IJNME 50(7):1587–1606](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.87)
- [Akgün et al. 2001 — ResearchGate full-text request](https://www.researchgate.net/publication/229456253_Fast_exact_linear_and_non-linear_structural_reanalysis_and_the_Sherman-Morrison-Woodbury_formulas)
- [Kirsch — Combined Approximations, a general reanalysis approach for structural optimization, Struct. Multidiscip. Optim.](https://link.springer.com/article/10.1007/s001580050141)
- [Kirsch — Reanalysis and sensitivity reanalysis by combined approximations, Struct. Multidiscip. Optim. 2009](https://link.springer.com/article/10.1007/s00158-009-0369-1)
- [Amir, Bendsoe & Sigmund (2009) — Approximate reanalysis in topology optimization, IJNME 78(12):1474–1491](https://onlinelibrary.wiley.com/doi/10.1002/nme.2536)
- [Amir & Sigmund (2012) — Efficient reanalysis techniques for robust topology optimization, CMAME 245–246:217–231](https://www.sciencedirect.com/science/article/abs/pii/S0045782512002289)
- [Davis & Hager — Algorithm 887: CHOLMOD supernodal sparse Cholesky factorization and update/downdate, TOMS 35(3), 2008](https://dl.acm.org/doi/abs/10.1145/1391989.1391995)
- [Davis & Hager — Multiple-Rank Modifications of a Sparse Cholesky Factorization, SIAM JMAA](https://people.engr.tamu.edu/davis/publications_files/MultipleRank_Modifications_of_a_Sparse_Cholesky_Factorization.pdf)
- [Davis & Hager — Dynamic Supernodes in Sparse Cholesky Update/Downdate, TOMS 2009](https://dl.acm.org/doi/abs/10.1145/1462173.1462176)
- [Eigen::CholmodSimplicialLDLT documentation](https://eigen.tuxfamily.org/dox/classEigen_1_1CholmodSimplicialLDLT.html)
- [Ren et al. (2020) — Structural Reanalysis Based on FRFs Using Sherman–Morrison–Woodbury Formula, Shock and Vibration](https://onlinelibrary.wiley.com/doi/10.1155/2020/8212730)
- [Jansen, Lombaert, Schevenels et al. (2014) — Topology optimization of fail-safe structures using a simplified local damage model, Struct. Multidiscip. Optim. 49:657–666](https://link.springer.com/article/10.1007/s00158-013-1001-y)
- [Zhou & Fleury (2016) — Fail-safe topology optimization, Struct. Multidiscip. Optim. 54(5):1225–1243](https://link.springer.com/article/10.1007/s00158-016-1507-1)
- [Adaptive topology optimization of fail-safe truss structures, Struct. Multidiscip. Optim. 2023](https://link.springer.com/article/10.1007/s00158-023-03585-x)
- [Integrated framework for optimization-based robust design to progressive collapse, Innov. Infrastr. Solut. 2025](https://link.springer.com/article/10.1007/s41062-025-02243-z)
- [Wilson & Habibullah (1987) — Static and Dynamic Analysis of Multi-Story Buildings Including P-Delta Effects, Earthquake Spectra 3(2):289–298 (GeoScienceWorld)](https://pubs.geoscienceworld.org/eeri/earthquake-spectra/article-abstract/3/2/289/584240/Static-and-Dynamic-Analysis-of-Multi-Story)
- [Wilson & Habibullah (1987) — Semantic Scholar record](https://www.semanticscholar.org/paper/Static-and-Dynamic-Analysis-of-Multi-Story-P-Delta-Wilson-Habibullah/0ae7142d90889c03e1a025dd352c053ad673141c)
- [RISA P-Delta Theoretical Background](https://risa.com/risahelp/adaptbuilder/Content/Analysis/2nd%20Order-PDelta-Analysis/Pdelta-Theoretical-Background.htm)
- [Bentley RAM P-Delta by Geometric Stiffness Method](https://docs.bentley.com/LiveContent/web/RAM%20Structural%20System%20Help-v4/en/GUID-98FA40E1-CC05-4897-885C-648934109607.html)
- [GSA Alternate Path Analysis Guidelines 2016 (PDF)](https://www.gsa.gov/system/files/Progressive_Collapse_2016.pdf)
