> 出處:研究輪 fan-out agent 查證報告(2026-06-10 補跑,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

現在我已收集到足夠的文獻根據,可以撰寫完整報告。

# WS-C P-Delta 二階分析的方法選型與 oracle

## 摘要

P-Delta 二階分析主要有三條路徑:(1) **pseudo-load 凍結分解迭代法**(Wilson & Habibullah 1987):以初次線彈性解的軸力生成等效側向「假載荷」,迭代疊加,等價幾何級數,收斂充要條件 P < P_cr;(2) **幾何勁度直接修正法** K_t = K_e + K_g:一步組裝切線勁度、重分解,無需迭代,但 K_g 以 Hermitian 多項式線性化,P/P_cr 高時有多項式截斷誤差;(3) **精確穩定函數法**:以 tan(kL)/tanh(kL) 型三角函式取代多項式,對單構件零誤差。OpenSees `PDeltaCrdTransf` 屬路徑 (2) 的結構層 pseudo-load 簡化版,採用小位移線性相容性假設,K_g 修正為 P/L 的剛體側移項,**不捕捉 P-δ(構件內撓曲)效應**,在 Newton-Raphson 迴圈內每步更新;拿它當 oracle 時必須知道它算的是 P-Δ(層間側移)近似。AISC B1 = C_m/(1 − P/P_e) 是單構件 P-δ 閉式放大器,為穩定函數精確解的一階展開,在高軸壓比 P/P_e > 0.6~0.7 時誤差可達 5–15%。

## 發現

### 1. Wilson & Habibullah 1987 pseudo-load / 幾何勁度一次分解法

**基本格式:**以線彈性解取得各層(或各構件)軸力 P_i,計算等效側向假載荷 ΔF_i = P_i · Δ_i / h_i(storey shear 型),疊加回原載荷重解;反覆直到位移增量收斂。因幾何勁度 K_g 只在組裝到整體後修正,初次分解 K_e 可凍結重用(factorize-once),每步只需 back-substitute。[LIT:https://pubs.geoscienceworld.org/eeri/earthquake-spectra/article-abstract/3/2/289/584240]

**收斂條件:**每次迭代位移增量比前一步多乘一個係數 r ≈ P/P_cr(對單層為穩定係數 θ),因此級數收斂充要條件為 r < 1,即 P < P_cr。θ → 1 時迭代發散;Wilson & Habibullah 的貢獻是利用「樓層質量不隨側向位移變動」的建築結構特性,把幾何勁度項一次性組裝入整體 K,使**多自由度框架只需一次分解而非每步重組**。[LIT:https://www.semanticscholar.org/paper/Static-and-Dynamic-Analysis-of-Multi-Story-P-Delta-Wilson-Habibullah/0ae7142d90889c03e1a025dd352c053ad673141c]

**迭代格式與幾何級數:**設第一步線彈性位移 Δ_0,pseudo-load 在第二步產生增量 Δ_1 = r·Δ_0,第三步 Δ_2 = r·Δ_1 = r²·Δ_0 …,總位移 = Δ_0(1 + r + r² + …) = Δ_0/(1−r) = Δ_0/(1−P/P_cr)。[THEORY]

**與每步重組 K_e+K_g 重分解法的精度差異:**若凍結分解版本採用初次線彈性軸力計算 K_g(即 K_g 不隨迭代更新),則相當於只做一次幾何勁度修正;若每步更新軸力再更新 K_g 重分解,二者差異在 P/P_cr < 0.3 時極小(< 0.1%),在 P/P_cr ≈ 0.8 時數值實驗顯示差距約 1–3%。幾何勁度矩陣本身已是 Hermitian 多項式線性化近似,在 P/P_cr = 0.9 時與精確穩定函數解的相對誤差約 **10%**。[LIT:https://pmc.ncbi.nlm.nih.gov/articles/PMC10719514/]

**「two-cycle iterative method」:**Chen & Lui 1991 一節定義「two-cycle = 先算線彈性得軸力→組 K_g→以 K_t = K_e + K_g 解第二遍」,本質上是幾何勁度法的兩遍流程,不是無限迭代;對於設計級精度(P/P_cr < 0.5 以下)誤差可控。[LIT:https://pmc.ncbi.nlm.nih.gov/articles/PMC10719514/]

### 2. 教科書方法分類 (McGuire, Gallagher & Ziemian 2e; Chen & Lui)

**McGuire et al. 2000「Matrix Structural Analysis」2nd ed.(免費 e-book at mastan2.com):**第 5 章(幾何非線性)→第 8 章(非線性平衡方程求解)。分類:(a) **直接迭代法(modified Newton)**:每步以初始或上一步切線勁度求解,不更新;(b) **Newton-Raphson**:每步以當前切線勁度重組重分解;(c) **弧長法(arc-length)**:跟蹤後屈曲路徑。幾何勁度線性化假設 = Hermitian 多項式形函式,元素級 K_g 為:

$$[K_g]_{e} = \frac{P}{L}\begin{bmatrix}0&0&0&0\\0&1&0&-1\\0&0&0&0\\0&-1&0&1\end{bmatrix}(\text{側移 DOF,平面}2\text{D})$$

加上端部轉角耦合項,完整 12×12 梁柱 K_g 見 McGuire §5.3。[LIT:麥書 Matrix Structural Analysis 2nd ed., McGraw-Wiley 2000, Ch.5–8; 免費下載 https://digitalcommons.bucknell.edu/books/7/]

**Chen & Lui 1987「Structural Stability: Theory and Implementation」:**分類相似,額外強調 stability functions s(kL), c(kL)(兩參數描述梁柱切線勁度)可取代 K_e + K_g 的分裂表達,對單構件等效精確解。Goto & Chen 1987 給出三維 stability function 形式。[LIT:W.F. Chen & E.M. Lui, Structural Stability: Theory and Implementation, Elsevier 1987; Scribd PDF https://www.scribd.com/document/440700188/1987-Chen-Lui-Structural-Stability-Theory-and-Implementation-pdf]

**直接迭代 vs Newton-Raphson:**直接迭代每步成本低但收斂慢(線性收斂,比率 ~ ρ(K_t⁻¹ K_g));N-R 每步重組切線勁度,收斂域大、二次收斂。對純 P-Delta(結構本質幾何非線性不強)工程實踐取直接迭代 2–5 步即可;N-R 在後屈曲或材料非線性時必要。[THEORY]

### 3. 精確 oracle:懸臂柱 beam-column 閉式解

**問題設定:**懸臂柱(固定端在 x=0,自由端在 x=L),軸壓力 P,自由端水平集中力 H。

**控制微分方程:**EI v'' + P v = H(L − x),令 k = √(P/EI)。

**精確解(tip deflection,自由端水平位移):**

$$\delta_{tip} = \frac{H}{P}\left(\frac{\tan(kL)}{kL} \cdot L - L\right) = \frac{HL^3}{3EI}\cdot\frac{3(\tan(kL) - kL)}{(kL)^3}$$

亦可改寫為:

$$\delta_{tip} = \frac{H}{Pk}\left(\frac{1 - \cos(kL)}{\cos(kL) \cdot \sin(kL)/kL - \cos(kL)}\right)$$

最緊湊形式(固定端-自由端邊界條件):

$$\delta_{tip} = \frac{H}{P}\left(\frac{\tan(kL)}{k} - L\right), \quad k = \sqrt{P/EI}$$

P → 0 時 tan(kL)/(k) → L + (kL)³/(3k) = L + PL³/(3EI),還原 H L³/(3EI)。[THEORY + LIT:https://pmc.ncbi.nlm.nih.gov/articles/PMC10719514/ (提供 ODE 解結構;完整閉式見任何「Beam-Column」章節,如 Timoshenko & Gere, Theory of Elastic Stability, 2nd ed., §1.3–1.4)]

**自由端固定端彎矩:**M_base = H·L + P·δ_tip = H·L / cos(kL)(精確).

**Timoshenko & Gere 引用:**Timoshenko S.P. & Gere J.M., Theory of Elastic Stability, 2nd ed., McGraw-Hill 1961, Ch.1,方程 (1-12)/(1-13)。[LIT:書目,網路無合法全文]

**AISC B1 = C_m/(1 − P/P_e) 的適用範圍與誤差:**
- B1 是穩定函數精確解對 P/P_e 的一階 Padé/Taylor 近似,端部彎矩 = B1 × M_0(一階彎矩)。
- 等端彎矩(最不利 C_m = 1)時 B1 最保守,一般情況 C_m < 1。
- 誤差來源:C_m 表達式不含 P/P_e 影響(LRFD Cm 僅含端彎矩比值 M_a/M_b)。
- 在 P/P_e < 0.6 時 B1 通常在精確解 ±5% 內;P/P_e > 0.7 雙曲率情況下 B1 可低估 P-δ 達 10–15%(LRFD Cm 低估端彎矩比的影響)。[LIT:https://www.sciencedirect.com/article/abs/pii/S0141029615005295;https://ascelibrary.org/doi/10.1061/(ASCE)0733-9445(1999)125:2(219)]
- AISC 360-16 Commentary C-C2.1 要求 B1 ≥ 1,並說明 B1 僅適用於無側移構件的 P-δ;**B2 才處理 P-Δ(層間側移)**。[LIT:AISC 360-16 Spec, Appendix 8 / Chapter C Commentary; https://www.aisc.org/globalassets/aisc/publications/standards/a360-16-spec-and-commentary_june-2018.pdf]

### 4. OpenSees PDeltaCrdTransf 的數學實質

**官方文檔說明:**「performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global coordinate system, considering second-order P-Delta effects.」採用小位移假設作相容性(線性相容),與 P-Δ 平衡項配套。[LIT:https://opensees.github.io/OpenSeesDocumentation/user/manual/model/geomTransf/PDelta.html]

**數學實質(Denavit & Hajjar 2013 NEU-CEE-2013-02 報告):**
- PDelta 轉換在 basic system → global 的雅可比 T 中保留 P/L 的一階幾何項,相當於在整體組裝後每個梁柱元素貢獻一個如下結構:剛度矩陣里加入 axial force × 端部側移 / L 的側移耦合項(等效 P/L 側移勁度)。
- **不包含 P-δ(構件中間撓曲)效應**,僅捕捉節點側移 Δ 引起的幾何效應(chord rotation)。
- 垂直方向的端部位移變化被忽略(線性相容性)。
- 對比 Corotational:Corotational 採用大位移相容性,追蹤 chord 轉動並更新局部座標,因此能捕捉幾何大位移(lateral drift > 10% 時有差別)且能正確處理拉力加勁效應;PDelta 在小-中等側移下(drift < 5–10%)與 Corotational 結果幾乎一致。[LIT:https://openseesdigital.com/2022/11/15/geometric-transformation/; Denavit & Hajjar 2013, Report NEU-CEE-2013-02, Northeastern Univ. https://repository.library.northeastern.edu/files/neu:376268]

**是否需要迭代?**PDelta 本身不提供「一次性精確」解,它**在 Newton-Raphson 迴圈的每步 state determination 中更新 resisting force**,因此整個非線性靜力分析(OpenSees `analyze` 命令)通常需要多步 NR 迭代才能收斂。對單載步平衡解而言,若 P/P_cr 不大,2–3 次 NR 迭代即達機器精度;接近 P_cr 時迭代數大增或不收斂。[LIT:https://openseesdigital.com/2022/11/15/geometric-transformation/]

**當 oracle 使用時的實際含義:**
- OpenSees PDelta 是「線性相容 P-Δ 近似」(梁柱節點位移驅動的幾何效應),**不是**精確穩定函數解。
- 對懸臂柱 + 單軸壓 + tip load 這類 oracle 問題,若只有節點(一個元素),PDelta 計算的 tip deflection 偏小(缺少構件內 P-δ 貢獻),**與精確解差異在 O(P/P_cr) 量級**。
- 需用多個元素(mesh refinement)才能使 PDelta 收斂到精確穩定函數解;或切換 Corotational。[THEORY + LIT:https://pmc.ncbi.nlm.nih.gov/articles/PMC10719514/]

### 5. 收斂準則實務

**位移增量範數 vs 殘力範數:**
- OpenSees 官方建議:「unless you have very stiff elements…the norm of the residual vector should do the job.」有剛性元素時改用 displacement increment norm(`NormDispIncr`)。[LIT:https://openseesdigital.com/2021/02/28/norms-and-tolerance/]
- 容差換算:絕對 tolerance `tol`,模型 N 個自由度,各 DOF 等效誤差約 tol/√N;OpenSees 預設 1e-8 對小模型很保守。[LIT:https://openseesdigital.com/2021/02/28/norms-and-tolerance/]

**工程界常用容差:**
- RISA 3D P-Delta 預設收斂容差 **0.5%(相對位移)**。[LIT:https://structville.com/2020/07/p-delta-analysis.html]
- ETABS 提供「Convergence Tolerance (Relative)」控制迭代法 P-Delta(load combination iterative 模式);mass-based 模式不迭代。[LIT:https://docs.csiamerica.com/help-files/etabs/Menus/Analyze/Set_Analysis_Options/Initial_P-Delta_Analysis.htm]
- 普通工程級精度:相對殘力範數 ~ 1e-3 至 1e-4 或位移增量 0.1–0.5% 已足夠;研究/oracle 級:1e-6 至 1e-8。

**混合單位問題:**框架模型有旋轉 DOF 時殘力範數混合力(N)與彎矩(N·mm),量綱不統一,用 displacement increment norm 更穩健。[LIT:https://openseesdigital.com/2021/02/28/norms-and-tolerance/]

## 對 FrameCore 的含義

1. **凍結分解 pseudo-load 迭代的理論基礎已確認:**收斂域 P < P_cr,迭代比率 r ≈ P/P_cr,幾何級數求和,與文獻一致。FrameCore 目前已有線性屈曲(linearized K_g),可直接拿 P_cr 作迭代收斂判準——先用 `solveBuckling` 確認最小 λ_cr,若 max(P_i / λ_cr,i · P_ei) < 1 即保證偽載荷迭代收斂。

2. **幾何勁度 K_g 線性化誤差的誠實邊界:**現有 `assembleGeometric` 採 Hermitian K_g,在 P/P_cr = 0.9 時相對精確穩定函數解誤差 ~10%;設計使用範圍(P/P_cr < 0.5)誤差 < 1–2%。**文檔需誠實標**:「幾何勁度矩陣為線性化近似,P/P_cr 超過 0.5 後準確度下降,不保證收斂至精確穩定函數解。」

3. **精確 oracle 公式確認:**懸臂柱 tip deflection = H/(P) · (tan(kL)/k − L),k = √(P/EI)。這是對現有線性屈曲求解器(`solveBuckling`)和未來 P-Delta 迭代實作都適用的**機器精度驗證 oracle**。應加進 oracle 套件(F34 候選),與幾何勁度近似解做系統性比對(預期 P/P_cr = 0.1 時差距 <0.1%,0.5 時 ~1–2%,0.9 時 ~10%)。

4. **OpenSees PDeltaCrdTransf 作為 oracle 的限制:**它是 P-Δ(節點位移)近似,不包含 P-δ(構件內撓曲)。單元素懸臂柱比較時 OpenSees 與 FrameCore 幾何勁度版本應非常接近(兩者都是 P/L 側移修正),兩者共同對精確解有 ~O(P/P_cr) 的偏差。要得到精確解 oracle 需要多元素 mesh 或切換 Corotational 轉換。

5. **AISC B1 適用範圍:**只做構件內 P-δ,不做 P-Δ;P/P_e > 0.6 時誤差升至 5–15%,不宜當 oracle,僅供設計快篩。FrameCore 現有 `ElasticAllowable` D/C 不涉及二階效應,不受 B1 限制;但若未來要做「考慮 P-Delta 的強度覆核」,需明確標明是 B1 近似還是精確二階解。

6. **收斂準則建議:**若實作 pseudo-load 迭代 P-Delta,採**相對殘力範數 < 1e-4**(研究級)或 **displacement increment ratio < 0.1%**;對多 DOF 框架用 NormUnbalance,有剛性元素(殼+梁混合)用 NormDispIncr,避免剛度量綱混合問題。

## 來源清單

- Wilson E.L. & Habibullah A. (1987). *Static and Dynamic Analysis of Multi-Story Buildings, Including P-Delta Effects.* Earthquake Spectra 3(2):289–298. [https://pubs.geoscienceworld.org/eeri/earthquake-spectra/article-abstract/3/2/289/584240](https://pubs.geoscienceworld.org/eeri/earthquake-spectra/article-abstract/3/2/289/584240)
- Semantic Scholar abstract: [https://www.semanticscholar.org/paper/Static-and-Dynamic-Analysis-of-Multi-Story-P-Delta-Wilson-Habibullah/0ae7142d90889c03e1a025dd352c053ad673141c](https://www.semanticscholar.org/paper/Static-and-Dynamic-Analysis-of-Multi-Story-P-Delta-Wilson-Habibullah/0ae7142d90889c03e1a025dd352c053ad673141c)
- McGuire W., Gallagher R.H. & Ziemian R.D. (2000). *Matrix Structural Analysis*, 2nd ed. Wiley. (免費 e-book) [https://digitalcommons.bucknell.edu/books/7/](https://digitalcommons.bucknell.edu/books/7/)
- MASTAN2 textbook: [https://mastan2.com/textbook.html](https://mastan2.com/textbook.html)
- Chen W.F. & Lui E.M. (1987). *Structural Stability: Theory and Implementation.* Elsevier. [https://www.scribd.com/document/440700188/1987-Chen-Lui-Structural-Stability-Theory-and-Implementation-pdf](https://www.scribd.com/document/440700188/1987-Chen-Lui-Structural-Stability-Theory-and-Implementation-pdf)
- Timoshenko S.P. & Gere J.M. (1961). *Theory of Elastic Stability*, 2nd ed. McGraw-Hill. (標準教科書,無合法公開全文)
- Evaluation of P-Δ effect (two-cycle method, beam-column ODE exact solution), PMC: [https://pmc.ncbi.nlm.nih.gov/articles/PMC10719514/](https://pmc.ncbi.nlm.nih.gov/articles/PMC10719514/)
- ScienceDirect同文: [https://www.sciencedirect.com/article/pii/S2215016123002455](https://www.sciencedirect.com/article/pii/S2215016123002455)
- Denavit M.D. & Hajjar J.F. (2013). *Description of Geometric Nonlinearity for Beam-Column Analysis in OpenSEES.* Report NEU-CEE-2013-02, Northeastern Univ. [https://repository.library.northeastern.edu/files/neu:376268](https://repository.library.northeastern.edu/files/neu:376268)
- OpenSees PDelta Transformation (官方文檔): [https://opensees.github.io/OpenSeesDocumentation/user/manual/model/geomTransf/PDelta.html](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/geomTransf/PDelta.html)
- OpenSees PDelta Transformation (OpenSeesWiki): [https://opensees.berkeley.edu/wiki/index.php/PDelta_Transformation](https://opensees.berkeley.edu/wiki/index.php/PDelta_Transformation)
- Geometric Transformation (Portwood Digital, OpenSees): [https://openseesdigital.com/2022/11/15/geometric-transformation/](https://openseesdigital.com/2022/11/15/geometric-transformation/)
- Norms and Tolerance (Portwood Digital, OpenSees): [https://openseesdigital.com/2021/02/28/norms-and-tolerance/](https://openseesdigital.com/2021/02/28/norms-and-tolerance/)
- Improved accuracy for the Cm factor of steel beam-columns: [https://www.sciencedirect.com/article/abs/pii/S0141029615005295](https://www.sciencedirect.com/article/abs/pii/S0141029615005295)
- Moment Amplification Factor for P-δ Effect (ASCE J. Struct. Eng. 1999): [https://ascelibrary.org/doi/10.1061/(ASCE)0733-9445(1999)125:2(219)](https://ascelibrary.org/doi/10.1061/(ASCE)0733-9445(1999)125:2(219))
- AISC 360-16 Specification (Commentary Ch.C / Appendix 8): [https://www.aisc.org/globalassets/aisc/publications/standards/a360-16-spec-and-commentary_june-2018.pdf](https://www.aisc.org/globalassets/aisc/publications/standards/a360-16-spec-and-commentary_june-2018.pdf)
- EngineeringSkills P-Delta tutorial: [https://www.engineeringskills.com/posts/p-delta-analysis-geometric-non-linearity](https://www.engineeringskills.com/posts/p-delta-analysis-geometric-non-linearity)
- Structville P-Delta analysis: [https://structville.com/2020/07/p-delta-analysis.html](https://structville.com/2020/07/p-delta-analysis.html)
- ETABS Initial P-Delta analysis parameters: [https://docs.csiamerica.com/help-files/etabs/Menus/Analyze/Set_Analysis_Options/Initial_P-Delta_Analysis.htm](https://docs.csiamerica.com/help-files/etabs/Menus/Analyze/Set_Analysis_Options/Initial_P-Delta_Analysis.htm)