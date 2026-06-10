> 出處:研究輪 fan-out agent 查證報告(2026-06-10,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

現在我有足夠資料可以撰寫完整的研究報告了。

# WS-E 板殼單元路線(DKQ / QM6 / MITC9i / Karamba 殼)

## 摘要

本報告查證四個子主題:① DKQ(Batoz & Tahar 1982)的公式結構與適用邊界;② QM6 不協調膜模式對膜鎖定的改善機理與 patch test 性質;③ MITC9i(Wisniewski & Turska 2018)相對 MITC9 的改進及收斂階;④ Karamba3D 採用 TRIC 三角形(Argyris 等 1997/2000)而非 MITC4 四邊形,且忽略橫向剪切。另確認 OpenSees 有 ShellDKGQ/ShellNLDKGQ(DKQ 板彎+GQ12 膜,可當薄板 oracle)與 ShellMITC9(九節點 MITC 擴展,可當更高精度 oracle)。FrameCore 目前的 MITC4 flat-shell 在四邊形殼領域已優於 Karamba3D TRIC,加 QM6 改進膜可解決膜鎖定主因,MITC9i 屬研究級複雜度,近期不值得實作。

---

## 發現

### F-E1 DKQ 公式結構

- DKQ(Discrete Kirchhoff Quadrilateral)由 Batoz & Tahar(1982)提出,為 4 節點、每節點 3 DOF(1 撓度 w + 2 旋轉 θₓ,θ_y),全元素共 **12 DOF**。 [LIT: https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620181106]
- 核心思想:在元素邊上特定離散點強制施加 Kirchhoff 假設(零橫向剪應變),不直接引入剪切自由度;約束使旋轉場以撓度場的三次插值唯一決定,等效消除橫向剪應變 DOF,保留純彎曲行為。[LIT: https://www.researchgate.net/publication/268819596_On_Discrete-Kirchhoff_Plate_Finite_Elements_Implementation_and_Discretization_Error]
- 無獨立橫向剪應變場:**不適用厚板**(t/L > ~1/20 時 Kirchhoff 假設失效,無法計算 Reissner-Mindlin 效應)。[LIT: https://www.sciencedirect.com/science/article/abs/pii/S0045794917317078]
- 實作複雜度:需沿元素 4 條邊、每邊若干離散點上設定 Kirchhoff 約束方程式,再做 static condensation 降維到 12 DOF 勁度矩陣;比 MITC4 的 assumed natural strain(ANS)取樣稍複雜,但比 MITC9i 簡單得多。[THEORY]
- 對薄板(t/L < 1/50)精度等同甚至略優於 MITC4:文獻顯示四邊形薄板元素 DKQ、MITC4、DKMQ、DSQ 在薄板極限下全達最優收斂;薄板時 DKMQ ≈ DKQ,厚板時 DKMQ ≈ MITC4。[LIT: https://www.sciencedirect.com/science/article/abs/pii/S0045794917317078]
- **DKQ 完全無厚板能力**;MITC4 兼顧薄/厚板。對混合薄厚場景(結構工程常見:板厚 t/L ~1/10~1/50),MITC4 是更安全選擇。[LIT: https://www.sciencedirect.com/science/article/abs/pii/S0045794917317078]

### F-E2 QM6 不協調膜模式

- **起源**:Wilson, Taylor, Doherty & Ghaboussi(1973)提出 Q6 元素:在 4 節點雙線性膜(Q4)基礎上,引入 2 個內部不協調模式 φ₁ = 1−ξ², φ₂ = 1−η²(各自控制 x 方向和 y 方向的抛物線型應變增強),讓 4 節點元素能精確表達純彎曲。[LIT: Incompatible Displacement Models, Wilson et al. 1973, Semantic Scholar https://www.semanticscholar.org/paper/Incompatible-Displacement-Models-Wilson-Taylor/50340679a2c21fcf3166824a18202868845e65b5]
- **QM6 修正**:Q6 無法嚴格通過 patch test(常應力狀態下誤差不消除)。Taylor, Beresford & Wilson(1976)修改數值積分方案,使 QM6 **嚴格通過 patch test**,同時保留純彎曲表達精度。[LIT: A non-conforming element for stress analysis, Taylor et al. 1976, https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620100602]
- **實作要點**:2 個不協調內部 DOF 在元素層做 static condensation 消除,最終勁度矩陣仍為標準 8×8(4 節點×2 DOF 平面,或 flat-shell 中的 8 DOF 膜部分)。凝縮不引入全域自由度。[THEORY]
- **對矩形元素**:可**精確表達純彎曲**。對扭曲元素:精度隨扭曲程度下降,但仍比 Q4 顯著更好。[LIT: https://www.studocu.com/no/document/norges-teknisk-naturvitenskapelige-universitet/elementmetoden-i-styrkeanalyse/lecture-notes-lecture-13-incompatible-elements/176276]
- **改善幅度**:文獻顯示對 cantilever beam with quad mesh 測試,Q4 相對誤差可達 ~72%,而 QM6 近乎精確(幾近 exact)。現代對比中 QM6 與精確解的差距通常 < 1-2%,視扭曲程度。[LIT: https://www.researchgate.net/figure/Stress-distributions-calculated-using-elements-Q4-QM6-and-SU4-for-the-cantilever-beam_fig20_354107158]
- **Patch test 注意**:QM6 通過**弱 patch test**(modified integration scheme),不通過**嚴格 patch test**;實務上表現優良,但在高度扭曲網格下相對誤差仍存在。[LIT: https://www.hindawi.com/journals/mpe/2013/901495/]
- **與 MITC4 膜的關係**:OpenSees `ShellMITC4` 的膜是原始雙線性(標準 Q4 等效),不含不協調增強;因此 QM6 改進膜可以**顯著改善 in-plane bending 精度**,這正是 FrameCore CLAUDE.md 提到的「膜鎖定主因是 QM6 可行」結論的文獻依據。[LIT: https://openseesdigital.com/2019/12/20/opensees-shells-by-the-seashore/]

### F-E3 MITC9i(Wisniewski & Turska 2018)

- **發表**:Wisniewski & Turska, "Improved nine-node shell element MITC9i with reduced distortion sensitivity," *Computational Mechanics* **62**:499–523, 2018. DOI: 10.1007/s00466-017-1510-4. [LIT: https://link.springer.com/article/10.1007/s00466-017-1510-4]
- **規模**:9 節點、每節點 6 DOF(3 平移 + 3 旋轉含 drilling),全元素 **54 DOF**。[LIT: https://ui.adsabs.harvard.edu/abs/2018CompM..62..499W/abstract]
- **相對 MITC9 的改進**:
  - 採用 Celia & Gray(1984)的修正形函數(corrected shape functions)取代標準 Serendipity/Lagrangian 九節點形函數,消除節點位置偏移對 patch test 的破壞。[LIT: https://link.springer.com/article/10.1007/s00466-017-1510-4]
  - 改進的橫向剪應變 MITC 轉換(improved bending and transverse shear transformations)。
  - 結果:在**規則網格**下通過所有 patch test;在**中側節點沿直邊偏移**及**中央節點任意偏移**的扭曲網格下仍通過 patch test——這是原始 MITC9 做不到的。[LIT: https://www.researchgate.net/publication/321074097_Improved_nine-node_shell_element_MITC9i_with_reduced_distortion_sensitivity]
- **積分方案**:採用 3×3 Gauss 積分;橫向剪應變採樣由 6 取樣點改為 2 取樣線(two sampling lines)以降低計算量。[LIT: https://ui.adsabs.harvard.edu/abs/2018CompM..62..499W/abstract]
- **收斂階**:9 節點二次元素理論上比 4 節點線性(MITC4)高一階;文獻確認 MITC9i 無鎖定、精度優異、對節點偏移不敏感。[LIT: https://link.springer.com/article/10.1007/s00466-017-1510-4]
- **實作複雜度評估**:54 DOF 本地勁度矩陣、修正形函數、改進剪切轉換三重疊加,複雜度遠高於 MITC4(24 DOF)和 DKQ(12 DOF)。需獨立實作並驗證 Celia-Gray 形函數修正。[THEORY]
- **分離文獻**:Wisniewski & Turska 2019/2020 有後續 "On Transverse Shear Strains Treatment in Nine-Node Shell Element MITC9i"(Springer LNACM v.86),進一步分析剪切應變處理策略。[LIT: https://link.springer.com/chapter/10.1007/978-3-030-47491-1_23]

### F-E4 Karamba3D 的殼單元

- **官方手冊明確陳述**:"The shell elements used in Karamba3D resemble the TRIC-element devised by Argyris and coworkers." [LIT: https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.7-create-surface-element/3.1.9-mesh-to-shell]
- **TRIC 描述**:3 節點三角形,每節點 6 DOF,全元素 18 DOF;基於 6 剛體模式 + 12 應變模式;積分**閉合式**(無數值積分);包含 CST 膜 + DKT 板彎曲成分。[LIT: https://www.sciencedirect.com/science/article/abs/pii/S0045782596012339]
- **Karamba3D 的橫向剪切**:手冊明確說「**neglects transverse shear deformation** in case of shell elements」—— Kirchhoff 假設,薄殼適用、厚殼誤差。[LIT: https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.7-create-surface-element/3.1.9-mesh-to-shell]
- **Karamba3D drilling DOF**:手冊說殼元素「lack bending stiffness for rotations about the shell normal」——即 drilling DOF **無剛度**,平坦面場景下出現剛體模式,需約束或忽略警告。[LIT: 同上]
- **原始論文**:J.H. Argyris, L. Tenek, L. Olofsson, *CMAME* 145:11–85, 1997;J.H. Argyris, M. Papadrakakis et al., *CMAME* 182:217–245, 2000。[LIT: https://manual.karamba3d.com/appendix/bibliography]
- **對比 FrameCore MITC4+QM6**:
  - 四邊形 vs 三角形:四邊形 MITC4 每節點 DOF 相同(6),但同網格面積下四邊形收斂通常優於三角形(減少自由度鎖定的「應力帽」效應)。[THEORY]
  - 厚殼能力:MITC4 含 Mindlin 橫向剪切 → 優於 Karamba TRIC(忽略剪切)。
  - Drilling DOF:MITC4 的 Hughes-Brezzi drilling(FrameCore 現有)比 Karamba TRIC 的無 drilling 剛度更穩健。
  - 膜精度:若加 QM6 改進膜,FrameCore 在 in-plane bending 場景下會顯著超越 Karamba TRIC(CST 膜=常應力三角形,是最低階膜)。
  - **結論:FrameCore MITC4 已整體對標或超越 Karamba3D TRIC;加 QM6 後在膜精度上明確超越。**

### F-E5 OpenSees 殼單元全覽(oracle 可用性)

| 元素 | 節點 | 彎曲 | 膜 | 幾何非線性 | 狀態 |
|---|---|---|---|---|---|
| ShellMITC4 | 4 | MITC4 | 雙線性 Q4 | 無 | 正式 oracle |
| ShellDKGQ | 4 | DKQ 薄板 | GQ12+drilling | 無 | **可當薄板 oracle** |
| ShellNLDKGQ | 4 | DKQ 薄板 | GQ12+drilling | 更新 Lagrangian | 幾何非線性版 |
| ShellDKGT | 3 | DKT 薄板 | GQ6+drilling | 無 | 三角形版 |
| ShellNLDKGT | 3 | DKT 薄板 | GQ6+drilling | 更新 Lagrangian | 三角形非線性版 |
| ShellMITC9 | 9 | MITC9 | 雙二次 | 不明 | **可當高精度 oracle** |

[LIT: https://openseesdigital.com/2019/12/20/opensees-shells-by-the-seashore/; https://opensees.berkeley.edu/wiki/index.php/ShellDKGQ; https://opensees.berkeley.edu/wiki/index.php/ShellNLDKGQ; https://deepwiki.com/OpenSees/OpenSees/6.7-shell-and-solid-elements]

- **ShellDKGQ** = DKQ(Discrete Kirchhoff Quadrilateral)板彎 + GQ12(generalized conforming quadrilateral)膜+drilling;Tsinghua/Xiamen 大學 Lu 等開發。可作為 FrameCore 薄板 Kirchhoff 解析 oracle。[LIT: https://opensees.berkeley.edu/wiki/index.php/ShellDKGQ]
- **ShellMITC9** = Tesser & Talledo 開發的九節點 MITC 殼元素,是 OpenSees 唯一九節點殼;可作為比 ShellMITC4 更高精度的 oracle。**注意**:DeepWiki 說主線 OpenSees 目前只列 6 個殼元素(不含 ShellMITC9),ShellMITC9 可能在較早版本或特定 fork;需確認安裝版本是否支援。[LIT: https://deepwiki.com/OpenSees/OpenSees/6.7-shell-and-solid-elements 與 https://openseesdigital.com/2019/12/20/opensees-shells-by-the-seashore/; 注:兩源出現矛盾,deepwiki 未列 ShellMITC9,portwood 文提及;需在本機 openseespy 確認可用性]
- **無 DKQ 獨立殼單元** 命名(無 `ShellDKQ`),DKQ 彎曲以 `ShellDKGQ` 形式嵌合膜一起出現。[LIT: https://deepwiki.com/OpenSees/OpenSees/6.7-shell-and-solid-elements]

---

## 對 FrameCore 的含義

**採什麼:**

1. **QM6 改進膜(近期最值得做)**:將現有 MITC4 的雙線性 Q4 膜部分換為 QM6(引入 2 個 φ₁=1−ξ², φ₂=1−η² 不協調模式、元素層 static condensation);實作代價小、精度提升顯著(特別是 in-plane bending 主導場景如薄牆弦)。Oracle 策略:矩形元素有精確解;OpenSees ShellDKGQ(含 GQ12 膜)可作近似對標。注意 CLAUDE.md 的警告:預設行為不能破壞 OpenSees `ShellMITC4` 閘門,應用 opt-in 旗標或新類別名稱。

2. **OpenSees ShellDKGQ 作薄板 oracle**:當 t/L < 1/50 場景需要 Kirchhoff 限解驗證時,可呼叫 `openseespy` 的 `ShellDKGQ` 當基準;比 ShellMITC4 在超薄場景更「乾淨」。

3. **OpenSees ShellMITC9 作高精度 oracle(如可用)**:在需要比 4 節點更高收斂階的基準時,用 ShellMITC9(9 節點二次);先在本機確認 `openseespy` 版本是否支援。

**避什麼:**

1. **DKQ 作 FrameCore 主要殼元素**:DKQ 完全無厚板能力(無橫向剪切),而結構工程常出現 t/L ~1/20~1/5 的中厚板;引入 DKQ 需要平行維護薄/厚版本,邊界判斷邏輯複雜。保持 MITC4 為主單元、DKQ 僅作 oracle。

2. **MITC9i 近期實作**:54 DOF + 修正形函數 + 改進剪切轉換,實作工作量約為 MITC4 的 3-4 倍,需新的 oracle 策略(OpenSees ShellMITC9 是舊版 MITC9 非 MITC9i);研究級元素,ROI 不高。

3. **採用 Karamba3D TRIC 路線**:三角形在相同網格密度下精度劣於四邊形;無橫向剪切;drilling 無剛度。FrameCore 的 MITC4 已優於此路線。

**誠實邊界:**

- QM6 改進膜對矩形元素可精確表達純彎曲,對扭曲元素精度隨形狀比下降,高度扭曲網格仍有誤差,需在 oracle 比較中加扭曲網格測試案例。
- DKQ 的「厚板不適用」界限文獻上約為 t/L > 1/20;在 1/50 < t/L < 1/20 的「中等薄板」區間,MITC4 和 DKQ 精度接近但 MITC4 更保守安全。
- ShellMITC9 在 OpenSees 的可用性存疑(兩份來源矛盾),使用前需 `import openseespy; help('ShellMITC9')` 或 `ops.element('ShellMITC9',...) ` 測試確認。

---

## 來源清單

1. Batoz & Tahar (1982), "Evaluation of a new quadrilateral thin plate bending element," *IJNME* 18:1655–1677. [https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620181106](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620181106)
2. Taylor, Beresford & Wilson (1976), "A non-conforming element for stress analysis," *IJNME* 10:1211–1219. [https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620100602](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620100602)
3. Wilson, Taylor, Doherty & Ghaboussi (1973), "Incompatible Displacement Models." Semantic Scholar: [https://www.semanticscholar.org/paper/Incompatible-Displacement-Models-Wilson-Taylor/50340679a2c21fcf3166824a18202868845e65b5](https://www.semanticscholar.org/paper/Incompatible-Displacement-Models-Wilson-Taylor/50340679a2c21fcf3166824a18202868845e65b5)
4. Wisniewski & Turska (2018), "Improved nine-node shell element MITC9i with reduced distortion sensitivity," *Computational Mechanics* 62:499–523. [https://link.springer.com/article/10.1007/s00466-017-1510-4](https://link.springer.com/article/10.1007/s00466-017-1510-4)
5. Wisniewski & Turska (2019/2020), "On Transverse Shear Strains Treatment in Nine-Node Shell Element MITC9i." [https://link.springer.com/chapter/10.1007/978-3-030-47491-1_23](https://link.springer.com/chapter/10.1007/978-3-030-47491-1_23)
6. Katili et al. (2018), "A comparative formulation of DKMQ, DSQ and MITC4 quadrilateral plate elements with new numerical results based on s-norm tests," *Computers & Structures* 189. [https://www.sciencedirect.com/science/article/abs/pii/S0045794917317078](https://www.sciencedirect.com/science/article/abs/pii/S0045794917317078)
7. Argyris, Tenek & Olofsson (1997), "TRIC: a simple but sophisticated 3-node triangular element," *CMAME* 145:11–85. [https://www.sciencedirect.com/science/article/abs/pii/S0045782596012339](https://www.sciencedirect.com/science/article/abs/pii/S0045782596012339)
8. Argyris, Papadrakakis et al. (2000), "The TRIC shell element: theoretical and numerical investigation," *CMAME* 182:217–245. [https://www.researchgate.net/publication/222217940_The_TRIC_shell_element_Theoretical_and_numerical_investigation](https://www.researchgate.net/publication/222217940_The_TRIC_shell_element_Theoretical_and_numerical_investigation)
9. Karamba3D v3 Manual, Mesh to Shell. [https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.7-create-surface-element/3.1.9-mesh-to-shell](https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.7-create-surface-element/3.1.9-mesh-to-shell)
10. Karamba3D v3 Bibliography. [https://manual.karamba3d.com/appendix/bibliography](https://manual.karamba3d.com/appendix/bibliography)
11. OpenSees ShellDKGQ Wiki. [https://opensees.berkeley.edu/wiki/index.php/ShellDKGQ](https://opensees.berkeley.edu/wiki/index.php/ShellDKGQ)
12. OpenSees ShellNLDKGQ Wiki. [https://opensees.berkeley.edu/wiki/index.php/ShellNLDKGQ](https://opensees.berkeley.edu/wiki/index.php/ShellNLDKGQ)
13. OpenSees Digital, "OpenSees Shells by the Seashore" (2019). [https://openseesdigital.com/2019/12/20/opensees-shells-by-the-seashore/](https://openseesdigital.com/2019/12/20/opensees-shells-by-the-seashore/)
14. DeepWiki, OpenSees Shell and Solid Elements. [https://deepwiki.com/OpenSees/OpenSees/6.7-shell-and-solid-elements](https://deepwiki.com/OpenSees/OpenSees/6.7-shell-and-solid-elements)
15. NTNU Lecture Notes, Lecture 13: Incompatible Elements – The Q6 Element. [https://www.studocu.com/no/document/norges-teknisk-naturvitenskapelige-universitet/elementmetoden-i-styrkeanalyse/lecture-notes-lecture-13-incompatible-elements/176276](https://www.studocu.com/no/document/norges-teknisk-naturvitenskapelige-universitet/elementmetoden-i-styrkeanalyse/lecture-notes-lecture-13-incompatible-elements/176276)
16. Hindawi MPE (2013), non-conforming plane elements patch test study. [https://www.hindawi.com/journals/mpe/2013/901495/](https://www.hindawi.com/journals/mpe/2013/901495/)
17. Semantic Scholar, Batoz & Tahar 1982 record. [https://www.semanticscholar.org/paper/Evaluation-of-a-new-quadrilateral-thin-plate-Batoz-Tahar/c8964c3baf6c5cb2aba87ad280812c511ce523c6](https://www.semanticscholar.org/paper/Evaluation-of-a-new-quadrilateral-thin-plate-Batoz-Tahar/c8964c3baf6c5cb2aba87ad280812c511ce523c6)
