> 出處:研究輪 fan-out agent 查證報告(2026-06-10,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

基於以上查證,現在整理完整的報告。

# WS-A1 Karamba3D v3 功能矩陣(官方手冊查證)

## 摘要

Karamba3D v3 是嵌入 Grasshopper/Rhino 的互動式 FEA 外掛,以 Timoshenko-Ehrenfest 梁單元 + TRIC 三角殼單元為核心,提供一到三階靜力分析、模態/屈曲分析、BESO 拓撲優化、EC3 截面優化等元件。**無材料非線性(塑鉸/pushover)**、**無反應譜**、**無時程動力分析**。免費試用版限制 ≤20 梁元素,殼元素上限約 50。大變形分析(LaDeform)與幾何非線性(Analyze Nonlinear WIP)均已提供但有明確精度限制。整體定位:設計探索工具而非 ETABS/SAP2000 級驗證工具。

---

## 發現

### 1. 分析類型

**1.1 Analyze(一階,Th.I)**
線性靜力,直接剛度法,處理全部載重工況與組合。輸出:最大節點位移(cm)、最大合力(kN)、內部變形能。適用所有元素類型(梁、桁架、殼、彈簧)。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.1-analyze.md]

**1.2 AnalyzeThII(二階,Th.II)**
考慮幾何非線性:軸壓力使系統軟化、軸拉力使系統硬化;以 N^II 修正幾何勁度矩陣,迭代至收斂。控制參數:RTol(相對位移增量容差)、MaxIter(預設50)。收斂失敗時元件標示橙色。可選 NoTenNII(僅取壓力 N^II)或 NoComNII(僅取拉力,用於薄膜防皺)。CombiNII 控制是否跨工況取最保守值。**屬 P-Δ 級別而非真大變形**。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.2-analyzethii.md]

**1.3 Analyze Nonlinear WIP(幾何非線性,標記 WIP)**
提供三種求解算法:Dynamic Relaxation、Newton-Raphson、Arc-Length。允許任意大位移,前提是應變仍小。**明確標注「殼的情況下可能在可接受時間內不收斂」**。收斂後結果數學上正確。**非材料非線性**。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.3-analyze-nonlinear-wip.md]

**1.4 Large Deformation Analysis(LaDeform,懸鏈/找形)**
純增量法:每步施加部分載重後更新幾何。**無內力輸出**——因為「這些性質會非常不準確」。主要用途是找形(仿 Gaudi/Isler 懸掛模型)。精度隨步數提升但有累積誤差。座標系在元素過垂直時可能翻轉。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.4-analyze-large-deformation.md]

**1.5 Buckling Modes(線性屈曲)**
輸入二階法的 N^II,解廣義特徵值問題 Ĉ·x + λ²·ĈG·x = 0,輸出挫屈乘數(ascending order)與挫屈形態。適用梁、桁架、殼。NModes 參數控制計算模態數,預設1。假設小位移至不穩定點。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.5-buckling-modes.md] [LIT:https://manual.karamba3d.com/appendix/a.4-background-information/a.4.5-natural-vibrations-eigen-modes-and-buckling.md]

**1.6 Natural Vibrations / Eigen Modes(模態分析)**
解廣義特徵值問題 Ĉ·x = ω²·M̂·x。梁用一致質量矩陣,桁架與殼用集中質量法。點質量只有平動慣量(無旋轉慣量)。輸出振型(最大分量歸一化為1)、模態質量、參與係數。NaturalVibrations 與 EigenModes 為兩個獨立元件(NaturalVibrations 含頻率輸出)。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.7-natural-vibrations.md] [LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.6-eigen-modes.md]

**1.7 時程 / 反應譜**
**無**——手冊明確回答 "I can't find any docs content about response spectrum or time history dynamic analysis in this site"。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.7-natural-vibrations.md?ask=response+spectrum+time+history+dynamic+analysis]

---

### 2. 元素類型

**2.1 梁(Beam)**
公式:**Timoshenko-Ehrenfest 梁理論**,考慮剪力變形;勁度矩陣使用 Helmut Rubin 的冪級數法計算,支援二階效應。可用 Line to Beam 或 Line to Truss 建立(桁架=僅軸力,等效於關閉梁的彎曲)。[LIT:https://manual.karamba3d.com/troubleshooting/4.3.-miscellaneous-problems/4.1.0-faq.md?ask=beam+element+formulation+Timoshenko]

**2.2 殼(Shell)**
公式:**TRIC 三角形元素**(Argyris & coworkers 1997/2000)。**任意四邊形自動拆成兩個三角形**。手冊明確:「Karamba3D neglects transverse shear deformation in case of shell elements」→ **薄殼假設(Kirchhoff-Love 等效),無 MITC 剪力鎖定補正**。Drilling DOF(殼法向旋轉)無彎曲剛度,屬數值處理層面。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.7-create-surface-element/3.1.9-mesh-to-shell.md] [THEORY:Argyris et al. 1997, Comput. Methods Appl. Mech. Engrg. 145:11–85]

**2.3 薄膜(Membrane)**
Mesh to Shell 加 Bending=False 等效於 Mesh to Membrane;純膜狀態,無彎曲剛度。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.7-create-surface-element/3.1.7.2-mesh-to-membrane.md]

**2.4 彈簧截面(Spring)**
可用 Spring Cross Sections 定義線彈簧元素。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.3-cross-section/3.3.3-spring-cross-sections.md]

**2.5 支承與 Release**
支承可設剛性或彈性彈簧(平動 Ct [kN/m]、旋轉 Cr [kNm/rad]),各方向獨立。樑端 joint/release:可全釋放或指定彈簧剛度(Ct-start/end, Cr-start/end),支援部分釋放。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.16-support.md] [LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.4-joint/3.3.6-beam-joints.md]

**2.6 Tension/Compression-Only**
由 Tension/Compression Eliminator 元件實作:**迭代到收斂**(不是一次性移除)。每步評估並以「soft-kill(極小剛度)」移除違反條件的元素,直至無變化或達 MaxIter。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.12-tension-compression-eliminator.md]

---

### 3. 優化元件

**3.1 Optimize Cross Section(截面優化)**
檢核標準:**Eurocode 3 (EN 1993-1-1)**,含:
- 受壓屈曲(buckling)
- 橫向扭轉屈曲(Lateral Torsional Buckling, LTB),屈曲長度自動計算或手動 BklLenLT 覆寫
- 力的疊加:EC3 Annex B Procedure 2,交互係數 C_my、C_mz、C_mLT 上限 0.9
迭代邏輯：兩層——外層位移迭代(滿足撓度限制)、內層利用率迭代(預設5次)。截面從 CrossSectionValues.bin 離散表中按高度升序選取第一個滿足的截面。**無 AISC 支援**。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.8-optimize-cross-section.md] [LIT:https://manual.karamba3d.com/appendix/a.4-background-information/a.4.6-approach-used-for-cross-section-optimization.md]

**3.2 BESO for Beams**
參數：TargetRatio(目標質量比)、MaxChangeIter(達目標所需迭代數)、MaxConvIter(達目標後額外收斂迭代)。敏感度依據：各力分量(軸力、剪、彎矩)的變形能密度加權(WTension, WCompr., WShear, WMoment)。BESOFac 啟用雙向(先多移再重啟活化)。MinDist 限制空間叢聚。WLimit 移除低於臨界利用率元素。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.9-beso-for-beams.md]

**3.3 BESO for Shells**
參數：TargetRatio、MaxIter。ER(演化率)公式：V_{i+1} = V_i·(1±ER)；預設值 ER = (1-TargetRatio)/MaxIter + AR_max/2 自動計算。敏感度過濾半徑 Rmin 預設 = √(totalArea/numberOfElements)·2(以 Rexp 作距離衰減)。**Soft-kill 減厚(KillThick 預設 0.00001m)而非直接移除元素**。Nhist/Conv 控制收斂判定。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.10-beso-for-shells.md]

**3.4 Optimize Reinforcement(RC 配筋優化)**
基於 Marti sandwich model,假設混凝土無拉力、無限壓力,鋼筋強度γ_s = 1.15(EC2)。分析仍為線彈性——使用彈性截面力換算配筋量,**無材料非線性**。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.11-optimize-reinforcement.md]

---

### 4. 載重

**通用載重**：重力(Gravity)、點力/點矩(Point-Load)、prescribed 位移(Point-Displacement,需配合支承)、初始應變(Initial Strain-Load,Eps0+Kappa0)、溫度載重(T+ΔT)、網格均布(Mesh Load Constant kN/m²)、網格變布(Mesh Load Variable)、預應力(Pretension)、不完美初始缺陷(Imperfection)。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.2-load/3.2.1-loads.md]

**梁元素載重**：集中力/矩(任意位置 t)、均布塊載(Block Load)、折線分布(Polylinear)、梯形(Trapezoidal)、缺口/不連續(Gap Load,用於影響線)、幾何不完美(Imperfection)。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.2-load/3.2.2-beam-loads.md]

**載重工況組合**：Load-Case-Combinator 支援字串規則展開(含 OR 運算、遞歸引用、正規表達式匹配)。**無類似 FrameCore `envelope()` 的包絡元件**——需外部程序疊加多組合結果。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.2-load/3.2.4-load-case-combinations/3.2.5.1-load-case-combinator.md]

---

### 5. 免費版限制

試用(temporary)授權上限：**≤20 梁元素**。殼元素隱含上限約 50(從「激活授權後可開啟超過 50 殼元素範例檔」反推)。試用授權需每次開啟 Rhino 執行 `karamba3dgetlicense` 從雲端取回，非永久安裝。商業版分 PRO/EDU/LAB 三層，均需付費。[LIT:https://manual.karamba3d.com/1-introduction/1.2-licenses.md] [LIT:https://manual.karamba3d.com/troubleshooting/4.3.-miscellaneous-problems/4.1.0-faq.md?ask=free+version+element+limit]

---

### 6. 材料非線性(塑鉸/塑性分析)

**無材料非線性。手冊明確回覆**："I cannot find information about that in the docs"(針對 plastic hinge + material nonlinearity 查詢)。具體確認：
- Optimize Cross Section 可設定塑性截面設計(Z vs S),但**分析仍假設線彈性**
- Beam-Joints 的轉動釋放/彈簧是機械鉸，不是塑性力矩-轉角模型
- Analyze Nonlinear WIP 僅處理幾何非線性
- 完全無 pushover、無 event-to-event 塑鉸、無纖維斷面
[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.8-optimize-cross-section.md?ask=plastic+hinge+material+nonlinear+plastic+analysis]

---

## 對 FrameCore 的含義

**可採的差異化方向(FrameCore 已有或優於 Karamba3D)**：

1. **Shell 公式**：Karamba3D 用 TRIC 三角形(忽略橫向剪力),FrameCore 用 MITC4 四邊形(有剪力鎖定補正,且 OpenSees ~1e-10)。對厚板/中厚殼精度明顯優勢;差距在 Karamba3D 四邊形須先拆三角,FrameCore 直接四邊形計算。這是可公開宣稱的定量差異。

2. **塑鉸 / 漸進倒塌**：FrameCore 已實作 event-to-event 塑鉸(w*=16Mp/L² ±2%)與 LSP 驅動器;Karamba3D 完全無。這是直接功能優勢。

3. **反應譜 / 模態疊加動力**：FrameCore 已有 SRSS/CQC + 模態疊加時域;Karamba3D 無。

4. **Factorize-once PreparedSystem 架構**：Karamba3D 手冊未提及類似「一次分解,多次 solve」的 API。FrameCore 的 `PreparedSystem` 對互動重解有明確速度優勢。

**需謹慎/避免過度宣稱的地方**：

5. **Analyze Nonlinear WIP 三算法**：Karamba3D 有 Newton-Raphson 與 Arc-Length,FrameCore 目前僅線性(no geometric NL solver)。若要競爭「大變形」場景,FrameCore 需補幾何非線性。

6. **EC3 截面優化(LTB,離散截面表)**：Karamba3D 已內建完整 EC3 流程+LTB;FrameCore 的 ElasticAllowable D/C 是彈性快篩,尚無 LTB 或標準截面表迭代。如定位「工程規範設計工具」需注意此差距。

7. **BESO 拓撲優化**：Karamba3D 有完整梁/殼 BESO;FrameCore 無,不在目前規劃。

**誠實邊界聲明建議**：
- FrameCore MITC4 vs Karamba3D TRIC:MITC4 在中等厚度板/殼精度更高(OpenSees ~1e-10),但兩者都是平面 facet 殼;曲面均隨網格細化收斂。
- FrameCore 線性分析套件完整度(8 階段含反應譜/模態疊加)在同類 Grasshopper 外掛中屬 unique。

---

## 來源清單

- Karamba3D v3 手冊首頁：https://manual.karamba3d.com
- Sitemap：https://manual.karamba3d.com/sitemap.md
- Analyze 元件：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.1-analyze.md
- AnalyzeThII 元件：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.2-analyzethii.md
- Analyze Nonlinear WIP：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.3-analyze-nonlinear-wip.md
- Large Deformation Analysis：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.4-analyze-large-deformation.md
- Buckling Modes：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.5-buckling-modes.md
- Eigen Modes：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.6-eigen-modes.md
- Natural Vibrations：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.7-natural-vibrations.md
- Optimize Cross Section：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.8-optimize-cross-section.md
- BESO for Beams：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.9-beso-for-beams.md
- BESO for Shells：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.10-beso-for-shells.md
- Optimize Reinforcement：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.11-optimize-reinforcement.md
- Tension/Compression Eliminator：https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.12-tension-compression-eliminator.md
- Mesh to Shell：https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.7-create-surface-element/3.1.9-mesh-to-shell.md
- Mesh to Membrane：https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.7-create-surface-element/3.1.7.2-mesh-to-membrane.md
- Beam Joints：https://manual.karamba3d.com/3-in-depth-component-reference/3.4-joint/3.3.6-beam-joints.md
- Support：https://manual.karamba3d.com/3-in-depth-component-reference/3.1-model/3.1.16-support.md
- Loads：https://manual.karamba3d.com/3-in-depth-component-reference/3.2-load/3.2.1-loads.md
- Beam Loads：https://manual.karamba3d.com/3-in-depth-component-reference/3.2-load/3.2.2-beam-loads.md
- Load Case Combinator：https://manual.karamba3d.com/3-in-depth-component-reference/3.2-load/3.2.4-load-case-combinations/3.2.5.1-load-case-combinator.md
- Cross Section Optimization Background：https://manual.karamba3d.com/appendix/a.4-background-information/a.4.6-approach-used-for-cross-section-optimization.md
- Natural Vibrations & Buckling Background：https://manual.karamba3d.com/appendix/a.4-background-information/a.4.5-natural-vibrations-eigen-modes-and-buckling.md
- Licenses：https://manual.karamba3d.com/1-introduction/1.2-licenses.md
- FAQ：https://manual.karamba3d.com/troubleshooting/4.3.-miscellaneous-problems/4.1.0-faq.md
- New in v3.1：https://manual.karamba3d.com/new-in-karamba3d-3.1.md
- Bibliography：https://manual.karamba3d.com/appendix/bibliography.md
- Argyris et al. 1997：Argyris, J.H., Tenek, L., and Olofsson, L. "Tric: a simple but sophisticated 3-node triangular element based on 6 rigid body and 12 straining modes." Comput. Methods Appl. Mech. Engrg. 145:11–85
- Argyris et al. 2000：Argyris, J.H., Papadrakakis, M., Apostolopoulou, C., and Koutsourelakis, S. "The tric shell element: theoretical and numerical investigation." Comput. Methods Appl. Mech. Engrg. 182:217–245
