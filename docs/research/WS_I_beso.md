> 出處:研究輪 fan-out agent 查證報告(2026-06-10 補跑,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

已蒐集足夠文獻。現在整合所有查到的資訊並撰寫報告。

# WS-I 拓撲優化(BESO for 桿系/殼)

## 摘要

BESO(Bi-directional Evolutionary Structural Optimization)是一種離散 0/1 heuristic 拓撲優化方法:以元素應變能密度(sensitivity number)排序,每次迭代移除(hard-kill)或軟化(soft-kill)低效元素,同時可補入高效元素;敏感度須跨迭代歷史平均以抑制振盪;收斂準則是最近 2N 迭代目標函數相對變化 < τ。對桿系的特化實作見 Karamba3D「BESO for Beams」,以截面力分量加權和為 sensitivity,支援 WTension/WCompr/WShear/WMoment 四權重。主要病態模式:孤立桿鏈、機構化——文獻已有連通性約束 BESO 先例。低體積率 compliance 暴增屬已知現象,文獻對策是 compliance-best 回退、早停或模擬退火接受概率。Michell 桁架提供已發布的 ground-structure 數值基準。Fail-safe TO(Jansen et al. 2014)是「逐元素移除倒塌約束」的直接先行技術。

## 發現

### 1. BESO 標準算法(Huang & Xie 2007/2010)

**Sensitivity number(元素應變能):**

```
α_i^e = (1/2) · u_i^T · K_i · u_i
```

即元素 i 在位移場 u_i 下的應變能;對線彈性最小 compliance 問題,此即節點位移對刪除元素 i 之目標函數一階近似。[LIT:Huang & Xie, *Convergent and mesh-independent solutions for BESO*, Finite Elements in Analysis and Design 43(14):1039–49, 2007; ResearchGate https://www.researchgate.net/publication/223369659]

**過濾(mesh-independency filter):**

節點敏感度 α_j^n 由周圍元素加權平均得出:

```
α_j^n = Σ_i (ω_i · α_i^e)
```

權重 ω_i 以元素中心到節點 j 的距離 r_{ij} 反比計算,保證 Σω_i = 1。然後元素從相鄰節點插值回元素敏感度。[LIT:同上]

**歷史平均(sensitivity averaging):**

```
ᾱ_i^k = (α_i^k + α_i^(k-1)) / 2
```

當前迭代 k 與前一迭代 k-1 的平均;消除離散 0/1 設計變數造成的跳躍振盪。多迭代版本可推廣為等權滑動窗,但兩步平均是最常見形式。[LIT:Huang & Xie 2007; confirmed in etd.lib.metu.edu.tr/upload/12619334/index.pdf]

**收斂準則:**

```
change = |Σ_{i=1}^{N} C_{k-i+1} - Σ_{i=1}^{N} C_{k-N-i+1}| / Σ_{i=1}^{N} C_{k-i+1} ≤ τ
```

比較最近 N 迭代與再前 N 迭代的 compliance 總和之相對差;Huang & Xie 2010 書中取 N = 5,τ = 0.1%(即 τ = 0.001)。此準則須「達到目標體積 **且** 滿足收斂條件」兩者同時成立。[LIT:Huang & Xie 2010, *Evolutionary Topology Optimization of Continuum Structures*, Wiley, Ch.3; Wiley Online Library https://onlinelibrary.wiley.com/doi/10.1002/9780470689486.ch3]

**Evolution Rate ER:**

體積更新規則:

```
V_{k+1} = V_k · (1 ± ER)
```

ER 典型值 0.02 ~ 0.05(即每步移除/加入 2%~5% 體積);Karamba3D 殼手冊 ER 預設可由使用者設定。[LIT:Karamba3D v3 manual, BESO for Shells § Advanced Settings, https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.10-beso-for-shells]

**ARmax(最大加入比例):**

ARmax 限制每迭代可加入之元素數上限為總元素數之一定比例,防止雙向操作在某步驟過度加入元素、造成回頭蕩盪。Karamba3D 殼手冊明確列為可調參數。[LIT:同上]

**Soft-kill vs Hard-kill:**

- Hard-kill:直接刪除元素(K 矩陣不含該元素);計算成本低,但可能造成奇異子矩陣(懸空自由度)與不連通拓撲。
- Soft-kill:以極小剛度 E_min(e.g., E_min = 10^{-9} · E_0)替換被刪元素;防止 K 奇異、允許後續重新加入;Karamba3D 殼手冊的 `KillThick = 0.00001 m` 即此做法。
- 文獻建議:如需機構偵測或懸空節點處理,soft-kill + 後置 hard-kill 或 hybrid 模式(ASME 2024 GPU BESO 論文)優於純 hard-kill。[LIT:Ghabraie, *An improved soft-kill BESO algorithm*, SMO 52(4):1215-1233, 2015; https://link.springer.com/article/10.1007/s00158-015-1268-2]

### 2. 桿系病態模式與連通性過濾

**病態模式分類:**

- **細長桿鏈(slender chain)**:串列同向桿、中間節點僅受單方向束制 → 數學上連通但幾何上退化。
- **孤立桿(isolated member)**:移除後形成沒有載重路徑的孤懸元素。
- **機構化(mechanism)**:移除後出現零能量模態,K 奇異。

這些病態在低體積率下最常見;文獻稱「hinge cancellation」技術(連串同向桿合併為一根)為清理手段之一。[LIT:Ohsaki, *Topology and Geometry Optimization of Trusses and Frames*, Kyoto U report, http://www.se-lab.archi.kyoto-u.ac.jp/ohsaki/pdf/b0201.pdf]

**連通性約束 BESO 先例:**

Zuo & Xie 2016 在 BESO 中加入連通性約束(Finite Elements in Analysis and Design):透過廣度優先搜尋(BFS/DFS 連通分量)偵測不連通元素群,移除後再補 sensitivity 重計算,使 hard-kill 結果恢復準確。論文指出「透過連通性約束可偵測並修正不準確的敏感度分析」。[LIT:Zuo & Xie, *A BESO algorithm with connectivity constraint*, FEA&D, 2016; ScienceDirect https://www.sciencedirect.com/science/article/abs/pii/S0168874X16305789]

FrameCore 現有 `analyzeConnectivity`(union-find)完全對應此先例的核心偵測能力。

**低體積率 compliance 暴增與停機策略:**

文獻已確認這是 ESO/BESO 的已知問題:

- 原因:volume-target 驅動迭代在目標體積附近可能意外刪除載重路徑關鍵桿,compliance 在最後幾步跳升一到兩個量級(實驗觀察:30% 體積時 52× 暴增)。
- 文獻做法 A(**compliance-best 回退**):記錄每迭代 (volume, compliance) 歷史,以 compliance 最小時的拓撲為最優設計輸出,而非以「到達目標體積」作為終止判斷——此為最常見的工程做法。[LIT:ESO/BESO 教科書普遍做法;明確討論見 Topology Optimization Part 2, FETraining, https://fetraining.net/topology-optimization-part-2/]
- 文獻做法 B(**SA 接受概率**):SA-BESO 引入模擬退火機制,當 compliance 上升時以 exp(-ΔC/T) 概率接受,跳出局部最優;代價是額外參數(溫度調度)。[LIT:Topology Optimization Based on SA-BESO, *Applied Sciences* 13(7):4566, 2023; https://www.mdpi.com/2076-3417/13/7/4566 — 此篇有抓到,但 403 擋住詳文]
- 文獻做法 C(**動態 ER**):DER-BESO 在接近目標體積時自動縮小 ER,減緩每步刪除量。[LIT:DER-BESO, SMO 2020; https://link.springer.com/article/10.1007/s00158-020-02588-2]
- **誠實標**:文獻中直接以「12×6 ground structure 砍到 30% compliance 暴增 52×」這類具體數字描述的論文查不到對應文獻;但 compliance 暴增現象本身是公認 BESO 已知問題,做法 A(回退最佳歷史拓撲)有充分文獻支持。[UNKNOWN:具體 52× 數字的文獻支撐]

### 3. SIMP vs BESO for 桿系的文獻比較

| 維度 | SIMP(桁架) | BESO(桿系) |
|------|-----------|-----------|
| 設計變數 | 連續截面積 A_i ∈ [A_min, A_max] | 離散 x_i ∈ {0,1}(或 soft-kill x_min) |
| 罰函數 | E(ρ) = ρ^p · E_0,p ≈ 3 | 無罰函數;直接按 sensitivity 排序刪加 |
| 敏感度 | ∂C/∂ρ_i 解析連續 | α_i^e = 1/2 u_i^T K_i u_i,離散近似 |
| 中間密度 | 允許(penalized) | 不允許(hard-kill)或以 E_min 替代(soft-kill) |
| 後處理 | 需閾值截取(threshold) | 直接得 0/1 拓撲 |
| 桁架適用性 | 天然適合(截面積即密度) | 需 ground structure 預定義;sensitivity 較粗糙 |

文獻結論:ground structure + LP(線性規劃)是桁架最優拓撲的 **全局最優**方法(計算費但數學嚴謹);BESO 是 heuristic 近似,對中大規模桁架計算快但無全局最優保證。SIMP 連續放鬆在 FEM 網格上更自然,桁架用 BESO 的優點主要是「輸出即離散拓撲、無後處理截取偽影」。[LIT:Review in Topology Optimization Applications on Engineering Structures, IntechOpen, https://www.intechopen.com/chapters/70489]

### 4. Michell 桁架 oracle 與 ground-structure 數值基準

**解析基準(Michell 1904):**

Michell 最小重量桁架理論:在應力約束 ±σ_allow、單點載重 P 下,最優結構體積 V* 滿足:

```
V* = P · δ* / σ_allow
```

其中 δ* 為虛功意義下的最優節點位移。對 2D 懸臂(半無限域)已有解析解族,由 Hencky 網格(等角斜交)描述。

**地基結構數值基準(Dorn, Gomory & Greenberg 1964 方法):**

Dorn et al. 1964 論文「Automatic design of optimal structures」建立了 ground structure LP 方法,是後續所有數值 Michell 近似的理論根基。[LIT:Dorn, Gomory, Greenberg 1964, *Journal de Mécanique* 3:25–52; cited in Ohsaki 2011 and GRAND/GRAND3 papers]

具體已發表數值基準:

- **GRAND (Paulino group 2014)**:ground structure LP 二維任意域;Michell 懸臂 20×10 節點格,σ_allow = 1,P = 1;文獻展示收斂至 Michell 解析值。[LIT:Paulino group 2014, Ground structure based topology optimization for arbitrary 2D domains, SMO; URL已 404]
- **GRAND3 (Paulino 2015)**:三維版本,同一問題框架。[LIT:https://paulino.princeton.edu/journal_papers/2015/SMO_15_GRAND3-GroundStructure.pdf]
- **Michell cantilever half-strip 解析值**:Lewinski et al. 系列論文給出 nondimensional compliance/mass 解析值;例如在特定長寬比下非量綱質量 M = 2.30259 ~ 3.35889(視 L/H 比);這些是 ground structure 數值解的收斂目標。[LIT:Lewinski et al., *Michell cantilevers constructed within a half strip*, SMO, https://link.springer.com/article/10.1007/s00158-010-0525-7 — 需訂閱]
- **99 行 Mathematica 基準**:Sokół 2011,SMO,「A 99 line code for discretized Michell truss optimization」;文中給出具體 compliance 值供對比。[LIT:https://link.springer.com/article/10.1007/s00158-010-0557-z]

### 5. Karamba3D「BESO for Beams」手冊摘錄

官方手冊 URL:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.9-beso-for-beams [LIT:Karamba3D v3 online manual]

**核心參數:**

| 參數 | 說明 |
|------|------|
| `TargetRatio` | 目標殘留質量 / 初始質量 |
| `MaxChangeIter` | 加/移除階段最大迭代數 |
| `MaxConvIter` | 收斂階段額外迭代數 |
| `WTension`, `WCompr`, `WShear`, `WMoment` | 截面力分量權重;加權後即為元素 sensitivity |
| `BESOFac` | 雙向因子:每迭代移除 (m+1)·n 根、再加回 m·n 根 |
| `MinDist` | 同一迭代中元素變更的最小空間距離(m),防空間聚集 |
| `WLimit` | 移除權重低於 WLimit × 平均權重 的元素(相對閾值) |
| `GroupIds` | 元素群組;整群同開關 |

**敏感度計算依據:**

「以截面力分量的變形能密度為基礎;各分量乘以對應權重後加總即為元素權重。」——為 Karamba3D 對桿系 BESO sensitivity 的官方定義。對應標準 BESO 的 α_i^e,但改為分力加權而非全應變能,使使用者可強調受拉/受壓/受彎等行為。

**Soft-kill 實作:**

BESO for Beams 以**降低截面剛度而非直接刪除**來實作「軟殺」,避免斷開連通性——Karamba3D 說明中明確指出此設計防止不連通結構出現。

**Karamba3D 與 ER/ARmax:**

BESO for Beams 不暴露 ER 與 ARmax 兩個連續域 BESO 的標準參數——改以 `BESOFac` 控制雙向加移除比例,`MaxChangeIter` 控制總迭代步數,兩者共同決定收斂速度。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.9-beso-for-beams]

### 6. Fail-Safe / Robustness TO 與「逐元素移除倒塌約束」

**Jansen et al. 2014(奠基論文):**

Jansen, M., Lombaert, G., Schevenels, M. et al., 「Topology optimization of fail-safe structures using a simplified local damage model」, *Structural and Multidisciplinary Optimization* 49(4):657–666, 2014. [LIT:https://link.springer.com/article/10.1007/s00158-013-1001-y]

核心方法:以「補丁(patch)內全部元素剛度歸零」模擬局部損傷;對所有可能補丁位置的損傷情境取 minimax(最差情境最小化 compliance);結果設計在任一位置被刪除後仍能傳力。

兩個已知缺點(文獻自述):(1) 需要幾乎等同於元素數量的 FEA 分析,計算極貴;(2) 損傷以連續補丁定義,而非離散構件——對桁架(離散構件)的適用性較差。

**Fail-safe truss TO(Stolpe, Svanberg 等系列):**

- Zhou & Rozvany 2001、Stolpe 2019「Fail-safe truss topology optimization」(*SMO* 58:1312–1334):直接對每根桿件定義損傷情境,以 LP/conic programming 求全局最優 fail-safe 桁架。[LIT:https://link.springer.com/article/10.1007/s00158-019-02295-7]
- Zhu et al. 2023「Adaptive topology optimization of fail-safe truss structures」(*SMO*):自適應加入主動損傷情境(damage-case adding),處理多達 16290 根桿件 + 16291 損傷情境至全局最優。[LIT:https://link.springer.com/article/10.1007/s00158-023-03585-x]

**與 FrameCore N2(漸進倒塌驅動器當 BESO 約束)的先行技術定位:**

FrameCore 的構想 N2 = 在 BESO 循環中,每步移除後跑 `runProgressiveCollapse` 評估倒塌裕度,以此作為 sensitivity 調整或可行性約束。文獻先行技術分析:

1. **概念最近的先例**:Fail-safe truss TO(Stolpe 2019, Zhu 2023)——直接構件移除損傷 × 靜力重分析,目標最小材料 + 滿足所有損傷情境下結構完整性,與 N2 幾乎等價,只差在 FrameCore 用 LSP 漸進驅動器(sequential linear analysis)而非 plastic LP。
2. **Jansen 2014**:連續域版本,方法論相似但適用於殼/板結構。
3. **N2 的差異化之處**:FrameCore 的 `runProgressiveCollapse` 是多步 sequential 移除(非單步單補丁),可評估**連鎖倒塌**而非單一損傷情境的靜力冗餘;這在文獻中更接近 Zhu 2023 的「progressive damage adding」框架,但後者是真彈塑性 LP,N2 是 sequential linear analysis——計算更快、但不是全局最優。[LIT:Stolpe 2019; Zhu 2023(同上)]
4. **誠實邊界**:N2 作為 BESO 嵌套約束,每迭代需跑一次完整的漸進倒塌驅動器,計算成本與 Jansen 2014 一樣昂貴(O(n_members) 次 factorize);FrameCore 有 factorize-once PreparedSystem 可部分緩解,但整體仍是 O(n_iter × n_members) 量級,需明確記錄。[THEORY]

## 對 FrameCore 的含義

**採什麼:**

1. **sensitivity number 公式直接可用**: `α_i = 0.5 * u_i^T * K_i * u_i` 是標準 BESO sensitivity,對 `Member` 元素可從 `recover` 後的局部位移 + 局部剛度矩陣計算;已有 `frame_cli.exe` 逐元素內力,可以內力重建應變能:
   `U_i = N²L/2EA + M_y²L/6EI_y + M_z²L/6EI_z + T²L/2GJ + (shear terms)`
   此即 Karamba3D「截面力分量加權」的解析版本,FrameCore 可直接用,亦可暴露各分力權重參數對應 W_tension/W_moment。

2. **歷史平均必須實作**: `ᾱ_i^k = (α_i^k + α_i^(k-1))/2` 是防止振盪最小代價,只需保存上一迭代 sensitivity 向量。

3. **收斂準則應採 2N 窗**: N=5, τ=0.001,同時要求「目標體積到達 **且** 收斂」,非先到即停。

4. **compliance-best 回退是必要實作**: 每迭代記錄 `(volume_ratio, compliance)`,輸出歷史最優拓撲而非最終拓撲——這是消除低體積率暴增最簡單有效的做法,有充分文獻支持。

5. **連通性約束利用現有 `analyzeConnectivity`**: FrameCore 的 union-find 連通分析完全可插入 BESO 迭代,移除後若出現懸空碎塊立即標記為不可選刪除對象,對應 Zuo & Xie 2016 先例。

6. **soft-kill 對 FrameCore 更安全**: 以 `member.active = false`(現有功能)實作 hard-kill 已可行,但若要避免機構化懸停導致 K 奇異,改為保留元素但縮比 EA/EI(例如 ×10^-6)是更穩健的 soft-kill;`PreparedSystem` 的 `pivotMargin` 可作為機構偵測守門。

7. **Michell 桁架做 oracle**: Sokół 2011 的 99 行代碼或 Dorn 1964 LP ground structure 對 2D 懸臂有精確數值,可作為 BESO 實作的對比基準(收斂到 Michell 解表示 sensitivity 正確)。

**避什麼:**

1. **不要以「到達目標體積」作為唯一停機條件**:文獻與實驗均顯示低體積率下 compliance 可暴增,需 compliance-best 回退機制。

2. **不要在無連通性檢查的情況下做 hard-kill**: 孤立桿、懸空節點會讓後續 `assembleAndFactor` 奇異(或機構→`singular=true`);FrameCore 的 `validateModel` 已能偵測懸空節點,但需在每步 BESO 移除後呼叫。

3. **不要宣稱全局最優**: BESO 是 heuristic;對桁架最優解的全局最優方法是 LP/SDP ground structure;FrameCore 的 BESO 實作應誠實標為「近似 heuristic,收斂至局部 Michell 近似」。

4. **N2(漸進倒塌約束 BESO)計算成本需明確標注**: 每迭代一次 `runProgressiveCollapse` = O(n_members) 次 `solveLoad`,對 100 根桿系每迭代約 100 次 sparse LDLᵀ 解;大型問題需快取或近似。

**誠實邊界:**

- BESO 的 sensitivity 基於線彈性假設:塑鉸/屈服後重分布不在範圍內。
- 連通性 + compliance-best 回退可消除大多數病態,但不能保證收斂至全局最優。
- Karamba3D 的 GroupIds 約束(整群同開關)在 FrameCore 中無現成對應;若要實作需額外維護 group-to-member 映射。
- N2 fail-safe 約束與 Jansen 2014 / Stolpe 2019 的本質差異:後者有數學最優性保證,N2 用 sequential linear analysis 近似倒塌,是工程近似。

## 來源清單

- Huang, X. & Xie, Y.M. (2007). Convergent and mesh-independent solutions for the bi-directional evolutionary structural optimization method. *Finite Elements in Analysis and Design*, 43(14), 1039–1049. https://www.researchgate.net/publication/223369659
- Huang, X. & Xie, Y.M. (2010). *Evolutionary Topology Optimization of Continuum Structures: Methods and Applications*. Wiley. Ch.3 BESO Method. https://onlinelibrary.wiley.com/doi/10.1002/9780470689486.ch3
- Karamba3D v3 Manual — BESO for Beams (§3.6.9). https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.9-beso-for-beams
- Karamba3D v3 Manual — BESO for Shells (§3.6.10). https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.10-beso-for-shells
- Zuo, Z.H. & Xie, Y.M. (2016). A Bi-directional Evolutionary Structural Optimisation algorithm with an added connectivity constraint. *Finite Elements in Analysis and Design*. https://www.sciencedirect.com/science/article/abs/pii/S0168874X16305789
- Ghabraie, K. (2015). An improved soft-kill BESO algorithm for optimal distribution of single or multiple material phases. *Structural and Multidisciplinary Optimization*, 52(4), 1215–1233. https://link.springer.com/article/10.1007/s00158-015-1268-2
- Jansen, M., Lombaert, G., Schevenels, M. et al. (2014). Topology optimization of fail-safe structures using a simplified local damage model. *Structural and Multidisciplinary Optimization*, 49(4), 657–666. https://link.springer.com/article/10.1007/s00158-013-1001-y
- Stolpe, M. (2019). Fail-safe truss topology optimization. *Structural and Multidisciplinary Optimization*, 58, 1312–1334. https://link.springer.com/article/10.1007/s00158-019-02295-7
- Zhu, J. et al. (2023). Adaptive topology optimization of fail-safe truss structures. *Structural and Multidisciplinary Optimization*. https://link.springer.com/article/10.1007/s00158-023-03585-x
- Topology Optimization Based on SA-BESO. *Applied Sciences* 13(7):4566, 2023. https://www.mdpi.com/2076-3417/13/7/4566
- Topology Optimization with Dynamic Evolution Rate (DER-BESO). *SMO* 2020. https://link.springer.com/article/10.1007/s00158-020-02588-2
- Sokół, T. (2011). A 99 line code for discretized Michell truss optimization written in Mathematica. *SMO*. https://link.springer.com/article/10.1007/s00158-010-0557-z
- Dorn, W.S., Gomory, R.E. & Greenberg, H.J. (1964). Automatic design of optimal structures. *Journal de Mécanique*, 3, 25–52. [書目,無公開 URL]
- Ohsaki, M. (2011). Topology and Geometry Optimization of Trusses and Frames. Kyoto U. http://www.se-lab.archi.kyoto-u.ac.jp/ohsaki/pdf/b0201.pdf
- BESO with connectivity constraint (Zuo & Xie 2016 abstract): https://www.sciencedirect.com/science/article/abs/pii/S0168874X16305789
- Topology Optimization Applications on Engineering Structures. IntechOpen. https://www.intechopen.com/chapters/70489