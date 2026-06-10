> 出處:研究輪 fan-out agent 查證報告(2026-06-10,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

現在我已收集到足夠的資料,可以撰寫完整的研究報告。

# WS-D Tension-Only 桿件迭代求解

## 摘要

Tension-only(僅受拉)桿件使結構勁度矩陣隨元素狀態變動,本質上是幾何非線性問題,必須以非線性分析處理。業界標準做法是 **active-set 不動點迭代**:逐步識別受壓桿件,將其停用(zero 或 negligible stiffness),重解,直到 active set 不再變動。這個迭代沒有一般性收斂保證——在對稱荷載下特定 X 型斜撐等構型會出現兩桿交替進出 active set 的振盪。從數學角度看,tension-only 問題等價於一個 **Linear Complementarity Problem (LCP)**:矩陣為 P-matrix(或對稱正定)時 LCP 有唯一解,對應的 active-set 迭代因此收斂;否則唯一解不保證存在。業界軟體(Karamba3D、RISA-3D、STAAD.Pro、Tekla Structural Designer)全部採迭代方案,並設 MaxIter 上限作保險,收斂判據一致為「兩次迭代間 active set 不再改變」。Karamba3D 採 soft-kill(零剛度不刪元素)且**不恢復**移除元素;RISA-3D 允許重新開啟;STAAD.Pro 預設 10 次迭代可用 `SET ITERLIM` 調整;Tekla 以 relaxation factor 處理振盪。Oracle 構想:X 型斜撐側向荷載時解析解乾淨,但若荷載恰好使雙對角線應力相等,小擾動便引發振盪。

---

## 發現

### 1. 標準做法 — active-set 迭代停用受壓桿

- 業界所有主流軟體(RISA-3D、STAAD.Pro、Tekla/SCIA、Karamba3D、ETABS/SAP2000)均採**重複線性解(iterative linear re-analysis)**:每輪識別受壓的 tension-only 桿→設為 inactive(零或極小勁度)→重解→直到 active set 穩定。[LIT:https://blog.risa.com/post/implementing-realistic-behavior-for-t-c-members-in-risa-3d]
- RISA-3D 明確描述:"the program removes the member from the model and updates the stiffness matrix. It then resolves the model and checks the updated status of T/C only members again." 判斷收斂 = "the status of all T/C only members stay unchanged during the iterations"。[LIT:https://blog.risa.com/post/implementing-realistic-behavior-for-t-c-members-in-risa-3d]
- RISA-3D 非收斂處理:顯示警告,並**強制所有問題桿件接受雙向力**後繼續計算(不是真解)。[LIT:https://blog.risa.com/post/implementing-realistic-behavior-for-t-c-members-in-risa-3d]
- STAAD.Pro 預設最多 **10 輪**,以 `SET ITERLIM i` 調高上限。文檔直言 "This method does not always converge and may become unstable."。[LIT:https://docs.bentley.com/LiveContent/web/STAAD.Pro%20Help-v18/en/STD_MEMBER_TENSION_COMPRESSION.html]
- Chan & Chui《Non-Linear Static and Cyclic Analysis of Steel Frames with Semi-Rigid Connections》(Elsevier, 2000)覆蓋非線性求解技術(Chapter 4)及塑鉸/半剛接非線性靜力分析(Chapter 6),但書籍本身主要聚焦半剛接框架的 Newton-Raphson,針對 tension-only active-set 的專章未能取得原文確認。[UNKNOWN — 書目確認:https://books.google.com/books/about/Non_linear_Static_and_Cyclic_Analysis_of.html?id=q3inlAEACAAJ]

### 2. 振盪/發散的處理

- **振盪根本原因**:當移除一根受壓桿後重解,原本受拉的另一根桿變受壓,下一輪又移除它,如此交替進出 active set 無法收斂。Tekla 文檔明確描述此現象 "when the model shuts the member off and re-runs the analysis, the would-be strain in the member is a tensile one and the program will turn the member back on again"。[LIT:https://support.tekla.com/article/why-will-the-analysis-not-complete-when-my-model-has-tension-only-bracing]
- **Relaxation factor(under-relaxation)**:Tekla Structural Designer 2017 起引入 relaxation factor — 允許 tension-only 桿件在迭代中承受「名義壓力」而不立即停用;收斂後這些桿件恢復零或拉力狀態。這是最常見的振盪緩解手段。[LIT:https://support.tekla.com/article/why-will-the-analysis-not-complete-when-my-model-has-tension-only-bracing]
- **單調移除(不恢復)**:Karamba3D 採 soft-kill(極小剛度 negligible stiffness),一旦元素被「消除」**不會恢復**;這等同強制單調縮小 active set,避免振盪但可能錯失真解(若某些桿件需先被壓再被放開才能收斂)。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.12-tension-compression-eliminator]
- **狀態哈希終止**:工業軟體文檔未見明確的「狀態哈希比對」術語,但 MaxIter 上限扮演類似角色——檢測到兩輪 active set 相同即停止,若無限循環則靠 MaxIter 截斷。[UNKNOWN — 狀態哈希作為正式策略未見於查到的文獻]

### 3. LCP 視角與收斂性質

- **Tension-only 問題可嚴格化為 LCP**:對於每根 tension-only 桿件,軸力 N ≥ 0,延伸量 δ ≥ 0,互補條件 N · δ = 0(要嘛受力要嘛縮短無反力)。以 K(active set)u = f 的框架表達,追蹤哪些桿件 active 即是求解 LCP(q, M),其中 M = reduced stiffness submatrix。[LIT:https://www.emerald.com/ec/article-abstract/42/9/3656/1298665/Tensegrity-structure-statics-with-slack-cables?redirectedFrom=fulltext][LIT:https://en.wikipedia.org/wiki/Linear_complementarity_problem]
- **唯一解條件**:LCP(q, M) 對所有 q 均有唯一解的充要條件是 M 為 **P-matrix**(所有主子式 > 0)。對稱正定矩陣是 P-matrix 的充分子集。對於良態支承的線彈性桁架/框架,完整 K 是對稱正定,但移除支承條件後的 reduced K 未必是 P-matrix(若存在機構則奇異)。[LIT:https://en.wikipedia.org/wiki/Linear_complementarity_problem][LIT:https://arxiv.org/abs/2508.05777]
- **Tin-Loi 等人**的互補性結構力學文獻(Archives of Computational Methods in Engineering, 2015)系統性整理了結構工程中非線性非光滑問題(含單側約束/tension-only)的 LCP/MPEC 框架,確認 LCP 視角的學術基礎。[LIT:https://link.springer.com/article/10.1007/s11831-015-9158-8]
- **收斂無一般性保證**:一般 tension-only active-set 迭代是 LCP 的 Gauss–Seidel 形式不動點迭代,**無一般性多項式時間收斂保證**,但若 M 是 P-matrix(線彈性框架充分支承情況下通常成立),則有限步終止。若結構因元素停用變為機構(M 奇異),LCP 可能無解或多解。[THEORY][LIT:https://en.wikipedia.org/wiki/Linear_complementarity_problem]

### 4. Karamba3D Tension/Compression Eliminator 原文

- 元件路徑:`3.5.12 Tension/Compression Eliminator`(原文 URL 已遷移至此)。
- 演算法:迭代,停止條件 = "no changes occur from one step to another" 或達 `MaxIter`。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.12-tension-compression-eliminator]
- Soft-kill:**保留元素結構但賦予 negligible stiffness**,而非從模型刪除。**不恢復**已移除元素。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.12-tension-compression-eliminator]
- `Compr` 參數:預設 `false` = 保留受拉元素(消除受壓);`true` = 反向。可由 `BeamInds` 指定篩選範圍。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.12-tension-compression-eliminator]
- 振盪處理:**文檔未提**。[UNKNOWN]

### 5. SAP2000/ETABS 業界基準

- 必須設定為 **Nonlinear Static 或 Nonlinear Time-History** 分析才能啟用 tension/compression limit 功能;線性分析完全不識別。[LIT:https://wiki.csiamerica.com/display/etabs/Tension-only+elements+in+ETABS]
- 設定方式:Assign > Frame/Line > Tension/Compression Limits,compression limit = 0(或負值)實現 tension-only。[LIT:https://wiki.csiamerica.com/display/etabs/Tension-only+elements+in+ETABS]
- 每個荷載工況單獨非線性求解 — "Each load case will conceivably have a different set of active members";線性疊加對 tension-only 完全無效,設計組合須逐一跑非線性。[LIT:https://resource.midasuser.com/en/blog/structure/analysis-and-design-method-of-structures-using-tension]
- CSI Analysis Reference Manual(2011-12)為官方最完整文檔,覆蓋 SAP2000/ETABS/SAFE/CSiBridge 的非線性求解控制。[LIT:https://docs.csiamerica.com/manuals/misc/CSI%20Analysis%20Reference%20Manual%202011-12.pdf]

### 6. Oracle 構想驗證 — X 型斜撐解析解

- **解析設定**:方形框架寬 L、高 L,兩根對角桿 EA 相同,水平側力 H 作用於頂部。假設兩桿為純軸力桿(鉸接端)、框架柱提供垂直約束、梁極剛。
- **解析解**:對角桿的方向角 θ = 45°,長度 = L√2。水平力由兩桿分擔;若兩桿均 active,由對稱性各承受水平分力 H/2 cosθ = H/(2cos45°) = H√2/2 的軸力,但實際一根受拉一根受壓。**正確解**:只有受拉桿 active,拉力 T = H/(cos45°) = H√2,壓桿停用後框架依然靜定(假設柱能提供垂直平衡)。[THEORY]
- **振盪構型**:若荷載為 **雙向對稱荷載(如豎向均佈)**,使兩根斜撐軸力相等且同號,active set 第一輪正確;但若荷載恰好在 **中性點**(如零側力+某個角度荷載使兩桿應力等量異號),任何微小數值誤差可使兩桿各被輪流停用 — 此為已知最容易觸發振盪的構型。更精確的已知振盪案例:X 型斜撐 + 側力大小介於兩桿「均受拉」與「一拉一壓」的轉換臨界值附近,由於非線性解不唯一(兩桿或一桿均是平衡態),迭代在兩態間震盪。[LIT:https://support.tekla.com/article/why-will-the-analysis-not-complete-when-my-model-has-tension-only-bracing]

---

## 對 FrameCore 的含義

1. **必須跑非線性迭代,不能線性一刀切**:FrameCore 現行 `assembleAndFactor` 是純線性;實作 tension-only 必須在外層加 active-set loop — 每輪設 `Member.active = false` for 受壓桿件,重新 `assembleAndFactor`(或增量更新),直到 active set 收斂。這與現有 `runProgressiveCollapse` 的 progressive removal 邏輯架構相近,可以共用。

2. **每個荷載工況獨立求解**:不能用線性疊加。`solve(combo)` 無法複用 single `PreparedSystem`;每個非線性工況必須從組合出發獨立收斂,或逐工況分別迭代。這衝擊 `factorize-once` 架構:可以 factorize-once per active-set-state 而不是 factorize-once for all combos。

3. **振盪處理要選策略**:
   - **Soft-kill + 不恢復**(Karamba3D 做法):最簡單實作,但若真解需要「桿件先被壓後恢復」則找不到。
   - **允許恢復 + relaxation factor**(Tekla 做法):更接近真解,但需調整 relaxation 參數,實作複雜。
   - **狀態哈希終止**:儲存每輪 active set 的 bitmask,若出現重複 state 立即截斷報警(cycle detected),是最便宜的振盪偵測。建議 FrameCore 採此作 guard。

4. **Oracle 設計**:X 型斜撐 + 純側向荷載解析解乾淨(拉桿 T = H√2,壓桿 T = 0),可直接做 F34 oracle。需補:荷載組合測試(D + L 對稱 → 兩桿均拉;W 側向 → 一拉一壓),以及振盪邊界構型(可主動觸發振盪後驗證偵測正確報告 CycleDetected 而非靜默 stale result)。

5. **LCP 正確性界限(誠實標)**:對充分支承的線彈性框架,減縮 K 在 active set 非空時是對稱正定(P-matrix),迭代有有限步收斂保證。若停用元素造成機構,LCP 可能無解或多解 — 此時應讓 `analyzeConnectivity`/`worstUtilization` 先行偵測機構後報告 `singular`,而非讓 active-set loop 無限跑。

6. **設計組合 N 個工況需 N 次迭代**:這是主要效能代價。若最終要與 Karamba3D 比較速度,benchmark 要計入「工況數 × 每工況迭代次數 × 每次重因子化」,而非只比單次求解。可以考慮 **incremental update**(移除一根桿件做 rank-1 downdate on LDLᵀ)降低重因子化代價,與 WS-B 的增量重分析可以共用技術路線。

---

## 來源清單

- [RISA-3D T/C Member Implementation](https://blog.risa.com/post/implementing-realistic-behavior-for-t-c-members-in-risa-3d)
- [Karamba3D Manual — Tension/Compression Eliminator](https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.12-tension-compression-eliminator)
- [STAAD.Pro Member Tension/Compression Specification (Bentley)](https://docs.bentley.com/LiveContent/web/STAAD.Pro%20Help-v18/en/STD_MEMBER_TENSION_COMPRESSION.html)
- [ETABS Tension-only elements — CSI Knowledge Base](https://wiki.csiamerica.com/display/etabs/Tension-only+elements+in+ETABS)
- [CSI Knowledge Base — Tension and Compression Limits Tips](https://wiki.csiamerica.com/display/kb/Tension+and+compression+limits+tips)
- [CSI Analysis Reference Manual 2011-12 (PDF)](https://docs.csiamerica.com/manuals/misc/CSI%20Analysis%20Reference%20Manual%202011-12.pdf)
- [Midas — Analysis and Design Method for Tension-Only Elements](https://resource.midasuser.com/en/blog/structure/analysis-and-design-method-of-structures-using-tension)
- [Tekla — Why will analysis not complete with tension-only bracing?](https://support.tekla.com/article/why-will-the-analysis-not-complete-when-my-model-has-tension-only-bracing)
- [SCIA Engineer — Model a tension-only element](https://www.scia.net/en/support/faq/scia-engineer/modelling/model-tension-only-element)
- [Wikipedia — Linear Complementarity Problem](https://en.wikipedia.org/wiki/Linear_complementarity_problem)
- [Emerald — Tensegrity Structure Statics with Slack Cables (LCP formulation)](https://www.emerald.com/ec/article-abstract/42/9/3656/1298665/Tensegrity-structure-statics-with-slack-cables?redirectedFrom=fulltext)
- [arXiv:2508.05777 — Existence and Uniqueness of LCP Solution in Contact Mechanics](https://arxiv.org/abs/2508.05777)
- [Springer — Complementarity Problems in Structural Engineering: An Overview (Tin-Loi)](https://link.springer.com/article/10.1007/s11831-015-9158-8)
- [Chan & Chui — Non-Linear Static and Cyclic Analysis of Steel Frames with Semi-Rigid Connections (Elsevier, 2000)](https://books.google.com/books/about/Non_linear_Static_and_Cyclic_Analysis_of.html?id=q3inlAEACAAJ)
- [RISA-3D Tension-Only Bracing Guide](https://blog.risa.com/post/how-to-model-tension-only-bracing-in-risa-3d)
- [Structural Academy — Tension and Compression Limits in Frame Elements](https://structuralacademy.com/article/en/Limites-De-Tracao)
