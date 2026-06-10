三份文件都已讀完。現在系統性交叉比對各個審查面向。

---

# 審查報告:novelty 措辭一致性與對標矩陣事實性

## 裁決

**FAIL** — 發現 3 個 MAJOR 級問題,需修正才能定稿。

---

## 發現(按嚴重度排序)

### [MAJOR] KARAMBA3D_ROADMAP.md §4:Karamba BESO 殼欄位措辭與 WS_A1 不一致

**位置**: 對標矩陣第「拓撲優化」列,Karamba3D 欄標 `✅ BESO 梁/殼`。

**問題**: WS_A1 §3.3 明確指出 Karamba BESO for Shells 採用 **soft-kill 減厚**而非直接移除元素(「Soft-kill 減厚(KillThick 預設 0.00001m)而非直接移除元素」),在概念上與 BESO(Bi-directional Evolutionary Structural Optimization 的二元 active/inactive 策略)存在本質差異。Roadmap §4 把它與「梁 BESO」並列為 `✅`,隱含兩者實作對等;但 §5 N2 卻明確寫「**我們=離散框架×序列線性倒塌×Stable/Collapsed 終點**」並指出 Karamba 殼 BESO 是「soft-kill 減厚+空間濾波」。矩陣欄位的 `✅` 沒有揭露這個重要的實作差異,讀者可能誤解 Karamba 殼 BESO 與梁 BESO 是同一類硬移除策略。

**建議修法**: 將該格改為 `🔶 BESO 梁(硬移除)/殼(soft-kill 減厚)`,並在腳注或備註欄補充「殼端以 soft-kill 減厚代替元素移除,見 WS_A1 §3.3」。

---

### [MAJOR] KARAMBA3D_ROADMAP.md §5 N2:roadmap 對 fail-safe TO 的措辭強度弱於 WS_N 明確警告

**位置**: §5 N2 欄「先行技術(誠實定位)」:「**勿自稱 fail-safe**」。

**問題**: WS_N §N2 結論對此給出了明確的理由和建議宣稱句(「以 GSA 替換荷載路徑(LSP 漸進倒塌模擬)為評估器的框架結構 BESO,在概念上類似 fail-safe 拓撲優化(Jansen 2014),但針對離散框架結構並使用序列線性倒塌分析而非靜態最壞情況柔度作為評估函式」)。但 §7 主線表格 S7 欄位寫「**S7 BESO + N2 韌性約束版**」,以及 §1 執行摘要第 2 點寫「N2 倒塌韌性 BESO」——兩處都未附任何 fail-safe TO 警告標記或與 WS_N 對齊的誠實定位標籤。

特別是 §7 的 S7 說明「**(+N2 韌性版,Karamba 無)**」,這是一個隱含優越性主張。WS_N 指出 Jansen 2014 fail-safe TO 已涵蓋「任意損傷後仍可承載」的概念,FrameCore N2 的差異在「離散框架 + LSP 評估器」,這個差異沒有在 S7 描述中出現,讀者看到「Karamba 無」會以為這是完全無先例的新算法。

**建議修法**: §7 S7 欄補註「(概念親緣 Jansen 2014 fail-safe TO,差異:離散框架+LSP 評估器,見 WS_N §N2)」;§1 執行摘要的「N2 倒塌韌性 BESO」後補一個誠實標籤如 `[LIT:Jansen2014,整合 novelty]`。

---

### [MAJOR] KARAMBA3D_ROADMAP.md §4:免費版殼元素上限敘述與 WS_A1 不精確對齊

**位置**: §4 對標矩陣備注(§2 摘要)以及 §9:「Karamba 免費版 ≤20 梁元素」,殼上限未提。

**問題**: WS_A1 §5 對殼上限的表述是:「殼元素隱含上限約 50(從『激活授權後可開啟超過 50 殼元素範例檔』**反推**)」。此推論依據的是反推(非官方直接聲明);WS_A1 並非明確查證到「殼 ≤50」,而是用反推並加上「約」字。Roadmap §2 把殼上限寫為確定性的「殼約 ≤50(WS_A2)」,引用標記指向 WS_A2 而非 WS_A1。但本輪只能讀到 WS_A1,WS_A1 §5 明確表示此為反推。若 WS_A2 有更強的直接依據,引用來源可能正確;但若 WS_A2 的結論也來自同一反推推論,則 roadmap 用確定性語氣記載了一個不確定的數字。

**建議修法**: 檢查 WS_A2 的殼上限根據;若仍是反推,將 §2 改為「殼元素試用上限**約** ≤50(反推,WS_A1/A2)」並補 `[PENDING:待官方確認]` 標記。

---

### [MINOR] KARAMBA3D_ROADMAP.md §4:影響線欄 Karamba 標記邏輯欠精確

**位置**: 對標矩陣「影響線/沉陷」列,Karamba3D 欄:「⛔/🔶(Gap Load 可做影響線)」。

**問題**: WS_A1 §4 梁元素載重確實提到 Gap Load 用於影響線,但 WS_A1 對此的描述是「Gap Load,用於影響線」的非正式說明,並非查證到 Karamba3D 有完整的 Müller-Breslau 影響線功能。Roadmap 的 `🔶` 標記可能過度肯定了 Karamba 對影響線的支援程度;而 FrameCore 的影響線是有 Müller-Breslau 互驗的完整實作。符號選擇上 `⛔/🔶` 格式不統一(矩陣其他格子只用單一符號),且這一行 FrameCore 欄位只寫 `✅` 沒有補說明,相對應的 Karamba 欄值被部分「肯定」,削弱了 FrameCore 的優勢對比清晰度。

**建議修法**: 改為 `⛔(無原生元件;Gap Load 可間接近似)`,或補一行腳注說明差異;FrameCore 欄補 `✅(Müller-Breslau 互驗)`。

---

### [MINOR] KARAMBA3D_ROADMAP.md §4:Karamba 二階欄標記未區分 P-Δ 與真大位移

**位置**: 對標矩陣「二階 P-Δ」列,Karamba3D 欄 `✅ ThII`。

**問題**: WS_A1 §1.2 明確指出 AnalyzeThII「**屬 P-Δ 級別而非真大變形**」;§1.3 Analyze Nonlinear WIP 才是真大位移,且「自標收斂限制」。Roadmap §4 把「二階 P-Δ」列 Karamba 標為 `✅ ThII`,下一列「幾何非線性(大位移)」標為 `🔶 WIP`——這個拆法是合理的,但 WS_A1 §1.2 另有一個細節:AnalyzeThII 有 `NoTenNII`/`NoComNII` 選項可單獨取壓力或拉力 N^II,用於薄膜防皺。這個細節在矩陣中未反映,對 Tension-only 欄位可能有連動(Karamba 的 Tension-only 能力與 ThII 有整合,不只是 Eliminator)。此問題較輕微,主要影響完整性。

**建議修法**: 在 Karamba ThII 格子後加備注或腳注說明 `NoTenNII/NoComNII` 選項與 Tension-only 欄的潛在重疊。

---

## 沒問題的部分

§5 先行技術定位欄(N1/N2/N3/N4)對 WS_N 結論的照抄基本忠實:N1 明確標「非數學新算法/工程整合」、N2 明確標「勿自稱 fail-safe」、N3 明確標「非新算法」、N4 標 `[UNKNOWN 級定位]`——每處均與 WS_N 的結論措辭一致,強度未超越;Tension/Compression Eliminator 的 `✅ 迭代 Eliminator` 標記、殼 TRIC 忽略橫剪的 `🔶` 標記、反應譜/時程 `⛔` 標記,與 WS_A1 結論均可逐一對應,未見矛盾;免費版梁元素「≤20」與 WS_A1 §5 直接對齊。