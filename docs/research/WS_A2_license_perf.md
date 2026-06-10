> 出處:研究輪 fan-out agent 查證報告(2026-06-10,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

충분한 정보를 수집했습니다. 보고서를 작성하겠습니다.

# WS-A2 Karamba3D 授權/商業模式/效能宣稱/已知限制

## 摘要

Karamba3D 是 Clemens Preisinger(與 Bollinger+Grohmann 合作)開發的參數化結構分析外掛,專為 Grasshopper/Rhino 環境設計。**閉源商業軟體**,核心求解器為 C++ native DLL(`karamba.dll`),外圍有 .NET 包裝層(`karambaCommon.dll`)供 C# scripting API。目前最新版 **3.1.60519(2025-05-20)**。授權分 PRO(€900/年)、EDU(€30/2年)、免費版(功能受限)。效能上對小/中型模型(初期設計階段)極快,但超過約 10,000 元素後顯著變慢,求解器型態推測為直接法(dense or sparse),具體 backend 廠商**未公開披露**。材料非線性正式**不支援**,幾何非線性為 WIP 狀態。

---

## 發現

### 1. 授權與價格

- **PRO License**：€900/年(Cloud 或 Network,1年訂閱制,無永久版);月租亦可(Cloud 1 Month,價格未列出)。[LIT:https://buy.karamba3d.com/products/pro-license]
- **EDU License**：€30/2年,2年後可免費續訂(需重新驗證學生身份)。含 PRO 全功能,僅非商業用途。[LIT:https://buy.karamba3d.com/products/edu-license]
- **免費版(Free)**：永久免費,非商業用途;beam/shell 元素數量**無限制**,但「其他功能」受限(手冊標示 limited,未詳列具體項目)。[LIT:https://manual-1-3.karamba3d.com/1-introduction/1.2-licenses]
- **試用版(舊 Full version 模式)**:限 ≤20 beam 元素、≤50 shell 元素。[LIT:https://manual-1-3.karamba3d.com/1-introduction/1.2-licenses]
- **大學 LAB License**:10 人×網路授權,無限時,非商業;教師可申請學期免費授權。[LIT:https://karamba3d.com/get-started/university/]
- **閉源確認**:核心求解器與所有 DLL 均**未開源**。GitHub 上 [karamba3d/K3D_tests](https://github.com/karamba3d/K3D_tests) 只是 C# 單元測試套件(呼叫 karambaCommon.dll NuGet 包),不含任何引擎源碼。[LIT:https://github.com/karamba3d/K3D_tests]

### 2. 技術棧與求解器 backend

- **三層架構**(由 Clemens Preisinger 本人在論壇確認):
  - `karamba.dll`:**C++ native DLL**,執行所有數值計算;「not a .NET assembly」,不可直接在 C# 專案引用。[LIT:https://discourse.mcneel.com/t/installing-karamba-dll-for-visual-studio/143546]
  - `karambaCommon.dll`:**.NET 包裝層**,提供幾何型別(vector/point/mesh)及 C# scripting API;獨立於 Rhino/Grasshopper。[LIT:https://scripting.karamba3d.com/1.-introduction/1.1-scripting-with-karamba3d]
  - `karamba.gha`:Grasshopper 外掛,將 GH 元件連接到 karambaCommon.dll。[LIT:https://scripting.karamba3d.com/1.-introduction/1.1-scripting-with-karamba3d]
- **C++ 求解器 backend 廠商/函式庫**:官方文檔及論壇從未公開披露是 Eigen/CHOLMOD/PARDISO/MKL 或自寫;搜遍手冊、scripting guide、GitHub 無任何字樣。**[UNKNOWN]**
- **稀疏或稠密**:官方文檔稱「計算時間隨 DOF 數以 n³ 成長(全連接結構)、帶狀結構以 ~0.5·n·n²_neigh 成長」,與直接法相符;30,000 節點以上用戶反映計算時間「急遽增加」,推測為直接稀疏或帶寬求解器。具體 sparse solver 型態 **[UNKNOWN]**。[LIT:https://manual.karamba3d.com/appendix/a.4-background-information/a.4.4-hints-on-reducing-computation-time]
- **殼元素型別**:Nature Architects 評測指出 Karamba3D 殼為「first-order elements」,超過一階精度的效果需用高階元素。手冊未明確命名 MITC4/DKT 等。**[UNKNOWN](官方從未公告殼元素形式名稱)**。[LIT:https://nature-architects.com/en/blog/1194/]
- **學術文獻**:Preisinger & Heimrath(2014)「Karamba—A Toolkit for Parametric Structural Design」,*Structural Engineering International* vol.24 pp.217-221 — 該文為工具介紹性論文,無 solver 實作細節。[LIT:https://www.semanticscholar.org/paper/Karamba%E2%80%94A-Toolkit-for-Parametric-Structural-Design-Preisinger-Heimrath/d6c6c120dd604c74ece02f41d2619e755534764b]

### 3. C# Scripting API

- **可用性**:karambaCommon.dll 提供完整 .NET API,可在 Grasshopper C# script component 內直接操作 Karamba3D 物件(Model、Element、Load 等)。[LIT:https://scripting.karamba3d.com]
- **官方手冊**:Scripting Guide 3.1.4 在 https://scripting.karamba3d.com,涵蓋迴圈、函式、物件導向程式概念、單元測試;範例同時提供 C#、IronPython、CPython。[LIT:https://scripting.karamba3d.com/1.-introduction/1.1-scripting-with-karamba3d]
- **GitHub 測試範例**:https://github.com/karamba3d/K3D_tests,C# 100%,示範如何以 API 建模/求解/讀結果。[LIT:https://github.com/karamba3d/K3D_tests]
- **限制**:karambaCommon.dll 不公開原始碼,API 文檔為 HTML reference(https://karamba3d.com/help/3-1-4),無法自訂求解器或替換演算法。[LIT:https://karamba3d.com/devs/]

### 4. 效能:官方說法與社群反饋

- **官方說法**:手冊稱計算時間取決於 DOF 數,全連接系統 O(n³),建議用桁架代替梁(「減半以上計算時間」)、多核心利用率、限制工況數、初期殼分析用粗網格。[LIT:https://manual.karamba3d.com/appendix/a.4-background-information/a.4.4-hints-on-reducing-computation-time]
- **官方宣稱模型規模上限**:「可有效處理約 10,000 元素的模型(取決於 CPU 與記憶體)」。[LIT:https://www.grasshopper3d.com/group/karamba3d/forum/topics/solving-speed]
- **社群反饋(McNeel Forum/Grasshopper)**:超過 30,000 節點計算時間「急遽增加」;多工況(Grasshopper 路徑設置不當可導致重複建多個模型)是常見陷阱。[LIT:https://www.grasshopper3d.com/group/karamba3d/forum/topics/solving-speed]
- **第三方評測(Nature Architects)**:線性靜力分析比 ANSYS 快 80–99%(0.004–0.014 秒 vs ANSYS 0.1–0.5 秒),適合初期設計參數研究;殼元素有**過剛現象**(扭轉剛度比實體元素高 58%),精度低於二階元素。[LIT:https://nature-architects.com/en/blog/1194/]

### 5. 已知技術限制

- **材料非線性**:正式**不支援**。手冊明確:「Karamba cannot calculate physically nonlinear phenomena like cracked concrete.」變通做法為 GH 迴圈迭代截面,非引擎內建。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.3-analyze-nonlinear-wip]
- **幾何非線性**:Analyze Nonlinear WIP 元件提供 Dynamic Relaxation / Newton-Raphson / Arc-Length 三種演算法,但標注 **WIP(work-in-progress)**;殼分析常有收斂問題,無一致切線剛度矩陣的二次收斂。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.3-analyze-nonlinear-wip]
- **二階理論(ThII)**:AnalyzeThII 提供幾何剛化效果,但 CombiNII=true 模式結果偏保守,收斂不保證,預設最多 50 次迭代。[LIT:https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.2-analyzethii]
- **殼元素精度**:第一階元素,曲面精度隨網格收斂但較慢;扭轉過剛。**[LIT:https://nature-architects.com/en/blog/1194/]**
- **平台限制**:PRO 授權僅限 PC(Windows);Cloud/Network Zoo 授權機制。[LIT:https://buy.karamba3d.com/products/pro-license]

### 6. 版本沿革

- 最新穩定版:**3.1.60519**,釋出日期 2025-05-20。[LIT:https://github.com/karamba3d/K3D_NightlyBuilds/releases]
- 3.1 主要新功能:ETABS(E2K)及 SAP2000(S2K)匯出支援、日本斷面資料庫、模糊搜尋截面名稱、Support Agent 元件、Curve-Curve Intersection 元件。[LIT:https://manual.karamba3d.com/new-in-karamba3d-3.1]
- 有 Nightly Build 管道(https://github.com/karamba3d/K3D_NightlyBuilds)持續更新。

---

## 對 FrameCore 的含義

**採什麼**
- FrameCore 技術棧(C++17/Eigen SimplicialLDLT)在求解器層面**優於或持平** Karamba3D:Eigen 稀疏求解器開源可驗證,Karamba 的 backend 完全黑盒。這是**「乾淨度」對標的核心優勢**。
- FrameCore 的 factorize-once PreparedSystem 架構天然支援多工況重用;Karamba3D 文檔提示大量工況是效能陷阱,FrameCore 在這一點架構上領先。
- karambaCommon.dll 的 C# API 設計可作為「FrameCore 如何暴露給 UE5 藍圖/C#」的設計參考:三層(native C++ solver / .NET wrapper / 應用層)是業界慣用模式。

**避什麼**
- 避免宣稱「任意規模都能跑」:Karamba3D 10,000 元素就是業界對輕量參數化引擎的公認門檻,FrameCore 若要超越需要有基線數據支撐。
- Karamba3D 殼元素扭轉過剛(+58%)是已知缺陷;FrameCore MITC4 同屬 flat-shell facet,曲面慢收斂源自膜鎖定,應在文檔裡誠實標,不要與 Karamba 的「first-order」行銷語言競爭,而要靠 Scordelis-Lo/pinched cylinder oracle 數據說話。
- 避免宣稱「材料非線性」—Karamba3D 也不宣稱,FrameCore 的 sequential linear + 塑鉸方案與其同屬「線性 + event-driven 近似」範疇,誠實標即可。

**誠實邊界**
- Karamba3D 的確切求解器 backend、稀疏圖型式、殼元素名稱均未公開,對標時無法做逐行比較;可比的只有功能矩陣與速度數據。
- 免費版「其他功能受限」的具體項目官方未列清單,若要對標免費版需自行測試。

---

## 來源清單

- [Karamba3D PRO License 購買頁](https://buy.karamba3d.com/products/pro-license)
- [Karamba3D EDU License 購買頁](https://buy.karamba3d.com/products/edu-license)
- [Karamba3D 授權手冊 v3](https://manual.karamba3d.com/1-introduction/1.2-licenses)
- [Karamba3D 授權手冊 v1.3.3(含免費版功能對照)](https://manual-1-3.karamba3d.com/1-introduction/1.2-licenses)
- [Karamba3D 大學授權頁](https://karamba3d.com/get-started/university/)
- [Scripting Guide 3.1.4 — 架構說明](https://scripting.karamba3d.com/1.-introduction/1.1-scripting-with-karamba3d)
- [Scripting Guide 首頁](https://scripting.karamba3d.com)
- [Karamba3D GitHub — K3D_tests(C# 測試套件)](https://github.com/karamba3d/K3D_tests)
- [Karamba3D GitHub — K3D_NightlyBuilds Releases(最新版本)](https://github.com/karamba3d/K3D_NightlyBuilds/releases)
- [Karamba3D Devs 頁(API 文檔連結)](https://karamba3d.com/devs/)
- [Karamba3D 效能優化手冊](https://manual.karamba3d.com/appendix/a.4-background-information/a.4.4-hints-on-reducing-computation-time)
- [Karamba3D AnalyzeThII 手冊](https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.2-analyzethii)
- [Karamba3D Analyze Nonlinear WIP 手冊](https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.3-analyze-nonlinear-wip)
- [Karamba3D What's New in 3.1](https://manual.karamba3d.com/new-in-karamba3d-3.1)
- [McNeel Forum — karamba.dll C++ 確認(Preisinger 本人回覆)](https://discourse.mcneel.com/t/installing-karamba-dll-for-visual-studio/143546)
- [Grasshopper Forum — Solving Speed 討論(30,000 節點瓶頸)](https://www.grasshopper3d.com/group/karamba3d/forum/topics/solving-speed)
- [Nature Architects — Karamba3D 效能評測(vs ANSYS,殼過剛)](https://nature-architects.com/en/blog/1194/)
- [Preisinger & Heimrath 2014, Structural Engineering International — Semantic Scholar](https://www.semanticscholar.org/paper/Karamba%E2%80%94A-Toolkit-for-Parametric-Structural-Design-Preisinger-Heimrath/d6c6c120dd604c74ece02f41d2619e755534764b)
