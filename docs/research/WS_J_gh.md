> 出處:研究輪 fan-out agent 查證報告(2026-06-10,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

現在我有足夠的資料可以撰寫完整的研究報告。

# WS-J Rhino 8 / Grasshopper 插件開發工具鏈與外部求解器橋接

## 摘要

Rhino 8 在 Windows 上同時支援 .NET 7.0(Core)與 .NET Framework 4.8 雙軌執行,推薦插件採多目標(multi-target)建置。GH_Component 子類的核心生命週期為 RegisterInputParams/RegisterOutputParams/SolveInstance 三件套;DataTree 透過 GH_Structure<T> 傳遞。官方以 McNeel RhinoVisualStudioExtensions dotnet 模板為標準起手式,.gha 是改名的 .NET assembly。Yak 是 CLI 套件管理工具,Food4Rhino 提供 UI 介面轉介至 Yak 發佈。外部 exe 橋接可行但有 UI 阻塞風險;最穩健的 async 模式為「背景執行緒 + RhinoApp.InvokeOnUiThread + GH_Document.ScheduleSolution」,或採用 Speckle 的 GrasshopperAsyncComponent 框架。Hops/Rhino.Compute 是 HTTP 版替代,適合雲端或長時運算但引入網路延遲。Karamba3D 的三層架構(C++ karamba.dll → .NET karambaCommon.dll → GH karamba.gha)是 FrameCore 橋接的最直接參考模型。C ABI P/Invoke 路徑技術可行,但結構 marshalling 成本與版本管理複雜度高於 CLI exe 橋接。

---

## 發現

### 1. Rhino 8 .NET 版本與插件目標框架

- **Rhino 8.0–8.19 預設 .NET 7.0;Rhino 8.20+ 升至 .NET 8.0。** McNeel 官方立場:「不在穩定版中間更換執行時期」,因此 net8.0 在 Rhino 8 週期內屬非官方(可跑但不受支援)。Rhino 9 才會正式切 net8+。[LIT:https://discourse.mcneel.com/t/net-core-8/200123]
- **推薦雙目標**:`<TargetFrameworks>net48;net7.0</TargetFrameworks>`,Windows 可選跑 .NET Framework 4.8(相容舊插件);Mac 只有 net7.0。[LIT:https://developer.rhino3d.com/guides/rhinocommon/moving-to-dotnet-core/]
- **NuGet 上的 Grasshopper SDK** 版本號如 `Grasshopper 8.17.25066.7001`,目前仍以 netstandard2.0/net48 為主要 TFM,供 Rhino 8 雙軌載入。[LIT:https://www.nuget.org/packages/Grasshopper/]
- **RhinoVisualStudioExtensions** 模板(McNeel 官方 VS 擴充)自動生成雙目標 .csproj,內含 Rhino 8 netfx/netcore 的 VS 除錯啟動配置與 Yak 打包任務。[LIT:https://github.com/mcneel/RhinoVisualStudioExtensions/releases]

### 2. GH_Component 子類生命週期

- **繼承 `GH_Component`**,覆寫:
  - `RegisterInputParams(GH_InputParamManager pManager)` — 宣告輸入埠(型別、名稱、nickname、description、可選性)
  - `RegisterOutputParams(GH_OutputParamManager pManager)` — 宣告輸出埠
  - `SolveInstance(IGH_DataAccess DA)` — 每次求解呼叫一次(每條 DataTree 路徑各一次,若 access 設為 item);讀取用 `DA.GetData / DA.GetDataList / DA.GetDataTree`,寫出用 `DA.SetData / DA.SetDataList / DA.SetDataTree`
- **DataTree** 在 SDK 側型別為 `GH_Structure<T>`,`T` 需實作 `IGH_Goo`;`DA.GetDataTree<T>(index, out GH_Structure<T> tree)` 取出整棵樹。[LIT:https://developer.rhino3d.com/api/grasshopper/html/M_Grasshopper_Kernel_IGH_DataAccess_GetDataTree__1.htm]
- **輸出在每次 SolveInstance 前自動清除**,不需手動 Clear。[LIT:https://developer.rhino3d.com/api/grasshopper/html/M_Grasshopper_Kernel_GH_Component_SolveInstance.htm]
- **.gha 檔即 .NET assembly 改名**:Grasshopper 在啟動時掃描 `%APPDATA%\Grasshopper\Libraries` 目錄下所有 .gha,用反射載入所有 `GH_Component` 子類。[LIT:https://developer.rhino3d.com/guides/grasshopper/your-first-component-windows/]

### 3. 除錯流程

- Visual Studio:按 F5 → VS 啟動 Rhinoceros → GH 自動載入 .gha → 在 SolveInstance 設中斷點即可攔截。專案範本已預設 `StartProgram = Rhino.exe`。[LIT:https://developer.rhino3d.com/guides/grasshopper/your-first-component-windows/]
- VS Code:使用 coreclr debugger(C# extension),不需要額外啟動專案。[LIT:https://developer.rhino3d.com/guides/rhinocommon/moving-to-dotnet-core/]

### 4. Yak 套件管理器與 Food4Rhino

- **Yak** = Rhino Package Manager 的 CLI 工具,位於 `C:\Program Files\Rhino 8\System\yak.exe`(Windows)。[LIT:https://developer.rhino3d.com/guides/yak/what-is-yak/]
- **標準指令**:
  - `yak spec` — 從 assembly 屬性自動產生 `manifest.yml` 骨架
  - `yak build` — 將目錄打包為 .yak(ZIP 改名)
  - `yak push` — 上傳至 McNeel package server
- **manifest.yml 格式**:
  ```yaml
  name: my-plugin
  version: 1.0.0
  authors:
  - Author Name
  description: 描述
  url: https://example.com
  icon: icon.png
  keywords: [structural, analysis]
  ```
- **Rhino 8 多目標套件**:manifest.yml 必須置於 framework 子目錄之外(頂層),內部可有 `net48/` 與 `net7.0/` 子目錄分別放對應 .gha。[LIT:https://developer.rhino3d.com/guides/yak/the-anatomy-of-a-package/]
- **Food4Rhino 關係**:Food4Rhino 是 McNeel 維護的插件展示平台;它提供 UI 按鈕「Push to Yak」,自動把已上傳的 .gha 推送至 Package Manager server。兩者互補:Food4Rhino 做展示/社群,Yak 做分發/自動更新。[LIT:https://www.food4rhino.com/en/faq]

### 5. 從 GH Component spawn 外部 exe 的實務

- **基本模式**:使用 `System.Diagnostics.Process`,設 `RedirectStandardInput = true`、`RedirectStandardOutput = true`、`UseShellExecute = false`,呼叫 `process.StandardInput.WriteLine(json)` 送資料,`process.StandardOutput.ReadToEnd()` 讀結果,最後 `process.WaitForExit()`。[LIT:https://learn.microsoft.com/en-us/dotnet/api/system.diagnostics.process.standardoutput?view=net-8.0]
- **阻塞問題**:SolveInstance 跑在 GH UI 執行緒;同步呼叫 `WaitForExit()` 會凍結整個 Rhino UI,使用者無法取消或拖動視窗。[LIT:https://discourse.mcneel.com/t/custom-component-wait-for-external-process-without-locking-gh-ui/71613]
- **推薦 async 模式**:
  1. SolveInstance 第一次進入 → 啟動背景 Task/Thread
  2. 背景完成後呼叫 `RhinoApp.InvokeOnUiThread(new Action(() => component.ExpireSolution(true)))` 或 `GH_Document.ScheduleSolution(delay_ms, callback)`
  3. SolveInstance 再次被呼叫時輸出快取結果
  - `ScheduleSolution(0, ...)` = 遞迴觸發(小心 stack overflow);`ScheduleSolution(1+, ...)` = 計時器觸發。[LIT:https://developer.rhino3d.com/api/grasshopper/html/Overload_Grasshopper_Kernel_GH_Document_ScheduleSolution.htm]
- **GrasshopperAsyncComponent 框架**:Speckle 開源的 `GH_AsyncComponent` 抽象基底類別,封裝「eager restart + cancellation token + DoWork on threadpool + InvokeOnUiThread expire」完整流程;已有 NuGet 套件 `GrasshopperAsyncComponent 2.0.1`。這是目前社群最成熟的 async 元件框架。[LIT:https://github.com/specklesystems/GrasshopperAsyncComponent]
- **手動觸發/debounce**:常見做法是在元件上加 Boolean Toggle 或 Button 輸入(只在上緣觸發一次);或繼承 `GH_Component` 覆寫 `CreateAttributes()` 自訂含按鈕的 attribute 類。社群廣泛採用此模式避免 slider 拖動過程反覆觸發重量計算。[LIT:https://discourse.mcneel.com/t/triggering-a-plugin-components-button/90478]
- **現有「外部 exe 求解器」模式的插件參考**:
  - **Alpaca4d (OpenSeesGH)**:Grasshopper 插件建在 OpenSees 之上,C# 94.7%、Tcl 4.5%,Tcl 部分暗示採用 OpenSees Tcl 腳本執行模式(可能是產生 Tcl 腳本 → 呼叫 OpenSees exe → 讀結果,但原始碼未公開確認)。[LIT:https://github.com/Alpaca4d/Alpaca4d][UNKNOWN:確切橋接模式(exe spawn vs embedded dll)未從公開文件確認]
  - **OSforGH**:另一個 OpenSees for Grasshopper,架構類似。[LIT:https://github.com/strdesigner/OSforGH]

### 6. Hops / Rhino.Compute 模式

- **Hops** 是 GH 內建元件(Rhino 7+),把目標 .gh 定義包裝成 HTTP POST 送至 Rhino.Compute 服務(或任何相容 REST endpoint);非同步求解,UI 不阻塞;支援平行多實例。[LIT:https://developer.rhino3d.com/guides/compute/hops-component/]
- **Rhino.Compute**:McNeel 開源的無頭 web server,在 ASP.NET Core 上跑 Rhino.Inside + Grasshopper,接受 JSON 格式定義/輸入 → 回傳結果。[LIT:https://developer.rhino3d.com/guides/compute/how-hops-works/]
- **對 FrameCore 橋接的比較**:

| 面向 | CLI exe 橋接(stdin/stdout) | Hops/HTTP | C ABI P/Invoke |
|---|---|---|---|
| 實作成本 | 低:寫 JSON 序列化 + Process | 中:需跑 compute server | 高:struct layout + DllImport |
| 延遲 | ~10–100ms 啟動(若常駐則近零) | HTTP round-trip:~1–50ms(本機) | 最低:零 IPC |
| UI 阻塞風險 | 有(需 async 包裝) | 無(Hops 原生 async) | 有(若同步呼叫) |
| 偵錯 | 簡單:exe 可獨立跑 | 複雜:需 server + client 同步 | 中:需 native debugger |
| 跨平台 | Windows only(FrameCore 本就 x64 Windows) | 理論上跨平台 | Windows x64 only |
| 部署 | 一個 exe 放旁邊 | 需安裝 Rhino.Compute | DLL 放旁邊 |

### 7. Karamba3D 核心元件流(UX 參考)

Karamba3D 三層架構:C++ `karamba.dll`(數值核心)→ .NET `karambaCommon.dll`(公開 API,獨立於 Rhino/GH)→ GH `karamba.gha`(元件 UI 層)。[LIT:https://scripting.karamba3d.com/1.-introduction/1.1-scripting-with-karamba3d]

**元件鏈(從上游到下游)**:

```
幾何輸入
  ↓
[Line to Beam] / [Mesh to Shell]    ← 幾何 → FE 元素轉換
  ↓
[Cross Section] + [Material Properties]  ← 截面/材料定義
  ↓
[Support]                           ← 邊界條件
  ↓
[Loads] / [Beam Loads] / [Prescribed Displacement]  ← 載重
  ↓
[Beam-Joints]                       ← 鉸接條件(可選)
  ↓
[Assemble Model]                    ← 把以上全部組裝成 IModel
  ↓
[Analyze] / [AnalyzeThII]           ← 求解(一階/二階)
  ↓
[Natural Vibrations] / [Buckling Modes]  ← 進階分析(可選)
  ↓
[Nodal Displacements] / [Beam Forces] / [Shell Forces]
[Reaction Forces] / [Utilization of Elements]  ← 結果讀取
```

- **Assemble Model** 是中心匯聚點:接收所有元素後輸出 `IModel` 物件,後續元件都以此物件為輸入。這個「單一 model 物件沿線傳遞」的模式讓參數化更新清晰。[LIT:https://grasshopperdocs.com/components/karamba3d/assembleModel.html]
- **Karamba3D 是商業閉源軟體**,提供 NuGet 套件 `KarambaCommon` 供第三方擴展,但核心求解器 karamba.dll 不開放原始碼。[LIT:https://manual.karamba3d.com/1-introduction/1.2-licenses]

### 8. C ABI P/Invoke 路徑評估

- **基本機制**:`[DllImport("framecore.dll", CallingConvention = CallingConvention.Cdecl)]` 聲明靜態 extern 方法,blittable 型別(int/float/double/指標)可直接傳遞;struct 需加 `[StructLayout(LayoutKind.Sequential)]`。[LIT:https://discourse.mcneel.com/t/help-with-invoking-c-dll-in-grasshopper-component/18818]
- **x64 限制**:Rhino 8 是純 x64,FrameCore 也是 x64,架構一致,無 32/64 混用問題。[THEORY]
- **工具鏈成本**:
  - 需維護一套 C 語言導出的 `extern "C"` wrapper 函式(不能直接導出 C++ class,因 ABI 不穩定)
  - 複雜資料結構(Eigen 矩陣、std::vector)無法直接 marshal,必須用 flat POD 陣列橋接
  - 版本升級時 struct layout 改動會靜默 crash,需嚴格版本鎖定
- **與 CLI exe 橋接比較**:P/Invoke 延遲最低(零 IPC),但開發/維護成本顯著高於 CLI exe 橋接;CLI exe 橋接可直接重用 FrameCore 現有的 `frame_cli.exe` 讀寫 JSON 的能力,且除錯極為方便。[THEORY]

---

## 對 FrameCore 的含義

**1. 目標框架選擇**:GH 插件專案採 `<TargetFrameworks>net48;net7.0</TargetFrameworks>` 雙目標。FrameCore standalone exe(`frame_cli.exe`)本身是純 C++,不受 .NET 版本影響;橋接層(C# GH component)採 net7.0 即可滿足 Rhino 8 全系列。

**2. 橋接模式推薦**:採 **CLI exe 橋接**(stdin/stdout JSON)而非 P/Invoke DLL,原因:
- FrameCore 已有 `frame_cli.exe` 且 JSON I/O 格式成熟
- 開發成本最低,C++ 側不需新增任何 extern "C" wrapper
- 除錯簡單(exe 可單獨跑,可加 --verbose)
- async 包裝用 `GrasshopperAsyncComponent` 框架或自實作「背景 Task + RhinoApp.InvokeOnUiThread + ScheduleSolution」模式,避免 UI 阻塞
- 唯一代價:每次 solve 有 exe 啟動延遲(~50ms);若需極低延遲可改為「持久化常駐 exe + named pipe」模式

**3. Karamba3D 元件鏈 UX 模板**:仿照其「Assemble Model → Analyze → Results」三階段設計 GH 元件組:
- `FrameCore.BuildModel`(輸入節點、梁、殼、截面、材料 → 輸出 `FrameModelGoo`)
- `FrameCore.Analyze`(輸入 FrameModelGoo + 載重工況 → 輸出 `FrameResultGoo`)
- `FrameCore.GetDisplacements` / `FrameCore.GetBeamForces` / `FrameCore.GetUtilization` 等結果讀取元件
- 中間物件用自訂 `GH_Goo<T>` 包裝,在 GH wire 上傳遞,不暴露 Eigen

**4. 發佈路徑**:先用 `yak spec` + `yak build` 製作 .yak,上傳 Food4Rhino 並按 Push to Yak,即可讓使用者在 Rhino 8 Package Manager 內一鍵安裝。manifest.yml 置於多目標套件頂層。

**5. 誠實邊界**:GH 插件內若採同步 exe 橋接,求解時 Rhino UI 會凍結;**必須**實作 async 模式才能給使用者可用的互動體驗。Hops/HTTP 模式適合需要雲端/離線計算的場景,但對本機 FrameCore 增加無謂的網路層成本,不推薦為首選。

---

## 來源清單

- [Rhino Developer Documentation](https://developer.rhino3d.com/)
- [Moving to .NET Core - Rhino Developer](https://developer.rhino3d.com/guides/rhinocommon/moving-to-dotnet-core/)
- [Your First Grasshopper Component (Windows)](https://developer.rhino3d.com/guides/grasshopper/your-first-component-windows/)
- [GH_Component API Reference](https://developer.rhino3d.com/api/grasshopper/html/T_Grasshopper_Kernel_GH_Component.htm)
- [IGH_DataAccess.GetDataTree Method](https://developer.rhino3d.com/api/grasshopper/html/M_Grasshopper_Kernel_IGH_DataAccess_GetDataTree__1.htm)
- [GH_Document.ScheduleSolution Method](https://developer.rhino3d.com/api/grasshopper/html/Overload_Grasshopper_Kernel_GH_Document_ScheduleSolution.htm)
- [Asynchronous Execution - Rhino Developer](https://developer.rhino3d.com/guides/scripting/advanced-async/)
- [RhinoVisualStudioExtensions Releases](https://github.com/mcneel/RhinoVisualStudioExtensions/releases)
- [What is Yak? - Package Manager Guide](https://developer.rhino3d.com/guides/yak/what-is-yak/)
- [Anatomy of a Yak Package](https://developer.rhino3d.com/guides/yak/the-anatomy-of-a-package/)
- [Creating a Rhino Plugin Package](https://developer.rhino3d.com/guides/yak/creating-a-rhino-plugin-package/)
- [Food4Rhino FAQ](https://www.food4rhino.com/en/faq)
- [The Hops Component - Rhino Developer](https://developer.rhino3d.com/guides/compute/hops-component/)
- [How Hops Works](https://developer.rhino3d.com/guides/compute/how-hops-works/)
- [Karamba3D Grasshopper Docs - Component List](https://grasshopperdocs.com/addons/karamba3d.html)
- [Karamba3D Assemble Model Component](https://grasshopperdocs.com/components/karamba3d/assembleModel.html)
- [Karamba3D Scripting Guide 3.1.4](https://scripting.karamba3d.com/1.-introduction/1.1-scripting-with-karamba3d)
- [Karamba3D Licensing](https://manual.karamba3d.com/1-introduction/1.2-licenses)
- [GrasshopperAsyncComponent (Speckle)](https://github.com/specklesystems/GrasshopperAsyncComponent)
- [NuGet: GrasshopperAsyncComponent 2.0.1](https://www.nuget.org/packages/GrasshopperAsyncComponent/)
- [Best practices for async in SolveInstance - McNeel Forum](https://discourse.mcneel.com/t/best-practices-to-run-asynchronous-code-inside-solveinstance/165737)
- [Custom component wait for external process - McNeel Forum](https://discourse.mcneel.com/t/custom-component-wait-for-external-process-without-locking-gh-ui/71613)
- [.NET Core 8 status for Rhino 8 - McNeel Forum](https://discourse.mcneel.com/t/net-core-8/200123)
- [How to Develop GH Components for Rhino 8 Using .NET 7 - McNeel Forum](https://discourse.mcneel.com/t/how-to-develop-grasshopper-components-for-rhino-8-using-net-7/165957)
- [Alpaca4d (OpenSeesGH) GitHub](https://github.com/Alpaca4d/Alpaca4d)
- [OSforGH (OpenSees for Grasshopper) GitHub](https://github.com/strdesigner/OSforGH)
- [P/Invoke with C++ DLL in GH - McNeel Forum](https://discourse.mcneel.com/t/help-with-invoking-c-dll-in-grasshopper-component/18818)
- [System.Diagnostics.Process.StandardOutput - Microsoft Learn](https://learn.microsoft.com/en-us/dotnet/api/system.diagnostics.process.standardoutput?view=net-8.0)
