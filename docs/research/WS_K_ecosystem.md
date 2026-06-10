> 出處:研究輪 fan-out agent 查證報告(2026-06-10,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

資料已足夠。整理輸出:

# WS-K Twinmotion 對接路徑 + 開源發佈 Pipeline

## 摘要

Twinmotion 2026.1(截至 2026-06 最新)支援 Datasmith Direct Link 與靜態 `.udatasmith` 檔案匯入兩條路徑。Rhino 6-8(含 Grasshopper)、Revit、SketchUp Pro 等 12 個應用程式支援 Direct Link。**方向是單向進入**:Twinmotion 是視覺化終端,UE 才是更下游的開發環境;TM→UE 有官方 Export to Datasmith 路徑,但 UE→TM 無官方反向路徑。Rhino Datasmith 傳輸幾何、材質、圖層、Named View、Attribute User Text 自訂屬性,但**不傳頂點色**(已知缺口,Epic 列入 backlog 未實作)。Datasmith SDK 隨 UE 源碼公開,第三方可自行寫 exporter。GitHub Actions 上用 `ilammy/msvc-dev-cmd` + MSVC 編 C++、`microsoft/setup-msbuild` 編 GH .gha 插件、`crashcloud/yak-publish` 推 Yak、`softprops/action-gh-release` 發佈 zip artifact 均有成熟先例。

---

## 發現

### 1. Twinmotion 匯入路徑官方現況

**Direct Link 支援的來源(含插件需求)**:

- 3ds Max、Archicad、BricsCAD(22+)、form-Z(10+)、Modo(17+)、Revit(2018-2024+)、**Rhino 6-8(含 Grasshopper)**、RIKCAD(11+)、SketchUp Pro(2020-2026)、SOLIDWORKS(2020-2025)、Vectorworks(2022+)、Allplan(2023.1+)
- 上述除 Allplan/BricsCAD/RIKCAD/Vectorworks/Revit(2024+) 內建外,其餘須另裝 plugin
- Navisworks(2019-2026) 僅支援靜態檔案匯入,無 Direct Link

**`.udatasmith` 靜態檔案匯入**:Twinmotion 自 2020.2 起原生吃 `.udatasmith`;使用者先在來源應用程式裝 Datasmith Exporter plugin 產出 `.udatasmith` + 材質資料夾,再在 TM 匯入,**無需 Direct Link 連線**。[LIT:https://dev.epicgames.com/documentation/en-us/twinmotion/datasmith-file-import-workflow-in-twinmotion]

**TM→UE 官方路徑**:Twinmotion File → Export to Datasmith → 產出 `.udatasmith` → UE 5.4+ 以 Datasmith 匯入;支援再匯入(Reimport)。[LIT:https://dev.epicgames.com/documentation/en-us/twinmotion/twinmotion-to-unreal-engine-workflow]

**UE→TM**:**無官方路徑**。方向是「設計工具→TM→(可選)UE」,不存在「UE 專案→TM」的反向匯出。[UNKNOWN 無官方文檔]

### 2. Rhino Direct Link 能帶什麼

| 資料類型 | 是否傳輸 | 備註 |
|---|---|---|
| 幾何 mesh/solid/surface | 是 | 每個獨立物件生成一個 Static Mesh |
| 材質(顏色、貼圖) | 是 | 轉為 Material Asset |
| **圖層結構** | 是 | 每個 Rhino 圖層對應一個 Actor(父節點) |
| Block Instance | 是 | 多重 instance 保留共用 Mesh |
| Named View | 是 | 轉為 CineCameraActor |
| Mesh Modifiers(Shut Lining/Displacement 等) | 是 | 明確文檔支援 |
| **Attribute User Text(自訂屬性)** | 是 | 作為 Datasmith Metadata 傳入 |
| **Rhino Display Color** | 否 | 明確排除:「ignores Rhino display colors completely」 |
| **頂點色(Vertex Colors)** | 否 | 已知 backlog;2021 Epic 確認「未做,待排期」;至今(2026-06)無更新 |
| Grasshopper 動態資料 | 間接 | 需先 Bake 到 Rhino 才能隨 Direct Link 更新 |

[LIT:https://dev.epicgames.com/documentation/en-us/unreal-engine/using-datasmith-with-rhino-in-unreal-engine]
[LIT:https://forums.unrealengine.com/t/rhino-datasmith-mesh-vertex-color/258622]

**「分析結果」可行路徑評估**:
- **頂點色路徑(最直觀)**:在 Grasshopper/Rhino 端把利用率烤成 Mesh Vertex Color → Datasmith Direct Link。**目前不通**,頂點色不傳。
- **材質/貼圖路徑**:在 GH 端把分析結果渲染成 UV 貼圖紋理 → 貼到 Rhino 材質 → Direct Link 傳材質貼圖。**可行**,但每次分析結果改變須重新產生貼圖並重賦材質。
- **Attribute User Text 路徑**:把 D/C 比值等純量寫入 Rhino 物件的 Attribute User Text → 傳到 TM 為 Datasmith Metadata。**可行傳數值**,但 TM 本身沒有「讀 Metadata → 驅動材質顏色」的機制;此路在 UE 端才有完整 Blueprint 驅動。
- **Alternative:glTF/OBJ + 頂點色**:GH 直接輸出帶頂點色的 glTF → TM 匯入 FBX/glTF。TM 2026.1 的 Interchange Import Pipeline 更新了 FBX/OBJ 匯入;glTF 頂點色支援需確認,**但此路不走 Datasmith Direct Link**,失去即時同步。[LIT:https://dev.epicgames.com/documentation/en-us/twinmotion/twinmotion-2026-1-preview-release-notes]

### 3. Datasmith SDK

Datasmith Export SDK 位於 UE 源碼 `Engine\Source\Programs\Enterprise\Datasmith\DatasmithSDK\`。
- **公開性**:包含在 Epic Games GitHub 的 UnrealEngine 倉庫中(需簽署 UE EULA 取得存取)。
- **授權**:UE EULA;**非 MIT/Apache 等寬鬆授權**;可用於開發第三方 exporter 但不得再發佈 SDK 本體。
- **能力**:DatasmithCore API(構建場景)+ DatasmithExporter API(寫入磁碟);支援 hierarchy、material、texture、collision、LOD、metadata、camera、自訂屬性結構。
- **第三方可自寫 exporter**:是,Epic 文檔明確鼓勵且以 Revit/SketchUp/3ds Max 為參考範例。[LIT:https://dev.epicgames.com/documentation/unreal-engine/datasmith-export-sdk-guidelines]
- **限制**:不得將 Datasmith Exporter plugin 二進位再發佈給他人(可以分享下載連結)。

### 4. GitHub Actions Windows MSVC 編 C++ 標準作法

**`ilammy/msvc-dev-cmd@v1`**(最常用):
```yaml
- uses: actions/checkout@v4
- uses: ilammy/msvc-dev-cmd@v1
  with:
    arch: amd64        # 可選:x64/x86/amd64_arm64
    # vsversion: 2022  # GitHub hosted runner 通常只有一個 VS,可省略
- name: Build
  run: |
    cmake -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release .
    nmake
  shell: cmd          # 重要:用 cmd 或 pwsh,避免 bash 覆蓋 link.exe 路徑
```
⚠️ `shell: bash` 會讓 GNU `link.exe` 遮蔽 MSVC `link.exe`,導致連結失敗。[LIT:https://github.com/ilammy/msvc-dev-cmd]

**`microsoft/setup-msbuild@v3`**:適合 `.sln`/`.vcxproj` 專案:
```yaml
- uses: microsoft/setup-msbuild@v3
- run: msbuild MyProject.sln /p:Configuration=Release /p:Platform=x64
```
[LIT:https://github.com/marketplace/actions/setup-msbuild]

**Release artifact zip 發佈**(`softprops/action-gh-release@v2`):
```yaml
- name: Package
  run: Compress-Archive -Path build\Release\*.exe,build\Release\*.dll -DestinationPath release.zip
- uses: softprops/action-gh-release@v2
  with:
    files: release.zip
    generate_release_notes: true
  env:
    GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
```
觸發條件通常設 `on: push: tags: ['v*']`。[LIT:https://github.com/softprops/action-gh-release]

### 5. Grasshopper Yak 自動發佈 GitHub Actions

有兩個成熟 action:

**`paramdigma/setup-yak`**:下載 yak.exe 到 runner PATH:
```yaml
- uses: paramdigma/setup-yak@1.0.0
  with:
    token: ${{ secrets.YAK_TOKEN }}
```
[LIT:https://github.com/Paramdigma/setup-yak]

**`crashcloud/yak-publish@main`**:一步完成 build manifest + push:
```yaml
- uses: crashcloud/yak-publish@main
  with:
    package-name: 'MyPlugin'
    token: ${{ secrets.YAK_TOKEN }}
    build-path: 'src/**/bin/**/'
    publish: 'production'   # 或 'test'
    platform: win
```
[LIT:https://github.com/marketplace/actions/deploy-yak]

YAK_TOKEN 取得方式:本機執行 `yak.exe --ci` 一次性登入後取得 token 存入 GitHub Secrets。Manifest (`manifest.yml`) 最少需要 `name`/`version`/`authors`/`description`,用 `yak spec` 指令自動產生骨架。[LIT:https://developer.rhino3d.com/guides/yak/creating-a-grasshopper-plugin-package/]

**GH .gha 插件 CI 構建範例**(已有公開先例):
```yaml
# windows-latest runner
- uses: microsoft/setup-msbuild@v3
- uses: NuGet/setup-nuget@v2
- run: msbuild /t:Restore /p:Configuration=Release
- run: msbuild /p:Configuration=Release
- uses: actions/upload-artifact@v4
  with:
    name: MyPlugin
    path: src/bin/MyPlugin.gha
```
[LIT:https://hiron.dev/GrasshopperCISample/en/build-with-github-actions.html]

### 6. 自動打包模組建議形態(對 FrameCore 的具體路徑)

建議的 one-shot zip 形態:
```
framecore-vX.Y.Z-windows-x64.zip
  ├── frame_cli.exe          (Standalone CLI)
  ├── FrameCore.dll          (若未來做 GH wrapper)
  ├── FrameCore.GH.gha       (若做 Grasshopper 元件)
  ├── examples/
  │   └── *.json / *.tcl
  └── README.txt
```

版本號注入:FrameCore 已有 `git SHA build banner`;可在 CI 用 `git describe --tags --abbrev=0` 取 semver tag,以 CMake `-DFRAMECORE_VERSION=...` 或 `#define` 注入,再用 `${{ github.ref_name }}` 作 artifact 命名 `framecore-${{ github.ref_name }}-windows-x64.zip`。

---

## 對 FrameCore 的含義

**採用**:
1. **`.udatasmith` 靜態匯入是最乾淨的對接路**:FrameCore 輸出變形幾何 + 分析結果 → 用 Datasmith SDK(或現有 Rhino plugin)產出 `.udatasmith` → TM 直接開檔;不需要 Direct Link 常駐連線。
2. **分析結果着色用「材質/貼圖路徑」而非頂點色路徑**:因為頂點色 Rhino→Datasmith 目前不通。具體做法是在 Grasshopper 端生成 UV 貼圖(D/C 色彩圖 png)後賦到材質,再 bake 到 Rhino 走 Direct Link/靜態匯出。
3. **Attribute User Text 傳數值到 TM/UE**:可把每個 Member 的 D/C 比值、安全係數寫入 Rhino 物件屬性,在 UE 端用 Blueprint 讀 Datasmith Metadata 驅動材質顏色 — 這才是完整管線。TM 側讀 Metadata 能力有限,主要在 UE 側用。
4. **發佈 pipeline**:一支 `.github/workflows/release.yml`,觸發條件 `push: tags: v*`;`ilammy/msvc-dev-cmd` + CMake NMake 編 `frame_cli.exe`(現有 `build.bat` 邏輯搬過來);`Compress-Archive` 打包;`softprops/action-gh-release` 上傳。若未來做 GH wrapper(.gha),加 `microsoft/setup-msbuild` + `crashcloud/yak-publish` 步驟推 Yak。

**避免**:
- 不要依賴頂點色通路傳分析著色到 TM(Epic 明確說不支援,無確切修復時間表)。
- 不要嘗試 UE 專案→TM 反向匯出(無官方支援)。
- GitHub Actions `shell: bash` 搭配 MSVC 有 `link.exe` 遮蔽風險,改用 `shell: cmd` 或 `shell: pwsh`。
- Datasmith SDK 非寬鬆開源,不得將其二進位直接打包再發佈;引用方式應在文檔中說明「需使用者自行裝 UE 源碼或 Epic 官方 plugin」。

**誠實邊界**:
- TM 側讀 Datasmith Metadata 的 UI/Blueprint 能力有限,「分析結果驅動著色」在 TM 端是二等公民,在 UE 端才是一等公民。
- Yak 發佈需 McNeel 帳號 + `yak.exe --ci` 取 token,**完全自動化需要先手動跑一次本機登入**。
- TM 2026.1 改用 Interchange Import Pipeline 處理 FBX/OBJ,glTF 頂點色支援狀態尚未確認,需實測。[UNKNOWN 2026.1 glTF 頂點色支援細節]

---

## 來源清單

- [Twinmotion Supported Design Applications for Datasmith](https://dev.epicgames.com/documentation/en-us/twinmotion/supported-design-applications-for-datasmith)
- [Datasmith Direct Link Workflow in Twinmotion](https://dev.epicgames.com/documentation/en-us/twinmotion/datasmith-direct-link-workflow-in-twinmotion)
- [Datasmith File Import Workflow in Twinmotion](https://dev.epicgames.com/documentation/en-us/twinmotion/datasmith-file-import-workflow-in-twinmotion)
- [Twinmotion to Unreal Engine Workflow](https://dev.epicgames.com/documentation/en-us/twinmotion/twinmotion-to-unreal-engine-workflow)
- [Using Datasmith with Rhino in Unreal Engine](https://dev.epicgames.com/documentation/en-us/unreal-engine/using-datasmith-with-rhino-in-unreal-engine)
- [Using Datasmith with Rhino - Twinmotion Docs](https://dev.epicgames.com/documentation/en-us/twinmotion/using-datasmith-with-rhino)
- [Rhino Datasmith mesh vertex color forum thread (Epic)](https://forums.unrealengine.com/t/rhino-datasmith-mesh-vertex-color/258622)
- [Datasmith Export SDK Guidelines - Unreal Engine 5.7](https://dev.epicgames.com/documentation/unreal-engine/datasmith-export-sdk-guidelines)
- [Twinmotion 2026.1 Preview Release Notes](https://dev.epicgames.com/documentation/en-us/twinmotion/twinmotion-2026-1-preview-release-notes)
- [ilammy/msvc-dev-cmd GitHub Action](https://github.com/ilammy/msvc-dev-cmd)
- [microsoft/setup-msbuild on GitHub Marketplace](https://github.com/marketplace/actions/setup-msbuild)
- [softprops/action-gh-release](https://github.com/softprops/action-gh-release)
- [crashcloud/yak-publish - Deploy Yak Action](https://github.com/marketplace/actions/deploy-yak)
- [paramdigma/setup-yak GitHub Action](https://github.com/Paramdigma/setup-yak)
- [McNeel - Creating a Grasshopper Plug-In Package (Yak)](https://developer.rhino3d.com/guides/yak/creating-a-grasshopper-plugin-package/)
- [GrasshopperCISample - Build with GitHub Actions](https://hiron.dev/GrasshopperCISample/en/build-with-github-actions.html)
- [Twinmotion 2025.2 Release Notes](https://dev.epicgames.com/documentation/en-us/twinmotion/twinmotion-2025-2-release-notes)
