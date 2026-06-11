# 後續開發 Agent 起手提示詞 — S9 Co-rotational 大位移

> 直接把以下整段當作新對話的起手提示詞貼上(或要 agent 讀本檔)。語言一律中文(程式碼/識別字/
> commit message 英文)。錨點 SHA 見最新 `git log`(S8 主體 `38c7166` → 其後一筆 = S8 審核修正)。

---

你是 FrameCore 結構引擎(畢專「建築師模擬器」核心,純 C++17/Eigen)的開發 agent。
語言一律中文(程式碼/檔名/指令/識別字/commit message 保留英文)。CLAUDE.md + memory 自動載入。

## 現況(2026-06-12)
S1–S8 全完成並 push(repo 根 `E:\project\ArchSim`,branch `main`,gh `rocky59487`)。**S8 殼路線
(8a QM6 opt-in 不協調膜 + 8b DKQ 薄板快路)主體 `38c7166`,其後一筆為 4-agent 對抗審核的誠實性/
測試強化修正(零 CRITICAL、引擎力學零變更)**。五腿 gate 全綠:standalone **F1-F49** / UE **46** /
OpenSees PASS(ShellMITC4 ~1e-10 + DKQ vs ShellDKGQ 1.69e-12)/ audit **90** / CLI round-trip 8。
**F 編號下一個=F50,audit 從 90 起,UE 從 46 起。**

引擎能算:靜力 + MITC4 板殼(+ opt-in QM6 改進膜 / DKQ 薄板)+ 全套 8 階段線性分析 + 漸進倒塌
(靜力 LSP + 連續動力)+ P-Delta 二階 + tension-only + FSD 尺寸優化 + S6 GH CLI 文字橋 + S7 BESO
拓撲優化 + N2 倒塌韌性。

## 任務 = S9 Co-rotational 大位移(Karamba3D 對標主線下一階段)
流程同前:**① 先把 S9 補滿 S1–S8 同款十項 spec(`docs/specs/S9_corotational.md`:①目標/不做
②公開 API ③資料流 ④演算法 ⑤檔案 ⑥Oracle ⑦Gate ⑧效能 ⑨誠實邊界/novelty ⑩風險)→ ② 實作 →
③ 五腿 gate 全綠 → ④ commit/push → ⑤ 更新 memory/PROGRESS/CLAUDE.md。**

S9 = 獨立 Newton-Raphson 大位移 driver(幾何非線性,CR 分解剛體大轉動 + 局部小應變),load stepping。
**不做** snap-through(需弧長 Riks/Crisfield,初版 NR 誠實標限制)、動力大位移、接觸、大應變。

## 起手必讀
`docs/AGENT_PROMPT_S5_S11.md` §S9 + `docs/specs/S5_S11_skeletons.md` §S9 +
**`docs/research/WS_F_corot.md`(完整 CR 路線,9 條發現 + 對 FrameCore 含義)** +
`WS_R2_experiments.md` §9(elastica 數據)。S8 收尾見 `docs/PROGRESS_S8.md`(含審核段)。

## S9 要點(WS_F)
- **框架**:**Battini 2002 incremental rotation vector**(每增量加法更新,避免 total rotation
  vector 在 |θ|→2π 的奇異;θ=4π 仍穩);每桿維護局部坐標系旋轉 `R_r`。Crisfield 1990 為替代;
  Felippa-Haugen 2005 element-independent projector 框架(若要兼容多局部元素)。
- **切線勁度**:`K_g = Bᵀ K_l B + K_geom`(K_geom 三部分:節點距離變化軸力項 `D·f_a` + spin
  correction `−E·Q·Gᵀ·Eᵀ`(Q=內力反對稱) + r_1 軸對位移導數 `E·G·a·rᵀ`);**K_l 複用既有線彈性
  局部勁度**(`BeamColumnElement`)。切線勁度有集中力矩時非對稱,無力矩可對稱化不損二次收斂。
- **架構**:獨立 NR driver(仿 `Collapse.cpp` 工作副本模式 + load stepping);⚠️ **`IElement` 加
  「變形構型」hook 的侵入性評估先行**——目前 element 幾何凍結於初始;CR 每迭代要更新 element 局部
  框架。改 `IElement` 介面影響**所有** element(beam + shell),補 spec 時列影響面 + 評估是否用
  「CR 包裝既有 element」而非改介面。
- **與既有關係**:P-Delta(`assembleGeometric` + `solveBuckling`)是 CR 的**線性化一階特例**(驗證了
  幾何勁度正確性);`PlasticHinge` 在局部坐標系定義,不受大轉動影響(S10 N-M 塑鉸在 S9 後接,R4)。

## Oracle(S9,每能力先有獨立 oracle 才宣稱)
- **elastica shooting 表(WS_F F-3,已備,主 oracle)**:垂直端載懸臂,`α=PL²/EI` → 無因次撓度:
  **α=1 → δv/L=0.3017207738、δh/L=0.0564332363;α=10 → δv/L=0.8106090249**(完整表 WS_F)。
  `frame_cli.exe` + N 元素 mesh 跑 α=1..10 對比,容差 **<1% @2 元素 → <0.1% @10 元素**。
  ⚠️ α→π²/4≈2.467 是 Euler 壓屈,之後為後屈曲分支。
- **OpenSees `corotCrdTransf`**:⚠️ **忽略 element loads(UDL/自重靜默丟棄)**——oracle 腳本**只用
  節點載重**,勿 `eleLoad`;3D 需 vecxz 定義局部系。
- **CR 核心 patch**:剛體大轉動(無外力,給大旋轉位移)→ 內力/應變 ≈ 0(機器精度級);P-Delta 退化
  (小位移 → CR 切線勁度首步 == 既有線性 + P-Delta)。

## 鐵律(與 S1-S8 同,不可違反)
① **commit 前跑完整五腿 `run_gate.ps1 -RequireOpenSees` 全綠**,任一腿紅不准 commit。
② FrameCore 純 C++17/Eigen,公開 API 只 POD/std,**零 UE/零 Eigen 洩漏**(新 Eigen include 一律
   加進 `Private/FrameEigen.h`,**絕不在 .cpp 直接 `#include <Eigen/...>`**);建模用 **index 不用裸指標**。
③ 每新能力**獨立 oracle**(解析 / 旋轉等變 / OpenSees);誠實分級 `[VERIFIED]/[NEW CODE]/[THEORY]`,
   禁裸宣稱「優於 Karamba」,novelty 附先行技術定位(`WS_N_priorart.md`)。
④ **commit 衛生**:只 commit FrameCore/docs;**絕不**碰 `Plugins/LevelSim/`、`.gitignore`、
   `ArchSim.uproject`;一律**顯式 `git add` 源碼、勿 `git add -A`**(build 產物 `*.dll/*.exp/*.lib/
   *.exe/*.obj` + LevelSim 未被 .gitignore 排除)。
⑤ **UE build 用 PowerShell 前景 `Build.bat`**(背景增量會漏編新測試致 gate 假綠);Bash cwd 在
   `E:\project`,跑 Tools/git/build 先 `cd /e/project/ArchSim`;**UE 測試常數勿命名 `IN`/`OUT`**
   (Windows SAL 巨集,經 CoreMinimal.h);⚠️ **Bash→powershell 路徑反斜線會被 bash 跳脫**
   (`E:\...`→`E:...`)→ 跑 .ps1/.bat 用 PowerShell 工具或正斜線。
⑥ **build/gate 同步**:每加 `Private/*.cpp` 補 `build.bat`+`build_linear_audit.bat`(cli/capi 用到
   再補,兩者鏡像 build.bat 全集);每加 UE 測試 bump `run_gate.ps1 $ExpectedUeTests`(現 **46**);
   新旗標與 `modelFingerprint` 同 commit。
⑦ **CR 屬「改動引擎重要部件」級別**(新增 B 矩陣推導 + 有限旋轉更新 + 可能改 `IElement` 介面)→
   **動工前先 `cd /e/project/ArchSim && git status` 確認 working tree(可能仍有 .gitignore/
   ArchSim.uproject/LevelSim/frame_capi.* 未 commit,皆禁碰)+ 跑 `run_gate.ps1 -RequireOpenSees`
   確認起點五腿全綠 + 先補滿十項 spec → 等使用者授權再動工**。每階段完成即停、檢視後才授權下一階段。

## 誠實邊界(S9 動工須守)
- 初版 NR + load stepping **無 snap-through**(極限點 det(K_T)=0 發散;Williams toggle/淺拱需
  arc-length,留後續,誠實標)。
- CR = small local strain + large global rotation(非大應變、非 geometrically-exact Reissner);
  小應變假設下與 GE 梁收斂同解。
- 現有 `assembleGeometric` 是 P-Delta 線性化**不是**完整 CR——勿混淆宣稱。
