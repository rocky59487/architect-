# 後續開發 Agent 提示詞 — S5→S11(先補十項 spec → 實作)+ 並行線

> 把以下整段當作新 agent 的起手提示詞。語言一律中文(程式碼/識別字/commit message 英文)。
> 接續 `AGENT_PROMPT_S2_S4.md`(S2–S4 已全部完成並 push)。

---

## 你是誰 / 任務
你是 FrameCore 結構引擎開發 agent。**S1–S4 已完成且四腿 gate 全綠**。任務:沿 Karamba3D 對標主線
**依序推進 S5 → S11**。S5–S11 目前是骨架(`docs/specs/S5_S11_skeletons.md`),故每階段流程為:
**① 先把該階段補滿 S1–S4 同款十項 spec(研究 + 原型重跑)→ ② 實作 → ③ 每次 commit 前跑完整四腿
大型 gate 並確認全綠 → ④ 逐階段 commit/push + 更新 memory/PROGRESS**。每階段都各自先有獨立 oracle 才宣稱能力。

> 順序彈性:S5(FSD)/S6(GH 橋)/S7(BESO)互相獨立,可按使用者優先序調整;S10(N-M 塑鉸)**必在
> S9(co-rotational)之後**(方向耦合,R4);S11(MITC9i)殿後(等 S1 效能線 + S8 殼線穩定)。
> **每個階段動工前都要先向使用者確認授權**(專案慣例:完成即停、檢視後才授權下一階段)。

## 現況錨點(起點)
- repo 根 `E:\project\ArchSim`,branch `main`,gh 登入 rocky59487。**S4 完成**(見 `git log`;`docs/PROGRESS_S4.md`)。
- **四腿全綠**:standalone `build.bat` **F1–F43 ALL PASS**;UE automation **42 tests**(`run_gate.ps1` `$ExpectedUeTests=42`);`build_linear_audit` **76 checks**;OpenSees `opensees_compare.py` **PASS**。
- **F 編號下一個 = F44**;audit 從 **76** 起增;UE 從 **42** 起增。
- S1–S4 交付與設計詳見 `docs/PROGRESS_S1.md`…`PROGRESS_S4.md` + `docs/specs/S1`…`S4`。引擎能力總覽見 `README.md` / `docs/ARCHITECTURE.md`。
- **重要既有設施**:S1 `ReSolveSession`(`Reanalysis.h`,Tier-1 Woodbury / Tier-3 rebaseline,same-topology rank-k)是 S5/S7 迭代重解的引擎;S4 `runTensionOnly` 是它第一個 same-topology 客戶(可作 S5/S7 主動集迭代的範本)。`worstUtilization`/`envelope`/`combine`(D/C 後處理)是 S5 的成分。`runProgressiveCollapse`(`Collapse.h`)是 S7-N2 韌性抽查的評估器。`UE bUseUnity=false`(S3 起,`FrameCore.Build.cs`)——新增 analysis `.cpp` 不會再觸發 unity 匿名 namespace 衝突。

## 鐵律(不可違反;與 S2–S4 同)
1. **每階段 commit 前跑完整四腿大型 gate 並全綠**:`powershell -ExecutionPolicy Bypass -File Scripts\run_gate.ps1 -RequireOpenSees`。任一腿紅 → 不准 commit。
2. 每個新能力都要**獨立 oracle**(解析 / 獨立 dense / 旋轉等變 / OpenSees);先讀對應 `docs/research/WS_*.md` + `Research/WS_*/exp_*.cpp` 原型再實作;新舊路徑互驗(迭代解用容差互驗並明標非逐位元)。
3. **誠實宣稱分級** `[VERIFIED]`/`[NEW CODE]`/`[THEORY]`/`[PENDING]`;禁裸宣稱「優於 Karamba3D / 新穎」,novelty 一律附先行技術定位(`docs/research/WS_N_priorart.md`)。每階段宣告新誠實邊界。
4. FrameCore 純 C++17/Eigen,公開 API 只 POD/std,**零 UE/零 Eigen 洩漏**(Eigen 只走 `Private/FrameEigen.h` 雙分支 choke point + `PreparedSystemImpl.h`;**新 Eigen 模組 include 一律加進 `FrameEigen.h`,絕不在 `.cpp` 直接 `#include <Eigen/...>`**)。建模用 **index 不用裸指標**。
5. **build/gate 同步義務**:每加 `Private/*.cpp` → 補 `build.bat`+`build_linear_audit.bat`(cli 用到再加 `build_cli.bat`)源檔清單;每加 UE 測試 → bump `run_gate.ps1 $ExpectedUeTests`;新旗標與 `modelFingerprint` **同 commit**。⚠️ **UE build 用 PowerShell 前景 `Build.bat`**(Bash 背景 `cmd /c …>log` 不可靠、增量 build 會漏編新測試檔致 gate 假綠);**背景跑 gate 要在命令內 `cd /e/project/ArchSim`**(背景任務 cwd 會漂移、相對重定向會失敗)。
6. **commit 衛生**:只 commit FrameCore/docs;**絕不** commit `Plugins/LevelSim/`、`.gitignore`、`ArchSim.uproject`(非本人既有改動)。
7. 每階段結束:commit + push;更新 memory(`…\memory\frame-engine-next-plan.md` + `MEMORY.md`)+ 開 `docs/PROGRESS_S{n}.md` + `E:\project\CLAUDE.md` 現況段。

---

## 各階段補 spec + 實作要點

### S5 — 尺寸優化(FSD,全應力設計)
- **目標 / 不做**:**做** `runSizeOptimization` 全應力 stress-ratio resizing `A ← max(A_min, A·σ/σ_a)`,多工況 `envelope` D/C 驅動,迭代重解吃 `ReSolveSession`(同 topology、只改截面值 → ⚠️ **改截面值會變 K_ff,超出 ReSolve same-topology 假設**,須走 rebaseline 或每次 `assembleAndFactor`;補 spec 時釐清:ReSolve 只省「移除/恢復」不省「改剛度值」,FSD 改 A 是改值 → 多半 fresh re-factor,效能段誠實標)。**不做** 真最佳化(KKT/MMA/SQP)、離散標準截面庫匹配(初版連續 A;離散表 opt-in 後處理 + 震盪守門)、LTB(`Section` 缺 warping/Cw,gap 清單見 WS_H)。
- **API 草**:`runSizeOptimization(model, combos, SizeOptOptions{maxIter, dcTol, Amin, sectionTable?}) → SizeOptResult{optimizedSections, dcHistory, converged, weightHistory}`。
- **演算法**:每迭代 = 對所有工況 solve(reuse factor 或 fresh)→ `worstUtilization` 取每桿跨工況最大 D/C → `A ← A·maxDC`(夾 `Amin`)→ 重組 → 收斂判據 `max|ΔA|/A < dcTol`(全應力)。離散表:就近匹配 + 狀態哈希震盪守門(**同 S4 cycle guard 模式**)。
- **Oracle**:10-bar truss(經典 Schmit/Berke)收斂到 **1593.16 lb**(`Research/WS_H_sizeopt/exp_size_opt.cpp` 已驗 24 迭代);全應力偏差 ≤1e-8;離散表震盪 → 守門有限終止。OpenSees 對照可選(WS_H 列 elasticBeamColumn 重設 A 的腳本)。
- **誠實邊界**:靜不定結構全應力 ≠ 全域最優(WS_H 文獻);無 LTB/局部挫屈;連續 A(離散匹配是 opt-in 後處理)。`[NEW CODE]`=stress-ratio 迭代 + 離散震盪守門。依據 `docs/research/WS_H_sizeopt.md`、`WS_R2_experiments.md §7`。

### S6 — Grasshopper 對接 MVP(CLI 橋 → daemon → C API)
- **目標 / 不做**:三步走 **J1 `frame_cli` 文字橋(MVP)→ J1.5 daemon 模式 CLI(常駐行程、stdin 多請求協議、持有 `PreparedSystem`/`ReSolveSession` 跨請求)→ J2 C API DLL(長期)**。**不做** 一步到位 C API(WS_R2 §10:行程開銷僅 6.7ms,~20k DOF 端到端 2.1s 可用 → J1 先夠;100k+ 受 factor 支配才催 daemon)。
- **內容**:`frame_cli` 補齊 `PDELTA`(✅ S3 已加)/`TONLY`(`MEMBER … tonly` token,S4 排到此)/`DYNC` 指令;GH C# 元件組(Assemble/Analyze/Buckling/PDelta/TensionOnly/DynCollapse 讀取器);`.gha` 建置(Rhino 8 .NET、Yak 發佈,WS_J 工具鏈);版本握手行(CLI 輸出 build SHA,沿用 git banner)。
- **Oracle / gate**:CLI 端到端整合測試(輸入→輸出黃金檔比對);**GH 元件不入四腿 gate**(獨立 repo/CI 議題,屬並行線生態)。daemon 模式加「同一連線多次 solve 結果 == 各自獨立 cli」一致性測。
- **誠實邊界**:GH 層在引擎之上(POD 邊界);daemon 狀態生命週期管理(stale model 偵測沿用 `modelFingerprint`)。依據 `WS_J_gh.md`、`WS_K_ecosystem.md`、`WS_R2 §10`。

### S7 — 拓撲優化 BESO + N2 倒塌韌性
- **目標 / 不做**:**做** hard-kill BESO(在 `Member.active`/`ShellQuad.active` 上零侵入,複用 S1 移除設施)+ **必修三項**:① 敏感度歷史平均(Huang & Xie,抑振盪)② **compliance 跳變停機準則**(⚠️ `Research/WS_I_beso/exp_beso_truss.cpp` 實測尾段 compliance **暴增 52×** → 沒停機會劣化)③ 機構守門用 **N1 capacitance 偵測**(`ReSolveSession` 移除後 capacitance 奇異 = 機構,零額外成本)。**N2 倒塌韌性約束版**:每 k 步以 `runProgressiveCollapse` 抽查關鍵桿移除情境,違反韌性 → 回滾 + 鎖該桿保護。**不做** soft-kill SIMP(連續密度,屬另一族)、應力約束 BESO(初版剛度/compliance 目標)。
- **API 草**:`runBESO(model, loads, BESOOptions{targetVolFrac, evolRate, filterRadius?, sensHistory, complianceJumpTol, redundancyCheckEvery?}) → BESOResult{activeHistory, complianceHistory, finalActive, stoppedReason}`。
- **Oracle**:Michell 桁架定性形態 + compliance **單調下降段**(停機前);N2:小框架「**無約束最優解 vs 韌性約束解**」可構造對比(約束解必通過倒塌抽查、無約束解必不通過)。
- **誠實邊界 / novelty**:LSP 級韌性評估(非真非線性倒塌);**勿自稱 fail-safe**——定位為「fail-safe TO 親緣,但離散框架 × 序列線性倒塌評估器 × Stable/Collapsed 終點的方法整合」(WS_N §N2)。`[NEW CODE]`=敏感度歷史平均 + compliance 跳變停機 + N2 抽查回滾。依據 `WS_I_beso.md`、`WS_N_priorart.md`、`WS_R2 §8`。

### S8 — 殼路線(8a QM6 opt-in 膜 → 8b DKQ 薄板快路)
- **目標 / 不做**:**做** QM6 不協調膜(opt-in,降 MITC4 平面 facet 的膜鎖定 = 曲面慢收斂主因,既有 P5 探索結論)+ DKQ 薄板快速路徑(非取代 MITC4,薄板專用)。**⚠️ 最大風險 / 必守**:**OpenSees `ShellMITC4` 是原始雙線性膜,改進膜必須 opt-in、預設保持 Q4**,否則破壞 OpenSees 殼閘門(現 ~1e-10 平板)。**不做** 完整曲殼幾何(facet 路線);殼幾何勁度(殼挫屈)另列。
- **Oracle**:膜 patch test 機器精度;QM6 vs Q4 在 Cook's membrane / 曲面 benchmark 的收斂對照(QM6 應更快);OpenSees 對照**僅在 opt-in 關閉(Q4)時**保持既有閘門。
- **誠實邊界**:facet 近似(曲面誤差隨網格收斂);QM6 不協調膜的 patch test 通過性需驗;DKQ 限薄板(厚板回 MITC4)。依據 `WS_E_shell.md` + 既有膜鎖定探索。

### S9 — Co-rotational 大位移(CR)
- **目標 / 不做**:**做** 獨立 Newton–Raphson 大位移 driver(Battini 2002 / Crisfield 路線,CR 分解剛體轉動 + 局部小應變),load stepping。**不做** snap-through(需弧長法,初版 NR 誠實標限制)、動力大位移、接觸。
- **架構**:獨立 NR driver(仿 `Collapse.cpp` 工作副本模式);`IElement` 加「變形構型」hook(目前 element 幾何凍結於初始;CR 需每迭代更新 element 局部框架)——**先做侵入性評估**(改 `IElement` 介面影響所有 element,補 spec 時列影響面)。
- **Oracle**:**elastica shooting 表已備**(`WS_R2 §9`:懸臂端點 α=PL²/EI=1 → δv/L=**0.3017207738**、δh/L 等,雙容差十位一致);OpenSees `corotCrdTransf` 對照。
- **誠實邊界**:NR 無 snap-through(弧長留後續);CR 小應變大轉動(非大應變);與塑鉸方向耦合(R4)→ **S10 必在 S9 後**。依據 `WS_F_corot.md`、`WS_R2 §9`。

### S10 — 材料非線性(N-M 互動塑鉸;選做,S9 後)
- **目標 / 不做**:**做** G1 路線 `Mp_eff = Mp · f(N/Np)`(AISC/EC3 軸力-彎矩互動式),銜接既有 event-to-event 塑鉸(觸發改 `|M| ≥ Mp_eff(N)`)。**不做** G2 集中彈簧塑性 / G3 纖維斷面(理由與成本見 WS_G);無卸載/反轉(誠實:仍 sequential linear,非真彈塑性)。
- **Oracle**:N-M 互動面解析點(短柱 N-M 包絡);OpenSees 對照(WS_G 列選項)。
- **誠實邊界**:單調載重(無卸載)、互動式近似(非纖維)、與 S9 CR 方向耦合需驗。依據 `WS_G_matnl.md`。

### S11 — MITC9i 高階殼(殿後,使用者指定)
- **目標**:Wisniewski & Turska 2018 的 9 節點插值協變 MITC9i。**先決 = 9 處引擎修改清單**(WS_E):`ShellQuad9` 型別 / 容器 / `SolveResult` 角點 4→9 / `modelFingerprint` / `validate` / `Connectivity` 邊 / element dispatch / build scripts / cli。動工前等 **S1 效能線 + S8 殼線穩定**。
- **Oracle**:patch(含扭曲網格)、vs MITC4 收斂階、文獻數字對表。依據 `WS_E_shell.md`。

---

## 並行線(不佔主線 F 編號;穿插或里程碑後)
- **C6–C8 可視化資料**(S5 後短衝刺,同餵 GH/UE):沿桿 BMD/SFD(Hermite 內插)/利用率場(複用 D/C)/贅餘度報告。POD 結果,獨立 oracle(內插值 vs 端點解析)。
- **生態系**(S6 後):Twinmotion 路徑(WS_K 三路線)、GitHub Actions CI(build+gate)、Yak/zip 自動打包(git SHA banner)。
- **UE 視覺層 U1–U21**(S6 後):S1+S2 使「玩家拆桿即時重算(ReSolve)+ 真動力倒塌回放(`DynCollapseFrame.u/v` 幀 + `FragmentCluster`→Chaos 帶初速碎塊)」成立;D/C 熱圖、振型、即時動態搖晃。**材料非線性渲染屬此層**(誠實標 stress visualization,非 material nonlinearity)。

## 研究定位項(不入主線;觸發條件見 `KARAMBA3D_ROADMAP.md §判定表`)
- AMG/SA-AMG、matrix-free EBE Krylov、全系統隱式暫態、Lanczos shift-invert(屈曲 -Kg indefinite):**R13 觸發**(>50 萬 DOF 互動需求成真才立項;研究輪實測 Eigen IC-PCG 在框架上**反而劣於 Jacobi**,通用預條件路線成本不可低估)。

## 一鍵指令
- **完整四腿 gate(commit 前必跑、必須全綠)**:`powershell -ExecutionPolicy Bypass -File E:\project\ArchSim\Scripts\run_gate.ps1 -RequireOpenSees`
- 快速 standalone:`Plugins\FrameSolver\Standalone\build.bat`(期望 `ALL PASS (failures=0)`)
- UE 模組 build(**PowerShell 前景**):`cmd /c "E:\project\UE_5.7\Engine\Build\BatchFiles\Build.bat ArchSimEditor Win64 Development -project=""E:\project\ArchSim\ArchSim.uproject"" -waitmutex"`
- 完整脈絡:`docs/KARAMBA3D_ROADMAP.md`、`docs/PROGRESS_S1`…`S4.md`、`docs/specs/S*.md`、`docs/research/WS_*.md`、`E:\project\CLAUDE.md`
