# S5–S11 階段骨架(研究結論依據;下輪逐階段升級為交接級 spec)

> 完成等級:骨架 = 目標/路線已由研究輪定案、oracle 構想成形;動工前該階段須先補滿
> S1–S4 同款十項規格。研究依據:`docs/research/WS_*.md` + `WS_R2_experiments.md`。

## S5 尺寸優化(FSD)
- 路線:stress-ratio resizing(`A←max(Amin, A·σ/σa)`)+ 多工況 envelope D/C;研究輪 10-bar
  收斂 1593.16 lb = 文獻經典解(WS_R2 §7)。複用 `worstUtilization`/`envelope`;迭代重解吃 ReSolve。
- API 草:`runSizeOptimization(model, combos, SizeOptOptions{maxIter, sectionTable?, dcTol}) → SizeOptResult{optimizedSections, dcHistory}`。
- Oracle:10-bar vs 文獻(WS_H 交叉引用);全應力偏差 ≤1e-8;離散截面表震盪守門(狀態哈希,同 S4 模式)。
- 誠實邊界:靜不定 FSD 非最優保證(WS_H 文獻);**初版無 LTB**(Section 缺 warping 等參數,
  gap 清單見 WS_H;Karamba OptiCroSec 有 LTB → 對標矩陣誠實標差距)。

## S6 GH 對接 MVP(CLI 橋 + daemon 模式)
- 研究結論(WS_J + WS_R2 §10):行程開銷僅 6.7ms;~20k DOF 端到端 2.1s 可用;100k+ 需跨解
  持久狀態 → **J1 CLI 橋為 MVP、J1.5 daemon 模式 CLI(常駐行程、stdin 協議多請求、
  持有 PreparedSystem/ReSolveSession)為第二步、J2 C API DLL 為長期**。
- 內容:GH C# 元件組(Assemble/Analyze/BucklingModes/PDelta/TensionOnly/DynCollapse 讀取器);
  frame_cli 加 `PDELTA`/`TONLY`/`DYNC` 指令與 `MEMBER ... tonly` token;`.gha` 建置(WS_J 工具鏈:
  Rhino 8 .NET、Yak 發佈);版本握手行(CLI 輸出 build SHA)。
- Oracle/gate:CLI 層端到端整合測試腳本(輸入→輸出黃金檔);GH 元件不入四腿 gate(獨立 repo/CI 議題)。

## S7 拓撲優化 BESO + N2 倒塌韌性
- 研究結論(WS_R2 §8):hard-kill BESO 在 active flag 上零侵入成立;**必修**:敏感度歷史平均
  (Huang & Xie)、compliance 跳變停機、機構守門改用 N1 capacitance 偵測(同機制零成本)。
- N2:倒塌韌性約束版(每 k 步以 `runProgressiveCollapse` 抽查關鍵桿移除情境,違反 → 回滾+鎖保護)。
  先行技術定位與措辭已定(WS_N §N2:類 fail-safe TO 但離散框架+LSP 評估器,勿用 fail-safe 一詞自稱)。
- Oracle:Michell 定性 + compliance 單調段;N2:小框架「無約束解 vs 韌性約束解」差異案例
  (約束解必須通過倒塌抽查、無約束解必須不通過 — 可構造)。

## S8 殼路線(7a QM6 opt-in 膜 → 7b DKQ)
- 依據 WS_E + 既有 P5 探索:QM6 = 最便宜膜精度增益(opt-in,預設 Q4 保 OpenSees 閘門);
  DKQ = 薄板快速路徑非取代 MITC4。Karamba 殼型式查證結論寫進對標矩陣。
- Oracle:patch test 機器精度、與 MITC4 收斂對照、OpenSees(QM6 須 opt-in 因 ShellMITC4 是原始膜)。

## S9 Co-rotational(CR)
- 依據 WS_F:Battini 2002 / Crisfield 路線;NR + load stepping;snap-through 需弧長(誠實標 NR 限制)。
- Oracle 數據已備:elastica shooting 表(WS_R2 §9,α=1→δv/L=0.3017207738 等,雙容差十位一致)
  + OpenSees corotCrdTransf;與塑鉸方向耦合(R4)→ S10 必在 S9 後。
- 架構:獨立 NR driver(仿 collapse 工作副本模式);`IElement` 加變形構型 hook 的侵入性評估先行。

## S10 材料非線性(N-M 互動塑鉸;選做)
- 依據 WS_G:建議 G1(`Mp_eff = Mp·f(N/Np)`,AISC/EC3 互動式)、排除 G2 彈簧/G3 纖維(理由與
  成本見 WS_G);無卸載(誠實:仍 sequential linear);與既有 event-to-event 塑鉸銜接點 = 觸發式
  `|M| ≥ Mp_eff(N)`。
- Oracle:互動面解析點;OpenSees 對照(WS_G 列選項)。

## S11 MITC9i(殿後,使用者備註指定)
- 依據 WS_E:Wisniewski & Turska 2018;先決條件 = 9 處引擎修改清單(ShellQuad9/容器/SolveResult
  角點 4→9/fingerprint/validate/Connectivity 邊/dispatch/build scripts/cli);等 S1 效能線+S8 殼線穩定。
- Oracle:patch(含扭曲網格)、vs MITC4 收斂階、文獻數字對表。

## 並行線(不佔主線編號)
- **生態系**(S6 後):Twinmotion 路徑(WS_K 三路線結論)、GitHub Actions CI(build+gate)、
  Yak/zip 自動打包(版本注入沿用 git SHA banner)。
- **C6–C8 可視化資料**(S3–S4 間短衝刺):BMD/SFD 沿桿、利用率場、贅餘度 → 同時餵 GH 與 UE。
- **UE 視覺層 U1–U21**(S6 後):S1+S2 使「玩家拆桿即時重算 + 真動力倒塌回放(幀快照+帶初速碎塊)」成立。

## 研究定位項(不入主線;觸發條件)
- AMG/SA-AMG、matrix-free EBE Krylov、全系統隱式暫態、Lanczos shift-invert:
  觸發條件與依據見 `KARAMBA3D_ROADMAP.md` §判定表 + R13(>50 萬 DOF 互動需求成真才立項;
  研究輪數據:Eigen IC-PCG 在框架上反而劣於 Jacobi,通用預條件路線不可低估成本)。
