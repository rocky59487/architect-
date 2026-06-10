# FrameCore 實作計畫(Implementation Plan)

> **這份是什麼**:把 `81639c1` 研究輪的全部產出(對標報告 + 獨創算法路線 + S1–S11 主線 + 實驗數據)
> 收斂成**一份可直接照著動工的開發計畫**。摘要級——每階段只給「目標 / 依賴 / 產出 / oracle+gate
> 摘要 / 驗收門檻 / 風險」;逐位元的 API 簽名、數值細節、檔案清單在各 spec,本檔以連結指回。
>
> **來源**:`docs/KARAMBA3D_ROADMAP.md`(主報告)、`docs/PERFORMANCE_BASELINE.md`(速度基線)、
> `docs/specs/S1..S4`(交接級 spec,逐位元就緒)、`docs/specs/S5_S11_skeletons.md`(骨架)、
> `docs/research/WS_*.md`(文獻查證)、`docs/research/WS_R2_experiments.md`(可重跑實驗數據)。
> **基準 commit**:`733833d`(研究輪未改引擎,只新增 docs/ 與 Research/ scratch)。

---

## 0. 怎麼用這份計畫

1. **讀順序**:本檔 §1(現況)→ §2(開發方法,每階段照辦)→ §3(依賴/排序)→ §4(逐階段)。
   要動工某階段時,先讀本檔該節拿到目標與驗收,再開對應 spec 拿逐位元細節。
2. **宣稱分級**(全檔沿用 roadmap §0 紀律):`[VERIFIED]` = 本輪 scratch 實驗 + oracle 過、可重跑;
   `[NEW CODE]` = spec 中無原型背書、實作時須新寫並由指定 oracle 把關;`[THEORY]` = 有推導(含外推);
   `[PENDING:…]` = 待實測。**「優於 Karamba3D / 新穎」禁止裸宣稱**,novelty 一律附先行技術定位(WS_N)。
3. **鐵則(不可違反,沿用 `CLAUDE.md`)**:
   - FrameCore 維持**純 C++17/Eigen**,公開 API 只用 POD/std,零 UE、零 Eigen 洩漏。
   - **每階段結束**:獨立 oracle + **四腿 gate 全綠**(見 §2)+ commit/push + 更新 memory 與本檔狀態欄。
   - 建模用 **index 不用裸指標**(`matIdx`/`secIdx`);`Research/` scratch **不入 gate**。
4. **本計畫不改變既有誠實邊界**(LSP ±30% 級、塑鉸非真彈塑性、D/C 非規範檢核、殼 facet…),
   只在每階段新增能力時隨之宣告新邊界(§9 彙整)。

---

## 1. 現況快照(基準 `733833d`)

**已有(主線起點)**:純 C++17/Eigen、零依賴、POD 公開 API;**靜力解早已稀疏**
(SimplicialLDLT/AMD)——舊 roadmap「換稀疏求解器」前提不成立,真正稠密殘留 = 屈曲 GEVP 與 modal
預設路徑;factorize-once `PreparedSystem` + fingerprint 防呆;MITC4 四邊形殼(OpenSees ~1e-10);
8 階段線性套件(組合/包絡/影響線/沉陷/模態/屈曲/反應譜/模態疊加動態);**漸進倒塌主線全套**
(移除/連通/碎塊 Chaos 交接/LSP 驅動器/殼 vM/event-to-event 塑鉸);四腿 gate(F1–F33、UE 34、audit 62、
OpenSees strict)。

**缺口(主線要補)**:幾何非線性、截面優化、生態系(零 DLL/序列化)、`build_perf.bat` 連結失敗
(R8 `[VERIFIED]`,S1 修)、動力倒塌、互動增量重解。

**主線一覽**(依序;✅=spec 就緒 🔶=骨架):

| 階段 | 一句話 | Spec | 獨創定位 |
|---|---|---|---|
| **S1** | ReSolve 三層重解階梯 + 稀疏屈曲 + 修 build_perf + 效能基線正式化 | ✅ [S1](specs/S1_resolve_ladder.md) | **N1 互動增量重解(主打)** |
| **S2** | N4 動量繼承連續動力倒塌(Ritz 基底 + 碎塊帶初速 Chaos 交接) | ✅ [S2](specs/S2_dynamic_collapse.md) | **N4 跨事件繼承(主打)** |
| **S3** | P-Delta 二階(N3 凍結分解 + Wilson 參考,雙路徑互鎖) | ✅ [S3](specs/S3_pdelta.md) | N3 架構整合 |
| **S4** | Tension-only 桿件(ReSolve 內迴圈) | ✅ [S4](specs/S4_tension_only.md) | — |
| **S5** | FSD 尺寸優化(stress-ratio;初版無 LTB,誠實標) | 🔶 [骨架](specs/S5_S11_skeletons.md) | — |
| **S6** | GH 對接 MVP:CLI 橋 → daemon 模式 → C API | 🔶 骨架 | — |
| **S7** | BESO 拓撲優化 + **N2 倒塌韌性約束版** | 🔶 骨架 | N2 方法整合 |
| **S8** | 殼:QM6 opt-in 膜 → DKQ 薄板快路 | 🔶 骨架 | — |
| **S9** | Co-rotational 幾何非線性(elastica oracle 已備) | 🔶 骨架 | — |
| **S10** | N-M 互動塑鉸(選做) | 🔶 骨架 | — |
| **S11** | MITC9i 高階殼(殿後) | 🔶 骨架 | — |

**並行線(不佔主線編號)**:C6–C8 可視化資料(S3–S4 間)、生態系(S6 後)、UE 視覺層 U1–U21(S6 後)。
**研究定位項(不入主線)**:AMG/SA-AMG、matrix-free EBE Krylov、全系統隱式暫態、Lanczos shift-invert
(觸發條件 = R13,§6)。

---

## 2. 開發方法(導入階梯——每階段一律照辦)

**四步導入階梯**(「不錯誤的前提下逐漸採用」的機制,roadmap §8):

1. **scratch 原型 + oracle**(本輪已完成 N1/N3/N4 + 屈曲/TO/FSD/BESO;見 [WS_R2](research/WS_R2_experiments.md))。
2. **進引擎**:opt-in 旗標 + **雙路徑互驗**(新舊算法同跑,audit check 鎖一致;Tier-2 類迭代解用容差互驗、文檔明標非逐位元)。
3. **gate 固化**:新增 F 編號 / UE 測試 / audit checks(各 spec ⑦ 節)。
4. **數據證明**:效能驗收(各 spec ⑧ 節)通過後**轉預設**,舊路徑降級 fallback + 永久 oracle。
   任何數值疑慮 → fallback 永遠正確(rebaseline / fresh factor / dense GEVP / 全基底)。

**四腿 gate**(每階段結束全綠才 commit;一鍵 `Scripts\run_gate.ps1 -RequireOpenSees`):

| 腿 | 內容 | 跑法 |
|---|---|---|
| ① standalone | FrameCore + F1..Fxx 解析/benchmark oracle | `Plugins\FrameSolver\Standalone\build.bat` → `ALL PASS (failures=0)` |
| ② UE automation | `FrameCore.*` headless 測試(`$ExpectedUeTests`) | `UnrealEditor-Cmd … Automation RunTests FrameCore`(讀 `Saved\Logs\ArchSim.log`) |
| ③ OpenSees strict | beam/shell/settlement/移除態/塑鉸態 逐位移比對 | `Tools\opensees_compare.py`(openseespy) |
| ④ linear_deep_audit | C++ standalone 深度審核(現 62 checks) | `Standalone\build_linear_audit.bat` |

**每加一個檔/測試的同步義務**(漏一個 = gate 假綠):
- 每加一個 `Private/*.cpp` → 補進 `build.bat` / `build_linear_audit.bat`(顯式源檔清單;`build_cli`/`build_perf` 視需要)。
- 每加一個 UE 測試 → bump `run_gate.ps1` 的 `$ExpectedUeTests`。
- 改 fingerprint 守門的旗標(如 S4 `tensionOnly`)→ **fingerprint 與旗標同 commit**,否則靜默 stale。

**效能驗收的機器綁定**(各階段 ⑧ 節通則):絕對毫秒出自研究輪單機(cl /O2 /MD、單緒),**驗收同機**
先重跑對應 `Research\bin\*.exe` 取當日基線,**倍率對倍率**比,不跨機比絕對時間;退步 >30% 視為驗收失敗
並回填 [PERFORMANCE_BASELINE.md](PERFORMANCE_BASELINE.md)。

---

## 3. 依賴圖 / 排序理由

```
S1 (ReSolve 地基) ──┬──► S2 (動力倒塌:事件重解吃 ReSolve)
                    ├──► S4 (TO 內迴圈 = ReSolve rank-6/翻轉)
                    └──► S3 (P-Delta 凍結路徑復用 ldlt0;與 N1 零摩擦)
S2 ─────────────────────► UE 視覺層(玩家拆桿即時重算 + 真動力倒塌回放)
S5/S7 ──(迭代重解)──► 吃 ReSolve
S6 (GH 橋) ───► 生態系並行線
S9 (CR) ───► S10 (N-M 塑鉸,必在 CR 後;R4 方向耦合)
```

- **S1 最先**:是所有迭代驅動器(S2/S4/S5/S7)的重解地基;依賴鏈源頭。
- **S2 次之**:獨創主打且資產齊備(FragmentCluster/模態疊加既有)——價值前置。
- **S3/S4**:小而 oracle 嚴,快速積累 gate。
- **S6 提前**:讓外界(GH)可用,催化生態。
- **S11 殿後**:使用者定調 + 9 處引擎修改成本高。
- **S10 必在 S9 後**:CR × 塑鉸方向耦合(風險 R4)。

---

## 4. 階段計畫(逐階段摘要)

> 格式:**目標 / 明確不做** → **依賴** → **產出(新檔·關鍵 API)** → **Oracle + Gate 成長** →
> **效能驗收** → **`[NEW CODE]` 風險點** → **spec 連結**。逐位元細節一律在 spec。

### S1 — ReSolve 三層重分析階梯 + 稀疏屈曲 ✅

- **目標**:`ReSolveSession`——桿/殼停用·恢復後的快速重解,三層自動階梯
  (Tier-1 Woodbury 精確低秩 → Tier-2 stale-LDLT 預條件 CG → Tier-3 全重分解;capacitance 奇異 = 機構偵測);
  `solveBuckling` 稀疏路徑(復用 `SparseEigsolver.h`,dense 留 fallback/oracle);**修 `build_perf.bat`**
  (R8:源檔清單補 `MITC4ShellElement.cpp`,獨立 bugfix 分開 commit);`PERFORMANCE_BASELINE.md` 正式化。
- **不做**:材料/截面值變更、支承 flag 變更走 Tier-3(fingerprint 已防呆);多執行緒;外部後端。
- **依賴**:無(主線地基)。
- **產出**:新 `Public/FrameCore/Reanalysis.h` + `Private/Reanalysis.cpp`(`ReSolveSession`、`ReanalysisOptions/Stats`);
  改 `BucklingAnalysis.cpp`(sparse overload,舊簽名不動)。
- **Oracle + Gate**:F34(tier-1,含 **UDL 桿等效載增量子案例**驗 `[NEW CODE]`)、F35(tier-2 容差)、
  F36(sparse buckling vs dense/Euler);UE 34→**37**;audit 62→**68**(+6:tier1 一致/恢復漂移/機構/tier2 容差/buckling sparse/全拉護衛)。
- **效能驗收**:單桿移除 ≥20×(研究輪 31×)、50 桿序列尾端 ≥12×(17.5×)、tier-2 160 桿 ≥5×(8.2×)、
  tower-M sparse buckling ≥10×(24.6×)。
- **`[NEW CODE]` 風險**:**F 增量記帳**(停用/恢復帶 UDL 桿/帶壓力殼時等效節點載同步離開·回到 F)——
  原型只測過節點載重,風險高於其餘 Tier-1 機制,**F34 的 UDL 子案例是唯一安全網**(漏測 = 帶 UDL 模型靜默錯解)。
  殼 rank-18 未在原型驗過 → 先加殼版 F34 子案例,不過則殼一律走 Tier-2/3。
- **spec**:[S1_resolve_ladder.md](specs/S1_resolve_ladder.md) ‖ 原型 `../Research/WS_N_incremental/exp_incremental_refactor.cpp`、`../Research/WS_B_solver/exp_sparse_buckling.cpp`

### S2 — N4 動量繼承連續動力倒塌 ✅

- **目標**:`runDynamicCollapse`——模態空間時間積分(Newmark β=1/4)、事件觸發元素停用、
  **跨事件狀態繼承**、**碎塊帶初速交接 Chaos**(`FragmentCluster` 加 `vel/angVel`);load-dependent Ritz 基底;模態 warm-start。
  把既有靜力 LSP 倒塌升級為連續動力,修掉「碎塊零初速」既有限制。
- **不做**:材料/幾何非線性動力、阻尼校準(吃使用者輸入)、Chaos 接手後雙向耦合(單向交接)。
- **依賴**:**S1**(事件重解吃 ReSolveSession;若 S2 先行,fresh assembleAndFactor 路徑亦可,介面留 hook)。
- **產出**:新 `Public/FrameCore/DynamicCollapse.h` + `Private/DynamicCollapse.cpp`;
  改 `Connectivity.h`(`FragmentCluster` + `Vec3 vel/angVel`,預設零 = 向後相容)、`SparseEigsolver.h`(`subspaceSmallest` 加 `X0` warm-start 末參,預設 nullptr 逐位元不變)。
- **Oracle + Gate**:F37(全基底跨事件等價 vs 全系統 Newmark,**含 detach 子案例**)、F38(動量帳 vs 閉式)、
  F39(雙質點解析);UE **37→40**(S1 後)或 34→37(S2 先行);audit **68→72**(+4)。
- **效能驗收**:事件間每步 O(basisSize)(目標 nf=20k、basis=30、每步 ≤0.5ms → 60fps 內 ≥30 步)`[PENDING]`;
  事件成本 ≤ 1.5× 一次 ReSolve。
- **`[NEW CODE]` 風險**:**load-dependent Ritz 生成**(Wilson 1985)——原型只比較了截斷誤差(純模態 m=40/108 仍 7.3%),
  生成代碼未原型化,須照 spec ④ 偽代碼新寫,**F37(b/c)把關**;truncationResidual 哨兵 > 門檻 → 自動擴基重投影,最終 fallback 全模態。
- **spec**:[S2_dynamic_collapse.md](specs/S2_dynamic_collapse.md) ‖ 原型 `../Research/WS_N_incremental/exp_dynamic_inherit.cpp`

### S3 — P-Delta 二階分析(雙路徑互鎖) ✅

- **目標**:`runPDelta` 兩條路徑——(A) 參考 = K_T = K_e + Kg 重組單次分解(Wilson 1987);
  (B) 凍結 = 重用 K_e 既有 LDLT 的 pseudo-load 迭代(預設,與 factorize-once/N1 零摩擦);雙路徑 audit 互鎖。
- **不做**:大位移(→S9 CR)、Kg 隨迭代更新軸力(N 凍結 = Th.II 慣例)、殼 Kg、動力 P-Delta。
- **依賴**:S1(凍結路徑復用 `ldlt0`;非硬依賴,可獨立)。
- **產出**:新 `Public/FrameCore/PDeltaAnalysis.h` + `Private/PDeltaAnalysis.cpp`;新 `Tools/pdelta_compare.py`。
- **Oracle + Gate**:F40(懸臂柱 f∈{0,0.3,0.95}:凍結 vs 參考 ≤1e-10、參考 vs 精確 ≤1e-3、P=0 逐位元=線性、
  f=1.05 兩路徑皆 `diverged`)、F41(OpenSees PDeltaCrdTransf);UE 40→**41**;audit 72→**75**(+3)。
- **效能驗收**:凍結路徑 f=0.5 迭代 ≤60(實測 46)且總時 < 參考路徑重分解(中型以上)。
- **`[NEW CODE]` 風險**:**保護式外推**(原型是近裸 Aitken,f=0.9 劣化 285→4742 = 本設計動機)——
  須照 spec ④ 加穩定窗 + 撤銷機制,F40 重測回填迭代數,**驗收底線:任何 f 不得劣於無外推版 1.2×**;
  **發散偵測器**(原型慢發散 f=1.05 未觸發)——須加 20 步滑動窗趨勢 + maxIter 雙保險,F40 旗標把關。
- **spec**:[S3_pdelta.md](specs/S3_pdelta.md) ‖ 原型 `../Research/WS_C_pdelta/exp_pdelta_convergence.cpp`

### S4 — Tension-only 桿件 ✅

- **目標**:`Member.tensionOnly` 旗標 + `runTensionOnly` 迭代驅動器(停用受壓 TO 桿、伸長判據再啟用、
  循環守門 + 單調 fallback 保證有限終止);內迴圈走 ReSolveSession(每翻轉 = rank-6 更新)。
- **不做**:compression-only(留欄位語意不實作)、初始鬆弛/預張力、纜索垂度非線性。
- **依賴**:S1(內迴圈 ReSolve)。
- **產出**:改 `Member.h`(+`tensionOnly`,**入 fingerprint**)、`FrameModel.cpp validate`、`FrameSolver.cpp modelFingerprint`
  (旗標與 fingerprint **同 commit**);新 `Public/FrameCore/TensionOnly.h` + `Private/TensionOnly.cpp`。
- **Oracle + Gate**:F42(X 斜撐:收斂 ≤3 迭代、**收斂解 vs 省略壓桿模型逐位元相等 ≤1e-15**、fingerprint 守門)、
  F43(循環守門有限終止 + 單調 fallback);UE 41→**42**;audit 75→**77**(+2)。
- **效能驗收**:n 根 TO 桿收斂 ≤ n+2 次重解;每次重解走 ReSolve tier-1 ≥10× vs fresh(中型以上)。
- **誠實負結果**:研究輪掃描**無 flip-flop 視窗**(物理 flip-flop 案例 `[PENDING:尋獲後升級]`),
  守門仍留(LCP 視角循環理論可能);收斂解可能順序相依(cycled 時偏保守,diagnostic 明標)。
- **spec**:[S4_tension_only.md](specs/S4_tension_only.md) ‖ 原型 `../Research/WS_D_tensiononly/exp_tension_only.cpp`

### S5 — FSD 尺寸優化 🔶(動工前先補滿 S1–S4 同款十項 spec)

- **目標**:stress-ratio resizing(`A←max(Amin, A·σ/σa)`)+ 多工況 envelope D/C;迭代重解吃 ReSolve;
  複用 `worstUtilization`/`envelope`。API 草:`runSizeOptimization(model, combos, SizeOptOptions) → SizeOptResult`。
- **依賴**:S1(迭代重解)。
- **Oracle**:10-bar truss vs 文獻經典解 1593.16 lb(`[VERIFIED]` 研究輪 24 迭代收斂);全應力偏差 ≤1e-8;
  離散截面表震盪守門(狀態哈希,同 S4 模式)。
- **誠實邊界**:靜不定 FSD **非最優保證**;**初版無 LTB**(Section 缺 warping 等參數;Karamba OptiCroSec 有 LTB → 對標矩陣誠實標差距)。
- **依據**:[S5–S11 骨架](specs/S5_S11_skeletons.md)、[WS_H](research/WS_H_sizeopt.md)、[WS_R2 §7](research/WS_R2_experiments.md)。

### S6 — GH 對接 MVP(CLI 橋 → daemon → C API)🔶

- **目標**:三步走——**J1** `frame_cli` 文字橋為 MVP(加 `PDELTA`/`TONLY`/`DYNC` 指令、`MEMBER … tonly` token、build SHA 握手行)→
  **J1.5** daemon 模式 CLI(常駐行程、stdin 協議多請求、持有 `PreparedSystem`/`ReSolveSession`)→ **J2** C API DLL(長期)。
  GH C# 元件組(Assemble/Analyze/Buckling/PDelta/TensionOnly/DynCollapse 讀取器);`.gha` 建置(Rhino 8 .NET、Yak 發佈)。
- **依賴**:S1–S4(暴露的驅動器);daemon 吃 factorize-once/ReSolve 紅利。
- **研究依據**:行程開銷僅 6.7ms、~20k DOF 端到端 2.06s 可用、100k+ 受 factor 支配 → 催生 daemon
  (`[VERIFIED]` [WS_R2 §10](research/WS_R2_experiments.md))。
- **Oracle/gate**:CLI 端到端整合測試(輸入→輸出黃金檔);**GH 元件不入四腿 gate**(獨立 repo/CI 議題)。
- **依據**:[WS_J](research/WS_J_gh.md)、[WS_K](research/WS_K_ecosystem.md)。

### S7 — BESO 拓撲優化 + N2 倒塌韌性 🔶

- **目標**:hard-kill BESO(在 `active` flag 上零侵入)+ **必修**:敏感度歷史平均(Huang & Xie)、
  compliance 跳變停機、機構守門改用 **N1 capacitance 偵測**(同機制零成本);**N2 韌性約束版**——
  每 k 步以 `runProgressiveCollapse` 抽查關鍵桿移除情境,違反 → 回滾 + 鎖保護。
- **依賴**:S1(capacitance 偵測)、既有 `runProgressiveCollapse`。
- **Oracle**:Michell 定性 + compliance 單調段;N2:小框架「無約束解 vs 韌性約束解」差異案例
  (約束解通過倒塌抽查、無約束解不通過)。研究輪 BESO 29 迭代到 vol 30%、終態 88 桿懸臂 Michell 式形態
  (`[VERIFIED]`;**尾段 compliance 暴增 52× → 必加停機準則**)。
- **獨創定位(WS_N §N2)**:概念親緣 fail-safe TO(Jansen 2014,連續介質/SIMP/靜態最壞柔度);
  **我們 = 離散框架 × 序列線性倒塌 × Stable/Collapsed 終點,方法整合 novelty,勿自稱 fail-safe**。
- **依據**:[WS_I](research/WS_I_beso.md)、[WS_N](research/WS_N_priorart.md)、[WS_R2 §8](research/WS_R2_experiments.md)。

### S8 — 殼路線(QM6 opt-in 膜 → DKQ)🔶

- **目標**:**7a** QM6 不協調膜(最便宜膜精度增益,**opt-in,預設 Q4 保 OpenSees 閘門**)→ **7b** DKQ 薄板快路(非取代 MITC4)。
- **Oracle**:patch test 機器精度、與 MITC4 收斂對照、OpenSees(QM6 須 opt-in 因 `ShellMITC4` 是原始雙線性膜)。
- **最大風險**:OpenSees `ShellMITC4` 是原始膜,改進膜**必須 opt-in、預設 Q4** 才不破壞 OpenSees 閘門。
- **依據**:[WS_E](research/WS_E_shell.md) + 既有 P5 探索(膜鎖定是 MITC4 曲面慢收斂主因)。

### S9 — Co-rotational 幾何非線性 🔶

- **目標**:獨立 NR driver(仿 collapse 工作副本模式)、Battini 2002 / Crisfield 路線、NR + load stepping;
  snap-through 需弧長(誠實標 NR 限制)。`IElement` 加變形構型 hook 的侵入性評估先行。
- **Oracle 數據已備**:elastica shooting 表(`[VERIFIED]` [WS_R2 §9](research/WS_R2_experiments.md):α=1→δv/L=0.3017207738、
  α=10→0.8106090249,雙容差十位一致)+ OpenSees corotCrdTransf。
- **依據**:[WS_F](research/WS_F_corot.md);與塑鉸方向耦合(R4)→ **S10 必在 S9 後**。

### S10 — N-M 互動塑鉸(選做)🔶

- **目標**:G1 路線 `Mp_eff = Mp·f(N/Np)`(AISC/EC3 互動式);與既有 event-to-event 塑鉸銜接 = 觸發式 `|M| ≥ Mp_eff(N)`。
  **排除** G2 彈簧 / G3 纖維(理由見 WS_G);無卸載(誠實:仍 sequential linear)。
- **依賴**:S9(CR;R4 方向耦合)。
- **Oracle**:互動面解析點 + OpenSees 對照。
- **依據**:[WS_G](research/WS_G_matnl.md)。

### S11 — MITC9i 高階殼(殿後)🔶

- **目標**:Wisniewski & Turska 2018;先決條件 = **9 處引擎修改**(ShellQuad9/容器/SolveResult 角點 4→9/
  fingerprint/validate/Connectivity 邊/dispatch/build scripts/cli)。等 S1 效能線 + S8 殼線穩定再動。
- **Oracle**:patch(含扭曲網格)、vs MITC4 收斂階、文獻數字對表。
- **依據**:[WS_E](research/WS_E_shell.md)(使用者定調殿後 + 修改成本)。

---

## 5. 並行線(不佔主線編號)

| 並行線 | 時機 | 內容 | 依據 |
|---|---|---|---|
| **C6–C8 可視化資料** | S3–S4 間(短衝刺) | BMD/SFD 沿桿(Hermite)、利用率場(複用 C3 D/C)、贅餘度報告 → 同餵 GH 與 UE | roadmap §7 |
| **生態系** | S6 後 | Twinmotion 路徑(WS_K 三路線)、GitHub Actions CI(build+gate)、Yak/zip 自動打包(git SHA banner) | [WS_K](research/WS_K_ecosystem.md) |
| **UE 視覺層 U1–U21** | S6 後 | 倒塌回放(吃 `DynCollapseFrame.u/v`)、碎塊帶初速 → Chaos、D/C 熱圖;**S1+S2 使「玩家拆桿即時重算 + 真動力倒塌回放」成立** | roadmap §7 |

---

## 6. 研究定位項(不入主線;觸發條件)

| 項 | 判定 | 觸發條件 / 依據 |
|---|---|---|
| AMG / SA-AMG | **研究定位** | 手寫 SA-AMG(6DOF/旋轉/剛體近核)數千行高風險;外部 lib 違零依賴;**實測 Eigen IC-PCG 在框架上劣於 Jacobi**(2,252 vs 1,248 迭代,框架剛度病態)→ 通用預條件非免費午餐 `[VERIFIED]` |
| matrix-free EBE Krylov | 研究定位 | Tier-2 隱式運算子是特例已入主線;全 EBE 觸發條件未到 |
| 全系統隱式暫態 | 研究定位 | N4 模態路線以遠低成本覆蓋遊戲需求 |
| Lanczos shift-invert | 研究定位 | subspace iteration 對少模數已夠(屈曲 24.6× `[VERIFIED]`) |

**R13 觸發條件(統一)**:「**>50 萬自由 DOF 的互動分析**」成為真實需求時立項,屆時重新評估
外部 lib(授權)vs 自寫(數千行)vs matrix-free EBE。在那之前 stale-LDLT 預條件(N1 Tier-2)覆蓋互動重解需求。
**直接法天花板**(`[VERIFIED]` 支撐此定位):186k DOF factor 522s/3.3GB、390k DOF 3,229s/10.4GB、
1M DOF 外推 ~9.3h/~39GB `[THEORY]` → 引擎價值區間 = 互動規模(≤~60k DOF)的重解速度 × 正確性文化。

---

## 7. Gate 成長軌跡(累計)

| 完成後 | standalone F | UE 測試 | linear_deep_audit | 備註 |
|---|---|---|---|---|
| 基準 `733833d` | F1–F33 | 34 | 62 | 起點 |
| **S1** | +F34–F36 | 37 | 68 | tier1/恢復漂移/機構/tier2/sparse buckling/全拉 |
| **S2** | +F37–F39 | 40 | 72 | 等價/動量/Ritz 殘差/能量帳 |
| **S3** | +F40–F41 | 41 | 75 | 雙路徑互驗/P=0 退化/發散旗標 |
| **S4** | +F42–F43 | 42 | 77 | 省略等價/fingerprint 守門 |
| S5+ | +F44… | 43… | 79… | 骨架,動工前定 |

> 數字是 spec 預測(各 spec ⑦ 節);實作時以實際新增為準並回填本表 + `run_gate.ps1`。
> UE 期望數基準取決於 S1/S2 先後(S2 spec ⑤ 已註)。

**實作狀態(2026-06-11,夜間無人監督)**:S1 進行中。已落地並 push:R8 `build_perf.bat` 修復(`0e2e500`)、
**稀疏屈曲** overload(`a91b171`,standalone **F34** + audit **63**;注意實際編號 F34 而非 spec 暫定 F36,
因 baseline 最高 logical fixture=F33)、`PERFORMANCE_BASELINE.md` 正式化 + 同機錨點。**ReSolveSession**
(三層階梯,本表的 F34/F35 預測項)尚未實作,已備妥逐位元交接筆記於 `docs/PROGRESS_S1.md`(建議拆
Tier-1+Tier-3 / Tier-2 兩個完整 commit)。UE automation + OpenSees 重 gate 依務實分層政策留待 S1 完成里程碑。
詳見 [PROGRESS_S1.md](PROGRESS_S1.md)。

---

## 8. 風險登記冊(濃縮自 roadmap §11)

| # | 風險 | 緩解 |
|---|---|---|
| R1 | −Kg 半正定(屈曲 GEVP) | 護衛 + 退 dense(`[VERIFIED]` 退化案例) |
| R2 | 直接法天花板 | 價值主張定位互動規模;百萬 DOF = 批次/離線域(R13) |
| R3 | CLI 吞吐(100k+ 143s) | J1.5 daemon 模式 |
| R4 | CR × 塑鉸方向耦合 | S10 在 S9 後 |
| R5 | 外部後端 vs 零依賴 | CHOLMOD 授權卡死(走零依賴 Woodbury);MKL 留 opt-in 議題 |
| R8 | `build_perf.bat` 連結失敗 | `[VERIFIED]`,**S1 修**(補 `MITC4ShellElement.cpp`) |
| R10 | N1 長序列漂移 | 實測 1.5e-15@R=600;仍設 R>2·maxRank rebaseline 哨兵 |
| R11 | Tier-2 確定性 | 同輸入同機可重現、與直接解差 tol 級(非逐位元) |
| R12 | N4 截斷誤差 | Ritz + 殘差哨兵 + 全基底 fallback |
| R13 | AMG 觸發 | >50 萬 DOF 互動需求成真才立項(§6) |

---

## 9. 誠實邊界(新能力一律隨之宣告)

既有邊界全部沿用(LSP ±30% 級、塑鉸非真彈塑性、D/C 非規範檢核、殼 facet/膜鎖定…)。主線新增:

- **S1 Tier-2**:容差級非逐位元(同輸入確定性可重現);`pivotMargin` 語意不延伸到 ladder(tier-1 報 capacitance pivot ratio,獨立指標)。
- **S2**:事件間線彈性模態空間;失效準則 = screening 級 D/C;截斷誤差由 `truncationResidual` 顯式報告;
  事件整步觸發(O(dt));Chaos 交接單向。
- **S3**:Th.II 線性化(N 凍結、小側移);f→1 收斂慢;不適用大位移/後挫屈;殼不貢獻 Kg。
- **S4**:軸力符號判據(無預張力/鬆弛長度);cycled 時收斂解順序相依且偏保守(diagnostic 明標)。
- **S5**:靜不定 FSD 非最優保證;初版無 LTB。 **S7**:BESO hard-kill 低體積過衝需停機準則。
- **通則**:1M DOF 數字是外推非實測;同機速度對標未做(`[PENDING:S6 後]`,免費版 ≤20 梁元素只能對微模型)。

**授權策略**:建議 **MIT**(Eigen MPL2 已設 `EIGEN_MPL2_ONLY` 允許靜態連結;CHOLMOD LGPL/GPL 不引入——
N1 走零依賴 Woodbury,劣勢變特色);對標賣點 = 「開源 + 可重跑 oracle 套件」vs 閉源商業。

---

## 10. 入口 / 交接

- **下一輪動工入口** = **S1**:讀本檔 §4 S1 → [specs/S1_resolve_ladder.md](specs/S1_resolve_ladder.md)(十項齊,含可移植原型路徑)→
  [WS_R2 §1–§3](research/WS_R2_experiments.md)。
- **鐵則不變**:四腿 gate 全綠才 commit;`Research/` 不入 gate;每階段結束更新 memory 與本檔狀態欄。
- **build/test 一鍵**(都在 `E:\project\ArchSim`):
  - standalone 閘門:`Plugins\FrameSolver\Standalone\build.bat`(期望 `ALL PASS (failures=0)`)。
  - 完整四腿:`powershell -ExecutionPolicy Bypass -File Scripts\run_gate.ps1 -RequireOpenSees`。
- **完整對標脈絡與判定依據**:[KARAMBA3D_ROADMAP.md](KARAMBA3D_ROADMAP.md)(§5 獨創路線、§6 技術判定、§8 導入階梯)。
