# S4 進度日誌 — Tension-only 桿件(`runTensionOnly`)

> 接續 `PROGRESS_S3.md`(S3 結案於 `2142e11`)。本階段把 S1 ReSolve 階梯接上第一個真正的
> same-topology 重複客戶:tension-only(纜索 / 細長 X 斜撐)主動集迭代——受壓桿停用、伸長
> 再啟用,**內迴圈每翻轉 = 一次 rank-6 Woodbury 更新**(非 fresh re-factor),收斂解與「把鬆弛
> 桿從 member 清單省略」逐位元相等。政策:每次 commit 前跑完整四腿 `run_gate.ps1 -RequireOpenSees` 全綠。

## 基準
- 起點 `2142e11`(S3)。working tree 既有雜項(`.gitignore`/`ArchSim.uproject`/`Plugins/LevelSim/`)**非本人改動,全程不碰**。
- baseline 四腿綠:standalone F1–F40 / UE 41 / OpenSees PASS / audit 74。
- 使用者「繼續」授權動工 S4。

## 交付內容(單一 S4 commit)
### 新增
- `Public/FrameCore/TensionOnly.h` — POD API(`TensionOnlyOptions`(maxIter/allowReactivation/axialTol/solve)、`TensionOnlyResult`(converged/cycled/iterations/finalState/slack))+ `runTensionOnly`。零 Eigen。
- `Private/TensionOnly.cpp` — 主動集 Jacobi 迭代驅動器,**內迴圈走 `ReSolveSession`**(`Reanalysis.h`,每翻轉 rank-6 Woodbury);`memberElongation`(走 `memberLocalAxes` + 解的節點位移,對 inactive 桿亦有效);FNV-1a 狀態哈希循環守門→單調(deactivate-only)fallback。**純 POD,不碰 Eigen**(ReSolveSession opaque + Vec3)。
- `Private/Tests/TensionOnlyTest.cpp` — UE 自動化 `FrameCore.TensionOnly.Eliminator`。

### 修改
- `Member.h`:`bool tensionOnly = false`(active 之後)。`FrameSolver.cpp` modelFingerprint:active hash 之後插 `tensionOnly`(兩檔同 commit,漏 = 靜默 stale)。
- `Standalone/main.cpp`:F42/F43(+ include)。`FrameTestFixtures.h`:`xBracedPortal`(X 斜撐 portal,兩斜撐 tension-only)。
- `linear_deep_audit.cpp`:`testTensionOnly()`(+2 checks,74→76)。
- `build.bat`/`build_linear_audit.bat`:源檔清單 +`TensionOnly.cpp`(`build_cli.bat` 不加——frame_cli 的 `tonly` token 排 S6)。`run_gate.ps1`:`$ExpectedUeTests` 41→42。

## 演算法 + 誠實定位
- **主動集不動點迭代**(Jacobi 全掃描):解 → 對每根 tensionOnly 桿:active 且 `N(endI) > axialTol`(壓,compression-positive)→ 停用;inactive 且 `(u_j−u_i)·x̂ > 0`(伸長)→ 再啟用(allowReactivation)。無翻轉 → converged。
- **內迴圈 = ReSolve**:`ReSolveSession::setMemberActive` 每翻轉一次 rank-6 Woodbury(Tier-1 exact);`solve()` = fresh assembleAndFactor+solveLoad 到 factorization round-off。**這是 S1 ReSolve 階梯的第一個 same-topology rank-k 客戶**(對照 S2/S3 事件重解走 fresh re-factor)。
- **有限終止**(WS_D 誠實定位:主動集迭代無普遍收斂保證,LCP 視角):FNV-1a 狀態哈希,任何狀態重現(確定性迭代 → 必進循環)→ `cycled`,自動以 `allowReactivation=false` 單調重跑(只停用永不恢復 → ≤ n_TO 步必終止);diagnostic 標順序相依保守性。
- **TO 旗標正交**:plain `assembleAndFactor`/`solveLoad` 忽略 tensionOnly(只 `runTensionOnly` 用),但入 fingerprint 堵 stale reuse。

## 數值證據(本輪實測)
| Oracle | 量測 | 門檻 | 實測 |
|---|---|---|---|
| **F42**(X 斜撐 portal,純橫向) | 收斂迭代數 | ≤3 | **2** |
| F42 | 鬆弛桿數(eliminate 受壓桿) | 1 | **1** |
| F42 | 存活斜撐受拉(N<0,compression-positive) | <0 | **N=−40044.86** |
| F42 | 收斂解 vs 省略壓桿模型(ReSolve round-off) | rel<1e-10 | **0.000000**(逐位元級) |
| F42 | 翻 tensionOnly → solveLoad 拒絕舊 factor | singular | **PASS** |
| **F43**(垂直壓潰) | 單調策略有限終止 | converged | **converged, 2 迭代** |
| F43 | 單調結果滿足 TO 符號約束 | — | **PASS** |
| F43 | 默認策略跨載重掃描全終止(無 hang) | — | **PASS** |
| **audit** 省略等價 | rel max\|u_TO − u_omitted\| | 1e-10 | **1.135e-15**(逐位元級) |
| audit fingerprint 守門 | 翻 tensionOnly → singular | — | **PASS** |

## 四腿 gate(commit 前)
`run_gate.ps1 -RequireOpenSees` → **GATE: PASS** — standalone **F1–F43** / UE **42 tests** / OpenSees **PASS** / deep audit **checks=76**。`bUseUnity=false`(S3 起)+ `TensionOnly.cpp` 純 POD(走 ReSolveSession opaque,零 Eigen),UE dual-build 乾淨。

## 誠實邊界
- 判據為軸力**符號**(無預張力 / 無鬆弛長度 / 無纜索垂度非線性 — 超出線性引擎)。
- `cycled` 解順序相依且偏保守(單調 fallback),diagnostic 明標;WS_D 掃描**未尋獲物理 flip-flop**(誠實負結果),故 F43 以單調有限終止 + 守門契約覆蓋,非物理循環案例。
- 內迴圈走 ReSolve(Woodbury)故收斂解 vs 省略模型 fresh solve 差 factorization round-off(實測達逐位元級 ~1e-15);active=false≡省略的逐位元性質本身由 F26(stage-1)釘死。
- **`validate` 未加 TO 特定檢查**:TO+release 正交合法、validate 為 pass/fail 無 warning 出口、TO+TrussPin 的 both-torsion 釋放已被既有 release 條件數守門擋;TO 旗標已入 fingerprint(核心守門)。
- **compression-only / 預張力 / OpenSees TO 比對**刻意排除(spec ①/⑤;frame_cli `tonly` token + OpenSees 對標排 S6 一起)。

## 下一步(待使用者授權)
F 編號下一個 = **F44**;audit 從 **76** 起增;UE 從 **42** 起增。S4 完成即停,回報後等使用者檢視再授權下一階段(S5 FSD 尺寸優化 / 或其他)。
