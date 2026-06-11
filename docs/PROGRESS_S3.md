# S3 進度日誌 — P-Delta 二階分析(`runPDelta`)

> 接續 `PROGRESS_S2.md`(S2 結案於 `7051485`)。本階段在 factorize-once 架構上加**線性化
> (Theory-II)P-Delta 二階分析**:幾何勁度 `Kg` 由一階軸力凍結組成,兩條互驗路徑——
> 凍結 pseudo-load 迭代(重用既有 LDLᵀ,零重分解)與 K_T 重組參考——抵達同一不動點。
> 政策(使用者 2026-06-11 定調):**每次 commit 前跑完整四腿 `run_gate.ps1 -RequireOpenSees` 並全綠**。

## 基準
- 起點 `7051485`。working tree 既有雜項(`.gitignore`/`ArchSim.uproject` 改動、`Plugins/LevelSim/`)**非本人改動,全程不碰**。
- baseline 四腿確認綠:standalone F1–F39 / UE 40 / OpenSees PASS / audit 71。
- 使用者授權後動工(S2「完成即停」閘已解除)。

## 交付內容(單一 S3 commit)
### 新增
- `Public/FrameCore/PDeltaAnalysis.h` — POD API(`PDeltaOptions`(maxIter/tolU/accelerate/refactorPath/solve)、`PDeltaResult`(converged/diverged/iterations/lastIncrement/finalState))+ `runPDelta`。零 Eigen。
- `Private/PDeltaAnalysis.cpp` — 雙路徑驅動器。資料流:`assembleAndFactor` → 一階 `solveLoad` → 凍結軸力 `Kg` → `reduceFF`;**固定外力 `F_ff = K_ff·u_lin_ff`**(從一階解反推,一次 sparse mat-vec,兩路徑共用故互驗不受重建捨入影響);**P=0 短路回線性逐位元**;凍結迭代 `u ← K_e⁻¹(F − Kg u)` 重用 `S.ldlt`;參考路徑 fresh LDLᵀ of `K_T`;`recoverState` scatter + `el->recover` 出二階桿端力。
- `Tools/pdelta_compare.py` — 獨立 P-Delta frac 掃描報告(沿用 opensees_compare 的 helper;不入 gate)。
- `Private/Tests/PDeltaTest.cpp` — UE 自動化 `FrameCore.PDelta.AmplificationOracle`。

### 修改
- `Standalone/main.cpp`:F40(+ include PDeltaAnalysis.h)。`FrameTestFixtures.h`:`pdeltaColumn`(垂直懸臂柱+軸壓+橫向載)。
- `linear_deep_audit.cpp`:`testPDelta()`(+3 checks,71→74)。
- `frame_cli.cpp`:`PDELTA path` 命令(輸出 `PDSTATUS conv div iters` + 二階 DISP/MF)。
- `Tools/opensees_compare.py`:P-Delta 區段(`geomTransf PDelta` + Newton vs 雙路徑,gate case @ P/Pcr=0.3)。
- `build.bat`/`build_linear_audit.bat`/`build_cli.bat`:源檔清單 +`PDeltaAnalysis.cpp`。`run_gate.ps1`:`$ExpectedUeTests` 40→41。
- **`FrameCore.Build.cs`:`bUseUnity = false`** — 見下節「UE unity build 修復」。

## 演算法 + [NEW CODE] 安全裝置
- **凍結 pseudo-load 迭代**(預設):固定點 `x_{k+1} = K_e⁻¹(F_ff − Kg_ff x_k)`,起點 `x_0 = u_lin_ff`。收斂率 = P/Pcr(P<Pcr 幾何收斂)。
- **保護式幾何外推**([NEW CODE]):Aitken `ρ = ⟨Δ_k,Δ_{k-1}⟩/⟨Δ_{k-1},Δ_{k-1}⟩`,**穩定窗** 0<ρ<0.95 且 |ρ_k−ρ_{k-1}|<0.2 才外推 `x += ρ/(1−ρ)Δ`;**撤銷**:外推後下一步 ‖Δ‖ 未降則回退裸迭代並停用 accelerate(本次求解內)——故任何 f 不劣於裸迭代。
- **發散偵測器**([NEW CODE]):20 步滑動窗 ‖Δ‖ 升過窗首且超初值 → diverged;maxIter 兜底。P<Pcr 幾何遞減永不誤觸,P>Pcr 遞增觸發。
- **參考路徑發散**:`K_T = K_e+Kg` 的 LDLᵀ `vectorD` 出現非正/近零 → 非正定 → diverged(鏡像 assembleAndFactor 機構偵測)。
- **誠實邊界**:Theory-II 線性化,軸力凍結於一階,小側移;殼不貢獻 Kg(同屈曲限制);非大位移(corotational 留 S9)。

## UE unity build 修復(`bUseUnity = false`)
- 症狀:加 PDeltaAnalysis.cpp/PDeltaTest.cpp 後 UE build 失敗,錯誤**全在既有檔**(`reduceFFsp` 在 DynamicCollapse.cpp/Reanalysis.cpp 重定義、`kPi` C4459 shadow-as-error)。
- 根因:各 analysis `.cpp` 在匿名 namespace 用**同名** TU-local helper(`reduceFF`×5 檔、`kPi`×多檔、`reduceFFsp`×2 檔)。UE unity/jumbo build 把數個 `.cpp` 併成一個 TU → TU-local 變衝突的 namespace-scope 符號;哪些檔同 blob 取決於檔數/順序,故**加新檔會引爆既有潛在衝突**(baseline 靠分組運氣而過)。
- 修法:對此純運算模組關 unity(一行 Build.cs)。根治現在+未來的同類衝突,helper 保持乾淨,build 與檔序無關,純運算模組編譯加速損失可忽略。屬 FrameCore 範圍(非禁區)。

## 數值證據(本輪實測)
| Oracle | 量測 | 門檻 | 實測 |
|---|---|---|---|
| **F40**(懸臂柱 P/Pcr=0.3) | 凍結 vs 參考 tip sway | 1e-10 | **1.62e-13** |
| F40(0.3) | 參考 vs beam-column 精確 | 1e-3 | **2.64e-07** |
| **F40**(P/Pcr=0.95) | 凍結 vs 參考 | 1e-10 | **1.45e-12** |
| F40(0.95) | 參考 vs 精確(8 元素離散) | 1e-3 | **3.72e-05** |
| **F40**(P=0) | 凍結 finalState vs 線性解 | 逐位元 | **bit-identical, 0 iters** |
| **F40**(P/Pcr=1.05) | 參考 K_T 非 PD / 凍結滑動窗 | diverged | **兩路徑 diverged**(凍結 20 步) |
| **audit** 雙路徑互鎖(0.5) | frozen vs reference | 1e-10 | **5.68e-14** |
| audit P=0 退化 | max\|u_pdelta − u_linear\| | 0 | **0.000**(逐位元) |
| audit 過 Pcr | 兩路徑皆 diverged | — | **PASS** |
| **OpenSees**(0.3) | frozen vs reference | 1e-9 | **0.000** |
| OpenSees(0.3) | ours vs beam-column 精確 | 1e-2 | **2.64e-07** |
| OpenSees(0.3) | OpenSees PDelta transf vs 精確 | 1e-2 | **1.37e-03**(線性化,忽略內部 P-δ) |
| OpenSees(0.3) | ours vs OpenSees | 1e-2 | **1.37e-03** |
| **sweep**(frac 0.1–0.9) | ours vs 精確 | 2e-3 | **2.3e-08 → 1.7e-05**(全 frac 極準) |
| sweep | 凍結迭代數 | — | 8 → 31(隨 f 增,如理論) |

## 四腿 gate(commit 前,全綠)
`run_gate.ps1 -RequireOpenSees` → **GATE: PASS** — standalone **F1–F40** / UE **41 tests** exit 0 / OpenSees **PASS**(含 P-Delta) / deep audit **PASS checks=74**。UE dual-build:`PDeltaAnalysis.cpp` 經 `FrameEigen.h` choke point(結構平行 BucklingAnalysis.cpp),零 Eigen 洩漏;`bUseUnity=false` 後無 unity 符號衝突。

## 誠實邊界
- Theory-II 線性化:軸力凍結於一階、小側移假設;f→1 收斂慢(迭代數 ~log tol/log f);不適用大位移/後挫屈(→S9 corotational);殼不貢獻 Kg(既有限制)。
- 我方 consistent 幾何勁度(含內部 P-δ)對解析精確解 ~1e-5;OpenSees `PDelta` geomTransf 是 P-large-Delta 線性化(忽略內部 P-δ、軸力 Newton 更新),~1e-3 級,故 gate case 設 f=0.3 + 寬容差,**緊 oracle 是解析 beam-column(F40),OpenSees 是跨工具一致性參考**。
- `finalState.reactions` 用線性 `K_e·u − F`(二階位移下的節點力殘差),非 K_T 平衡(幾何項是虛擬應力勁化算子非外力);桿端力是線性元素內力 evaluated at 二階位移 = 二階內力。

## 下一步(待使用者授權)
F 編號下一個 = **F42**;audit 從 **74** 起增;UE 從 **41** 起增。S3 完成即停,回報後等使用者檢視再授權下一階段(S4 tension-only ReSolve 翻轉 / 或其他)。
