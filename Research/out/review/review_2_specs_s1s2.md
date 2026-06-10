現在我有足夠資料進行全面審查。讓我整理所有發現。

# 審查報告:S1/S2 Spec 可實作性(面向 2)

## 裁決

**FAIL** — 發現 6 項 MAJOR 與數項 MINOR 問題。問題集中在:S1 的 F 增量機制是新代碼(非移植)卻以「原型已驗」語氣呈現;S2 的 Ritz 基底生成演算法描述不完整到無法實作;兩份 spec 對 build .bat 的修改要求與真實清單不一致;fingerprint 互動未指定;若干 oracle 門檻的機器判定條件模糊。

---

## 發現(按嚴重度排序)

### S1

- **[MAJOR] S1 §④「F 同步」自述為新代碼但隱藏在「Tier-1 原型已驗」區塊內**
  `S1_resolve_ladder.md` §④ 第 4 段明確寫:「研究原型未測 UDL — 實作時 F 增量為新代碼,oracle F34 必須含一根帶 UDL 的桿」。但 §④ 的標題是「Tier-1(原型已驗,直接移植)」,且 §① 頭部聲明「原型代碼可大段移植」。
  閱讀原型 `exp_incremental_refactor.cpp` 確認:`L.init(ps.impl.get(), nodalLoadVector(model))` 只傳入節點等效載向量,整份原型無任何 `addEquivalentNodalLoads` 或 UDL 貢獻追蹤。亦即 F 增量機制完全不在原型裡。
  **問題**:新來的實作 agent 讀到「原型已驗、直接移植」,極可能跳過對 F 增量的嚴格 oracle 設計、低估這段代碼的風險。
  **建議修法**:在 §④ 中把「F 同步」分離為獨立小節,明確標記 [PENDING] 或 [NEW CODE];在 §① 的「做」清單加「F 增量(UDL 桿等效載追蹤)(新代碼,需獨立 oracle)」;oracle F34 的 UDL 子案例要求語氣從「必須含」改為更強的邊界說明(說明測不到會漏什麼)。

- **[MAJOR] S1 稀疏屈曲的 `subspaceSmallest` 介面呼叫簽名與真實引擎 header 不符**
  S1 §④ 指定:
  ```
  subspaceSmallest(Kff, S.ldlt, −Kg_ff, nev=1)
  ```
  但 `SparseEigsolver.h` 的實際簽名是:
  ```cpp
  bool subspaceSmallest(const SpMat& A, const LDLTSolver& Ainv, const SpMat& B,
                        int nev, VecX& lambda, MatX& vec,
                        int maxIter = 300, real tol = 1e-11)
  ```
  Spec 漏掉了 `VecX& lambda` 與 `MatX& vec` 兩個 out-parameter,且沒說明回傳 `bool`。實作 agent 按 spec 拼 call site 會得編譯錯誤。
  此外 spec 說「`nf > 500`(或新 overload 旗標)時走稀疏路徑」,但真實 `BucklingAnalysis.cpp` 現在做的是稠密路徑(整個 nf×nf dense solve),沒有任何旗標介面。Spec 沒有說明如何在既有 `solveBuckling` 中加入 overload 旗標而不破壞現有 POD 邊界(`BucklingResult` 標頭沒有選項結構)。
  **建議修法**:把 `subspaceSmallest` 完整呼叫範例補入 spec;明確說明是新增 `solveBuckling(PreparedSystem, FrameModel, BucklingOptions)` overload 還是修改既有簽名,以及 `BucklingOptions` 的 header 路徑。

- **[MAJOR] S1 build .bat 修改清單遺漏 `build_linear_audit.bat`,且對 `build_perf.bat` 的補丁說明部分錯誤**
  S1 §⑤ 寫「4 支 build .bat 全部加 `Reanalysis.cpp`」。真實存在的 4 支是:
  `build.bat`、`build_cli.bat`、`build_linear_audit.bat`、`build_perf.bat`。
  - `build_cli.bat` 只列 `ModalAnalysis.cpp` 之前的核心模組(沒有 `Collapse.cpp` 等),是刻意的精簡清單。若 `Reanalysis.cpp` 不需要 CLI,就不應加;spec 沒有說明豁免 `build_cli.bat` 的理由。
  - S1 §③「R8:build_perf.bat 補 MITC4ShellElement.cpp」是有效發現:現有 `build_perf.bat` 確實只到 `BeamColumnElement.cpp` 就跳到 `FrameSolver.cpp`(缺 `MITC4ShellElement.cpp`)。但 S1 §⑤ 同時說「build_perf 補 MITC4ShellElement.cpp(R8)」——這只補了 shell,若 `frame_perf.cpp` 沒用到 `Reanalysis.cpp`,那 `build_perf.bat` 不必加 `Reanalysis.cpp`。Spec 把 R8 修復與 S1 新增混在同一欄,容易讓 agent 以為只要補 R8 就好。
  **建議修法**:在 §⑤ 逐支列出是否要加 `Reanalysis.cpp` 及理由;把 R8 修復(獨立 bugfix)與 S1 新增分成兩欄。

- **[MAJOR] S1 fingerprint 互動未指定**
  `PreparedSystemImpl.h` 的 `fingerprint` 是結構/幾何/UDL 的 hash,`solveLoad` 依賴它拒絕 stale factor。`ReSolveSession` 建構時拷貝 `base` 再 `assembleAndFactor`,建出自己的 `PreparedSystem`。Spec 整份沒有說:
  (a) `ReSolveSession::solve()` 回傳的 `SolveResult` 的 `pivotMargin` 要填來自 baseline LDLT 還是 Tier-1 capacitance 比值(兩者語意不同,§⑨ 只說「文檔分開講」但沒指定填哪個);
  (b) `rebaseline()` 之後新 `PreparedSystem` 的 fingerprint 是用「修改後活躍集」生成的,但 spec 沒說明消費者拿到的 `SolveResult` 應如何讓外部的 `solveLoad` 重用(或不重用)場景保持一致;
  (c) `ReSolveSession` 的 `valid()` 應於何時 false —— baseline singular 時,還是 validate() 失敗時,還是兩者?
  **建議修法**:在 §② 的 API 註解或新增 §③ data-flow 補充段中明確說明這三點。

- **[MINOR] S1 §⑧ 效能驗收門檻的基線場景未綁定機器/編譯旗標**
  「XXL 單桿移除 ≥20×」等數字來自研究原型實測,但 spec 沒說基線機器(研究輪的 `frame_perf.exe` 在什麼 CPU 上跑)、編譯旗標(原型用 `-O2 /MD`,引擎端同)、或重現指令。退步 30% 的警戒線如果在不同機器上跑根本不可比。
  **建議修法**:加「上述數字在 [機型 or 相對指標] 上取得;實作時先在同機器重跑 `frame_perf.exe` 取基線,再與 ReSolveSession 比較」。

- **[MINOR] S1 §③ data-flow 的 F 增量記帳描述缺少「殼等效載」處理路徑**
  `setShellActive` 的 F 增量只說「隨停用/恢復離開/回到 F」但殼的等效載(`ShellPressure` 的等效節點載)與梁 UDL 不同(殼使用 `assembleEquivalentLoads`),原型完全沒有殼相關代碼。Spec 的「殼 rank ≤18/片」只在 API 說明欄提及,§③ 的 F 增量記帳只舉梁 UDL 例。
  **建議修法**:在 §③ 補一行說明殼等效載的追蹤方式,或明確說「殼壓力貢獻留 oracle F34 sub-case 確認」。

---

### S2

- **[MAJOR] S2 §④ Ritz 向量生成演算法描述不完整,無法直接實作**
  S2 §④ 描述 load-dependent Ritz 生成:
  > x₁=K⁻¹f, xᵢ₊₁=K⁻¹Mxᵢ, 逐支 M-正交化(Gram-Schmidt ×2 pass), 投影小特徵問題對角化

  缺少以下關鍵決策,導致實作 agent 無法不研究就動工:
  (a) **正規化**:每支 Ritz 向量在 M-正交化後是 M-正規化還是 K-正規化?未說明。
  (b) **退化處理**:若某支 xᵢ 被前面向量 span 完全覆蓋(norm 接近零),如何處理?
  (c) **「投影小特徵問題對角化」的矩陣形式**:是 Φᵀ K Φ q = ω² Φᵀ M Φ q 嗎?投影後的矩陣用哪種求解器(dense GeneralizedSelfAdjointEigenSolver)?
  (d) **事件後 Ritz 種子**:「以該事件的不平衡力向量當 Ritz 種子」——不平衡力向量是 F' - K' u(t_e) 嗎?還是 ΔF(僅新增載)?若事件後 K' 是由 ReSolveSession 提供,如何取得 K'?
  (e) **basisSize** 是否包含靜力修正向量?若採 mode-acceleration 修正,Ritz 生成與修正項的交互未說明。
  原型 `exp_dynamic_inherit.cpp` 的 `partEquivalence()` 使用 `denseModes()` 直接全特徵分解,並沒有 Ritz 生成代碼;`partWarmStart()` 使用的是 `subspaceSmallestX0`(特徵模態 warm-start),不是 Ritz 向量。**S2 spec 宣稱 load-dependent Ritz 是研究輪驗證的結論,但原型裡沒有 Ritz 生成的 C++ 實作**,只有數值結果比較(`[trunc]` 表格)。
  **建議修法**:補一段「Ritz 生成偽代碼」或完整 C++ 原型片段;或明確說「Ritz 生成代碼尚未原型化,實作 agent 需參考 Wilson(1985)教科書流程後加 oracle F37 對照」,並把 spec 狀態從「實作就緒」降為「演算法決策已定,代碼未原型化」。

- **[MAJOR] S2 §② `DynCollapseHistory` 依賴 `CollapseOutcome` 與 `CollapseHingeEvent` 但未說明 #include 路徑**
  `DynCollapseHistory` 使用 `CollapseOutcome::Invalid` 和 `CollapseHingeEvent`,這兩個型別定義在 `Collapse.h`。S2 spec 說新增的 header 是 `Public/FrameCore/DynamicCollapse.h`,但沒有說是否 `#include "FrameCore/Collapse.h"`。考慮到 PIMPL/POD 邊界規定(鐵則 1),若 `DynamicCollapse.h` 引入 `Collapse.h`,會連帶引入 `Connectivity.h`(間接)與 `ISectionStrength.h` 等一批依賴,需確認不會造成循環 include。
  **建議修法**:在 §⑤ 的「新增 DynamicCollapse.h」欄位補一行「#include 依賴:Collapse.h, FrameModel.h, SolveOptions.h」,或改為在 spec 的 API 代碼塊顯示完整 include 清單。

- **[MAJOR] S2 §③ data-flow 中「模態 warm-start」路徑的 SparseEigsolver 呼叫未說明**
  data-flow 第 3 步:「新基底 Φ'(Ritz 重生成;模態路徑 warm-start 以 Φ 起始)」。
  `SparseEigsolver.h`(`subspaceSmallest`)使用固定隨機種子(無 X0 參數)。Warm-start 需要 `X0` 參數,但那是 `exp_dynamic_inherit.cpp` 裡的研究副本 `subspaceSmallestX0`,**不在引擎 header 裡**。S2 spec 若要實作 warm-start,必須在 `SparseEigsolver.h` 加入 `X0` 參數(或提供重載),但 spec 沒說要修改這個 header。
  **建議修法**:在 §⑤ 加入「修改 `SparseEigsolver.h`:加 `const MatX* X0 = nullptr` 參數(研究副本 `subspaceSmallestX0` 已驗,直接移植)」。

- **[MINOR] S2 §⑥ F37 oracle 門檻「relMax ≤ 1e-8」條件含糊**
  「全基底 + useRitzVectors=false → vs 全系統 Newmark relMax ≤ 1e-8(實測 1.97e-12)」。原型只跑了一次事件(step 240 of 600),使用的是 `tower(2,1,3)`,nf=108。F37 oracle 要求「含 detach 案例驗」(§⑩ 風險欄),但 §⑥ 的 F37 描述沒有明確說 detach 案例的預期 relMax 門檻是否相同(full detach 後 nf 縮小,投影殘差語意不同)。
  **建議修法**:在 F37 中分兩個 sub-case:有 detach 與無 detach,分別標出容差;或說明「含 detach 的 sub-case 採用同一 1e-8 門檻(detach 後的較小系統)」。

- **[MINOR] S2 §② `FragmentCluster` 的 `vel`/`angVel` 欄位加入方式與 `Connectivity.h` 真實宣告潛在衝突**
  S2 §② 說「Connectivity.h 修改:加 `Vec3 vel; Vec3 angVel;`,預設零 = 完全向後相容」。現有 `Connectivity.h` 的 `FragmentCluster` 是純 POD 結構(有 `real inertia[6]`,無建構子),加 `Vec3` 成員需要確認 `Vec3` 預設建構子(有:`Vec3() = default` with field-initializer `{0,0,0}`)不違反 trivially-constructible 假設。
  更重要的問題:現有 `Connectivity.h` 的 header comment 明確說「碎塊零初速交接(靜力引擎不估分離速度)」;加 `vel`/`angVel` 後,靜力驅動器(`runProgressiveCollapse`)仍輸出 0,但這個欄位現在暗示動力版本才有非零值。**現有的 `Collapse.h` 使用 `std::vector<FragmentCluster> detached` 在 `CollapseStep` 裡**,若 `FragmentCluster` 長大了 6 個 `real` 欄位,每個 step 的 `detached` 記憶體佔用略增,且現有的 OpenSees 比對 script (`opensees_compare.py`) 若有讀 fragment cluster 的序列化,需確認沒有 size 假設。
  **建議修法**:在 §⑤ 補充說明「靜力路徑 FragmentCluster 的 vel/angVel 維持零值;Collapse.h/opensees_compare.py 不需修改,但需確認沒有 sizeof(FragmentCluster) 的硬編碼」。

- **[MINOR] S2 $ExpectedUeTests 37→40 與 S1 34→37 的累加順序假設有問題**
  S1 §⑦ 說「34→37」(加 3 UE 測試),S2 §⑦ 說「37→40」(加 3 UE 測試),且 S2 §⑤ 說「run_gate.ps1 UE 37→40」。這隱含 S2 依賴 S1 先完成(否則基準是 34 不是 37)。但 S2 §⑤ 說「依賴:S1」—— S2 先行時提供了 fresh 路徑。若 S2 單獨實作(S1 未完成),`run_gate.ps1` 的 `$ExpectedUeTests` 應從 34 改到 37(跳過 S1 的 3 個測試)。
  **建議修法**:在 S2 §⑦ 補「若 S1 未完成,本輪 UE gate 數從 34→37(跳過 S1 的 3 個);若 S1 已完成則 37→40」。

---

## 沒問題的部分

S1 §② `ReSolveSession` 的 API 簽名(PIMPL 結構、zero-Eigen 公開邊界、FRAMECORE_API 標記、move-only 語意)與引擎既有慣例(`PreparedSystem` PIMPL 模式)完全相容;S1 Tier-1/Tier-2 的數值演算法主體(Woodbury 公式、stale-LDLT PCG、`StalePrecond` wrapper)在原型 `exp_incremental_refactor.cpp` 裡有完整可移植的 C++ 代碼,機制偵測邏輯(`capacitance pivot ratio < 1e-10`)也清楚;S2 §④ Newmark β=1/4 γ=1/2 的模態空間積分公式與原型一致,動量抽取(`T_trans`/`T_rot` 投影)有原型代碼支撐且數值主張邊界(own-axis 截面極慣量差異)誠實標出;兩份 spec 的誠實邊界(§⑨)措辭符合專案鐵則,先行技術定位清楚、無裸宣稱新穎性。