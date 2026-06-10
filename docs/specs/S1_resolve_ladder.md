# S1 交接規格 — ReSolve 三層重分析階梯 + 稀疏屈曲(實作就緒)

> 狀態:研究輪已以 scratch 原型驗證全部數值主張(`docs/research/WS_R2_experiments.md` §1–§3)。
> 本 spec 供下一輪實作 agent 直接動工;原型代碼 `Research/WS_N_incremental/exp_incremental_refactor.cpp`、
> `Research/WS_B_solver/exp_sparse_buckling.cpp` 可大段移植。

## ① 目標 / 明確不做

**做**:
1. `ReSolveSession`:拓撲增量(桿/殼 停用/恢復)後的快速重解,三層自動階梯
   (Tier-1 Woodbury 精確低秩 → Tier-2 stale-LDLT 預條件 CG → Tier-3 全重分解)。
2. **F 增量記帳(等效載追蹤)— `[NEW CODE]`**:停用/恢復帶 UDL 桿(或帶壓力殼)時其等效
   節點載同步離開/回到 F。原型只測過節點載重,此段無原型背書,風險級別高於其餘 Tier-1
   機制,oracle F34 的 UDL 子案例是它唯一的安全網(漏測會造成帶 UDL 模型的靜默錯解)。
3. `solveBuckling` 稀疏路徑(復用 `SparseEigsolver.h`),dense 為 fallback 與 oracle。
4. 修復 `build_perf.bat`(R8:源檔清單補 `MITC4ShellElement.cpp`,已實錘 LNK2001;
   獨立 bugfix,與本階段新增無耦合)。
5. `docs/PERFORMANCE_BASELINE.md` 正式化(研究輪數據為初版)。

**不做**(本階段邊界):材料/截面值變更走 Tier-3(fingerprint 已防呆);支承 flag 變更走 Tier-3;
prescribed 值與節點載重本就免重分解(既有 solveLoad 語意);多執行緒;外部後端(R13 觸發條件未到)。

## ② 公開 API(完整簽名)

```cpp
// Public/FrameCore/Reanalysis.h(新檔;POD/std 邊界,零 Eigen 洩漏,PIMPL)
namespace frame {

struct ReanalysisOptions {
    int  maxRank      = 96;      // Tier-1 累積 rank 上限(~16 桿);超過 → Tier-2
    real pcgTol       = 1e-10;   // Tier-2 相對殘差(Eigen CG 語意 ||Ax-b||/||b||)
    int  pcgMaxIter   = 500;
    bool allowTier2   = true;    // false: 超過 maxRank 直接 Tier-3
    real mechPivotTol = 1e-10;   // capacitance |piv|min/|piv|max < 此值 -> 機構
    SolveOptions solve;          // baseline assembleAndFactor 透傳
};

struct ReanalysisStats {
    int  tier        = 0;        // 實際走的層(1/2/3;0=無增量直接 baseline)
    int  rank        = 0;        // 當前累積 rank
    int  pcgIters    = 0;
    real relResidual = 0;        // tier2 殘差;tier1 填 0
    bool refactored  = false;    // 本次觸發了 rebaseline
    bool mechanism   = false;    // tier1 capacitance 奇異(= 移除後成機構)
};

class FRAMECORE_API ReSolveSession {
public:
    explicit ReSolveSession(const FrameModel& base, const ReanalysisOptions& opts = {});
    ~ReSolveSession();
    ReSolveSession(ReSolveSession&&) noexcept;
    ReSolveSession& operator=(ReSolveSession&&) noexcept;
    ReSolveSession(const ReSolveSession&) = delete;
    ReSolveSession& operator=(const ReSolveSession&) = delete;

    bool valid() const;                          // baseline 是否成功(非機構、validate 過)
    const std::string& diagnostic() const;

    // 增量:回傳 false = 未知 id。與現態相同 = no-op(回傳 true)。
    bool setMemberActive(MemberId id, bool active);
    bool setShellActive(int shellId, bool active);   // 殼 rank ≤18/片(24-6 剛體)

    // 解「目前活躍集」。結果與 fresh assembleAndFactor+solveLoad 的契約見 ⑥。
    // 機構時 SolveResult.singular = true(tier1 由 capacitance 判定;tier3 由 LDLT 判定)。
    SolveResult solve(ReanalysisStats* stats = nullptr);

    void rebaseline();                           // 強制 Tier-3:以目前活躍集重建 baseline
};

} // namespace frame
```

**語意釘死(審查補強)**:
- `valid()` = base 拷貝 `validate()` 通過 **且** baseline `assembleAndFactor` 非 singular;
  兩者任一失敗 → false + `diagnostic()` 給原因。
- `SolveResult.pivotMargin`:Tier-0/1/2 填 **baseline 構型**的 LDLT 值(語意=基準構型健康度,
  文檔明標);Tier-3 後填新分解值。Tier-1 的 capacitance pivot ratio **只**經
  `ReanalysisStats`(獨立指標,不混入 pivotMargin)。
- Session **不暴露**內部 `PreparedSystem`;外部既有 `solveLoad`+fingerprint 防呆機制
  與 session 無交互(rebaseline 純屬內部事務)。

## ③ 資料流

```
建構:base 拷貝(同 collapse driver 慣例,呼叫者模型永不被改)
  → assembleAndFactor → baseline Impl{K0, ldlt0, fmap, elems, F0(=節點載重+活躍元素等效載)}
setMemberActive/setShellActive:
  → deltaSet[(elemKind, idx, ±)] 記帳
  → F 增量記帳 [NEW CODE]:該元素 addEquivalentNodalLoads 的貢獻隨停用/恢復離開/回到 F
    (桿=UDL 等效載;殼=ShellPressure 等效載——兩者都走 IElement::addEquivalentNodalLoads
     同一介面,但原型皆未覆蓋;F34 須各含一個子案例)
solve():
  R = Σ rank(deltaSet)        // 桿 ≤6、殼 ≤18(逐元素特徵分解快取)
  if R == 0:        u = ldlt0.solve(F_ff)                        [tier 0]
  elif R ≤ maxRank: Tier-1(W/s/Z/M 增量同步 → C=diag(1/s)+M → FullPivLU
                    → pivRatio < mechPivotTol ? mechanism : u = u0' − Z C⁻¹ Wᵀ u0')
  elif allowTier2:  Tier-2(K' = K0ff + W·diag(s)·Wᵀ 顯式組裝(稀疏三重積,實測 ~6ms@XXL)
                    → CG<stale precond=ldlt0.solve> guess=上次解 → 收斂? : Tier-3)
  else:             Tier-3(fresh assembleAndFactor(work) 接管為新 baseline,清空 W/Z/deltaSet)
  → recover:以「目前活躍集」的 prepared elements 對 u 做桿端力/殼內力回復 → SolveResult
```

## ④ 演算法(數值細節)

**Tier-1(原型已驗,直接移植)**:
- 元素全域剛度 12×12(殼 24×24)由 `BeamColumnElement`/`MITC4ShellElement` `prepare+assemble` 取得
  (與引擎逐位元同源);drop 固定 DOF 後 `SelfAdjointEigenSolver`,保留 λ > 1e-9·λmax(梁 rank 恰 6)。
- `W`(nf×R,密存;列 nnz ≤24)、`Z = ldlt0.solve(W)`(密,nf×maxRank×8B = 14MB@nf=18.7k)、
  `M = WᵀZ` 增量成長(新塊 B1 = WnewᵀZold;對稱性免重算舊塊)。
- `C = diag(1/s) + M`,`FullPivLU`;pivot ratio = min|diag U|/max|diag U|。
  **機構判定**:K' 奇異 ⇔ C 奇異(行列式引理);實測健康 ~1e-8 vs 機構精確 0,門檻 1e-10。
- 恢復(active=false→true)= 追加 +λ 更新;精確抵銷實測漂移 1.46e-15@R=600。
- (F 增量見 ①-2 與 ③:`[NEW CODE]`,不在「原型已驗」範圍內。)

**Tier-2(原型已驗)**:
- `Eigen::ConjugateGradient<SpMat, Lower|Upper, StalePrecond>`;`StalePrecond::solve(b) = ldlt0.solve(b)`
  (wrapper 原型在 exp_incremental_refactor.cpp,~15 行,直接移植)。
- 實測:移除 160 桿(rank 960)18 次迭代、relErr 2.9e-12@tol 1e-10;Jacobi 對照 1273 迭代。
- 確定性:CG 無隨機,同輸入同機可重現;與直接解差 ≤ tol 級(**非逐位元** — 文檔/audit 都用容差)。

**Tier-3 / rebaseline 政策**:R > 2·maxRank 連續累積、Tier-2 不收斂、或任何疑慮 → rebaseline。
正確性永遠可由 rebaseline 兜底。

**稀疏屈曲**:介面決策 = **新增 overload,舊簽名不動**:
```cpp
// Public/FrameCore/BucklingAnalysis.h(修改:加 POD options + overload;舊單參版委派預設值)
struct BucklingOptions { int denseThreshold = 500;   // nf 超過走稀疏;<=0 = 強制稀疏
                         int nev = 1; int maxIter = 300; real tol = 1e-11; };
FRAMECORE_API BucklingResult solveBuckling(const PreparedSystem&, const FrameModel&,
                                           const BucklingOptions& opts);
```
稀疏路徑完整呼叫(注意 out-params 與 bool 回傳,原型 `exp_sparse_buckling.cpp` 同款):
```cpp
VecX lambda; MatX vec;
const bool ok = subspaceSmallest(Kff, S.ldlt, negKgff, opts.nev, lambda, vec,
                                 opts.maxIter, opts.tol);
if (!ok) { /* dense fallback(既有路徑,不動) */ }
else     { R.criticalFactor = lambda(0); /* mode 由 vec.col(0) scatter 回 6N */ }
```
`gtrips.empty()`(全拉)→ 既有診斷不變。殘差驗收:||(−Kg)φ − γKφ|| / (γ||Kφ||) 記入
(實測 1e-7~1e-10)。

## ⑤ 檔案清單(含建置同步)

| 動作 | 檔案 |
|---|---|
| 新增 | `Source/FrameCore/Public/FrameCore/Reanalysis.h` |
| 新增 | `Source/FrameCore/Private/Reanalysis.cpp` |
| 修改 | `Source/FrameCore/Private/BucklingAnalysis.cpp`(sparse 路徑 + fallback) |
| 修改 | `Standalone/main.cpp`(F34–F36 fixtures) |
| 修改 | `Standalone/linear_deep_audit.cpp`(+6 checks,見 ⑦) |
| 修改 | build scripts 逐支:`build.bat` +`Reanalysis.cpp`(gate 主建置);`build_linear_audit.bat` +`Reanalysis.cpp`(audit 用到 ladder 互驗);`build_cli.bat` **不加**(frame_cli 不用 ladder;若 S6 加 daemon 指令再補)— 豁免理由要寫在 .bat 註解;`build_perf.bat` 視 frame_perf 是否加 ladder 計時段決定(預設不加) |
| 修改 | **`build_perf.bat` 補 `MITC4ShellElement.cpp`(R8 獨立 bugfix,與上一列分開 commit 訊息點名)** |
| 修改 | `Scripts/run_gate.ps1` `$ExpectedUeTests` 34→37 |
| 新增 | UE 測試 `FrameCoreReanalysisTest.cpp`(×2)、buckling sparse 測試(×1) |
| 修改 | `README.md` / `docs/ARCHITECTURE.md`(能力與誠實邊界段) |

## ⑥ Oracle(數學 + 門檻)

- **F34 ReSolve tier-1**:塔(3,2,4)+ 一根帶 UDL 的桿;序列移除 8 桿(含該 UDL 桿)+ 全恢復:
  每步 vs fresh `assembleAndFactor+solveLoad` relMax ≤ **1e-10**(實測 1e-13 級);恢復後 vs baseline ≤ 1e-12;
  機構案例:2 元素鏈移除底元素 → `mechanism=true` 且 fresh `singular=true`;portal 移除梁 → 兩者皆穩定。
- **F35 ReSolve tier-2**:塔(3,2,4)移除 20 桿(rank>96)→ tier=2、relErr vs fresh ≤ **1e-8**、迭代數 ≤ 100。
- **F36 sparse buckling**:column8(vs dense ≤1e-6 且 vs Euler 解析 ≤1e-4)、tower-S、tower-M(vs dense ≤1e-6);
  全拉退化 → 與 dense 同診斷。
- OpenSees:既有「移除態」逐位移比對場景改走 ReSolveSession 路徑重跑一次(同容差)。

## ⑦ Gate 影響

standalone F34–F36;UE `FrameCore.Reanalysis.LadderAgreesFresh`、`FrameCore.Reanalysis.MechanismDetection`、
`FrameCore.Buckling.SparseAgreesDense`(34→37);audit +6:tier1 一致性、恢復漂移、機構偵測、tier2 容差、
buckling sparse vs dense、全拉護衛(62→68)。

## ⑧ 效能驗收(記入 PERFORMANCE_BASELINE.md)

XXL 單桿移除 ≥20×(研究輪 31×);50 桿序列尾端 ≥12×(實測 17.5×);tier-2 160 桿 ≥5×(實測 8.2×);
tower-M sparse buckling ≥10×(實測 24.6×)。退步超過 30% 視為驗收失敗。
**機器綁定**:上述絕對時間出自研究輪單機(cl /O2 /MD,單緒);驗收時**同機**先重跑
`Research\bin\exp_incremental_refactor.exe` 取當日基線,倍率對倍率比,不跨機比絕對毫秒。

## ⑨ 誠實邊界(文檔措辭)

- Tier-1/2 = 「同構型拓撲增量的精確/容差級重解」:節點集、支承、材料截面值不變(變更自動 Tier-3)。
- Tier-2 結果與直接解差 ≤ pcgTol 級,非逐位元;同輸入確定性可重現。
- `pivotMargin` 語意不延伸到 ladder:tier-1 報 capacitance pivot ratio(獨立指標,文檔分開講)。
- 先行技術(WS_N):SMW 重分析=Akgün et al. 2001、rank-k downdate=Davis & Hager;
  FrameCore 的貢獻=「三層自動階梯 × 倒塌/優化驅動器 × 互動編輯」的零依賴工程整合,非數學新算法。

## ⑩ 風險 / fallback

| 風險 | 緩解 |
|---|---|
| 長序列漂移 | 實測 1.5e-15@R=600 無虞;仍設 R>2·maxRank 強制 rebaseline + audit 殘差哨兵 |
| Z 記憶體 | nf×maxRank×8B 封頂(maxRank 可調);MEGA 96 rank = 47MB |
| 殼 rank-18 未在原型驗過 | 實作時先加殼版 F34 子案例;不過 → 殼一律 Tier-2/3(階梯仍成立) |
| Tier-2 不收斂(K' 與 K0 差太遠) | 自動 Tier-3;統計 stats.refactored 供調 maxRank |
| 任何數值疑慮 | rebaseline 永遠正確(等價 fresh 路徑) |
