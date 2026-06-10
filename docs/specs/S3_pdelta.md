# S3 交接規格 — P-Delta 二階分析(N3 凍結分解 + 參考路徑雙路互驗;實作就緒)

> 研究輪驗證:凍結分解 pseudo-load 迭代與重組 K_T 參考解同一不動點(≤5.4e-13)、
> 對 beam-column 精確解 2.6e-7(f=0.3, 8 元素)、迭代數精準符合幾何級數、f>1 發散(物理正確)。
> 原型 `Research/WS_C_pdelta/exp_pdelta_convergence.cpp`;數據 WS_R2 §4;文獻定位 WS_N §N3
> (= Wilson & Habibullah 1987 / 業界 nodal-shear 法的架構整合,非新算法 — 文檔照此措辭)。

## ① 目標 / 不做

**做**:`runPDelta` 驅動器,兩條路徑:(A) 參考路徑 = K_T = K_e + Kg 重組單次分解
(Wilson 1987 形);(B) 凍結路徑 = 重用 K_e 既有 LDLT 的 pseudo-load 迭代(預設,
與 factorize-once/N1 零摩擦)。雙路徑 audit 互鎖。
**不做**:大位移(CR 是 S9)、Kg 隨迭代更新軸力(N 凍結於一階解 = Th.II 慣例,誠實標)、
殼幾何勁度(既有限制沿用)、動力 P-Delta。

## ② 公開 API

```cpp
// Public/FrameCore/PDeltaAnalysis.h(新檔)
namespace frame {
struct PDeltaOptions {
    int  maxIter    = 200;
    real tolU       = 1e-12;    // ||Δu||/||u|| 收斂
    bool accelerate = true;     // 保護式幾何外推(見 ④;裸 Aitken 禁用)
    bool refactorPath = false;  // true = 參考路徑(K_T 單次分解)
    SolveOptions solve;
};
struct PDeltaResult {
    bool converged = false;
    bool diverged  = false;     // P 超過臨界(或不穩定):凍結路徑增長 / 參考路徑 K_T 非正定
    int  iterations = 0;
    real lastIncrement = 0;     // 最後 ||Δu||/||u||
    SolveResult finalState;     // recover 過的完整結果(含桿端力)
};
FRAMECORE_API PDeltaResult runPDelta(const FrameModel& model, const PDeltaOptions& opts = {});
} // namespace frame
```

## ③ 資料流

```
solve 線性 → axial N(凍結) → assembleGeometric → Kg_ff
refactorPath: ldlt(K_ff + Kg_ff) → u(一次);K_T 非正定/奇異 → diverged
凍結路徑:   u ← ldlt0.solve(F_ff − Kg_ff·u) 迭代(ldlt0 = 既有 PreparedSystem 分解)
            (收斂域 P < Pcr;收斂率 = P/Pcr,文檔附迭代數表)
→ recover(最終 u)→ PDeltaResult
```

## ④ 演算法細節(研究輪教訓內建;誠實分級)

- **保護式外推 — `[NEW CODE]` 設計,非原型已驗**:原型實作的是近裸 Aitken(僅 0<ρ<0.999
  上限,無穩定窗、無撤銷),其 f=0.9 劣化 285→4742、f=0.95 改善 569→57 兩個數字**都屬
  近裸版**——劣化案例正是本保護設計的動機。保護式規格:
  ρ_m = ⟨Δ_m,Δ_{m-1}⟩/⟨Δ_{m-1},Δ_{m-1}⟩;僅當 0<ρ<0.95 且 |ρ_m−ρ_{m-1}|<0.2 才外推
  x ← x + ρ/(1−ρ)Δ;外推後下一步 ||Δ|| 未降 → 撤銷該外推並停用 accelerate(本次求解內)。
  實作時以 F40 重測保護版迭代數並回填本表(`[PENDING:S3 實測]`);驗收底線=任何 f 下
  不得劣於無外推版 1.2×。
- **發散偵測 — `[NEW CODE]` 設計**:原型的單步增長比(>4×)在 f=1.05 慢發散(每迭代 ~×1.05)
  下**未觸發**,以 maxIter 兜底(`div=0, conv=0`)——此即新設計動機。規格:20 步滑動窗
  ||Δ|| 趨勢上升且超過初值 → diverged;maxIter 兜底(到頂未收斂 → converged=false
  並回報 lastIncrement)。
- P=0 退化:第一步即收斂、結果逐位元 = 線性解(F40 鎖)。

## ⑤ 檔案

新 `Public/FrameCore/PDeltaAnalysis.h`、`Private/PDeltaAnalysis.cpp`;`main.cpp`(F40–F41);
`linear_deep_audit.cpp` +3;4 支 .bat;`run_gate.ps1` UE +1(40→41);`Tools/pdelta_compare.py` 新增;
README/ARCHITECTURE 能力段。

## ⑥ Oracle

- **F40**:懸臂柱(8 元素)f∈{0, 0.3, 0.95}:凍結 vs 參考 ≤1e-10;參考 vs 精確
  δ=H(tan(kL)−kL)/(Pk) ≤1e-3(實測 2.6e-7@0.3、3.7e-5@0.95);P=0 逐位元=線性;
  f=1.05:參考路徑 K_T 非正定 → diverged=true;凍結路徑 → **diverged=true(新偵測器
  須抓到;原型以 maxIter 兜底 conv=0 div=0,故此旗標是對新偵測器的驗收,非原型重現)**。
- **F41/OpenSees**:`Tools/pdelta_compare.py` 以 PDeltaCrdTransf + elasticBeamColumn 比對
  (其近似語意先由 WS_C 文獻報告釘死,再定容差;預設 1e-6 起跳,實測後收緊)。

## ⑦ Gate

F40–F41;UE `FrameCore.PDelta.AmplificationOracle`;audit +3(雙路徑互驗、P=0 退化、發散旗標)(72→75)。

## ⑧ 效能驗收

凍結路徑 f=0.5 迭代 ≤60(實測 46)且總時間 < 參考路徑重分解(中型以上模型);記入 baseline。

## ⑨ 誠實邊界

Th.II 線性化:N 凍結於一階解、小側移假設;f→1 收斂慢(迭代 ~log tol/log f);
不適用大位移/後挫屈(→S9 CR);殼不貢獻 Kg(既有限制)。措辭:「實作 Wilson & Habibullah (1987)
線性化 P-Delta 與業界 pseudo-load 迭代法,於 factorize-once 架構零重分解整合」。

## ⑩ 風險/fallback

外推失誤 → 自動停用降級裸迭代(仍收斂,只慢);f 接近 1 → 迭代數爆 → maxIter 回報未收斂
+ 文檔指引用參考路徑;參考路徑 K_T 接近奇異 → LDLT 診斷回報 diverged。
