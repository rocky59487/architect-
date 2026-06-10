# S4 交接規格 — Tension-only 桿件(實作就緒)

> 研究輪驗證:X 斜撐迭代 2 步收斂,收斂解與「該桿直接從 member 清單省略」**逐位元相等**;
> 振盪掃描無 flip-flop(誠實負結果,守門仍留);原型 `Research/WS_D_tensiononly/exp_tension_only.cpp`;
> 數據 WS_R2 §5;文獻 WS_D(Karamba Eliminator 行為查證、LCP 視角)。

## ① 目標 / 不做

**做**:`Member.tensionOnly` 旗標 + `runTensionOnly` 迭代驅動器(停用受壓 TO 桿、
伸長判據再啟用、循環守門、單調 fallback);內迴圈走 ReSolveSession(每翻轉 = rank-6 更新)。
**不做**:compression-only(對稱擴充留欄位語意但不實作)、初始鬆弛/預張力(未來欄位)、
纜索垂度非線性(超出線性引擎)。

## ② 公開 API

```cpp
// Public/FrameCore/Member.h(修改):
    bool tensionOnly = false;   // true = 僅受拉:驅動器迭代中受壓即停用。入 fingerprint
                                // (改旗標 = 改驅動行為語意,堵 stale reuse 歧義;成本零)。

// Public/FrameCore/TensionOnly.h(新檔)
namespace frame {
struct TensionOnlyOptions {
    int  maxIter = 32;
    bool allowReactivation = true;   // false = 單調停用(保守、必然有限終止)
    real axialTol = 0.0;             // N > axialTol(壓)才停用;0 = 機器級
    SolveOptions solve;
};
struct TensionOnlyResult {
    bool converged = false;
    bool cycled    = false;          // 觸發循環守門(已自動降級單調重跑,結果為單調解)
    int  iterations = 0;
    SolveResult finalState;
    std::vector<MemberId> slack;     // 收斂時被停用的 tension-only 桿
};
FRAMECORE_API TensionOnlyResult runTensionOnly(const FrameModel& model,
                                               const TensionOnlyOptions& opts = {});
} // namespace frame
```

## ③ 資料流 / ④ 演算法

```
工作副本 → ReSolveSession(base)
迭代(Jacobi 全掃描):
  解 → 對每根 tensionOnly 桿:
    active 且 N(endI) > axialTol(壓,compression-positive)→ 停用
    inactive 且軸向伸長 (u_j−u_i)·x̂ > 0(memberLocalAxes)→ 再啟用(allowReactivation)
  無翻轉 → converged
  (state,flips) 轉移哈希重現 → cycled → 自動以 allowReactivation=false 單調重跑
    (單調版每步只停用、永不恢復 → ≤ n_TO 步必終止;diagnostic 標註順序相依的保守性)
```
收斂性質誠實定位(WS_D):主動集不動點迭代無普遍收斂保證(LCP 視角);
驅動器以「哈希守門 + 單調 fallback」保證**有限終止**。

## ⑤ 檔案

改 `Member.h`(+tensionOnly)、`FrameModel.cpp` validate(TO 桿須 active 初始 true 不強制,
但 TO+release 全放鬆組合警示)、`FrameSolver.cpp` modelFingerprint:**目前 member 迴圈
hash 到 `active` 為止、無 tensionOnly——在 `active` 之後插一行
`h = fpMix(h, mem.tensionOnly ? 1ull : 0ull);`(兩檔必須同 commit,漏 fingerprint = 靜默 stale)**;
新 `Public/FrameCore/TensionOnly.h`、`Private/TensionOnly.cpp`;F42–F43;audit +2;
4 支 .bat;UE +1(41→42);frame_cli 的 `MEMBER` 行尾 `tonly` token 排 S6 一起。

## ⑥ Oracle

- **F42**:X 斜撐 portal:收斂 ≤3 迭代;拉桿 N<0;**收斂解 vs 省略壓桿之模型逐位元相等**
  (≤1e-15;實測 0.0);fingerprint:翻 tensionOnly 旗標後 solveLoad 拒絕舊 factor。
- **F43**:循環守門:以合成翻轉序列直接驅動內部狀態機驗有限終止;單調 fallback 結果
  = 該順序下的有效平衡(驗 N 符號全符合 TO 約束)。物理 flip-flop 案例 [PENDING:掃描未尋獲,
  尋獲後升級為真實案例]。

## ⑦ Gate

F42–F43;UE `FrameCore.TensionOnly.Eliminator`;audit +2(省略等價、fingerprint 守門)(75→77)。

## ⑧ 效能驗收

n 根 TO 桿收斂 ≤ n+2 次重解;每次重解走 ReSolve tier-1(rank-6/翻轉)≥10× vs fresh(中型以上)。

## ⑨ 誠實邊界

軸力符號判據(無預張力/鬆弛長度);收斂解可能順序相依(cycled 時)且偏保守 — diagnostic 明標;
Karamba「Tension/Compression Eliminator」行為對照(WS_D 查證結論)寫進文檔對標矩陣。

## ⑩ 風險/fallback

循環(理論可能)→ 守門+單調降級(必終止);TO 桿全停用後機構 → ReSolve mechanism /
LDLT singular 正常診斷回報;與塑鉸/倒塌驅動器的組合語意(TO 桿在 collapse driver 內)
留 S2.1/S7 整合議題,本階段 runTensionOnly 為獨立驅動器。
