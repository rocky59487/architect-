# S2 交接規格 — N4 動量繼承連續動力倒塌(實作就緒*;Ritz 生成除外,見 ④)

> 狀態:核心數學已以 scratch 原型驗證(`docs/research/WS_R2_experiments.md` §6):
> 全基底跨事件繼承 vs 全系統 Newmark = 1.97e-12;動量帳對碎塊閉式 0 差;截斷誤差數據
> 直接驅動本 spec 的基底選型。原型 `Research/WS_N_incremental/exp_dynamic_inherit.cpp`。
> **誠實分級**:Newmark/投影/動量抽取 = 原型已驗可移植;**load-dependent Ritz 生成 =
> `[NEW CODE]`**(選型依據是截斷誤差數據,生成代碼未原型化,④ 給偽代碼與決策,F37 把關)。
> 本階段把既有靜力 LSP 倒塌(`runProgressiveCollapse`)升級為連續動力模擬,並修掉
> 「碎塊零初速交接」的既有誠實限制。獨創定位(WS_N):模態疊加/投影=教科書;
> 「跨拓撲事件的狀態繼承 × 遊戲物理動量交接」整合無直接先例(誠實措辭見 ⑨)。

## ① 目標 / 明確不做

**做**:`runDynamicCollapse` 驅動器(模態空間時間積分、事件觸發元素停用、跨事件狀態繼承、
碎塊帶初速交接);`FragmentCluster` 加 `vel`/`angVel`;load-dependent Ritz 基底;模態 warm-start。
**不做**:材料非線性動力、幾何非線性(大位移)、阻尼參數校準(吃使用者輸入)、
土壤/基礎、Chaos 接手後的雙向耦合(單向交接)。

## ② 公開 API

```cpp
// Public/FrameCore/DynamicCollapse.h(新檔)
namespace frame {

struct DynCollapseOptions {
    real dt            = 1e-3;    // s(引擎單位制 N-mm-tonne-s)
    real maxTime       = 10.0;
    int  basisSize     = 30;      // Ritz/模態基底數
    bool useRitzVectors = true;   // false = 純特徵模態(對照/除錯路徑)
    real rayleighAlpha = 0.0;     // C = aM + bK(模態空間保持對角)
    real rayleighBeta  = 0.0;
    real removeThreshold = 1.0;   // ElasticAllowable D/C 篩(沿用靜力驅動器語意)
    int  screenEvery   = 5;       // 每 N 步做一次 D/C 篩(recover 成本控制)
    real quietKineticRatio = 1e-6;// Stable:動能 < ratio × max 歷史動能,持續 1 個 T1
    int  maxEvents     = 64;
    std::vector<MemberId> initialRemovals;     // t=0 情境(GSA 突移)
    std::vector<int>      initialShellRemovals;
    int  frameStride   = 10;      // 每 N 步存一幀 (u,v) 快照
    SolveOptions solve;
};

struct DynCollapseEvent {
    real t = 0;
    FailMode mode = FailMode::None;
    std::vector<MemberId> removedMembers;
    std::vector<int>      removedShells;
    std::vector<CollapseHingeEvent> formedHinges;     // 塑鉸模式留 S2.1 擴充位
    std::vector<FragmentCluster> detached;            // ★ 含 vel/angVel(見下)
    real truncationResidual = 0;   // ||u − Φq||/||u|| 投影殘差(誠實指標)
    real energyBefore = 0, energyAfter = 0;           // 事件能量帳
};

struct DynCollapseFrame { real t = 0; std::vector<real> u, v; };   // 6N 快照(UE 回放)

struct DynCollapseHistory {
    CollapseOutcome outcome = CollapseOutcome::Invalid;
    std::string diagnostic;
    std::vector<DynCollapseEvent> events;
    std::vector<DynCollapseFrame> frames;
};

FRAMECORE_API DynCollapseHistory runDynamicCollapse(const FrameModel& model,
                                                    const DynCollapseOptions& opts = {});
} // namespace frame
```

```cpp
// Public/FrameCore/Connectivity.h(修改:FragmentCluster 加欄位,POD、預設零 = 完全向後相容)
    Vec3 vel;       // 分離瞬間質心線速度(mm/s, global)。靜力驅動器仍輸出 0(舊行為)。
    Vec3 angVel;    // 分離瞬間角速度(rad/s, global)= I⁻¹ L_com。
```

## ③ 資料流

```
靜平衡 u0(既有 solve;當前載重 = caller 先 compose/addSelfWeight,同靜力驅動器契約)
→ 基底 Φ0(Ritz:見 ④;特徵模態路徑用 SparseEigsolver + warm-start)
→ 模態座標初始化(q 由 u0 投影;q̇=0;從靜平衡起步 → 後續事件提供激振)
→ 時間迴圈(Newmark β=1/4 γ=1/2 per mode):
    每 screenEvery 步:u=Φq(+靜力修正) → recover → D/C 篩(checkSection/checkShellSurface 峰值)
    觸發(D/C>threshold 取最危,決定性 tie-break 同靜力驅動器):
      t_e 整步觸發 → 事件:
        1. 停用元素(工作副本;K'、M'、F' 經 ReSolveSession 或 fresh)
        2. analyzeConnectivity → detached:對每個 fragment
             p = Σ T_trans,kᵀ M_frag v(t_e)、L = Σ T_rot,kᵀ M_frag v(t_e)
             vel = p/mass、angVel = I_cluster⁻¹ L(I = 既有閉式慣量張量)
           → fragment 整塊停用 + 節點 pin(同靜力驅動器),帶初速交接 Chaos
        3. 新基底 Φ'(Ritz 重生成;模態路徑 warm-start 以 Φ 起始)
        4. 狀態繼承:q' = Φ'ᵀ M' u(t_e)、q̇' = Φ'ᵀ M' v(t_e);殘差 r=||u−Φ'q'||/||u|| 記錄
        5. 能量帳:E 前/後記入 event
→ 終止:Stable(動能準則)/ Collapsed(新構型機構=ReSolve mechanism 或 LDLT singular)/ maxTime
```

## ④ 演算法(數值細節)

- **Newmark**(原型直接移植):β=1/4、γ=1/2,無條件穩定+線性系統保能(實測漂移 6.6e-11);
  模態空間 k̂ᵢ = ωᵢ² + 2ξᵢωᵢ·(γ/βdt)… Rayleigh 阻尼下 ξᵢ=(α/ωᵢ+βωᵢ)/2 對角保持,公式照教科書。
- **基底選型(研究輪結論,關鍵)**:純特徵模態截斷誤差大(m=5/10/20 均 ~39%、m=40/108 仍
  7.3%,靜力修正僅砍半)→ **預設 load-dependent Ritz vectors**(Wilson 1985)。
  **`[NEW CODE]`:生成代碼未原型化**(原型只比較了截斷誤差,Ritz 生成須照下列偽代碼新寫,
  正確性由 F37 全基底等價 + F39 對照把關):
  ```
  // 輸入:K(ldlt 已備)、M、種子向量 g;輸出:Φ(nf×m,M-正交規範)、w2(m)
  x1 = ldlt.solve(g);  x1 /= sqrt(x1ᵀ M x1)                  // M-正規化
  for i = 2..m:
      xi = ldlt.solve(M x_{i-1})
      xi -= Σ_j (xjᵀ M xi) xj   (兩遍 Gram-Schmidt)           // M-正交化 ×2
      β = sqrt(xiᵀ M xi);  if β < 1e-8·‖x1‖_M → 退化:換隨機向量補一支重來(計數防呆)
      xi /= β
  // 投影:Kr = XᵀKX (m×m), Mr = XᵀMX (=I 數值上;仍顯式算)
  GeneralizedSelfAdjointEigenSolver(Kr, Mr) → (w2, Q);  Φ = X·Q   // Φ M-正交規範
  ```
  種子 g:初始相 = 當前載重 F;**事件後 = 殘差 r = F′ − K′·u(t_e)**(不平衡力,
  天然涵蓋局部破壞激發內容);若 ‖r‖≈0 退回 F。靜力修正與 Ritz 並用時 basisSize 不含
  修正向量(修正項按 ③ 公式另計)。
  對照路徑 useRitzVectors=false(純特徵模態)留作 audit 互驗(全基底時兩者等價)。
- **動量抽取**(原型已驗,逐式移植):T_trans = 平移剛體向量、T_rot = 對 com 旋轉剛體向量
  (平移分量 e×(r−r_com)、轉動 DOF 分量 e);p、L 由 fragment 元素一致質量矩陣投影。
  實測:質量/橫向角動量 vs `FragmentCluster` 閉式 = 0 差;own-axis 含截面極慣量項(細長桿
  可忽略;spec 決議:**p、L 用 FE 一致質量抽取(能量一致),慣量張量用 cluster 閉式**)。
- **事件定位**:整步觸發(O(dt) 時間誤差,誠實標);dt 由 caller 控,文檔給 T1/100 建議值。
- **Collapsed 判定**:事件後新構型 factor 奇異(機構)→ Collapsed;沿用靜力驅動器
  「不區分局部/全域機構」的誠實邊界。

## ⑤ 檔案清單

| 動作 | 檔案 |
|---|---|
| 新增 | `Public/FrameCore/DynamicCollapse.h`(include 依賴:`Collapse.h`(CollapseOutcome/CollapseHingeEvent)、`Connectivity.h`(FragmentCluster)、`FrameModel.h`、`SolveOptions.h`——皆 POD 公開 header,無循環;`Collapse.h` 已含 `Connectivity.h`)、`Private/DynamicCollapse.cpp` |
| 修改 | `Public/FrameCore/Connectivity.h`(FragmentCluster + `Vec3 vel; Vec3 angVel;`,預設零建構不變 POD 性;**靜力驅動器 `runProgressiveCollapse` 維持輸出零值(舊語意)**;header 的「零初速」註解改寫為「靜力路徑零初速;動力路徑(S2)回填」;grep 確認無 `sizeof(FragmentCluster)` 硬編碼) |
| 修改 | `Private/SparseEigsolver.h`(`subspaceSmallest` 加 `const MatX* X0 = nullptr` 末參數;研究副本 `subspaceSmallestX0` 已驗,直接移植;預設 nullptr = 行為逐位元不變) |
| 修改 | `Standalone/main.cpp`(F37–F39)、`linear_deep_audit.cpp`(+4 checks) |
| 修改 | build scripts:`build.bat`/`build_linear_audit.bat` +`DynamicCollapse.cpp`(cli/perf 不加,S6 再議);`run_gate.ps1` UE 期望數 +3(**基準取決於 S1 是否已完成:S1 已完成 37→40;S2 先行則 34→37**) |
| 新增 | UE 測試 ×3(等價、動量、終止狀態) |
| 修改 | `README.md`/`ARCHITECTURE.md`(動力倒塌能力+誠實邊界);`Hinge.h` 註解不動 |

依賴:S1(ReSolveSession 供事件重解;若 S2 先行,fresh assembleAndFactor 路徑亦可,介面留 hook)。

## ⑥ Oracle

- **F37 跨事件等價**:塔(2,1,3) nf=108,basisSize=nf(全基底):
  (a) useRitzVectors=false(純模態)vs 全系統 Newmark 參考(測試內建小型直接積分器,
  原型移植)逐時刻 relMax ≤ **1e-8**(實測 1.97e-12);(b) useRitzVectors=true 同案例同門檻
  (全基底時 Ritz 與模態 span 相同 → 把關 Ritz 生成 `[NEW CODE]` 正確性);
  (c) **含 detach 子案例**(鏈塔,事件切斷上段):同 1e-8 門檻,參考解在 pin 清理後的
  同一縮減系統上積分;事件前能量漂移 ≤ 1e-9。
- **F38 動量帳**:鏈分離場景:(a) 剛體測試運動下 p、L vs 閉式恰等(≤1e-12);
  (b) 動力分離瞬間 fragment p 非零且方向與 v(t_e) 一致;(c) `vel=p/m`、`angVel=I⁻¹L` 欄位回填正確。
- **F39 單自由度解析**:雙質點彈簧鏈斷開 → 留存部自由振動振幅/相位 vs 解析(≤1e-10);
  Ritz 路徑同案例 ≤1e-8。
- OpenSees(選做,strict 不擋):小框架 Newmark 同參數事件前段比對。

## ⑦ Gate 影響

F37–F39;UE `FrameCore.DynCollapse.{InheritanceEquivalence,MomentumHandoff,Outcomes}`(37→40);
audit +4(等價、動量閉合、Ritz vs 全基底殘差上界、能量帳一致)(68→72)。

## ⑧ 效能驗收

事件間每步成本 O(basisSize)(目標:nf=20k、basis=30 時每步 ≤0.5ms → 60fps 預算內 ≥30 步);
事件成本(基底重建+投影)≤ 1.5× 一次 ReSolve;數據記入 PERFORMANCE_BASELINE。

## ⑨ 誠實邊界(文檔措辭)

- 事件間線彈性模態空間;失效準則 = screening 級 D/C(同靜力驅動器,非規範檢核)。
- 截斷誤差由 `truncationResidual` 顯式報告;基底不足時誤差可觀(研究輪:純模態 m=40/108 → 7.3%)
  — Ritz 緩解但不消除;全基底 = 精確(小模型可驗)。
- 事件整步觸發(O(dt));Chaos 交接單向、碎塊離開後不回饋。
- 先行技術:模態疊加/Ritz 向量/投影皆教科書(Wilson);**宣稱措辭**:「跨拓撲事件的
  模態狀態繼承與動量保存碎塊交接,在開源參數化結構工具中未見先例([UNKNOWN] 級定位,
  WS_N 查證);數學構件均為已知方法的工程整合。」

## ⑩ 風險 / fallback

| 風險 | 緩解 |
|---|---|
| Ritz 基底仍漏高頻內容 | truncationResidual 哨兵 > 門檻 → 自動擴基(+k 支)重投影;最終 fallback = 全模態(小模型)|
| 連鎖事件密集(每幾步一次) | 事件率哨兵:> 上限改建議靜力 LSP 驅動器(文檔指引)|
| pin 清理與投影互動 | 投影在 pin 之後的新 fmap 上做(原型同序);F37 含 detach 案例驗 |
| dt 過大漏失效尖峰 | screenEvery×dt 文檔警語 + 峰值取樣選項留擴充位 |
| 能量「不守恆」質疑 | 事件能量帳顯式記錄(物理上移除元素本就改變能量;不做隱藏修正)|
