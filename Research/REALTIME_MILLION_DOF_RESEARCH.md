# 即時百萬 DOF 結構求解 — 研究方向:壓低 factor,單一 direct 算法

> **開新視窗用這份無縫接續。** 前一個對話太長。**方向已修正(見「方向轉變」)**:不走 HP 迭代/seed 路線,改**正面壓低 LDLᵀ factor**,讓單一 sparse-direct 算法 scale 到百萬 DOF,廢掉 HP 雙車道。語言一律中文。

## 一句話目標
即時(互動級 ~16–100ms / 重解)求解百萬(~10⁶)DOF 線彈性,精度 `‖u−u_exact‖/‖u_exact‖ ≤ 1e-9`,**用單一 direct(Cholesky/LDLᵀ)算法**(factor 一次 + 每幀 backsub),不要 HP 雙車道。

## 方向轉變(關鍵 — 使用者洞察,推翻前一版 handoff)
前一版提議用 matrix-free PCG + AMG 繞過 LDLT factor wall。**那是錯的優先級。** 使用者指出:
- **HP 依賴 seed = 重傷**:out-of-subspace / 單次 / 任意 RHS 慘敗 40–700x(benchmark 實證),與自家 LDLT 形成「窄賽道自我競爭」的雙車道,維護負擔大。
- **direct 的 backsub 本來就 ~17x 勝 PCG**(WS_B 反覆驗證的硬事實);direct 唯一痛點是 **factor 貴**。
- **正解是壓低 factor,不是繞過它**:若 factor 能壓到百萬 DOF 可行,則「factor 一次(結構固定,setup 攤銷)+ 每幀 backsub(~ms)」全面勝出 —— 無 seed、無 PCG、無 AMG、無雙車道,**單一乾淨算法**。HP 連同其 seed 重傷整個退役。

## factor wall 很可能是「假牆」(誠實技術判斷,待驗)
FrameCore 現用 **`Eigen::SimplicialLDLT`** — simplicial(逐列、無 supernode、不吃 BLAS3 dense kernel),naive 實作。WS_B 實測 factor 48k DOF 已 20.7s,**但這很可能是 Eigen simplicial 太弱,不是 factor 的本質下限**:
- **fill-reducing ordering**(AMD / METIS / nested dissection):3D 問題用 nested dissection 把 fill 壓到近最優;Eigen 預設 AMD,METIS ND 對 3D 通常更好 —— **單這步可能就推牆數倍**。
- **supernodal / multifrontal**(CHOLMOD / PARDISO / MUMPS):factor 組織成 dense supernode blocks 走 BLAS3,典型比 simplicial 快 **10–50x** + 多執行緒 + out-of-core。
→ 換好 ordering + supernodal,factor wall 可能從 48k 推到百萬級(setup 秒–分鐘、記憶體 GB 級,**對固定結構 factor-once 可接受**)。

## 三方法(達成「壓低 factor → 單一 direct」)
1. **Fill-reducing ordering**:Eigen 現用 ordering vs METIS nested dissection 對 3D 彈性的 fill / factor 時間影響。
2. **Supernodal / multifrontal factor**:(a) 工業庫 CHOLMOD(supernodal)/ PARDISO / MUMPS — 快但**新依賴**(授權 + dual-build + UE 整合,鐵則張力,先評估);(b) 自實作 supernodal Cholesky — 無依賴但大工作。
3. **Factor-once + backsub-many**:固定結構 factor 一次(PreparedSystem 既有設計),每幀只 backsub(O(fill) ~ms)。量百萬 DOF backsub 是否即時。

## 載體:先只做可行性評估 + 小原型(使用者要求)
不投入大工程。先 go/no-go:壓低 factor 能否讓單一 direct scale 到百萬 DOF + 即時 backsub。**成功 → 廢 HP 雙車道,引擎統一單一算法(這才是真目標)。** research-only,評估階段不進五腿 gate、不改 default solve,現有 LDLT / 解析解為 oracle。

## 可行性評估要回答的(go/no-go)
1. **ordering 的威力**:同一 3D 彈性網格,Eigen 預設 vs METIS nested dissection,factor 時間 + fill(記憶體)差多少?(可能單步推牆數倍)
2. **simplicial vs supernodal**:Eigen SimplicialLDLT vs CHOLMOD supernodal(或自實作)在 ~10萬 DOF 的 factor 時間,外推百萬。快幾倍?
3. **百萬 DOF factor 可行性**:最佳 ordering + supernodal 下,百萬 DOF 3D factor setup(秒?分鐘?)+ fill 記憶體(GB?fit RAM?)。
4. **每幀 backsub 即時性**:百萬 DOF backsub 時間(~ms?)→ 達互動 / 幀率?
5. **依賴取捨**:CHOLMOD/PARDISO/MUMPS 的授權 + 純 C++17/Eigen 雙車道 + UE dual-build + 五腿 gate 整合代價 vs 自實作 supernodal 的工作量。
6. **壓低後 HP 是否完全多餘**:確認「factor-once + backsub」在固定結構全面勝 HP → 可廢雙車道(結構**改變**的增量更新是 ReSolve/S1 領域,非 HP seed)。
7. **factor 的本質下限**:3D nested dissection 的 O(n²) flops / O(n^4/3) fill 在百萬 DOF 是否真可接受,還是有硬天花板(若有 → 才回頭考慮迭代,但 HP seed 重傷仍在)。

## 現狀 / 路徑 / 工具(交接)
- branch `research/hpfem-solver-v1`,HEAD = 本檔最新 commit。**WS_B seeded HpSession 已完成但窄、未發布、本方向下大機率退役**(保留作對照 oracle / 結構變動的 ReSolve 增量場景)。
- LDLT 在 `Source/FrameCore/Private/FrameEigen.h`:`using LDLTSolver = Eigen::SimplicialLDLT<SpMat>`;`Private/FrameSolver.cpp` 的 assembleAndFactor/solveLoad;`Private/PreparedSystemImpl.h` 持 K/fmap/ldlt(factor-once 既有設計)。
- benchmark 範本:`Standalone/hp_bench.cpp`(改成量 factor / backsub vs ordering / solver)。
- 背景:`Research/WS_B_solver/HPFEM_RESEARCH_NOTES.md`(factor wall scale story、reused-backsub ~17x 勝 PCG 的硬事實);記憶 `hpfem-solver-research.md`(durable 教訓 ①–⑨ + 環境踩雷)。

## 鐵則
research-only(評估階段不進五腿 gate、不改 default solve);現有 LDLT / 解析解為正確性 oracle;**誠實:報實測 factor/backsub 時間 + 記憶體,別宣稱未驗證的 scaling**(上一輪最大教訓);過 gate 才 commit;顯式 `git add` 不 `-A`;不碰 `.gitignore`/`ArchSim.uproject`/`Plugins/LevelSim/`/`Research/WS_C`/`Research/WS_N`/build 產物。**新依賴(CHOLMOD/METIS/...)先在評估報告明列授權 + dual-build + 鐵則張力,不擅自引入。** 環境:bash cwd 在 `/e/project` → `git -C /e/project/ArchSim`;UE build 引號/PowerShell 正斜線見記憶。

## 建議第一步
小原型:可調規模 3D structured 彈性 grid(~1萬→10萬 DOF),量:(a) Eigen SimplicialLDLT factor 時間 + 記憶體 vs 規模 + ordering(natural / AMD / 若可接 METIS);(b) backsub 時間 vs 規模;(c) 外推百萬 DOF 的 factor setup / 記憶體 / backsub。若 ordering 單步就大幅推牆 → 可能不需工業庫。據此寫 **go/no-go:單一 direct 能否 scale 到百萬 DOF 即時** → 能則 HP 退役、引擎統一。**先評估,別造大輪子。** 用 Workflow / 平行 agent + 對抗式查核維持誠實。
