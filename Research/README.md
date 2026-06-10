# Research 目錄 — 研究輪 scratch 實驗(非引擎)

**本目錄全部是研究實驗 scratch 代碼,不是 FrameCore 引擎的一部分。**

- 不修改 `Source\FrameCore\Public\`、`Source\FrameCore\Private\`、`Standalone\` 下任何既有檔案。
- 不入驗證閘門(`Scripts\run_gate.ps1` 不包含本目錄)。
- 實驗 exe 透過 `#include` 引擎 Public + Private headers 並連結引擎 obj 編成獨立執行檔
  (Private 內省是研究特權:正式功能進引擎時走公開 API + 導入階梯)。
- 一鍵建置:`build_research.bat`(全部)/ `build_research.bat -skipcore`(引擎 obj 已在時跳過重編)。
- 產物:`bin\exp_*.exe`;實驗輸出記錄在 `out\*.txt`(可重現:重跑同 exe 同參數)。

| 實驗 | 對應研究線 | 問題 |
|---|---|---|
| `WS_N_incremental\exp_incremental_refactor.cpp` | N1 增量重分析 | Woodbury rank-≤6 更新解元素移除/恢復,精度/漂移/加速比/機構偵測 |
| `WS_B_solver\exp_sparse_buckling.cpp` | WS-B | subspace iteration 稀疏屈曲 vs dense GEVP,含退化護衛 |
| `WS_B_solver\exp_million_dof.cpp` | WS-B | 百萬自由 DOF 真跑(時間/記憶體;OOM 也是數據) |
| `WS_B_solver\exp_solver_compare.cpp` | WS-B | Eigen 內建 solver 比較(LDLT/LLT/SparseLU/CG) |
| `WS_C_pdelta\exp_pdelta_convergence.cpp` | WS-C + N3 | 凍結分解 P-Delta vs 重分解參考版 vs 精確解 |
| `WS_D_tensiononly\exp_tension_only.cpp` | WS-D | tension-only 迭代收斂/振盪防護 |
| `WS_H_sizeopt\exp_size_opt.cpp` | WS-H | 10 桿桁架 FSD 收斂 |
| `WS_I_beso\exp_beso_truss.cpp` | WS-I + N2 | ground structure BESO(Michell 定性) |
| `WS_N_incremental\exp_dynamic_inherit.cpp` | N4 | 跨事件模態繼承等價性、碎塊動量帳、modal warm-start |
| `WS_F_corot\elastica_check.py` | WS-F | elastica 精確解數值表(純 Python) |
| `WS_J_gh\cli_throughput.py` | WS-J | frame_cli 端到端吞吐(行程開銷/1.6k/19.5k/107k DOF) |

研究結論寫在 `..\docs\research\WS_*.md` 與 `..\docs\KARAMBA3D_ROADMAP.md`。
