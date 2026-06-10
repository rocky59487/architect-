# FrameCore 效能基線(初版:2026-06-10 研究輪)

> 量測代碼:`Research/`(scratch,不入 gate);原始輸出 `Research/out/*.txt`;單機
> (Windows 11,cl /O2,Eigen 3.4.0 SimplicialLDLT/AMD,單執行緒)。所有數字 `[VERIFIED]`
> 除非標注外推。S1 起本檔轉正式:每階段效能驗收(spec ⑧ 節)結果回填於此,
> 退步 >30% 視為驗收失敗。

## 1. 靜力直接解規模階梯(`exp_million_dof`,塔式 benchmark 家族)

| 自由 DOF | nnz(K) | assembleAndFactor | solveLoad | factor 後記憶體 |
|---|---|---|---|---|
| 18,720 | 0.82M | 1.55 s | 17.9 ms | 120 MiB |
| 61,560 | 2.67M | 30.6 s | 108 ms | 585 MiB |
| 186,000 | 7.99M | 522 s | 805 ms | 3,257 MiB |
| 389,664 | 16.66M | **3,229 s(53.8 min)** | 3,749 ms | 10,357 MiB |
| 1,007,964(模型實測,分解未跑) | — | **未實跑**;4 點擬合 ~O(n^2.46) 外推 ≈ 9.3 h、~39 GB `[THEORY:外推]` | — | — |

註:186k/390k 兩階段與其他研究作業並行執行(CPU 餘裕充足、全程 CPU-bound 無 swap,
peak≈WS 證實),計時受並行影響估 <10%;1M 在分解前截斷(模型建構/統計為實測)。

讀法:factor 隨 n 超線性(~n^2.4,3D 框架 AMD 排序填入支配);**互動價值區間 ≤ ~60k DOF**;
重解(solveLoad)始終毫秒級 = factorize-once 架構的紅利;百萬 DOF 屬批次/離線域(R13 觸發條件)。
量測變動註:同一 nf=18,720 模型在 `exp_incremental_refactor` 獨立執行中 factor=1,669ms
(vs 本表 1,546ms),±8% 屬執行間 OS 調度變動,跨文檔引用時以各自原始輸出為準。

## 2. ReSolve 階梯(N1;`exp_incremental_refactor`,XXL nf=18,720)

| 操作 | 時間 | vs fresh(1.4–1.8s 量測變動,中位 ~1.6s) | 精度 |
|---|---|---|---|
| Tier-1 單桿移除 | 54.5 ms | **31×** | 7.7e-14 |
| Tier-1 第 50 桿(R=300) | 82 ms | 17.5× | 5.7e-13 |
| Tier-1 移除+恢復 ×50(R=600) | 67.5 ms/解 | — | 漂移 1.46e-15 |
| Tier-2 批量 16 桿(PCG 12 迭代) | 155 ms | 10.3× | 2.1e-12 |
| Tier-2 批量 160 桿(18 迭代) | 221 ms | 8.2× | 2.9e-12 |
| 純回代(重解下限) | 10.9 ms | — | — |

## 3. 特徵值問題

| 問題 | dense | sparse(subspace) | 加速 | 精度 |
|---|---|---|---|---|
| 屈曲 nf=1,440 | 2,109 ms | 85.6 ms | **24.6×** | 2.6e-14 |
| 屈曲 nf=18,720 | 不可行(~5.3 GB) | 5,193 ms | — | 殘差 1.2e-7 |
| 模態 warm-start(nev=10, nf=288) | — | cold 10 → warm 7 迭代 | 1.4× | λ 一致 7.5e-14 |

## 4. 直接 solver 橫評(同一 K_ff;`exp_solver_compare`)

XXL(nf=18,720):SimplicialLDLT 1,364ms ≈ SimplicialLLT 1,403ms < SparseLU 3,492ms(慢 2.6×;
小模型差距更大:tower-M 27.8ms vs 5.4ms ≈ 5.1×);
CG+Jacobi 1,248 迭代/423ms(單發);**CG+IncompleteCholesky 2,252 迭代/2,402ms(劣於 Jacobi)**
→ 通用 IC 預條件對框架剛度病態無效;引擎預設(SimplicialLDLT)維持正確選擇。

## 5. CLI 端到端(GH 橋參考;`cli_throughput.py`)

| 模型 | 端到端中位 | 備註 |
|---|---|---|
| 行程開銷(12 DOF) | 6.7 ms | spawn+parse 下限 |
| 1,620 DOF | 26 ms | 互動無感 |
| 19,500 DOF | 2.06 s | 人工觸發可用 |
| 107,502 DOF | 142.7 s | factor 支配 → daemon/C API 域 |

## 6. 動力(N4;`exp_dynamic_inherit`,nf=108)

全基底模態 Newmark 與全系統 Newmark 等價 1.97e-12;截斷誤差:m=5/10/20 均 ~39%、
m=40 才降至 7.3%(nf=108;mode-acceleration 修正約砍半至 19.5%/5.1%)→ S2 spec 改推
load-dependent Ritz;事件投影 O(basis·nf);warm-start 省 30% 特徵迭代。每步成本與基底數線性 — S2 驗收目標:nf=20k、basis=30、
每步 ≤0.5ms(60fps 預算內 ≥30 步)`[PENDING:S2]`。
