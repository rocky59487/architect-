# WS-R2 實驗彙編 — 研究輪 scratch 數據(2026-06-10)

> 所有實驗代碼在 `Research\`(非引擎、不入 gate),建置 `Research\build_research.bat`,原始輸出在 `Research\out\*.txt`。
> 機器:Windows 11、VS2026 preview cl /O2、Eigen 3.4.0(UE ThirdParty)。本檔所有數字皆 `[VERIFIED]`(可重跑)。

## 1. N1 ReSolve 階梯(`exp_incremental_refactor.exe`,out: `incremental_xxl.txt`)

模型:XXL 塔(3,250 節點、9,816 桿、nf=18,720)。Baseline:factor 1,669ms、engine solveLoad 18.9ms、純回代 10.9ms(同模型在 scale_ladder 獨立執行 factor=1,546ms,±8% 執行間變動)。

**Tier-1 Woodbury(逐桿移除,rank 6/桿)**:
- 單桿移除:ladder 54.5ms vs fresh 1,687ms = **31.0×**;relErr vs fresh = 7.7e-14。
- 50 桿連續移除:relErr 由 7.7e-14 緩慢增長至 5.7e-13(仍機器精度量級,未發散);速度 31×→17.5×(R=300 時)優雅退化。
- 移除 50 + 恢復 50(累積 R=600 更新):**對 baseline 漂移 1.46e-15** = 機器精度,精確抵銷成立。
- 機構偵測:鏈狀懸臂移除底元素 → capacitance pivot ratio = **0.000e+00**,與 fresh `singular=1` 一致;portal 移除梁 → piv=1.7e-3 判穩定,relErr vs fresh 8.2e-14。健康 piv ~1e-8、機構 = 精確 0,分離清晰(門檻 1e-10)。
- 誠實註記:ladder 計時只含解 u(不含桿端力 recover);fresh 含全套。純回代 10.9ms 供公平讀數。

**Tier-2 stale-LDLT 預條件 PCG(整批移除,tol 1e-10)**:
| 移除桿數 | rank | stale 迭代 | stale 總時間 | fresh | 加速 | relErr | Jacobi 對照迭代 |
|---|---|---|---|---|---|---|---|
| 16 | 96 | 12 | 155ms | 1,589ms | 10.3× | 2.1e-12 | 1,218 |
| 40 | 240 | 14 | 176ms | 1,588ms | 9.0× | 1.9e-12 | 1,221 |
| 80 | 480 | 16 | 197ms | 1,689ms | 8.6× | 2.0e-12 | 1,246 |
| 160 | 960 | 18 | 221ms | 1,803ms | 8.2× | 2.9e-12 | 1,273 |

結論:stale-LDLT 預條件把迭代數壓到 Jacobi 的 1/71–1/101(16 桿 ~1/101、160 桿 ~1/71,移除越多比值越小);移除 160 桿(1.6% 元素)仍只要 18 次迭代。三層階梯(≤16 桿走 Tier-1、更多走 Tier-2、巨變走重分解)成立。

## 2. 稀疏屈曲(`exp_sparse_buckling.exe`,out: `sparse_buckling.txt`)

`subspaceSmallest(A=K_ff, Ainv=既有 ldlt, B=−Kg_ff)` 之最小 λ = criticalFactor,vs 引擎 dense GEVP:

| 案例 | nf | λ sparse vs dense | dense | sparse | dense 記憶體估 |
|---|---|---|---|---|---|
| column8(Euler 解析) | 48 | 7.0e-15(vs 解析 2.1e-6=離散誤差) | 0.2ms | 0.2ms | — |
| tower-S | 288 | 2.7e-14 | 20ms | 5.2ms | 1 MiB |
| tower-M | 1,440 | 2.6e-14 | **2,109ms** | **85.6ms(24.6×)** | 32 MiB |
| tower-L | 18,720 | (dense 不可行) | — | 5,193ms,殘差 1.2e-7 | **5,347 MiB** |
| 全拉退化 | — | Kg=0 → 護衛 λ=+inf;引擎 dense 路徑正確回 "no positive buckling eigenvalue" | | | |

結論:S1 可直接用既有 `SparseEigsolver.h` 替換 dense GEVP;tower-L 殘差 1.2e-7(特徵值收斂 tol 1e-11、向量殘差較鬆)→ spec 設殘差驗收欄位。

## 3. 規模基線(`exp_million_dof.exe` 階梯,out: `scale_ladder.txt`)+ solver 比較

| nf | nnz(K) | factor | solve | factor 後記憶體 |
|---|---|---|---|---|
| 18,720 | 0.82M | 1.55s | 17.9ms | 120 MiB |
| 61,560 | 2.67M | 30.6s | 108ms | 585 MiB |
| 186,000 | 7.99M | **522s** | 805ms | 3,257 MiB |
| 389,664 | 16.66M | **3,229s(53.8 min)** | 3,749ms | 10,357 MiB |
| 1,007,964(nf 實測,分解未跑) | — | **未實跑**:4 點擬合 ~O(n^2.46) 外推 ≈ 9.3 小時、~39GB [THEORY:外推] → 直接法天花板,觸發 R13 AMG/matrix-free 條件論述 | | |

`exp_solver_compare.exe`(同一 K_ff/F_ff;out: `solver_compare.txt`):
- SimplicialLDLT ≈ SimplicialLLT(引擎預設已是最佳直接法);SparseLU 慢 2.5×。
- XXL:CG+Jacobi 1,248 迭代/423ms(單發尚可,但 factorize-once 重解 10ms 碾壓);**CG+IncompleteCholesky 反而 2,252 迭代/2,402ms**(Eigen IC 對框架剛度比例病態無效)→ 支撐「stale-LDLT 預條件」而非通用 IC/AMG 的設計;AMG 觸發條件見 R13。

R8 實錘:`build_perf.bat` 連結失敗(5 個 `MITC4ShellElement` LNK2001;源檔清單缺該 .cpp)→ 既有 frame_perf.exe 為舊產物,S1 修。

## 4. P-Delta / N3(`exp_pdelta_convergence.exe`,out: `pdelta.txt`)

懸臂柱 nElem=8、H=1kN、P=f·Pcr;精確解 δ=H(tan(kL)−kL)/(Pk):

| f=P/Pcr | 凍結迭代數(理論) | Aitken | ref vs 精確 | frozen vs ref |
|---|---|---|---|---|
| 0 | 1 | 1 | 1.6e-14 | 逐位 0 |
| 0.3 | 27(26.8) | 15 | 2.6e-7 | 1.5e-13 |
| 0.5 | 46 | 21 | 1.0e-6 | 5.1e-14 |
| 0.9 | 285(306) | **4,742(劣化!)** | 1.7e-5 | 4.7e-13 |
| 0.95 | 569(628) | **57** | 3.7e-5 | 5.4e-13 |
| 1.05 | 發散(物理正確) | 發散 | ref=−53.18 vs 精確 −53.18(後臨界一致) | — |

- 凍結分解(N3)與重組 K_T 參考版同一不動點(≤5.4e-13);迭代數精準符合幾何級數 log(tol)/log(f)。
- 網格 8→16 元素:離散誤差 2.6e-7 → 1.7e-8(~16×)。
- **Aitken 教訓**:f=0.95 砍 10×,但 f=0.9 反而劣化 → S3 spec 用保護式外推(ρ 穩定才外推+失敗回退),裸 Aitken 不可入引擎。
- 發散偵測器在慢發散(增長率 ~1.05/迭代)下未觸發 → spec 用「迭代上限+殘差趨勢窗」雙保險。

## 5. Tension-only(`exp_tension_only.exe`,out: `tension_only.txt`)

- caseA(側載 X 斜撐):2 迭代收斂,壓桿停用、拉桿 N=−40,045(拉力✓);**收斂解 vs「該桿直接從 member 清單省略」= 0.000e+00 逐位元相等**。
- caseB(豎載+側載掃描 H=100..3e4):全部 ≤3 迭代收斂,**無 flip-flop 視窗(誠實負結果)**;H=1k..3k 區間有「停用→再啟用→穩定」路徑(再啟用機制被運動)。狀態哈希守門留設計(LCP 視角下循環理論上可能)。

## 6. N4 動量繼承(`exp_dynamic_inherit.exe`,out: `dynamic_inherit.txt`)

塔(2,1,3) nf=108、階躍載重、C=0、dt=T1/200、第 240 步移除一斜撐:
- **全基底模態繼承 vs 全系統 Newmark:1.97e-12**(數學等價性成立);事件前能量守恆漂移 6.6e-11(Newmark 平均加速度保能);事件釋放能量 ΔE=−4.69 可量化(物理數據)。
- **截斷誤差(關鍵發現)**:m=5/10/20 → 39% 誤差、m=40 → 7.3%;靜力修正(mode-acceleration)只砍半 → **S2 spec 改推 load-dependent Ritz vectors(Wilson)取代純特徵模態**。
- **動量帳全精確**:碎塊質量 FE vs `FragmentCluster` 閉式 relDiff=0;線動量平移測試 0;橫向角動量 L_FE vs I_cluster·ω relDiff=**0.000e+00**;own-axis 例外:FE 帶截面極慣量項 205 t·mm² vs 桿模型 0(細長桿可忽略,spec 註記)。
- warm-start:nev=10、nf=288:cold 10 迭代 → warm 7(省 30%,中等效益,特徵值一致 7.5e-14)。

## 7. FSD 尺寸優化(`exp_size_opt.exe`,out: `size_opt.txt`)

經典 10-bar truss(stress-only):24 迭代收斂 **1593.1632 lb**、斷面 [7.94, 0.1, 8.06, 3.94, 0.1, 0.1, 5.74, 5.57, 5.57, 0.1] in² = 文獻經典解(交叉引用見 WS_H);非 min-gauge 桿全應力偏差 3.4e-11。註:尾段 maxBendShare 飆高是 min-gauge 桿 N≈0 的除零假象;受載桿彎曲份額 2–6%(桁架近似成立)。

## 8. BESO(`exp_beso_truss.exe`,out: `beso.txt`、`beso_topology.txt`)

12×6 ground structure(306 桿)→ vol 30%:29 迭代完成(終態 88 桿,圖:`beso_topology.png`,呈懸臂 Michell 式上下弦+斜腹桿定性形態);連通守門擋不住彎曲機構 → **逐桿 factor 試移除**(小模型可負擔;引擎版將改用 N1 capacitance 奇異偵測 = 同一機制零額外成本)。**尾段 compliance 暴增 52×**(hard-kill 低體積過衝)→ spec 需:敏感度歷史平均(Huang & Xie 標準)+ compliance 跳變停機。

## 9. Elastica(`elastica_check.py`,out: `elastica.txt`)

Shooting 法表(雙容差 1e-10/1e-12 十位數一致;線性極限檢查 1.1e-9):α=PL²/EI=1 → δv/L=0.3017207738;α=10 → δv/L=0.8106090249(文獻交叉驗證見 WS_F)。= 未來 S9 CR 的 oracle 數據表。

## 10. CLI 吞吐(`cli_throughput.py`,out: `cli_throughput.txt`)

| 模型 | DOF | 端到端中位 |
|---|---|---|
| trivial(行程開銷) | 12 | **6.7ms** |
| tower | 1,620 | 26ms |
| tower | 19,500 | 2.06s |
| tower | 107,502 | 142.7s |

結論:CLI 橋對 GH 人工觸發在 ~20k DOF 內完全可用;100k+ 受 factor 支配 → 催生 **J1.5 daemon 模式 CLI**(常駐行程、同模型多解、吃 PreparedSystem/ReSolve 紅利)或 C API 長期路。
