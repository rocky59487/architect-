# FrameCore × Karamba3D — 對標研究 + 獨創算法路線 + 開發主線(2026-06-10 研究輪定稿)

> 本檔 = 研究輪主報告。細節:`docs/research/WS_*.md`(文獻查證)、`docs/research/WS_R2_experiments.md`
> (實驗數據,全部可重跑)、`docs/specs/S1..S4_*.md`(交接級實作規格)、`docs/specs/S5_S11_skeletons.md`、
> `docs/PERFORMANCE_BASELINE.md`(速度基線)。基準 commit:`733833d`。

## 0. 宣稱紀律(全文適用)

`[VERIFIED]` = 本輪 scratch 實驗 + oracle 通過,可重跑;`[LIT:…]` = 有來源;`[THEORY]` = 有推導
(含標明的外推);`[PENDING:…]` = 待驗;`[UNKNOWN]` = 查無/不知;`[NEW CODE]` = spec 中
無原型背書、實作時須新寫並由指定 oracle 把關的段落(S1 F 增量、S2 Ritz 生成、S3 保護式
外推與發散偵測)。
「優於 Karamba3D」「新穎」類主張禁止裸宣稱;novelty 一律附先行技術定位(WS_N)。
**審核狀態**:本報告與 specs/WS_R2/baseline 已過五面向對抗式查核(宣稱紀律/S1-S2 與
S3-S4 可實作性/novelty 一致性/跨文檔數據一致性;零 CRITICAL,全部 MAJOR/MINOR 已修;
審查原文 `Research/out/review/`;抽查 8 個 `[VERIFIED]` 數字全數與原始輸出吻合)。

## 1. 執行摘要

1. **對標事實成立但要重框**:Karamba3D 是設計探索工具(閉源商業、免費版 ≤20 梁元素
   `[LIT:WS_A1 §5]`),其分析核心在多個維度**落後 FrameCore 既有能力**:無材料非線性/塑鉸/
   倒塌、無反應譜/時程、殼用忽略橫剪的 TRIC 三角形(四邊形拆兩半)`[LIT:WS_A1]`。
   FrameCore 落後處:幾何非線性(Karamba 有 ThII+NR/弧長 WIP)、EC3 截面優化(含 LTB)、
   BESO、GH 生態(它的主場)。
2. **獨創路線定案**:「互動式結構重分析 × 連續動力倒塌」——N1 ReSolve 三層階梯
   (單桿移除 31×、批量 stale-PCG 8–10×、機構=capacitance 奇異 `[VERIFIED]`)+
   N4 動量繼承動力倒塌(跨事件繼承 1.97e-12、碎塊動量帳 0 差 `[VERIFIED]`)+
   N3 凍結分解 P-Delta + N2 倒塌韌性 BESO(概念親緣 fail-safe TO,Jansen 2014;整合 novelty)。
   先行技術全部查證、措辭已定(§5)。
3. **使用者新增六項判定完成**(§6):PCG 入主線(=N1 Tier-2)、模態疊加+動量繼承入主線(=N4/S2)、
   AMG/SA-AMG/matrix-free/全隱式暫態 = 研究定位+觸發條件(R13;實測 Eigen IC-PCG 在框架上
   劣於 Jacobi `[VERIFIED]`,通用預條件不是免費午餐)。
4. **直接法天花板實測**:186k DOF factor 522s/3.3GB、390k DOF factor 3,229s/10.4GB
   `[VERIFIED]`;1M DOF(nf=1,007,964 實測)4 點擬合外推 ~9.3h/~39GB `[THEORY:外推]`
   → 引擎的價值主張是「互動規模(≤~60k DOF)上的重解速度與正確性文化」,不是百萬 DOF 競賽。
5. **主線 S1–S11 定稿**(§7),S1–S4 已有交接級 spec(新 agent 可直接動工)。

## 2. Karamba3D 現況(官方手冊查證,詳 WS_A1/WS_A2)

- 核心:Timoshenko 梁(Rubin 冪級數)、TRIC 三角殼(**忽略橫向剪變形**;四邊形自動拆兩三角)。
- 分析:Th.I、Th.II(迭代 P-Δ 級)、Analyze Nonlinear **WIP**(NR/弧長/動力鬆弛,殼可能不收斂,
  官方自標)、LaDeform 找形(**不輸出內力**)、線性屈曲、模態。**無**反應譜、時程、材料非線性、
  塑鉸、pushover、漸進倒塌。
- 優化:OptiCroSec(EC3 + LTB + 離散截面表;無 AISC)、BESO for Beams/Shells(殼用 soft-kill
  減厚+空間濾波)、RC 配筋(Marti sandwich,分析仍線彈性)。
- Tension/Compression Eliminator:**迭代式**,soft-kill(極小剛度)至收斂。
- 載重組合:字串規則展開;**無 envelope 元件**。
- 授權:閉源商業(PRO/EDU/LAB),試用 ≤20 梁元素、殼約 ≤50(WS_A2)。

## 3. FrameCore 現況快照(基準 733833d)

純 C++17/Eigen、零依賴、POD 公開 API;稀疏 SimplicialLDLT(**靜力解早已稀疏**——roadmap
「換稀疏求解器」前提不成立,真正稠密殘留 = 屈曲 GEVP 與 modal 預設路徑);factorize-once
`PreparedSystem` + fingerprint 防呆;MITC4 四邊形殼(OpenSees ~1e-10);8 階段線性套件
(組合/包絡/影響線/沉陷/模態/屈曲/反應譜/模態疊加動態);**漸進倒塌主線全套**(移除/連通/
碎塊 Chaos 交接/LSP 驅動器/殼 vM/event-to-event 塑鉸);四腿 gate(F1–F33、UE 34、audit 62、
OpenSees strict)。已知缺口:幾何非線性、優化、生態系(零 DLL/序列化;`frame_cli` 文字格式
是現成交換格式)、`build_perf.bat` 連結失敗(R8 `[VERIFIED]`,S1 修)。

## 4. 對標矩陣(✅=已有 ⛔=無 🔶=部分/限制;詳 WS_A1)

| 能力 | FrameCore 今天 | Karamba3D v3 | 主線後 |
|---|---|---|---|
| 梁(Timoshenko) | ✅(+release/塑鉸) | ✅ | — |
| 殼 | ✅ MITC4 四邊形(剪變形✓) | 🔶 TRIC 三角、忽略橫剪 | S8 QM6/DKQ、S11 MITC9i |
| 載重組合/包絡 | ✅ combine+envelope | 🔶 組合有、無 envelope | — |
| 影響線/沉陷 | ✅ | ⛔/🔶(Gap Load 可做影響線) | — |
| 模態 | ✅(opt-in 稀疏) | ✅ | — |
| 線性屈曲 | ✅(dense;S1 轉稀疏) | ✅ | S1 |
| 反應譜/模態疊加時程 | ✅ SRSS/CQC+Newmark | ⛔ | — |
| 二階 P-Δ | ⛔ | ✅ ThII(NoTenNII/NoComNII 選項與 tension-only 行為部分重疊) | S3(雙路徑) |
| 幾何非線性(大位移) | ⛔ | 🔶 WIP(自標收斂限制) | S9 CR |
| 材料非線性/塑鉸 | ✅ event-to-event(誠實標 sequential linear) | ⛔ | S10 N-M 互動 |
| 漸進倒塌 | ✅ LSP 驅動器+碎塊交接 | ⛔ | S2 升級動力版 |
| **動力倒塌+動量交接** | ⛔(本輪驗證) | ⛔ | **S2(獨創主打)** |
| **互動增量重解** | 🔶 factorize-once | ⛔(未見等價 API) | **S1 ReSolve(獨創主打)** |
| Tension/Compression-only | ⛔ | ✅ 迭代 Eliminator | S4 |
| 截面優化 | ⛔(D/C 快篩有) | ✅ EC3+LTB | S5(FSD;**初版無 LTB,誠實標**) |
| 拓撲優化 | ⛔ | 🔶 BESO 梁(硬移除)/殼(soft-kill 減厚,非元素移除,WS_A1 §3.3) | S7(+N2 韌性約束版;概念親緣 fail-safe TO(Jansen 2014),差異=離散框架+LSP 評估器,WS_N §N2) |
| GH 生態 | ⛔ | ✅(主場) | S6 CLI 橋→daemon→C API |
| 開源 | ✅(目標 MIT) | ⛔ 閉源商業 | — |
| 驗證文化 | ✅ 四腿 gate+OpenSees | `[UNKNOWN]`(無公開 oracle 套件) | gate 持續成長 |

速度對標:同機實測延後到 S6 後(免費版 ≤20 梁元素,只能對微模型,屆時誠實設計實驗);
本輪建立絕對基線(`PERFORMANCE_BASELINE.md`)。

## 5. 獨創算法路線(N-track;先行技術定位 = WS_N,實驗 = WS_R2)

**路線宣言**:「互動式結構重分析 × 連續動力倒塌」。每個新算法都有舊路徑當精確 oracle
(雙路徑互驗),「不錯誤的前提下逐漸採用」由導入階梯(§8)機器驗證。

| | 內容 | 實驗結果 `[VERIFIED]` | 先行技術(誠實定位) | 落地 |
|---|---|---|---|---|
| **N1** | ReSolve 三層階梯:Woodbury rank-≤6/桿 精確重解 → stale-LDLT 預條件 PCG → 全重分解;capacitance 奇異=機構偵測 | 單桿 31×/relErr 7.7e-14;50 桿序列無漂移;移除+恢復 R=600 漂移 1.5e-15;批量 160 桿 18 迭代 8.2×(Jacobi 對照 1273) | SMW 重分析=Akgün et al. 2001;rank-k downdate=Davis & Hager(CHOLMOD,授權不相容 MIT);**我們=零依賴三層自動階梯×倒塌/優化驅動器×互動編輯的工程整合,非數學新算法** | S1(spec 就緒) |
| **N2** | 倒塌韌性 BESO:LSP 倒塌驅動器當約束評估器 | BESO 基礎設施驗證(29 迭代到 vol 30%;過衝/守門教訓入 spec) | fail-safe TO=Jansen 2014(連續介質/SIMP/靜態最壞柔度);**我們=離散框架×序列線性倒塌×Stable/Collapsed 終點,方法整合 novelty;勿自稱 fail-safe** | S7 |
| **N3** | 凍結分解 P-Delta(pseudo-load 迭代於既有 LDLT) | 與重組 K_T 參考同不動點 ≤5.4e-13;迭代數=幾何級數理論;保護式外推教訓(裸 Aitken f=0.9 劣化) | = Wilson & Habibullah 1987 / 業界 nodal-shear 法;**我們=與 factorize-once/N1 零摩擦的架構整合,非新算法** | S3(spec 就緒) |
| **N4** | 動量繼承動力倒塌:模態空間跨事件狀態繼承 + 碎塊帶初速 Chaos 交接 | 全基底等價 1.97e-12;動量帳 0 差;能量帳可量化;截斷誤差大(m=40/108→7.3%)→ **spec 改推 load-dependent Ritz** | 模態疊加/Ritz/投影=教科書(Wilson);跨拓撲事件繼承×遊戲物理動量交接的整合「未見直接先例」`[UNKNOWN 級定位,WS_N]` | S2(spec 就緒) |

## 6. v3 新增技術項判定(使用者點名+自判)

| 項 | 判定 | 一句話依據 |
|---|---|---|
| PCG | **入主線**(N1 Tier-2) | stale-LDLT 預條件 12–18 迭代 vs Jacobi 1218–1273(nb=16..160)`[VERIFIED]`;Eigen 內建零依賴 |
| 隱式/matrix-free Krylov | Tier-2 隱式運算子=特例已入;全 EBE 研究定位 | 觸發條件未到(§R13) |
| AMG / SA-AMG | **研究定位,不入主線** | 手寫 SA-AMG(6DOF/節點/旋轉 DOF/剛體近核)數千行高風險;外部 lib 違零依賴;實測 Eigen IC-PCG 反而劣於 Jacobi(框架剛度病態)→ 通用預條件成本被低估;觸發=「>50 萬 DOF 互動」需求成真 |
| 模態疊加 | 已存在 → 延伸入主線(N4) | `solveModalStepResponse` 既有;新工作=跨事件繼承 |
| 動量繼承 | **入主線**(N4 核心) | 修掉「碎塊零初速」既有限制;動量帳 oracle 乾淨 `[VERIFIED]` |
| 模態 warm-start(自判補) | 入主線(N4 內) | cold 10→warm 7 迭代 `[VERIFIED]`(中等效益,誠實標) |
| Lanczos shift-invert(自判補) | 研究定位 | subspace iteration 對少模數已夠(屈曲 24.6× `[VERIFIED]`) |
| 全系統隱式暫態(自判補) | 研究定位 | N4 模態路線以遠低成本覆蓋遊戲需求 |

## 7. 開發主線(每階段:獨立 oracle + 四腿 gate + commit/push;導入階梯 §8 適用)

| 階段 | 內容 | Spec | Gate 成長(預測) |
|---|---|---|---|
| **S1** | ReSolve 階梯 + 稀疏屈曲 + 修 build_perf + 效能基線 | `specs/S1`(就緒) | F34–36、UE 37、audit 68 |
| **S2** | N4 動量繼承動力倒塌(Ritz 基底、FragmentCluster+vel/angVel) | `specs/S2`(就緒) | F37–39、UE 40、audit 72 |
| **S3** | P-Delta 雙路徑(N3 凍結 + Wilson 參考) | `specs/S3`(就緒) | F40–41、UE 41、audit 75 |
| **S4** | Tension-only(ReSolve 內迴圈) | `specs/S4`(就緒) | F42–43、UE 42、audit 77 |
| S5 | FSD 尺寸優化(無 LTB 誠實標) | 骨架 | F44+ |
| S6 | GH MVP:CLI 橋 → daemon 模式(J1.5)+ .gha/Yak | 骨架 | CLI 黃金檔測試 |
| S7 | BESO + **N2 韌性約束版** | 骨架 | F4x ×2 |
| S8 | 殼:QM6 opt-in 膜、DKQ | 骨架(P5 數據) | patch+OpenSees |
| S9 | Co-rotational(elastica oracle 表已備) | 骨架 | F4x ×2 |
| S10 | N-M 互動塑鉸(選做) | 骨架 | F4x |
| S11 | MITC9i(殿後,使用者備註) | 骨架 | F4x ×2 |

並行線:生態系(Twinmotion/CI/打包;WS_K)於 S6 後;C6–C8 可視化於 S3–S4 間;UE 視覺層
U1–U21 於 S6 後(S1+S2 = 「玩家拆桿即時重算+真動力倒塌回放」成立)。
排序理由:S1 是所有迭代驅動器的地基(依賴鏈);S2 是獨創主打且資產齊備(價值前置);
S3/S4 小而 oracle 嚴(快速積累);S6 提前讓外界可用;S11 殿後(使用者定調+9 處修改成本)。

## 8. 導入階梯(「不錯誤的前提下逐漸採用」的機制)

① scratch 原型 + oracle(本輪完成:N1/N3/N4 + 屈曲/TO/FSD/BESO)
② 進引擎:**opt-in 旗標 + 雙路徑互驗**(新舊算法同跑,audit check 鎖一致;Tier-2 類迭代解
   用容差互驗並文檔明標非逐位元)
③ gate 固化(F 編號 / UE 測試 / audit checks,各 spec ⑦ 節)
④ 數據證明(效能驗收欄,各 spec ⑧ 節)後**轉預設**;舊路徑降級 fallback + 永久 oracle。
任何數值疑慮 → fallback 路徑永遠正確(rebaseline / fresh factor / dense GEVP / 全基底)。

## 9. 授權與開源策略

建議 **MIT**:Eigen MPL2(`EIGEN_MPL2_ONLY` 已設)允許不改源靜態連結;CHOLMOD(LGPL/GPL)
不可引入(N1 因此走零依賴 Woodbury 而非 cholmod_updown——劣勢變特色);對標賣點 =
「開源 + 可重跑 oracle 套件」vs 閉源商業(Karamba 免費版 ≤20 梁元素、殼約 ≤50——後者由官方範例檔行為反推,WS_A1 §5)。`[LIT:WS_A2、WS_N]`

## 10. 誠實邊界(對外文檔一律隨能力宣告)

既有邊界全部沿用(LSP ±30% 級、塑鉸非真彈塑性、D/C 非規範檢核、殼 facet…)。本輪新增:
Tier-2 容差級非逐位元;N4 截斷誤差顯式報告(`truncationResidual`)、事件整步觸發、單向 Chaos
交接;P-Delta = Th.II 線性化(N 凍結);FSD 靜不定非最優保證、初版無 LTB;BESO hard-kill
低體積過衝需停機準則;1M DOF 數字是外推非實測;同機速度對標未做(`[PENDING:S6 後]`)。

## 11. 風險登記冊

R1 −Kg 半正定(護衛+退 dense,`[VERIFIED]` 退化案例);R2 直接法天花板(186k=522s、390k=3,229s/10.4GB 實測、
1M 外推 ~9.3h/~39GB → 價值主張定位互動規模);R3 CLI 吞吐(6.7ms 開銷/19.5k 2.06s/107k 143s
`[VERIFIED]` → J1.5 daemon);R4 CR×塑鉸方向耦合(S10 在 S9 後);R5 外部後端 vs 零依賴
(CHOLMOD 授權卡死;MKL 留 opt-in 議題);R6 MIT(§9);R7 宣稱紀律(§0);R8 build_perf
連結失敗 `[VERIFIED]`(S1 修);R9 novelty 措辭已由 WS_N 定稿(§5 第四欄照抄);
R10 N1 漂移(實測機器精度,仍設 rebaseline 哨兵);R11 Tier-2 確定性(同輸入可重現、
與直接解差 tol 級);R12 N4 截斷(Ritz+殘差哨兵+全基底 fallback);**R13 AMG 觸發條件**:
「>50 萬自由 DOF 的互動分析」成為真實需求時立項,且屆時重新評估外部 lib(授權)vs 自寫
(數千行)vs matrix-free EBE——在那之前 stale-LDLT 預條件覆蓋互動重解需求。

## 12. 下一輪交接(給實作 agent 的入口)

讀順序:本檔 §5–§8 → `specs/S1_resolve_ladder.md`(十項齊,含可移植原型路徑)→
`research/WS_R2_experiments.md` §1–§3。鐵則不變:四腿 gate 全綠才 commit;
`Research/` 不入 gate;每階段結束更新 memory 與本檔狀態欄。
