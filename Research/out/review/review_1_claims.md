現在我有完整資料進行審查。以下是完整審查報告。

---

# 審查報告:宣稱紀律與過度宣稱

## 裁決

**PASS(有若干 MINOR/MAJOR 問題,無 CRITICAL,可在下輪文檔修訂中一起修)**

---

## 發現(按嚴重度排序)

---

### [MAJOR] PERFORMANCE_BASELINE.md §2:Tier-1 baseline factor 時間數字與實驗輸出不一致

**位置**:PERFORMANCE_BASELINE.md 第 25 行,表格 Tier-1 基準欄

**問題**:文檔寫「vs fresh(1.6–1.8s)」,但 `incremental_xxl.txt` 第 1 行原始輸出顯示 `factorMs=1668.8`(即 1.669s),各步的 `freshMs` 實測值範圍在 1385ms 到 1803ms 之間,並非穩定在「1.6–1.8s」區間——最低有 1385ms(1.4s 量級)、最高達 1803ms。若用來算加速倍率,分母不穩定。同一文件的 `[baseline]` 行寫 `factorMs=1668.8`,而正文 §2 表頭寫「vs fresh(1.6–1.8s)」略微誇大了最低端。

**影響**:表中 step=50 的 fresh=1435ms 對應加速比 `1435/82=17.5×` 是正確的,但 step=1 用 1687ms 算 31×也正確。問題在於「1.6–1.8s」範圍漏掉了 fresh 在較小 R 時的較短計時區間(如 step=23 freshMs=1467ms),雖在 1.4–1.8s 更寬範圍內,文件說法仍基本成立,但邊界略緊。

**建議修法**:將「vs fresh(1.6–1.8s)」改為「vs fresh(1.4–1.8s,timer variance)」或直接寫「vs fresh ~1.6s 中位」並加 note 說明有 ±15% 量測變動。

---

### [MAJOR] WS_R2_experiments.md §1:Tier-2 baseline factor 時間與 PERFORMANCE_BASELINE.md 不一致

**位置**:WS_R2 第 8 行(baseline row)vs PERFORMANCE_BASELINE.md §2 表格

**問題**:
- WS_R2 §1 寫 `Baseline:factor 1,669ms`
- PERFORMANCE_BASELINE.md §2 表中 Tier-1 行寫「vs fresh(1.6–1.8s)」
- 原始輸出 `incremental_xxl.txt` 第 1 行:`factorMs=1668.8`(即 1.669s,對應 WS_R2 的 1,669ms)
- 但 `scale_ladder.txt` 同模型(nf=18,720)記錄的是 `factorMs=1546`(1.55s)

**具體衝突**:同一個 nf=18,720 XXL 塔,`scale_ladder.txt` 測得 factor=1546ms,`incremental_xxl.txt` 測得 1669ms。WS_R2 用 1,669ms 當 baseline,PERFORMANCE_BASELINE.md §1 表格卻用 1.55s。這是兩次不同執行的量測差異,但兩份文件交叉引用時讓 factor 時間出現 8% 差距。

**建議修法**:在 PERFORMANCE_BASELINE.md §1 腳注加一句:「同一 nf=18,720 模型在 exp_incremental_refactor 獨立執行中 factor=1669ms;差異屬正常 OS 調度 ±10%」,消除讀者疑惑。

---

### [MAJOR] KARAMBA3D_ROADMAP.md §6:Tier-2「Jacobi 對照 1273」與實驗輸出不完全一致

**位置**:ROADMAP §6 PCG 判定行:「stale-LDLT 預條件 12–18 迭代 vs Jacobi ~1250 `[VERIFIED]`」

**問題**:原始輸出 `incremental_xxl.txt` 顯示 Jacobi 迭代數為 1218、1221、1246、1273 四個值,對應 nb=16/40/80/160 桿移除。文件引用了 1273(最大值),WS_R2 的表格也正確列出了全部四個值。但 ROADMAP 正文中寫「Jacobi ~1250」與 WS_R2 中「Jacobi 對照迭代 1,273」(表格右欄,只顯示 160 桿那行)數字有差異:1250 是近似數、1273 是實際 nb=160 那行的最大值。

**分析**:僅 ROADMAP 正文用了「~1250」的近似值,WS_R2 的原始表格正確列出了所有四個數字。近似值本身沒有錯(1250 是 1218–1273 的合理均值近似),但讀者可能拿 ROADMAP 的「1273」(執行摘要 §1 提到)與正文 §6 的「~1250」比對,造成前後不一致感。

**建議修法**:ROADMAP §6 PCG 行改為「vs Jacobi 1218–1273(nb=16..160)」或「vs Jacobi 1200–1280」,與 WS_R2 表格一致。

---

### [MAJOR] PERFORMANCE_BASELINE.md §4:SparseLU「慢 2.5×」計算來源不明

**位置**:PERFORMANCE_BASELINE.md 第 42–43 行:「SparseLU 慢 2.5×」

**問題**:原始輸出 `solver_compare.txt` 顯示:
- XXL(nf=18,720):`SimplicialLDLT factorMs=1364.4`、`SparseLU factorMs=3491.7`
- 計算倍率:3491.7/1364.4 = **2.56×**
- tower-M(nf=1440):`SimplicialLDLT factorMs=5.4`、`SparseLU factorMs=27.8`
- 計算倍率:27.8/5.4 = **5.1×**

文件只說「慢 2.5×」,但這是 XXL 模型的數字,而在 tower-M(小模型)SparseLU 慢了 5.1×。「2.5×」本身對 XXL 是正確的,但沒有指出規模依賴性,對讀者可能誤導為一般規律。

**建議修法**:加一句「(XXL 模型;小模型差距更大,tower-M 約 5×)」。

---

### [MINOR] KARAMBA3D_ROADMAP.md §1 執行摘要:「批量 160 桿 18 迭代 8.2×(Jacobi 對照 1273)」的括號來源

**位置**:ROADMAP §1 第 2 點 N1 描述行

**問題**:括號內「Jacobi 對照 1273」只取了 nb=160 那行的最大值,沒說明這是最壞 case。建議加「(nb=160 最壞 case)」或改為「Jacobi 1200–1280+」。這是次要問題,不影響正確性。

---

### [MINOR] PERFORMANCE_BASELINE.md §6 動力(N4):截斷誤差描述不精確

**位置**:PERFORMANCE_BASELINE.md 第 58 行:「截斷誤差大(m=40/108→7.3%)」

**問題**:原始輸出 `dynamic_inherit.txt` 顯示:
- m=5: errPlain=3.886e-01(38.9%)、errModeAccel=1.945e-01(19.5%)
- m=10: errPlain=3.885e-01(38.9%)、errModeAccel=1.950e-01(19.5%)
- m=20: errPlain=3.892e-01(38.9%)、errModeAccel=1.947e-01(19.5%)
- m=40: errPlain=7.300e-02(7.3%)、errModeAccel=5.056e-02(5.1%)

ROADMAP §10 誠實邊界寫「N4 截斷誤差顯式報告」,WS_R2 §6 寫「m=5/10/20 → 39% 誤差、m=40 → 7.3%」,但 PERFORMANCE_BASELINE.md §6 只說「截斷誤差大(m=40/108→7.3%)」,讓讀者誤以為 m=40 已經夠好。實際上 m=5/10/20 全都是 39%,代表截斷問題在低 m 時非常嚴重,7.3% 只在 m=40 時才達到,且 mode-acceleration 修正只砍半(5.1%)。

PERFORMANCE_BASELINE.md 也沒有提到 m=5/10/20 都是 39% 的關鍵事實——只有 ROADMAP 誠實邊界提到了,這使得讀者若只看 BASELINE 文件可能低估截斷風險。

**建議修法**:PERFORMANCE_BASELINE.md §6 改為「截斷誤差:m=5/10/20 均 ~39%、m=40 降至 7.3%(full-basis nf=108);mode-acceleration 約砍半;→ spec 改推 load-dependent Ritz」,並保持與 WS_R2 §6 一致的完整描述。

---

### [MINOR] WS_R2_experiments.md §1:Tier-1「5.7e-13(無漂移問題)」的括號措辭可能誤導

**位置**:WS_R2 §1 第 12 行:「50 桿連續移除:relErr 全程 7.7e-14 → 5.7e-13(無漂移問題)」

**問題**:原始輸出確認:step=50 的 relErr=5.725e-13,而 step=1 是 7.707e-14。誤差確實從 step=1 到 step=50 增長了約 7.4 倍。括號「無漂移問題」的說法屬實(仍在機器精度量級),但單就括號本身讀起來像「誤差沒有增長」,而實際上有穩定但可觀的增長趨勢。這是輕微的措辭問題。

**建議修法**:改為「50 桿連續移除:relErr 全程 7.7e-14 → 5.7e-13(仍機器精度量級,未發散)」。

---

### [MINOR] PERFORMANCE_BASELINE.md §5:CLI「107,502 DOF」數字

**位置**:PERFORMANCE_BASELINE.md 第 53 行:「107,502 DOF / 142.7s」

**問題**:原始輸出 `cli_throughput.txt` 顯示 `dof=107502 medianMs=142746.3`(即 142.746s)。PERFORMANCE_BASELINE.md 寫「142.7s」,正確(四捨五入)。但 WS_R2 §10 寫「142.7s」也正確。這個數字本身沒問題,僅備注 `incremental_xxl.txt` 的 baseline factor=1668ms 是兩次不同執行的量測,一致性已在上面 [MAJOR] 中說明。

此條目不需修正,純粹是確認數字正確的記錄。

---

### 五項宣稱紀律專項核對結果

**(a) 「優於 Karamba」「獨創/未見先例」主張標記與依據**

- ROADMAP §5 N-track 表:每個算法均有「先行技術(誠實定位)」欄,措辭謹慎。
  - N1:明確說「非數學新算法」,工程整合 novelty。`[VERIFIED]` 附實驗數字。**合格**。
  - N2:明確說「方法整合 novelty;勿自稱 fail-safe」。**合格**。
  - N3:明確說「與 factorize-once 的架構整合,非新算法」。**合格**。
  - N4:「未見直接先例 `[UNKNOWN 級定位,WS_N]`」——謹慎用 UNKNOWN 級定位,不裸宣稱。**合格**。

**(b) [VERIFIED] 數字是否對應實驗輸出(抽查 8 個數字)**

| 文件中的宣稱 | 原始輸出數字 | 吻合? |
|---|---|---|
| Tier-1 單桿 31× | `incremental_xxl.txt` step=1: speedup=31.0 | ✅ 完全吻合 |
| relErr 7.7e-14 | step=1: relErr=7.707e-14 | ✅ 完全吻合 |
| 漂移 1.5e-15(或 1.46e-15) | restore行: relErrVsBaseline=1.457e-15 | ✅ 完全吻合 |
| Tier-2 批量 160 桿 18 迭代 8.2× | nb=160: staleIters=18, speedup=8.2 | ✅ 完全吻合 |
| 屈曲 nf=1440 dense 2109ms sparse 85.6ms 24.6× | sparse_buckling.txt: tDenseMs=2109.0, tSparseMs=85.6 | ✅ 完全吻合 |
| 全基底等價 1.97e-12 | dynamic_inherit.txt: fullBasisRelErr=1.965e-12 | ✅ 吻合(文件寫 1.97e-12,原始為 1.965e-12,四捨五入正確) |
| FSD 24 迭代 1593.1632 lb | size_opt.txt: iter=24 weightLb=1593.1632 | ✅ 完全吻合 |
| BESO 29 迭代 vol 30% | beso.txt: iter=29 volFrac=0.2984 | ✅ 吻合(0.2984 ≈ 30%) |

所有抽查數字均與原始輸出匹配,無虛構數字。

**(c) 外推/未實跑的數字是否清楚標注**

- 1M DOF:ROADMAP §1 寫「~9–10h/30GB+ `[THEORY:外推]`」,PERFORMANCE_BASELINE.md §1 在表格中寫「**未實跑**;~O(n^2.4) 外推 `[THEORY:外推]`」。標注清楚。✅
- `scale_ladder.txt` 顯示 nf=389,664 那行在輸出截斷前只有 `scale-begin` 記錄(18:56:40 後就截斷了),推測 389k 的結果可能也未完成——但文件在 PERFORMANCE_BASELINE.md 表格中「389,664」那行已標「(回填中:scale_ladder.txt)」,誠實標注。✅

**(d) 有沒有把「研究 agent 查證」當絕對事實而缺 URL 的情況**

- ROADMAP §2 對 Karamba 的描述標明「`[LIT:WS_A1 §5]`」「`[LIT:WS_A1]`」,但 WS_A1 本身是否有具體 URL 需另查文獻文件。以下屬待查:「Karamba 免費版 ≤20 梁元素」(`[LIT:WS_A2]`)、「殼 ≤50」(`[LIT:WS_A2]`)——這些數字是否有官方手冊截圖或 URL 支撐,在本次審查範圍的三份文件內看不到 WS_A1/WS_A2 的 URL 展開。**這不算宣稱紀律違反**(文件已正確用 `[LIT:WS_A2]` 標記),但若 WS_A1/WS_A2 文件本身沒有具體 URL,則這些 Karamba 規格宣稱只靠研究 agent 查證而無原始 URL 可溯源,偏向 `[PENDING]` 而非 `[LIT]` 等級。此條目屬輕微風險,建議確認 WS_A1/WS_A2 有具體 URL。

**(e) 誠實邊界段是否涵蓋所有已知限制**

逐項核對 ROADMAP §10:

| 限制 | 文件中有無覆蓋 |
|---|---|
| 模態截斷 39% | WS_R2 §6 有:「m=5/10/20 → 39% 誤差」,ROADMAP §10 有「N4 截斷誤差顯式報告」。但 PERFORMANCE_BASELINE.md §6 描述不完整([MINOR] 已提) |
| BESO 過衝 52× | WS_R2 §8:「尾段 compliance 暴增 52×」,ROADMAP §5 N2 行有「過衝/守門教訓」,ROADMAP §10 有「BESO hard-kill 低體積過衝需停機準則」。✅ 涵蓋 |
| Aitken 劣化(f=0.9 時 4742 迭代) | WS_R2 §4、ROADMAP §5 N3 行:「保護式外推教訓(裸 Aitken f=0.9 劣化)」,ROADMAP §10 有對應。✅ 涵蓋 |
| 1M DOF 未實跑 | ROADMAP §10 有「1M DOF 數字是外推非實測」。✅ 涵蓋 |
| 同機速度對標未做 | ROADMAP §10 有「同機速度對標未做([PENDING:S6 後])」。✅ 涵蓋 |

誠實邊界段涵蓋所有要求項目,無遺漏。

---

## 沒問題的部分

[VERIFIED] 標記的 8 個抽查數字全部與 `Research/out/*.txt` 原始輸出逐位元或四捨五入匹配,無虛構數字;novelty 宣稱全部附先行技術定位且措辭謹慎(N1/N2/N3 不自稱新算法,N4 用 UNKNOWN 級定位);外推數字均清楚標注 `[THEORY:外推]`;誠實邊界段涵蓋全部五項要求限制;[PENDING] 標注用法正確。