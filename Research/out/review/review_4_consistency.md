所有檢查已完成。現在輸出審查報告。

---

# 審查報告:跨文檔數據一致性

## 裁決

**PASS(有 3 個 MINOR 問題,無 CRITICAL/MAJOR)**

## 發現(按嚴重度排序)

- **[MINOR]** `docs/PERFORMANCE_BASELINE.md` §2 vs `docs/research/WS_R2_experiments.md` §1 + `docs/specs/S1_resolve_ladder.md` §④:
  同一數字(restore 後對 baseline 的漂移)在不同文檔間使用不同有效位數。
  原始輸出 `incremental_xxl.txt` = `relErrVsBaseline=1.457e-15`。
  WS_R2 和 S1 spec 均寫 **1.46e-15**(3 位有效數字);PERFORMANCE_BASELINE 寫 **1.5e-15**(2 位有效數字)。
  兩者皆為合法捨入,但跨文檔不一致,讀者可能誤以為是不同實驗的數字。
  建議:統一為 `1.46e-15`(3 sig fig,與 WS_R2/S1 保持一致)。

- **[MINOR]** `docs/KARAMBA3D_ROADMAP.md` §11 風險登記冊 R3 vs `docs/research/WS_R2_experiments.md` §10 / `docs/PERFORMANCE_BASELINE.md` §5:
  CLI 19.5k DOF 端到端時間。原始輸出 `cli_throughput.txt` = `medianMs=2059.4ms = 2.06s`。
  WS_R2 和 PERFORMANCE_BASELINE 均寫 **2.06s**;ROADMAP R3 寫 **2.1s**。
  ROADMAP 對 2.059 直接捨入到 2.1s,稍為激進,三份文檔結論不統一。
  建議:ROADMAP R3 改寫 `2.06s` 或標注 `(~2.1s)`。

- **[MINOR]** `docs/research/WS_R2_experiments.md` §2 結論文字:
  "stale-LDLT 預條件把迭代數壓到 Jacobi 的 ~1/100"。
  原始輸出 `incremental_xxl.txt` 四個批量案例的實際迭代比值為:
  nb=16 → 1218/12 = **101×**、nb=40 → 1221/14 = **87×**、nb=80 → 1246/16 = **78×**、nb=160 → 1273/18 = **71×**。
  只有 nb=16 嚴格符合 ~100×;nb=160 僅 71×,與 "~1/100" 有 30% 偏差。
  在「最典型的互動移除量級(16–40 桿)」下成立,但作為廣泛宣稱有過度簡化之嫌。
  建議:改為 "迭代數約為 Jacobi 的 1/70–1/100(移除桿數越多比值越小,16 桿時 ~1/100、160 桿時 ~1/71)"。

## 沒問題的部分

所有其他核心數字(31×/7.7e-14/24.6×/522s/1.97e-12/0差動量/6.7ms/1593.1632/12–18 迭代/P-Delta 各迭代數/能量漂移 6.6e-11/截斷誤差 39%/7.3%/稀疏屈曲記憶體與殘差/Tier-2 四行表格所有欄位/規模階梯 nnz 與 factor 時間)在各文檔間全部一致,且與原始輸出逐項吻合(精度均屬合理捨入);表格無錯行、單位無誤、6 個抽驗倍率計算均正確。