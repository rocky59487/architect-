現在我有足夠的資訊進行完整審查。以下是審查報告:

---

# 審查報告:S3/S4 spec + S5-S11 骨架可實作性

## 裁決

**PASS(有條件)** — 未發現必須修才能定稿的 CRITICAL 問題;但存在 2 個 MAJOR 與若干 MINOR 問題需在升格為交接級 spec 前修訂。

---

## 發現(按嚴重度排序)

### [MAJOR] S3_pdelta.md §④:Aitken 加速器的「保護式」描述與原型實驗實際實作不一致

**位置**:S3_pdelta.md §④「保護式外推」段落 + `exp_pdelta_convergence.cpp` 第 38-47 行

**問題**:spec 描述的保護式觸發條件為「僅當 `0<ρ<0.95` 且 `|ρ_m−ρ_{m-1}|<0.2` 才外推,外推後下一步 `||Δ||` 未降 → 撤銷並停用」。但原型程式碼(`frozenIterate` with `aitken=true`)的實際邏輯是:當 `rho > 0 && rho < 0.999` 就執行外推,**沒有「穩定性條件 `|ρ_m−ρ_{m-1}|<0.2`」的檢查,也沒有「外推後未降則撤銷並停用」的回退機制**。spec 描述的是更謹慎的「保護式 Aitken」,但原型實驗實作的是更接近裸 Aitken(只多一條 `rho < 0.999` 上限)。

實驗數據亦佐證:f=0.9 時裸平均版 285 步 → Aitken 版竟爆炸到 4742 步(spec 第 51 行說「f=0.9 實測裸版劣化 285→4742」),但 spec 的說法是「裸 Aitken」才會這樣劣化,而原型 Aitken 版正是 4742 步,即原型中的「Aitken」版就是裸版,與 spec 聲稱的「保護式」不同。

**影響**:實作 PDeltaAnalysis.cpp 時若照 spec 的保護式邏輯寫,將與研究輪的數值基礎(凍結 vs 參考的收斂對比、itersAitken 數值)不符,且 f=0.9 的 4742 步這個「壞例子」在實際保護式版本中可能不再重現,使 F40 oracle 的教訓數字失去支撐。

**建議**:兩種修法擇一:(a) 修正原型程式碼使其真正實作保護式邏輯並重跑實驗更新數據,再讓 spec 描述與新數據一致;(b) 修正 spec 描述與實際原型邏輯一致(裸 Aitken,條件 `0<ρ<0.999`),並重新措辭「f=0.9 Aitken 版 4742 步」為「裸 Aitken 在此點劣化」而非聲稱已有保護。選(a)較佳,因為保護式邏輯確實更好。

---

### [MAJOR] S3_pdelta.md §④:發散偵測描述與原型實驗邏輯不符

**位置**:S3_pdelta.md §④「發散偵測」段落 vs `exp_pdelta_convergence.cpp` 第 35-37 行,以及 `pdelta.txt` f=1.05 那行

**問題**:spec 描述的發散偵測為「20 步窗 `||Δ||` 趨勢上升且超過初值 → diverged」。但原型程式碼的實際邏輯是單步檢測:`dn > 4.0 * normPrev && dn / xn > 1.0` — 只看**相鄰兩步**增長超過 4 倍,而非 20 步滑動窗口趨勢。

f=1.05 的實驗結果(`itersPlain=5000 conv=0 div=0`)顯示凍結路徑跑到 maxIter=5000 後回報 `converged=false`(非 `diverged=true`)。也就是說:原型平坦版本的發散偵測在 f=1.05 這個案例完全沒觸發(div=0),而是靠 maxIter 兜底。**spec 聲稱 f>1 → `diverged=true`(兩路徑)的 F40 oracle 承諾,與現有原型實驗數據不符 — 平坦版 div=0。**

**影響**:F40 oracle 的 `f=1.05 → diverged=true(兩路徑)` 承諾若按原型邏輯無法保證,需要確認實作版 PDeltaAnalysis.cpp 是否會做真正的發散偵測使此 oracle 成立,還是需修正 spec 承諾。

**建議**:明確標示「原型實驗凍結路徑在 f=1.05 以 maxIter 兜底而非 div=true」,F40 oracle 把 f=1.05 凍結路徑的期望由 `diverged=true` 修正為 `converged=false`(或在實作中補實際發散偵測邏輯並說明與原型差異)。

---

### [MINOR] S4_tension_only.md §⑤:fingerprint 擴充說明正確但缺現狀比對說明

**位置**:S4_tension_only.md §⑤「改 `FrameSolver.cpp` `modelFingerprint`(+tensionOnly hash)」

**問題**:現行 `FrameSolver.cpp` 的 `modelFingerprint`(第 48-57 行)對每個 member 已 hash `id/i/j/matIdx/secIdx/refVec/release[12]/active`,**但尚未包含 `tensionOnly` 欄位**,因為 `Member.h` 現在也沒有此欄位。spec 的修改計劃本身是正確的(加 `tensionOnly` 入 hash),但 spec 未明說「當前 fingerprint 沒有 tensionOnly,實作 S4 時必須同步加」,若實作者沒注意到兩個檔案都要動,可能只加了 `Member.h` 欄位而漏改 fingerprint。

**建議**:在 §⑤ 明確標注「`FrameSolver.cpp` modelFingerprint 的 member 迴圈內需在 `active` hash 後加一行 `h = fpMix(h, mem.tensionOnly ? 1ull : 0ull);`」,避免漏改。

---

### [MINOR] S3_pdelta.md §①:精確解容差數字「2.6e-7」與實驗數據一致但 f=0.95 數字未引用

**位置**:S3_pdelta.md §⑥ oracle F40 描述 vs pdelta.txt

**問題**:spec 的 F40 oracle 提到「f=0.3 精確 2.6e-7、f=0.95 精確 3.7e-5」。對照 pdelta.txt:f=0.30 `relRefVsExact=2.639e-07`(一致),f=0.95 `relRefVsExact=3.716e-05`(一致,spec 寫 3.7e-5 略微四捨五入,可接受)。數字本身正確,但 **spec 沒提 f=0.5 的 1.026e-06 與 f=0.7 的 3.353e-06 和 f=0.9 的 1.666e-05**。oracle 只標兩個點不算錯,但 [VERIFIED] 的嚴格性下,建議 F40 引用的是 oracle 承諾用的 f 值而非全部掃描點;現在已足夠。此項純粹紀錄。

---

### [MINOR] S4_tension_only.md §⑥:F42 「收斂 ≤3 迭代」與實驗數據輕微出入

**位置**:S4_tension_only.md §⑥ F42 oracle vs tension_only.txt

**問題**:spec 寫「F42:X 斜撐 portal 收斂 ≤3 迭代」,實驗數據 `[caseA] iters=2`。數據支撐「≤3」的承諾(2 < 3),字面上沒錯。但 spec 同時承諾「收斂解 vs 省略壓桿之模型逐位元相等(≤1e-15;實測 0.0)」— 實驗 `relErrVsOmitted=0.000e+00` 完全支撐。此點標為 PASS,附記供參考。

---

### [MINOR] S5_S11_skeletons.md S5:「FSD 靜不定非最優保證」標記正確,但 10-bar 收斂值可驗

**位置**:S5_S11_skeletons.md §S5

**問題**:spec 標「10-bar vs 文獻 1593.16 lb(WS_R2 §7)」且 `size_opt.txt` 應有對應數據。size_opt.txt 存在於 `Research/out/`,但骨架未引用實際數字,僅說「WS_H 交叉引用」。骨架等級文件這樣寫可接受,但若升格交接 spec 時應補。此 MINOR 是骨架完成度的正常狀態,非缺失。

---

### [MINOR] S5_S11_skeletons.md S9:「elastica shooting 表,α=1→δv/L=0.3017207738 等,雙容差十位一致」

**位置**:S5_S11_skeletons.md §S9 oracle 數據描述

**問題**:這組數字宣稱「雙容差十位一致」— 十位數字的精度宣稱。此項需有獨立計算或文獻可查可重現,目前 `elastica.txt` 有對應實驗結果。若 `elastica.txt` 已有數據支撐則標 [VERIFIED] 是可行的,但骨架未明示來源是實驗結果還是文獻引用。升格 spec 前需補清來源。

---

### 不需修改、確認正確的部分

- **S4 fingerprint 設計方向正確**:S4 要求 `tensionOnly` 入 fingerprint,現有 `modelFingerprint` 的 member 迴圈結構(第 48-57 行)已有 `active` hash 的位置,要插入 `tensionOnly` 只需一行,無架構衝突。
- **S3 雙路徑架構**:「凍結路徑重用 LDLT + 參考路徑重分解」在 `PreparedSystem::Impl` 結構(`S.ldlt`/`S.K`/`S.elems` 均公開在 impl 內)下完全可讀取,零架構阻力。
- **S4 循環守門描述**:哈希守門 + 單調 fallback 設計,與原型程式碼邏輯(`seen.insert(stateBits*4+nextBits)`)一致,僅注意原型的 hash key 是「(前態, 後態)轉移對」而非只有「後態」,spec 描述「state-hash」未細說但不矛盾。
- **S5-S11 骨架整體評估**:均正確定位為「研究結論定案、oracle 構想成形、動工前補全十項 spec」級別,沒有超出已驗數據的過度宣稱。S7 N2 的「類 fail-safe TO 但離散框架+LSP 評估器,勿用 fail-safe 一詞」誠實邊界措辭到位。S9 彈性線與 OpenSees corotCrdTransf 雙對照的思路正確。S10 N-M 互動「仍 sequential linear;無卸載」誠實標清楚。
- **S3/S4 的 `[PENDING]` 標記使用規範**:S4 F43 明確標 `[PENDING:掃描未尋獲,尋獲後升級為真實案例]`。這是正確的誠實標記用法。
- **性能驗收承諾(S3 §⑧)**:「凍結路徑 f=0.5 迭代 ≤60(實測 46)」— pdelta.txt 顯示 f=0.5 `itersPlain=46`,完全吻合 [VERIFIED]。

---

## 沒問題的部分(一句話)

S3/S4 的架構設計(雙路徑互驗、凍結路徑重用 LDLT、fingerprint 擴充計劃、循環守門+單調 fallback)均與既有引擎及實驗數據相符,S5-S11 骨架誠實標記到位、無過度宣稱,但 S3 的「保護式 Aitken」描述與原型實作邏輯以及 f=1.05 diverged 承諾與實驗數據有兩處需修正。