# S1 進度日誌(無人監督夜間執行 — 2026-06-10 23:35 起)

> 睡前授權:照 `docs/IMPLEMENTATION_PLAN.md` 從 S1 起盡量往後推、不半成品、務實分層 gate、持續 commit/push。
> 本檔記錄每個交付點的決策與假設,供使用者醒來檢視。gate 政策:**standalone(build.bat + build_linear_audit.bat)綠就 commit**;UE automation + OpenSees 重 gate 留到 **S1 整體完成的里程碑** 跑一次。

## 基準
- 起點 commit `81639c1`(研究輪)。working tree 既有雜項(`.gitignore`/`ArchSim.uproject` 改動、`Plugins/LevelSim/`)**非本人改動,全程不碰**。
- baseline standalone gate 確認綠:`ALL PASS (failures=0)`,F1..F34(原 F1..F33;F33=event-to-event 塑鉸驅動,與 F32 共用 `[F32]` printf 表頭)。

## 已完成交付點

### 1. `docs/IMPLEMENTATION_PLAN.md` 入庫 — commit `63b80bc` ✅
單獨 commit,嚴守 commit 衛生(未掃入雜項)。

### 2. R8:修 `build_perf.bat` 連結失敗 — commit `0e2e500` ✅
- 源檔清單補 `MITC4ShellElement.cpp`(FrameSolver 的 assembleAndFactor 即使純梁柱模型也引用殼 dispatch 符號 → LNK2001)。
- 順手補齊 `vswhere` 目錄上 PATH 的一致性修復(build.bat/build_cli/build_linear_audit 都有,build_perf 漏了 → vcvars 裸 vswhere 警告)。
- 驗證:`frame_perf.exe` 連結乾淨無警告。獨立 bugfix,與 S1 新增無耦合。

### 3. 稀疏屈曲 overload + F34 oracle — commit `a91b171` ✅
- **`BucklingOptions`**(POD:`denseThreshold=500`/`nev=1`/`maxIter=300`/`tol=1e-11`)+ `solveBuckling` 三參 overload;舊兩參版委派 `BucklingOptions{}`。
- `BucklingAnalysis.cpp` 重構:共用步驟(參考解→軸力→Kg)→ 依 `nf` vs `denseThreshold` 分派 **稀疏**(`subspaceSmallest(Kff, S.ldlt, negKgff, nev)`,復用既有 LDLT,`lambda(0)`=criticalFactor)或 **稠密**(GES,**bit-identical** 保 F23 不變);稀疏失敗→稠密 fallback(永遠正確)。
- **決策/假設(誠實標)**:
  - **F 編號**:稀疏屈曲=**F34**(實作順序;spec ⑥ 暫定 F36,但 baseline 最高 logical fixture=F33,故依 PLAN §7「以實際新增為準」順序編號;ReSolve tier1/tier2 接 F35/F36)。
  - **殘差閘**:在 spec ④ 之外加一道寬鬆 pencil 殘差閘(`‖(−Kg)φ−γKφ‖/(γ‖Kφ‖) < 1e-3` 才採信;研究輪實測 1e-7~1e-10,故對良態模型永不誤觸,只擋「收斂但不準」→ 退稠密)。屬保守安全網,文檔已標。
  - **退化護衛**:全拉(`gtrips.empty()`)→ 稠密路徑,沿用既有「no compression」診斷;強制稀疏在全拉時亦退稠密(`lambda≤0` 或 Kg 空),不偽造正屈曲因子。
  - **F34 用矩形斷面**(`Rectangular(60,100)`):避免方斷面 y/z 簡併模態,讓最低屈曲值單一、稀疏迭代乾淨收斂。
- **驗證**:
  - standalone F34:pinned-pinned n=10/fixed-free n=10/pinned-pinned n=24,**sparse==dense rel = 2e-14 / 2e-13 / 1e-12**(非零→確認走稀疏非 fallback)、sparse==Euler rel = 1.35e-5 / 8.4e-7 / 4.1e-7、全拉護衛三檢全過。`ALL PASS (failures=0)`。
  - linear_deep_audit +1 check「sparse path agrees with dense」rel=5.48e-13 → **checks 62→63 PASS**。
  - 原型參照:`Research/WS_B_solver/exp_sparse_buckling.cpp`(同款 reduceFF + subspaceSmallest 映射)。
- **無新 `.cpp`** → build 腳本源檔清單不需改。

### 4. ReSolveSession Tier-1 + Tier-3(commit A,本次)✅
- 新 `Public/FrameCore/Reanalysis.h`(POD `ReanalysisOptions`/`ReanalysisStats` + PIMPL `ReSolveSession`,零 Eigen 洩漏)+ `Private/Reanalysis.cpp`。
- **Tier-1 Woodbury**(移植原型 `exp_incremental_refactor.cpp`):暫態 `BeamColumnElement`/`MITC4ShellElement` prepare+assemble 抽 bit-consistent 元素 K → 降 free + `SelfAdjointEigenSolver` 取正模態低秩 → `K_cur_ff = K0_ff + W diag(s) Wᵀ`;solve = `u0ff − Z·C⁻¹·(Wᵀ u0ff)`,**C 奇異(pivRatio < mechPivotTol)= 機構**(從因子判定,不靠連通)。
- **Tier-3 rebaseline**:rank > maxRank → fresh `assembleAndFactor`(永遠正確);`fmap` 由支承決定、對元素啟用不變,故重建後 reduced RHS 仍有效。
- **F 增量 [NEW CODE]**:每次 solve 從「當前 active 集」重建 `F_cur`(nodal + active 元素 `addEquivalentNodalLoads`)→ 停用桿/殼的 UDL/壓力等效載自動離開、恢復自動回來;另組 `K_cur`(不分解)供 reactions + prescribed term;recover 用當前 active 元素。
- **驗證**:standalone **F35**:portal+UDL 桿移除 reRel **2.43e-12**、restore drift 1.3e-16、column 移除 5e-12、2 桿鏈移底→**機構偵測**、clamped plate 殼 facet 移除 reRel **1.70e-12**(殼 F 增量也對;殼 rank=6 因 clamped plate 只中心節點自由)。audit **+3**(Tier-1==fresh 4e-14 / restore 2.6e-16 / 機構)→ **63→66**。`build.bat` + `build_linear_audit.bat` 補 `Reanalysis.cpp` 源檔。
- **誠實標**:Tier-1 精確但 vs fresh ~1e-13(浮點路徑差,非逐位元);每次 solve 重建元素 + 組 K_cur(省 factorization、不省組裝,speedup 比 u-only 原型 31× 略低、同量級,精確數待 commit B/PERFORMANCE 量測);Tier-2 PCG 未做(rank 超 maxRank 走 Tier-3,正確但較頻繁重分解)。

## 待補(S1 里程碑重 gate 一次處理)
- **UE automation**:加 `FrameCore.Buckling.SparseAgreesDense`(+ 之後 ReSolve 的 `FrameCore.Reanalysis.LadderAgreesFresh`/`MechanismDetection`)→ bump `run_gate.ps1` `$ExpectedUeTests`(34→實際數)→ 跑 headless UE 測試。
- **OpenSees strict**:既有「移除態」逐位移場景改走 ReSolveSession 重跑一次(同容差)。
- 理由:務實分層 gate 政策(額度受限),重 gate 集中里程碑跑。

## 下一個交付點:ReSolveSession **Tier-2 PCG**(commit B;commit A=Tier-1+Tier-3+F 增量已完成,見「已完成 4」)

> ✅ **commit A 已完成**(Tier-1 Woodbury + Tier-3 rebaseline + F 增量 + recover + 機構 + F35;standalone 全綠、audit 66)。
> **剩 commit B = Tier-2 stale-LDLT PCG**(純效能中間層,可選):rank 落在 (maxRank, 2·maxRank] 時改用 `Eigen::ConjugateGradient<SpMat, Lower|Upper, StalePrecond>`(precond = baseline ldlt,guess = 上次解),取代直接跳 Tier-3。原型 `Research/WS_N_incremental/exp_incremental_refactor.cpp` 的 `StalePrecond`(~15 行)+ `batchUpdateSparse` 可直接移植;**容差互驗(非逐位元)**,F 編號 F36、audit +1。下面的引擎接點筆記對 Tier-2 仍適用。

**已確認的引擎接點(寫 Reanalysis.cpp 直接用)**:
- 型別:`MemberId=int`、`gdof(nodeIdx,d)=6*nodeIdx+d`(`FrameTypes.h`,public,`DOF_PER_NODE=6`)。
- baseline:`assembleAndFactor(base, opts)` → `const auto& S = *ps.impl`;拿 `S.K`(全域 N×N 稀疏)、`S.fmap`(free map,-1=約束)、`S.nf`、`S.ldlt`(**K_ff 的 LDLT**,可直接當 Woodbury 的 K0^-1)、`S.elems`(**只含 active** 元素)。
- **元素 K 抽取**(每根停用/恢復的 member/shell):建暫態 `BeamColumnElement(memberIdx)`(12 DOF)或 `MITC4ShellElement(shellIdx)`(24 DOF)→ `prepare(model,opts,why)` → `assemble(trips)` 拿全域三元組(與 baseline **逐位元同源**,因同一份元素碼)。從 trips 收集相異全域 DOF → 組 ne×ne 密集 Ke。
- **降到 free + 低秩**:keep=fmap≥0 的 local idx,rid=其 fmap 值;`SelfAdjointEigenSolver(Ks)` 取 λ>1e-9·λmax 的正模態 → W 欄(nf 空間),s=sign·λ(−1 移除/+1 恢復)。梁 rank≤6、殼≤18。
- **Woodbury**(直接移植原型 `Ladder`,`exp_incremental_refactor.cpp` L54-140):`Z=ldlt.solve(W)`;`Mc=Wᵀ Z` 增量長(新塊 `B1=Wnᵀ Z`,對稱免重算);solve:`C=Mc+diag(1/s)`→`FullPivLU`→`pivRatio=min/max|diag U|`,`<mechPivotTol(1e-10)`=機構;`u_ff=u0ff − Z·(C⁻¹·(Wᵀ u0ff))`。
- **F 增量 [NEW CODE](原型未驗,F35 須把關)**:baseline `Fff = reduceVec(nodal + Σ active elems addEquivalentNodalLoads(F))`。停用 member/shell 時:暫態元素 `addEquivalentNodalLoads(Fe=Zero(N))` → reduceVec → **從 Fff 減掉**(該桿 UDL/該殼壓力的等效載隨元素離開);恢復則加回。`u0ff = ldlt.solve(Fff_current)` **每次 solve 重算**(F 變了)。⚠️ 純節點載時 ΔF=0(原型情境);UDL/殼壓力才有 ΔF——這正是風險點。
- **recover**:`u_ff` scatter 回全域 `u(N)`(free 填 uf、約束填 prescribed)→ 對**目前 active 集**建暫態元素 prepare → `recover(u,R)` 出桿端力/殼內力;先 stamp 所有 id(鏡像 `FrameSolver.cpp` L266-273,inactive 列留 id+0 力)。`R.pivotMargin` 填 baseline 構型值(語意=基準健康度,文檔標)。
- **Tier-3 rebaseline**:`assembleAndFactor(work)`(work=base 拷貝+當前 active flags)→ 換新 baseline,清空 W/Z/Mc/deltaSet/Fff 重建。永遠正確(=fresh 路徑)。
- **dispatch**:R=Σrank;R==0→tier0(`ldlt.solve(Fff)`);R≤maxRank→Tier-1;否則 (commit A) →Tier-3 rebaseline / (commit B) →Tier-2 PCG 否則 Tier-3。
- 需自寫小工具(research_common 不在引擎):`reduceVec(Fglobal,fmap,nf)`、`scatterVec(uff,fmap,N)`(`reduceFF` 已在 BucklingAnalysis.cpp,可複用或內聯)。
- **build 同步**:新 `Private/Reanalysis.cpp` 要補進 `build.bat` + `build_linear_audit.bat` 源檔清單(顯式,顯式豁免 build_cli/build_perf 並註明)。
- **API**:`Public/FrameCore/Reanalysis.h`(POD `ReanalysisOptions{maxRank=96,pcgTol=1e-10,pcgMaxIter=500,allowTier2,mechPivotTol=1e-10,SolveOptions solve}`、`ReanalysisStats{tier,rank,pcgIters,relResidual,refactored,mechanism}`、PIMPL `ReSolveSession`,簽名見 `docs/specs/S1_resolve_ladder.md` ②)。零 Eigen 洩漏。
- **F35 oracle(spec ⑥;新編號 F35,非 spec 暫定)**:塔(3,2,4)序列移除 8 桿(**含一根帶 UDL 的桿**)每步 vs fresh `assembleAndFactor+solveLoad` relMax≤1e-10;**+ 殼壓力子案例**(F 增量殼版,殼 rank≤18 未在原型驗過→不過則殼走 Tier-3);全恢復 vs baseline≤1e-12;機構:2 元素鏈移底→mech=true 且 fresh singular;portal 移梁→兩者皆穩。audit +(tier1 一致/恢復漂移/機構/[Tier-2 容差待 commit B])。
- ⚠️ 殼 rank-18 風險:先加殼版 F35 子案例;不過 → 殼一律走 Tier-2/3(階梯仍成立),文檔誠實標。
