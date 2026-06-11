# S8 進度 — 殼路線(8a QM6 opt-in 不協調膜 + 8b DKQ 薄板快路)

> 接續 `PROGRESS_S7.md`(S7 BESO `ec157a0`)。規格 `docs/specs/S8_shell.md`,研究依據
> `docs/research/WS_E_shell.md` + 膜鎖定 P5 探索(plan `moonlit-wishing-charm.md`)。

## 摘要

S8 在既有 MITC4 flat-shell facet 上加**兩條 `SolveOptions` opt-in 分支**,**共用** 24-DOF facet
框架 / 組裝 / 3D 旋轉 / 質量 / recover 骨幹,**預設兩旗標 false = 原始 MITC4 逐位元不變** →
OpenSees `ShellMITC4` 平板閘門(~1e-10)不破。**不新增 element 類別 / 不改 dispatch / SolveResult /
modelFingerprint / 不新增 `.cpp`**(全寫進 `MITC4ShellElement.cpp`,四 build 腳本源檔清單免動)。

- **8a `useIncompatibleMembrane`** = QM6 不協調膜(Wilson 1973 + Taylor 1976):2 bubble 模式
  φ₁=1−ξ², φ₂=1−η²(4 內部 DOF)+ 中心 Jacobian + static condensation,根治 in-plane 膜鎖定。
- **8b `useDKQPlate`** = DKQ 離散 Kirchhoff 薄板(Batoz & Tahar 1982):8 節點 serendipity 旋轉場 +
  沿邊 Kirchhoff 約束消邊中 DOF → 12-DOF 純彎曲薄板,**無橫向剪**(Qx=Qy=0),薄板專用(t/L<~1/20)。

**五腿全綠**:standalone **F1-F49**(F48 QM6 / F49 DKQ)/ UE **46**(+2:`Shell.IncompatibleMembrane`
/ `Shell.DKQPlate`)/ OpenSees PASS(ShellMITC4 預設 ~1e-10 + **DKQ vs ShellDKGQ 1.69e-12**)/ audit
**90 checks**(+8)/ CLI round-trip(frame_cli `OPT` 擴充 2 token)。

## 8a QM6 不協調膜

**數學**:位移增廣 u=ΣNᵢuᵢ+P₁a₁+P₂a₂, v=ΣNᵢvᵢ+P₁a₃+P₂a₄。`B_inc`(3×4)用**元素中心常數
Jacobian J₀**(Taylor 1976 修正)→ ∫B_inc=0 → 常應力狀態 bubble 零能 → **畸變網格 patch test 通過**。
增廣勁度 static condensation `K* = Kc − Kca Kaa⁻¹ Kcaᵀ`(Kc 重用 `membraneK` 避免公式漂移)。

**⚠️ 關鍵發現(drilling 耦合)**:初版只把 bubble 耦合到**膜應變**,QM6 細長懸臂只解鎖到 −38%
(P5 純膜 −1.1%)。根因:**Hughes-Brezzi drilling 用協調 u,v(不含 bubble)算旋轉 ω**;純彎曲需
`v~x²`(bubble 提供)才能讓 drilling 應變 `thz+½u,y−½v,x` 歸零,但 drilling 看不到 bubble →
殘留偽 drilling 能量**把膜 bubble 釋放的鎖又鎖回去**。修正:加 `BdrillIncQM6`,讓 bubble **同時**
耦合膜 + drilling(`Kca += γ·Bdᵀ·Bdi`)→ QM6 解鎖到 **−0.9%**(吻合 P5)。Q4 對 P5 純膜逐項吻合
(−75.5% vs −75.6%)佐證膜公式正確。

**oracle F48**:(a) 畸變膜 patch 機器精度(eNxx=0);(b) 細長 4×1 in-plane 懸臂 Q4 −75.5%(鎖死)
→ QM6 −0.9% vs Euler-Bernoulli;(c) Cook's skew membrane 自洽收斂(QM6 N=4 −4.2% 收斂快於 Q4
−27.5%,+ Q4 單調趨近 QM6 收斂值的相容性論證,refC≈25.0 vs 文獻 23.96 差 ~4% 為 drilling penalty +
角點取樣慣例,非缺陷)。

**recover**:中心膜 recover **無需改**——`B_inc(0,0)=0`(∂P₁/∂ξ=−2ξ|₀=0、∂P₂/∂η=−2η|₀=0),
且 `ShellElementForces.N` 是中心值,故膜應變 = `Bm(0,0)·dm` 自動受益於 QM6 改善的位移解,bubble
中心無貢獻(per-corner 膜輸出才需凝聚 bubble,未來加再處理)。

## 8b DKQ 離散 Kirchhoff 薄板

**數學**:`serendipity8`(8 節點值+導數)+ `dkqEdgeCoeffs`(Batoz a,b,c,d,e)+ `dkqH`(Hx,Hy 12-vec,
Batoz (w,θx,θy) DOF)+ `BdkqMine`(曲率 3×12,Jacobian 轉 x,y + DOF remap)+ `plateK_DKQ`(2×2
Gauss,無 shear)。

**⚠️ 關鍵發現(DOF 符號 remap)**:Batoz native (w,θx,θy) → 我的 (w,bx,by) 的對映:
`bx=θy, by=−θx`。**初版符號寫反**:**撓度收斂完美(0.05% to Kirchhoff)但 patch moment 大錯
(11.5 vs 2.747)**。原因:`K = Bᵀ·Db·B` 對 remap 整體符號**不敏感**(平方消去),故撓度對 ≠ remap
完全對;`recover` 的 `M = Db·B·dp` 才暴露符號。constant-curvature patch test 正確仲裁出
`B.col(bx)=+Bb.col(θy), B.col(by)=−Bb.col(θx)`,修正後 patch 機器精度。

**oracle F49**:(a) 畸變常曲率 patch 機器精度(regular + parallelogram,eMxx=0);(b) SS 薄板
(t/a=1/200)DKQ 收斂 Kirchhoff 0.00406qa⁴/D(N16 0.05%,無 Mindlin overshoot);(c) clamped 薄板
0.00126qa⁴/D(N16 1.52%)。**OpenSees ShellDKGQ 對照 = 1.69e-12**(flat 薄板純彎曲,膜不參與,只比
共用 Batoz DKQ 板;機器精度級 = 第三方獨立元素鐵證)。

**recover**:DKQ 分支板矩用 `BdkqMine`,**Qx=Qy=0**(Kirchhoff 無橫向剪,刻意非遺漏);per-corner
同 MITC4 機制。

## 檔案改動

- `Public/FrameCore/SolveOptions.h`:+`useIncompatibleMembrane` / `useDKQPlate`(2 bool,預設 false)。
- `Private/MITC4ShellElement.h`:+`useQM6_` / `useDKQ_`(prepare 設定,recover 分支用)。
- `Private/MITC4ShellElement.cpp`:+型別別名 Mat3x4/Mat12x4/Mat4/Row4;+`BincQM6` / `BdrillIncQM6` /
  `membraneK_QM6`(8a);+`serendipity8` / `dkqEdgeCoeffs` / `dkqH` / `BdkqMine` / `plateK_DKQ`(8b);
  prepare 讀 opts 選分支;recover DKQ 分支(板矩 + Qx=Qy=0);+include `SolveOptions.h`。
- standalone `main.cpp`:F48 / F49;`FrameTestFixtures.h`:+`slenderMembraneCantilever` / `cooksMembrane`。
- UE `Private/Tests/ShellS8Test.cpp`(新,+2 測試);`run_gate.ps1` `$ExpectedUeTests` 44→46。
- `linear_deep_audit.cpp`:+`testShellS8`(8 checks,82→90)。
- `frame_cli_core.cpp`:`OPT` 指令 +2 向後相容 token(im/dk);`Tools/opensees_compare.py`:
  `run_frame_cli_shell` / `run_opensees_shell` +`dkq` 參數 + DKQ vs ShellDKGQ leg。
- **不動**:dispatch(`FrameSolver.cpp`)/ `Reanalysis.cpp` / `Shell.h` / `SolveResult.h` /
  `FrameModel.cpp`(fingerprint/validate)/ 四 build .bat 源檔清單(無新 .cpp)。

## 誠實邊界 / novelty

- **QM6 = Wilson 1973 + Taylor 1976,DKQ = Batoz-Tahar 1982**:皆成熟文獻元素,FrameCore 為其**實作**,
  **非新方法**。`[VERIFIED]` = 文獻方法 oracle 化;`[NEW CODE]` 僅限「opt-in 旗標整合 + bubble 同時
  耦合膜/drilling 的工程整合」+「DKQ 接 24-DOF facet 的 DOF remap」。
- **QM6 只幫膜主導 / in-plane bending**;純板彎曲(MITC4 assumed-shear 已良好)+ flat-facet O(1/N²)
  幾何 floor 不受益。QM6 過**弱** patch test(modified integration),高度扭曲殘差仍存。
- **DKQ 完全無厚板能力**(無橫向剪;t/L>~1/20 必回 MITC4);`Qx=Qy=0` 刻意;薄板快路 / Kirchhoff
  oracle,非主殼元素(WS_E「避什麼」)。
- **預設安全**:兩旗標 false 保 OpenSees ShellMITC4 ~1e-10;改進膜/板第三方對標各走 ShellDKGQ(DKQ)/
  文獻收斂(QM6 無直接 OpenSees 膜對應)。對標 Karamba3D TRIC 須附先行技術定位(WS_E F-E4)。

## ⚠️ 踩雷(durable)

1. **QM6 必須讓 bubble 同時耦合膜 + drilling**,否則細長彎曲只解鎖一半(drilling 用協調 ω 鎖回)。
2. **DKQ DOF 符號 remap 靠 patch test 仲裁**:K 對 remap 整體符號不敏感(平方),撓度對 ≠ remap 對;
   constant-curvature patch 的 recover moment 才暴露符號。`bx=+θy, by=−θx`。
3. **Bash→powershell 路徑**:反斜線被 bash 跳脫(`E:\...` → `E:...`),跑 .ps1 gate 用 PowerShell 工具
   或正斜線。
4. **UE build 用 PowerShell 前景 Build.bat**(背景增量漏編新測試致假綠);run_gate 前先 Build.bat。
5. **UE 測試常數勿用 IN/OUT**(SAL 巨集,經 CoreMinimal.h)。
6. **frame_cli `OPT` 擴充用 `if (ss>>tok)`** 保向後相容(舊 `OPT er ut pt` 不破)。
7. **Cook's membrane 絕對值依賴設定**(載重分配/量測點/drilling),用自洽收斂(細網格自參考)+ 相容性
   論證,勿硬綁文獻 23.96。

## 下一步

S9 Co-rotational(WS_F,Battini 2002 / NR + load stepping;elastica shooting 表已備)→ S10 N-M 互動
塑鉸(必在 S9 後,R4 方向耦合)→ S11 MITC9i(殿後)。提示詞 `docs/AGENT_PROMPT_S5_S11.md`。
DKQ 後續:per-corner 膜 QM6 輸出、變厚度、DKQ 殼幾何勁度(薄板挫屈)。
