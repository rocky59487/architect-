# S10 進度 — N-M 互動塑鉸(軸力-彎矩互動降伏)

> 接續 S9 家族收尾(`67db412`)。規格 `docs/specs/S10_nm_interaction.md`。**S10 = 材料非線性(塑鉸)線收尾**:
> 把 stage-4b event-to-event 塑鉸的固定 `Mp` 換成軸力降伏的 `Mp_eff(N)`。

## 摘要

S10 在既有崩塌驅動器 `runProgressiveCollapse` 加 opt-in 軸力-彎矩互動(`CollapseOptions::nmInteraction`):
hinge-capable 桿的塑性彎矩在軸力下降為 `Mp_eff(N) = Mp·max(0,1−(N/Ny)²)`(`Mp=fy·Z`、`Ny=fy·A`),
觸發判據與 residual 從 `Mp` 改為 `Mp_eff(N)`,端對齊軸力。**預設 off → 餵 N=0 → 逐位元同 stage-4b**。

- **新增**:`Public/FrameCore/NMInteraction.h`(inline `reducedPlasticMoment`,純 POD/零 Eigen)、
  `CollapseOptions::nmInteraction`(bool, 預設 false)。
- **改動**:`Collapse.cpp` 候選迴圈(`#include NMInteraction.h` + Mp→Mp_eff(N) 替換)。
- **不新增 .cpp**(header-only)→ 四 build 腳本免動;**不入 fingerprint**(純後處理旗標,不改 K);
  不新增 element/dispatch/SolveResult 欄位。

**五腿全綠**:standalone **F1-F54** / UE **50**(+1 `Collapse.NMInteraction`)/ OpenSees PASS / audit **104**(+3)/
CLI round-trip ALL PASS。`$ExpectedUeTests` 49→50;audit 101→104。

## 力學 / 演算法

- **降伏式** `Mp_eff(N)=Mp·(1−(N/Ny)²)`(負截 0):
  - **矩形實心 = 一次原理精確**:中性軸偏移,中央帶 `2c=N/(b·fy)` 承軸力,殘餘 `M_N=Mp−N²/(4·b·fy)`
    ≡ `Mp(1−(N/Ny)²)`(EC3 §6.2.9 / 經典塑性);
  - **圓形 = 保守**(真包絡更飽滿,鉸不晚於實際);
  - **單軸 |N|**:平方項使拉/壓對稱、符號無關;**無 My-Mz 雙軸耦合**。
- **觸發**(`Collapse.cpp`):`Mpk=reducedPlasticMoment(Mp0[k],N_end,Ny)`,dof 4/5 用 `endI.N`、10/11 用 `endJ.N`;
  `ratio=|M|/Mpk`、`ratio≥1` 成鉸、residual `±Mpk`。`Mpk≤0`(軸力 squash)→ 跳鉸候選,脆性軸力篩主導。
- **凍結**:residual 存 `PlasticHinge.Mp`,成鉸端 recover `M=0`(F32 契約)→ 不重估、不重複觸發 = 凍結自然成立。
- **OFF 逐位元**:N=0 → `reducedPlasticMoment(Mp,0,Ny)=Mp·(1−0)=Mp`(FP 恆等);`if(!(Mpk>0))continue` 在 OFF
  永不觸發(`Mp0[k]=fy·Z>0`,hingeCapable 已 gate)→ 行為與 stage-4b 完全一致。

## Oracle(實測)

### F54(standalone,12 checks 全綠)
- **公式**:`Mp_eff(0)=Mp`/`(Ny/2)=0.75Mp`/`(Ny)=0`/`|N|` 對稱/`2Ny` 截 0 — **rel 0**(tol 1e-15)。
- **矩形精確**:vs `Mp−N²/(4bfy)` 中性軸偏移 — **rel 0**(tol 1e-12)。
- **驅動器括弧**(靜定 X 懸臂,基端 `wL²/2`、軸力 P 皆 stiffness-independent;`P=Ny/2→Mp_eff=0.75Mp`):
  `0.99w*`(on)Stable 無鉸 → `1.01w*`(on)基端鉸→機構→Collapsed → 觸發比 **1.01**(rel 0)、凍結 residual
  **0.75Mp**(rel 0);**同載 off → Stable**(`|M|/Mp=0.7575<1`)= 軸力互動決定崩塌。

### audit testNMInteraction(+3,101→104)
- (a) 矩形 `Mp_eff` == 中性軸一次原理 rel 0;
- (b) `nmInteraction=false` 對**軸載**模型 == stage-4b 逐位元(outcome/steps/maxDC/triggerRatio 全等)= strict no-op;
- (c) 軸力互動決定性:同載 on=Collapsed、off=Stable。

### UE `FrameCore.Collapse.NMInteraction`(+1,49→50)
F54 懸臂括弧鏡像(公式點檢 + on 1.01w* Collapsed / off 同載 Stable + 凍結 residual)。

## 誠實聲明

- **非新方法**:G1 軸力降伏塑性鉸 = 教科書 sequential-linear(矩形包絡 EC3 §6.2.9);`[NEW CODE]` 僅「固定 Mp →
  Mp_eff(N) 替換 + 凍結」。**不自稱真彈塑性**(無卸載/反轉/N-M 切線)、**不自稱 pushover/纖維斷面**(WS_G 排除 G2/G3)。
- 邊界承襲崩塌驅動器:LSP 級 sequential linear、文獻 ±30%;N-M 使彎曲鉸載更貼真實(保守端)。

## 踩雷 / durable

- ⚠️ **上一輪 `Collapse.h` 的 `nmInteraction` 旗標+註解已落地但 `Collapse.cpp` 觸發判據沒落地、`NMInteraction.h`
  沒建**(Edit 顯示成功 ≠ 真落地)→ 開工先 `git status`/grep 核對(本輪即此情形)。
- ⚠️ 成鉸端 recover `M=0`(非 Mp)是鉸不重複觸發的關鍵機制(F32 line 1496 契約);故 S10 **不需**「跳過已成鉸端」守門。
- ⚠️ `(N/Ny)²` 平方使 `endI.N` 壓正/拉負符號無關 → 不必煩惱端力符號慣例。
- ⚠️ UE 測試:`std::vector` 用 `.size()`/`.empty()`(非 TArray `.Num()`);新 .cpp 由 UE Adaptive Build 自動編
  (排除出 unity);`$ExpectedUeTests` 計數守門可擋新測試漏編的假綠。
- ⚠️ 矩形 `Section::Rectangular(b,d)`:`b`=寬(local z)、`d`=深(local y);`Zz=b·d²/4`、`A=b·d`。

## 下一步

S9 家族 + S10 完成 → **co-rotational 大位移線 + 塑鉸材料非線性線皆收尾**。主線選項(待授權):S11 MITC9i
高階殼(殿後)、可視化資料線(C6/C7/C8)、UE5 視覺層。詳見 `docs/IMPLEMENTATION_PLAN.md` / `docs/KARAMBA3D_ROADMAP.md`。
