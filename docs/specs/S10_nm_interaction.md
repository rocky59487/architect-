# S10 規格 — N-M 互動塑鉸(軸力-彎矩互動降伏)

> 接續 S9 家族(S9/S9b/S9c co-rotational,收尾 `67db412`)。研究依據 `docs/AGENT_PROMPT_S5_S11.md` §S10、
> `docs/specs/S5_S11_skeletons.md` §S10(WS_G 路線 G1)。**S10 必在 S9 後**(R4:塑鉸方向耦合;本實作的塑鉸
> 在局部端力座標定義,與大轉動正交,故 S9 完成後可獨立銜接)。

## 目標 / 不做(誠實邊界)

- **做**:G1 路線 `Mp_eff(N) = Mp · max(0, 1−(N/Ny)²)`,`Mp=fy·Z`、`Ny=fy·A`(軸力-彎矩互動式),銜接既有
  event-to-event 塑鉸驅動器(stage 4b):觸發判據與 residual 從固定 `Mp` 改為 `Mp_eff(N)`。
- **不做**:G2 集中彈簧塑性 / G3 纖維斷面(理由與成本見 WS_G);無卸載/反轉(仍 sequential linear,**非真彈塑性**);
  無 My-Mz 雙軸彎矩耦合(各彎曲軸用同一軸力項各自降);無 N-M 切線耦合。

## 力學 / 演算法

- **降伏式**:`Mp_eff(N) = Mp·(1−(N/Ny)²)`,負值截 0(`reducedPlasticMoment`,`NMInteraction.h`,header-only inline)。
  - **矩形實心截面 = 一次原理精確**:塑性中性軸偏移,中央深度 `2c=N/(b·fy)` 帶承軸力,殘餘彎矩
    `M_N = Mp − N²/(4·b·fy)`,代數恆等於 `Mp(1−(N/Ny)²)`(教科書塑性 N-M 包絡,EC3 EN1993-1-1 §6.2.9 實心矩形)。
  - **圓形實心截面 = 保守**:真包絡較此拋物線飽滿 → 預測之鉸不晚於實際(崩塌偏安全)。
  - **AISC H1.1** 是更保守的雙線性**設計**檢核,非此包絡 — 不混為一談;此為截面塑性容量,非規範驗收比。
  - **單軸**:`(N/Ny)²` 平方項 → 拉/壓同 `|N|=Ny` 降伏、且 N 符號天然無關。
- **觸發判據改動**(`Collapse.cpp`,候選迴圈):每端/每軸 `Mpk = reducedPlasticMoment(Mp0[k], N_end, Ny)`,
  端對齊軸力(dof 4/5 用 `endI.N`、dof 10/11 用 `endJ.N`);`ratio = |M|/Mpk`,`ratio≥1` 成鉸,residual = `±Mpk`。
  `Mpk≤0`(軸力達 squash)→ 跳過該鉸候選,改由脆性軸力篩主導(誠實:全軸力降伏是軸力破壞非延性彎曲鉸)。
- **凍結於形成**:residual 存入 `PlasticHinge.Mp`(stage 4a 既有機制),成鉸後該端 recover 出 `M=0`(F32 契約)→
  後續步驟不再重估、不重複觸發。capacity 凍結是 sequential-linear 的自然結果。

## API / 架構

- **新增**:`CollapseOptions::nmInteraction`(`bool`,預設 `false`)。`Public/FrameCore/NMInteraction.h`(inline
  `reducedPlasticMoment(Mp,N,Ny)`,純 POD/零 Eigen)。
- **預設 `false` = 餵 `N=0` → `Mp_eff==Mp` 逐位元同 stage-4b**;不新增 element / dispatch / SolveResult 欄位 /
  **不入 `modelFingerprint`**(純驅動器後處理旗標,不改 K);**不新增 `.cpp`**(header-only)→ 四 build 腳本免動。
- 只在 `plasticHinges==true` 時有意義(non-hinge-capable 桿、脆性殼維持原樣)。

## Oracle(F54,實測全綠)

1. **公式**:`Mp_eff(0)=Mp`、`(Ny/2)=0.75Mp`、`(Ny)=0`、`|N|` 符號對稱、`2Ny` 截 0 — 機器精度(rel=0)。
2. **矩形一次原理精確**:`reducedPlasticMoment` vs `Mp−N²/(4bfy)` 中性軸偏移 — rel 0(tol 1e-12)。
3. **驅動器括弧法**:靜定 X 懸臂(基端 `Mz=wL²/2`、軸力=端載 `P`,皆 stiffness-independent)。`P=Ny/2 → Mp_eff=0.75Mp`;
   - `0.99 w*`(N-M on)→ Stable 無鉸;`1.01 w*`(N-M on)→ 基端鉸 → 機構 → Collapsed;
   - **同載 N-M off → Stable**(`|M|/Mp=0.7575<1`)→ 證明軸力互動是崩塌的決定因素;
   - 觸發比 = 1.01(rel 0)、凍結 residual `|Mp|=Mp_eff(P)=0.75Mp`(rel 0)。
- **audit(+3,checks 101→104)**:(a) 矩形精確 rel 0、(b) `nmInteraction=false` 對軸載模型**逐位元同 stage-4b**
  (strict no-op)、(c) 軸力互動決定性(同載 on=Collapsed / off=Stable)。
- **UE `FrameCore.Collapse.NMInteraction`**(+1,`$ExpectedUeTests` 49→50):F54 懸臂括弧鏡像。

## 五腿 gate

standalone **F1-F54** / UE **50** / OpenSees PASS / audit **104** / CLI round-trip ALL PASS。

## 誠實聲明 / novelty

- **非新方法**:G1 軸力降伏的塑性鉸是教科書 sequential-linear(矩形包絡 = EC3 §6.2.9 / 經典塑性理論);
  `[NEW CODE]` 僅 = 把既有 stage-4b 驅動器的固定 `Mp` 換成 `Mp_eff(N)` 並凍結。**不自稱真彈塑性**(無卸載/反轉/
  N-M 切線);**不自稱 pushover/纖維斷面**(刻意排除,WS_G G2/G3)。
- 與既有「崩塌驅動器是 LSP 級 sequential linear,文獻 ±30%」邊界一致;N-M 互動使彎曲鉸載更貼近真實(偏保守端)。

## 下一步

S9 家族 + S10 完成 → **co-rotational 大位移線 + 材料非線性(塑鉸)線皆收尾**。剩餘主線選項(待授權):
S11 MITC9i 高階殼(殿後,9 處引擎修改先決)、可視化資料線(C6 BMD/SFD、C7 利用率場、C8 贅餘度)、UE5 視覺層。
詳見 `docs/IMPLEMENTATION_PLAN.md`、`docs/KARAMBA3D_ROADMAP.md`。
