# WS-F2 — 3D 通用 Co-rotational 梁的逐行可實作公式(供 S9b)

> 出處:S9 研究輪 fan-out agent(2026-06-12),錨定 **OpenSees `CorotCrdTransf3d.cpp` 開源逐行** +
> Crisfield 1990 (CMAME 81) + Battini 2002 (KTH 博論)。S9 v1 只做**平面** CR(WS_F);本檔保存 **3D
> 通用** CR 的完整公式,供 S9b 直接動工。標記:**[OS]**=OpenSees 逐行可查、**[CR90]/[BAT02]**=理論來源、
> **[整理]**=本報告符號統一。

## 來源
- [OS] `https://github.com/zhuminjie/OpenSeesLite/blob/master/SRC/coordTransformation/CorotCrdTransf3d.cpp`
- [CR90] Crisfield 1990, *CMAME* 81:131-150;[BAT02] Battini 2002, KTH `diva2:9068`。

## 0. 符號
節點 I,J;當前座標 `x_I,x_J`;`L0/Ln` 初始/當前長;`R_I,R_J∈SO(3)` 節點旋轉矩陣;元素框架
`E=[e1|e2|e3]`(列向量);全域 12-DOF `q=[u_I(3) θ_I(3) u_J(3) θ_J(3)]`(與 `localStiffness12` 同序)。

## 1. 旋轉更新(spatial incremental — 避 2π 奇異)
`R_node ← exp(skew(Δθ))·R_node`(左乘空間增量;Δθ = NR 步的旋轉 DOF 修正)。
Rodrigues:`exp(skew(θ)) = I + (sinθ/θ)skew(θ) + ((1−cosθ)/θ²)skew(θ)²`。θ→0 級數:`sinθ/θ≈1−θ²/6`、
`(1−cosθ)/θ²≈½−θ²/24`(ε<1e-7 切換)。log:`θ=arccos((tr R−1)/2)`、`log R=θ/(2sinθ)(R−Rᵀ)`(θ→0
取 ½;θ→π 用四元數)。**為何不用 total rotation vector**:其切映射 `T_s(Ψ)` 在 `|Ψ|→2kπ` 奇異 [BAT02 §2.2]。

## 2. 元素框架 E [OS update() L574-644]
- `e1 = (x_J+u_J − x_I−u_I)/Ln`。
- 平均三元組 `Rbar` = `R_I→R_J` 測地中點:`dR=R_J R_Iᵀ` → 四元數半角 `w_half=w(dR)/2` → `Rbar=R(w_half)R_I`。
- `r1,r2,r3 = Rbar` 各列;`e2 = r2 − ((r2·e1)/(1+r1·e1))(e1+r1)`、`e3 = r3 − ((r3·e1)/(1+r1·e1))(e1+r1)`
  (Gram-Schmidt 對 e1;分母 `1+r1·e1≈0` = 翻轉 180° 需分支)。

## 3. 局部變形 [OS update() L632-643]
`u_bar = Ln − L0`;6 個局部節點旋轉(扭轉@I、繞 e3/e2 彎@I、扭轉@J、繞 e3/e2 彎@J):
`ul[k] = asin( ½(r_α^T S(e_β) e_γ − r_β^T S(e_γ) e_α) )`(α,β,γ 按各 DOF 循環;`r`=R_I/R_J 列,`S`=skew)。
**剛體零內力**:剛體 R_g 下 `R_{I,J}=R_g R_0`、`E=R_g E_0`,叉積中 R_g 抵消 → `ul=0,u_bar=0` [CR90 §2]。

## 4. 內力 f_int(12×1)[OS getGlobalResistingForce() L1303-1352]
`pb` = natural 力(7×1:6 端力矩/扭矩 + 軸力 N)= `K_local·[ul;u_bar]`。
`f_int = Tᵀ·Tpᵀ·pb`(虛功一致;T 見 §5)。軸力沿 e1、端彎矩經框架轉全域、端彎矩生等效橫向剪分配兩節點。

## 5. 轉換矩陣 T(7×12)[OS compTransfMatrixBasicGlobal() L817-1046]
`A = (I − e1 e1ᵀ)/Ln`;`getLMatrix(ri)`(L1430-1515):
`L1 = (ri·e1/2)A + A ri (e1+r1)ᵀ/2`、`L2 = skew(ri)/2 − (ri·e1/4)skew(r1) − skew(ri)e1(e1+r1)ᵀ/4`、
`L(12×3)=[L1;L2;−L1;L2]`。T 各行(扭/彎@I、扭/彎@J、軸)由 `A·r{I,J}{2,3}`、`Lr2/Lr3` 與 `/(2cosθ)` 組裝
(軸力行 = `[−e1ᵀ,0,e1ᵀ,0]`)。

## 6. 一致切線 Kt [OS getGlobalStiffMatrix() L1355-1596]
`Kt = Tᵀ Kl T + Ksigma1 + Ksigma2 + Ksigma3 + Σ_j m_j tan(ul_j) T_jᵀT_j`,`m_j=pl_j/(2cos ul_j)`。
- `Ksigma1 = N·[[A,0,−A,0],[0,0,0,0],[−A,0,A,0],[0,0,0,0]]`(軸力幾何,**必加**;對稱)。
- `Ksigma3`(彎矩-軸耦合,Lr2/Lr3·m)、`Ksigma2`(spin 修正,`getKs2Matrix` L1830-1849,非對稱來源)。
- 無集中外力矩(conservative)時 Kt 對稱;OpenSees 保留非對稱(LU),亦可 `(Kt+Ktᵀ)/2`(小步長誤差可忽略)。
- **務實**:Ksigma3/Ksigma2 對中小轉角影響小;NR 收斂解只取決於 `f_int`(殘差),切線只影響收斂速度 →
  S9b 可先 `TᵀKlT + Ksigma1` 求位移正確,再補完整 Kgeom 提收斂速度。

## 7. S9b oracle 增量
3D 任意軸剛體旋轉 patch(`f_int≈0` 機器精度)、空間(非平面)elastica、OpenSees `corotCrdTransf` 對照
(⚠️ **忽略 element loads → 只用節點力、勿 eleLoad**;3D 需 vecxz)。平面退化 == S9 v1 逐位元(回歸保護)。
**S10(N-M 塑鉸)在 S9b 後**(R4 方向耦合)。

## 8. 陷阱表
e2/e3 `1+r1·e1` 奇異(翻轉 180°,四元數備援)、`asin` 出域(小步長 |ul|<π/2)、Rodrigues θ→0/π(級數/
四元數)、`1/(2cosθ)` θ→π/2(步長控制)、Kt 非對稱(LU 或對稱化)、更新順序(先 R 後 u)。
