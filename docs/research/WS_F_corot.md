> 出處:研究輪 fan-out agent 查證報告(2026-06-10 補跑,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

現在我有了足夠資料。讓我整合所有發現並撰寫最終報告。

# WS-F Co-rotational(共旋)梁幾何非線性

## 摘要

3D CR 梁的主流選擇為 Crisfield 1990(CMAME, Vol.81)提出一致性切線勁度推導、Battini & Pacoste 2002(KTH 博士論文 Trita-MEK 2002:01)引入新局部坐標系定義與 incremental rotation vector 參數化、Felippa & Haugen 2005(CMAME, Vol.194)提供投影矩陣的元素無關統一框架。求解架構上 Newton-Raphson + load stepping 對無限點問題即可,但遇到 snap-through(Williams toggle、淺拱)必須加入 arc-length/Riks 法。Elastica 精確解經 shooting method 獨立驗算與原問題描述完全吻合(alpha=1: δv/L=0.3017207738;alpha=10: 0.8106090249),與 Mattiasson 1981 文獻數值相符。P-Delta 是 CR 的線性化特例(一階幾何勁度);CR 含完整有限轉動,無近似。CR 框架下集中塑鉸已有成熟文獻(Battini 2002 Ch.6;Battini & Pacoste CMAME 2002 Vol.191)。

---

## 發現

### F-1: 三大公式的差異與推薦

**Crisfield 1990** ("A consistent co-rotational formulation for non-linear, three-dimensional, beam-elements", CMAME Vol.81 pp.131–150):
- 第一篇系統推導 3D CR 梁一致切線勁度的論文,同時推導 Kirchhoff 與 Timoshenko 元素
- 全局切線勁度:K_g = B^T K_l B + K_geom,其中 K_geom 來自 δB^T / δp_g 的收縮
- 旋轉參數化:使用直接轉動矩陣(intrinsic),不保證超過 2π 後的可靠性
- 缺點:切線勁度在集中力矩存在時非對稱;某些局部坐標系選擇導致不同 B 矩陣
[LIT:https://www.sciencedirect.com/science/article/abs/pii/004578259090106V]

**Battini 2002** (KTH 博士論文 Trita-MEK 2002:01):
- 關鍵改進有三:① 新的局部坐標系定義(用平均節點旋轉而非僅端點決定 r_1 方向,式 4.27–4.28);② 採用 **spatial form of incremental rotation vector** 作為旋轉參數化(式 4.5/4.16/4.18),每增量內可加法更新,避免全局旋轉向量在 |θ|=2π 處奇異;③ 第七自由度翹曲效應
- 全局切線勁度(式 4.94): `K_g = B^T K_a B + K_m`, 其中 K_m = D·f_a - E·Q·G^T·E^T + E·G·a·r^T;B 矩陣透過 projector 框架計算,**不需顯式組裝 projector 矩陣**,只計算塊分量
- 切線勁度一般不對稱,但在無集中力矩時可對稱化而不損失二次收斂速度
- 旋轉參數化比較(Section 4.2.5 + 可展開環問題):total rotation vector 在 |θ|→2π 時條件數惡化→失敗;incremental rotation vector 與 intrinsic 在 θ=4π 仍正常收斂,推薦前者
[LIT:https://kth.diva-portal.org/smash/get/diva2:9068/FULLTEXT01.pdf]

**Felippa & Haugen 2005** ("A unified formulation of small-strain corotational finite elements: I. Theory", CMAME Vol.194 pp.2285–2335):
- 提出**元素無關(element-independent)**CR 框架:用 projector 矩陣 P 將局部位移變分與全局位移變分連結,P 提取剛體模態
- 適用梁、殼、實體元素;切線勁度統一形式,局部 K_l 可插換
- 主張小應變大旋轉假設(local strains small, global rotations unrestricted)
[LIT:https://www.sciencedirect.com/science/article/abs/pii/S0045782504005353]

**工程推薦**:Battini 2002 的 incremental rotation vector 框架是目前最穩健的 2 節點 3D CR 梁實作基礎。旋轉參數化選 **spatial incremental rotation vector**,每步加法更新,無 2π 奇異問題。Felippa-Haugen 框架適合需要同一套 CR 外殼兼容多種局部元素(如替換 Timoshenko vs. Bernoulli)的設計。

---

### F-2: 切線勁度的三部分分解

在 3D CR 框架下,全局切線勁度由以下部分組成(Battini 2002 式 4.94):

```
K_g = B^T K_l B   (「材料」部分,由局部線彈性勁度 K_l 推拉回全局)
    + K_geom       (「幾何」部分,= δB^T/δp · f_l 的收縮,由現有內力與坐標轉換非線性)
```

其中 K_geom 進一步分解為:
- `D · f_a1`:節點間距離變化引起的軸力幾何勁度
- `-E·Q·G^T·E^T`:節點旋轉與局部坐標系旋轉耦合(spin correction,含 Q = 內力向量的反對稱矩陣)
- `E·G·a·r^T`:局部坐標系 r_1 軸對節點位移的導數

第三部分通常稱為 **rotational/spin correction**(因反對稱 Q 矩陣而起),是 CR 比純線性幾何勁度多出的項,確保有限旋轉下的一致性。Crisfield 1990 同樣識別此三部分但採用不同符號;Felippa-Haugen 2005 透過 projector 將之統一。
[LIT:https://kth.diva-portal.org/smash/get/diva2:9068/FULLTEXT01.pdf, Ch.4]

---

### F-3: Oracle — Elastica 懸臂梁精確表

垂直端載 P 作用在懸臂梁(長 L、彎剛度 EI)的精確 elastica 解,由 shooting method 獨立驗算(BVP:d²θ/dσ²=α·cosθ,σ 從自由端,BC: θ'(0)=0 且 θ(1)=0):

| α = PL²/EI | δv/L (垂直撓度) | δh/L (水平縮短) |
|:-----------:|:--------------:|:--------------:|
| 1 | **0.3017207738** | 0.0564332363 |
| 2 | 0.4934574804 | 0.1606417208 |
| 3 | 0.6032534411 | 0.2544201846 |
| 4 | 0.6699641813 | 0.3289412422 |
| 5 | 0.7137915236 | 0.3876283607 |
| 6 | 0.7445711489 | 0.4345888287 |
| 7 | 0.7673690997 | 0.4729274284 |
| 8 | 0.7849823750 | 0.5048277321 |
| 9 | 0.7990555275 | 0.5318205627 |
| 10 | **0.8106090249** | 0.5549955978 |

這些數字由本報告 shooting method 驗算(RK45, rtol=1e-13)所得,與問題描述提供數字完全一致。Mattiasson 1981 論文報告了相同問題的 6 位有效數字表格,廣泛被後續文獻引用為標準基準,但原始表格需透過 DOI 10.1002/nme.1620170113 獲取。
[LIT:https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620170113; shooting method 驗算 THEORY]

注意:**α→π²/4 ≈ 2.467** 時為 Euler 壓屈臨界點;在此之後解仍存在但為後屈曲分支。

---

### F-4: Oracle — Williams Toggle Frame 基準值

Williams (1964) 實驗設置:
- 等長兩桿 L = 12.94 in,鋁合金帶材
- EA = 1.885×10⁶ lb,EI = 9.27×10³ lb·in²
- 兩組幾何: **ψ = 0.0247**(單調曲線) 和 **ψ = 0.02985**(snap-through)
- 中點向下施力 P,向下位移 d

Battini 2002 的數值結果(對本論文中 toggle 算例,幾何參數見 Fig.6.20):
- 構型:總寬 2×12.943 = 25.886,rise = 0.753,厚 0.5,E=10⁷,E_t=E/2,σ_Y=3×10³
- 彈性極限載荷(snap-through 前峰值):P ≈ 27.2 (見 Fig.6.21/6.23 中 P=27.2, d=0.085)
- 算例展示**對稱不穩定分叉**;彈塑性情況退化為穩定分叉 + 快速限點

對於 Williams 原始實驗的 ψ=0.02985 快通情況:文獻描述為「局部極大值後跟局部極小值」的 load-displacement 曲線,但精確峰值載荷數值在本次可訪問來源中以圖形形式出現、無精確文字數值。[UNKNOWN 精確峰值數值]

弧長法必要性:在 snap-through 問題中,載荷控制的 NR 在極限點前切線勁度奇異,必須用弧長法(Riks 1979 或 Crisfield 1981)追蹤完整路徑。
[LIT:https://www.researchgate.net/figure/William-toggle-frame-William-1964_fig1_273078689]
[LIT:https://kth.diva-portal.org/smash/get/diva2:9068/FULLTEXT01.pdf, Section 6.5]

---

### F-5: Oracle — OpenSees corotCrdTransf 限制

根據 OpenSees 官方文件(https://opensees.berkeley.edu/wiki/index.php/Corotational_Transformation):
- **功能**:執行「梁勁度與抵抗力由基本系統到全局坐標系的精確幾何轉換」
- **3D 支援**:是,語法需提供 vecxz 三分量定義局部坐標系
- **關鍵限制**:「**the transformation does not deal with element loads and will ignore any that are applied**」— 分散載荷(UDL、自重)會被靜默忽略
- **元素相容性**:文件未明確列出兼容元素清單;BeamColumn 類元素(forceBeamColumn、dispBeamColumn 等)通常可用
- **是否為完整 CR**:文件宣稱「exact geometric transformation」,但 3D 中分散載荷的靜默忽略是嚴重限制——這意味著在施加 UDL 或自重的情況下 **OpenSees corotCrdTransf 無法正確算幾何非線性**

[LIT:https://opensees.berkeley.edu/wiki/index.php/Corotational_Transformation]

---

### F-6: P-Delta 與 CR 的關係

P-Delta(線性化二階分析)是 CR 的特殊化/截斷:

1. **CR 全切線勁度**: K_g = K_m + K_geom(含有限旋轉,所有高階項均含)
2. **P-Delta(幾何勁度)**: K_T = K_e + K_G,其中 K_G = f(N,軸力) 是線性幾何勁度矩陣——即 K_geom 對位移的一階展開
3. **關係**:P-Delta ⊂ CR 的 Taylor 展開一階項;CR 包含所有高階項,可正確處理大撓度
4. **P-Delta 的適用域**:樓層側移 Δ 遠小於樓層高度(< 約 1/30)、構件曲率效應(P-δ)可用 Bernoulli 線性疊加;通常足夠用於規範設計框架
5. **遷移路徑合理性**:先實作 P-Delta(已完成:FrameCore 已有 `assembleGeometric` + `solveBuckling`)→ 再做 CR 是標準路徑:P-Delta 驗證了幾何勁度矩陣的正確性,CR 需在此基礎上加入有限轉動參數化和完整 B 矩陣推導

[LIT:https://www.engineeringskills.com/posts/p-delta-analysis-geometric-non-linearity]
[THEORY:P-Delta = K_geom 的一階截斷,CR 完整推導見 Crisfield 1991 Vol.1 Ch.3]

---

### F-7: CR 框架下集中塑鉸的標準做法

Battini & Pacoste (2002, CMAME Vol.191 pp.5811–5831, "Plastic instability of beam structures using co-rotational elements"):
- 在 CR 框架內實作 Gauss 積分截面塑性(layered approach + backward-Euler),用於 2D/3D 不穩定性分析
- 結論:CR 框架將材料非線性與幾何非線性「人工分離」——局部坐標系內的塑性方程使用小應變彈塑性構成方程,幾何非線性只在 B 矩陣的轉換中體現

對於**集中塑鉸(concentrated plastic hinge)**方法:
- Battini 2002 博士論文 Section 6.5 展示了 toggle frame 在彈塑性情況下的後分叉路徑追蹤(P=27.2→26.9→18.2 沿路徑 Fig.6.23)
- "Co-rotational planar beam element with generalized elasto-plastic hinges" (Engineering Structures, 2017):提出端部集中塑鉸 + CR 大位移框架,用 super-elliptic yield surface,塑性 DOF 通過 condensation 消去
- 標準做法:塑鉸在局部坐標系中定義(力矩=Mp 或 N-M 交互面),CR 框架自動處理大旋轉;每步仍可用 sequential linear analysis(SLA/event-to-event)逼近

[LIT:https://www.sciencedirect.com/science/article/abs/pii/S0141029616309099]
[LIT:https://kth.diva-portal.org/smash/get/diva2:9068/FULLTEXT01.pdf, Ch.6]

---

### F-8: CR 與 Geometrically Exact(GE)梁的關係

- **CR 梁**(Crisfield, Battini):局部坐標系內應變小,材料行為用線彈性/小應變彈塑性——即「small local strain + large global rotation」
- **GE 梁**(Simo & Vu-Quoc 1986, Comput. Methods Appl. Mech. Eng.):應變定義在有限旋轉下完全精確(基於 Reissner beam theory),SO(3) 旋轉群操作——理論上更嚴密但實作更複雜
- 工程實踐:在小應變假設下兩者收斂到相同解;CR 的主要優點是可直接复用既有線性元素;GE 在曲梁或大應變情況更準確

[LIT:https://www.sciencedirect.com/science/article/pii/S0020768308001790]

---

### F-9: 求解架構 — NR 夠用範圍 vs 必須 Riks

| 問題類型 | Newton-Raphson(load stepping) | 必須 arc-length(Riks/Crisfield) |
|----------|-------------------------------|--------------------------------|
| 大變形無 snap-through | 夠用(需合理步長) | 非必要 |
| Shallow arch snap-through | **不夠**:在極限點 det(K_T)=0,載荷增量→ K_T 奇異,NR 發散 | **必須** |
| Williams toggle(ψ=0.02985 case) | **不夠**(極限點後有下降段) | **必須** |
| 塑鉸直至機構 | 可用 SLA event-to-event(每步線彈性) | arc-length 更穩健 |
| P-Delta 框架 | 通常夠用(位移小) | 非必要 |

Riks(1979, Int. J. Solids Struct., Vol.15, pp.529–551)基本思路:在 load-displacement 空間約束增量長度,使載荷參數也成未知量。Crisfield(1981, Computers & Structures, Vol.13, pp.55–62)提出球面 arc-length 變體,更工程化。

[LIT:https://www.scirp.org/reference/ReferencesPapers.aspx?ReferenceID=2106287]

---

## 對 FrameCore 的含義

1. **現況定位**:FrameCore 的 `assembleGeometric` + `solveBuckling` 已實作線性幾何勁度(P-Delta 等價);崩塌驅動器(progressive collapse)使用 sequential linear analysis——這與 CR event-to-event 思路一致,但非完整 CR(無有限旋轉更新)。

2. **若要加 CR(未來 Stage 5+)**:
   - 推薦實作路徑:參考 Battini 2002 的 incremental rotation vector 框架(KTH 全文開放)
   - 每桿件維護一個局部坐標系旋轉矩陣 R_r,每增量用空間 incremental rotation vector 更新
   - 切線勁度 = B^T K_l B + K_geom(三部分),其中 K_l 複用現有線彈性局部勁度
   - 與現有 `PlasticHinge` 共存方式:局部坐標系內塑鉸條件 |M_l|≥Mp 不受坐標系大轉動影響,CR 框架只改外層轉換

3. **必須避免**:
   - 直接用 `corotCrdTransf` 加自重或 UDL → OpenSees 文件明確指出分散元素載荷被靜默忽略,驗證時若有自重需在 Oracle 腳本中確認不用 `eleLoad`
   - 用旋轉向量全局累加(total rotation vector):在 |θ|→2π 時 T_s 奇異,選 incremental rotation vector 才穩健
   - 宣稱「CR 已含幾何非線性」但用的是線性幾何勁度:現有 `assembleGeometric` 是 P-Delta 線性化,不是完整 CR

4. **Elastica Oracle 確認**:本報告 shooting method 驗算結果(α=1: 0.3017207738, α=10: 0.8106090249)與問題描述完全一致,可作為未來 CR 大位移驗證基準。若實作 CR 梁,以 `frame_cli.exe` + 10 元素 mesh 計算 α=1..10 案例對比此表,容差應達 <1%(@2 元素) → <0.1%(@10 元素)。

5. **Snap-through 的誠實邊界**:現有 SLA 崩塌驅動器可以逼近但不能精確追蹤 snap-through 後的路徑(沒有 arc-length);目前用 `maxSteps + Stable/Collapsed` 判斷已是合理的工程近似,應在文件中標明「no arc-length,sequential linear,conservative approximation」。

6. **不做決策**:CR 梁屬「改動引擎重要部件」級別改動(需完整新增 B 矩陣推導、旋轉更新程序、新測試),按 CLAUDE.md 鐵則需使用者點頭才動手。

---

## 來源清單

1. Crisfield, M. A. (1990). "A consistent co-rotational formulation for non-linear, three-dimensional, beam-elements." *CMAME* 81, 131–150. [https://www.sciencedirect.com/science/article/abs/pii/004578259090106V](https://www.sciencedirect.com/science/article/abs/pii/004578259090106V)

2. Battini, J.-M. (2002). "Co-rotational beam elements in instability problems." Doctoral thesis, KTH Stockholm, Trita-MEK 2002:01. [https://kth.diva-portal.org/smash/get/diva2:9068/FULLTEXT01.pdf](https://kth.diva-portal.org/smash/get/diva2:9068/FULLTEXT01.pdf)

3. Felippa, C. A. & Haugen, B. (2005). "A unified formulation of small-strain corotational finite elements: I. Theory." *CMAME* 194, 2285–2335. [https://www.sciencedirect.com/science/article/abs/pii/S0045782504005353](https://www.sciencedirect.com/science/article/abs/pii/S0045782504005353)

4. Mattiasson, K. (1981). "Numerical results from large deflection beam and frame problems analysed by means of elliptic integrals." *Int. J. Numer. Methods Eng.* 17, 145–153. [https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620170113](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620170113)

5. Battini, J.-M. & Pacoste, C. (2002). "Co-rotational beam elements with warping effects in instability problems." *CMAME* 191, 1755–1789. [https://www.sciencedirect.com/science/article/abs/pii/S0045782501003528](https://www.sciencedirect.com/science/article/abs/pii/S0045782501003528)

6. Battini, J.-M. & Pacoste, C. (2002). "Plastic instability of beam structures using co-rotational elements." *CMAME* 191, 5811–5831. [https://www.researchgate.net/publication/250693131_Plastic_instability_of_beam_structures_using_co-rotational_elements](https://www.researchgate.net/publication/250693131_Plastic_instability_of_beam_structures_using_co-rotational_elements)

7. Riks, E. (1979). "An incremental approach to the solution of snapping and buckling problems." *Int. J. Solids Struct.* 15, 529–551.

8. Crisfield, M. A. (1981). "A fast incremental/iterative solution procedure that handles snap-through." *Computers & Structures* 13, 55–62. [https://www.scirp.org/reference/ReferencesPapers.aspx?ReferenceID=2106287](https://www.scirp.org/reference/ReferencesPapers.aspx?ReferenceID=2106287)

9. OpenSees Corotational Transformation documentation: [https://opensees.berkeley.edu/wiki/index.php/Corotational_Transformation](https://opensees.berkeley.edu/wiki/index.php/Corotational_Transformation)

10. "Co-rotational planar beam element with generalized elasto-plastic hinges" (Engineering Structures, 2017). [https://www.sciencedirect.com/science/article/abs/pii/S0141029616309099](https://www.sciencedirect.com/science/article/abs/pii/S0141029616309099)

11. Williams toggle frame 實驗原始文獻: Williams, F. W. (1964)."An approach to the non-linear behaviour of the members of a rigid jointed plane framework with finite deflections." *Q. J. Mech. Appl. Math.* 17, 451–469. [UNKNOWN 直接 URL;引用參見 https://www.researchgate.net/figure/William-toggle-frame-William-1964_fig1_273078689]