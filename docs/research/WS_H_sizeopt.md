> 出處:研究輪 fan-out agent 查證報告(2026-06-10 補跑,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

現在我已掌握所有核心資料,可以整合撰寫報告。

# WS-H 尺寸優化(fully-stressed design 與規範檢核)

## 摘要

Fully-Stressed Design(FSD)是以 stress-ratio resizing 迭代使每根構件在至少一種工況下達到應力限值的啟發式優化手法。對靜定結構已有嚴格最優性證明,對靜不定結構不保證最優且理論上可能不收斂,但實務中收斂快速、結果接近最優。經典 10-bar truss 基準(僅應力約束)文獻最優解 1593.2 lb,我方 FSD 實驗收斂 1593.16 lb 與文獻吻合。Karamba3D v3 的 Optimize Cross Section 以 EC3 EN 1993-1-1(含 buckling、LTB)對截面表逐一掃描選段,預設 5 次 ULS 迭代,以「不再改動截面」為收斂準則,最大化利用率靶值(MaxUtil 預設 1.0)。離散截面選取的震盪問題有充分文獻記載,最簡單且保守的緩解策略是永遠選比連續最優面積大一號的截面(向上取整)。EC3 LTB 完整檢核需要 Iw(翹曲常數)、It(扭轉常數)、Iz,目前 FrameCore `Section` 只有 A/Iy/Iz/J/Asy/Asz/Zy/Zz,缺 Iw 與重心至剪力中心距離 zg,無法做完整 EC3 6.3.2 檢核。多工況 FSD 標準做法是取每根構件跨全部工況的最大 D/C 作為 resizing 依據。

---

## 發現

### 1. FSD 收斂性質

**靜定結構**:內力與截面積無關,stress-ratio 一步迭代即達到全局最優解,可嚴格證明 FSD 解 = 最輕解。數學依據:目標函式 W = Σ ρ_i L_i A_i 在固定力 F_i 下最小化,令每根 A_i = |F_i| / σ_allow 直接得到最優。  
[LIT: Haftka R.T., Gürdal Z., *Elements of Structural Optimization*, Kluwer/Springer, 1992, §2.1 pp. 17–20; §7.2 pp. 236–238]  
[LIT: Schmit L.A., "Structural Design by Systematic Synthesis," *Proc. 2nd Conf. Electronic Computation*, ASCE, 1960, pp. 105–132]

**靜不定結構**:截面積改變時剛度改變,內力重分布,FSD 不保證最優解,存在多個不同的 fully-stressed 固定點,部分固定點在 stress-ratio 迭代下為不穩定點(發散)。  
[LIT: Haftka & Gürdal 1992, §7.2–7.3]  
[LIT: Mueller C., Liu J., *Fully Stressed Design of Frame Structures and Multiple Load Paths*, ASCE J. Struct. Eng., 2002, 128(6): 806–814, https://ascelibrary.org/doi/10.1061/(ASCE)0733-9445(2002)128:6(806)]  
[LIT: Burns S.A. et al., *Multiple Fully Stressed Structural Designs and the Stress Ratio Method*, Comput.-Aided Civ. Infrastruct. Eng., 1995, 10(4): 285–295, https://onlinelibrary.wiley.com/doi/10.1111/j.1467-8667.1995.tb00386.x]

**全局收斂定理(限 truss)**:Mäkinen K. et al. 1997 年證明:對受應力約束的 truss sizing「perturbation」版 stress-ratio 方法,迭代序列永遠收斂到全局最優。該結論限於 truss(軸力為主)且需加小擾動;一般框架不在此保證範圍內。  
[LIT: Mäkinen K. et al., "Global convergence of the stress ratio method for truss sizing," *Struct. Multidiscip. Optim.*, 1997, 14(1): 55–58, https://link.springer.com/article/10.1007/BF01742935]

### 2. 10-bar truss 基準(應力僅約束)

問題定義:2 個 bay 各 360 in、E = 10^7 psi、σ_allow = ±25 000 psi、P = 100 kip × 2(節點 2/4)、ρ = 0.1 lb/in³、A_min = 0.1 in²、無位移約束。

**文獻最優解(Design Problem 1)**:  
最優重量 = **1593.2 lb**  
最優截面(in²):[7.94, 0.10, 8.06, 3.94, 0.10, 0.10, 5.74, 5.57, 5.57, 0.10]  
桿 2/5/6/10 貼下界,桿 1/3/7/8/9 完全受力。  
[LIT: Haftka R.T., Gürdal Z., 1992, pp. 238 & 244 — 由 NASTRAN CoFE 頁面引用確認]  
[LIT: NASTRAN CoFE Example, https://vtpasquale.github.io/NASTRAN_CoFE/2._Examples/b._Optimization/1._10-Bar_Truss_Sizing/]

我方 FSD 實驗值 1593.16 lb 與文獻 1593.2 lb 差距 0.04 lb(<0.003%),在 Amin bound 容忍範圍內,視為完全吻合。文獻常見引用「1593.18」[UNKNOWN — 查無對應單一文獻數字,可能為不同收斂容差版本],1593.2 是 Haftka & Gürdal 原著數字。

加入位移約束(±2 in)的 Design Problem 2 最優解 = 5060.85 lb,截面大幅增加,屬不同問題。

### 3. 離散截面選取的震盪問題與緩解

連續最優解向下取整可能違反約束,向上取整保守但收斂;交替在兩個截面間跳躍是標準震盪模式。

已記錄的緩解策略:
- **永遠向上取整(round-up)**:保守、不震盪、不保證最輕但一定可行。最簡單實作。
- **阻尼/移動限制(damping / move limits)**:限制每步只升降一個截面等級,減少大跨跳躍。
- **偽離散法(pseudo-discrete rounding)**:先做連續優化,再逐步固定離散變量並重優化餘下部分。[LIT: Thanedar P.B. et al., "A pseudo-discrete rounding method for structural optimization," *Struct. Optim.*, 1992, https://www.researchgate.net/publication/226284044]
- **Karamba3D 做法**:截面族依面積(重量)排序;每輪選「第一個夠用的截面」等效於 round-up;MaxUtil 閥值讓截面在「足夠用」時停下,不強求完全滿利用。

[LIT: Erbatur F. & Al-Hussainy M.M., "Oscillation problem in element optimization with discrete sections and control," *Comput. Struct.*, 1989, 30(3): 537–543, https://www.sciencedirect.com/science/article/abs/pii/004579498990240X]

### 4. Karamba3D Optimize Cross Section 手冊細節

規範:**EC3 EN 1993-1-1**(考慮 buckling + LTB),互動係數依 Annex B 計算(k_yy, k_yz);material safety factors γ_M0 = 1.0(非屈曲主導)、γ_M1 = 1.1(屈曲主導);SwayFrame = true 時 C_my/C_mz 限制不低於 0.9。

迭代邏輯(ULS):
1. 用當前截面計算各構件截面力
2. 從截面族頭部掃描,選第一個滿足 MaxUtil 閾值的截面
3. 若無任何修改或已達 Util Iter(預設 5)次則停止;否則重複

靜不定結構因截面剛度互相耦合需要多輪迭代。多工況:所有 LCasesUtil 指定的工況同時評估,以最大利用率決定選段。

截面族組織:截面表有三種排序版本:依高度(`CrossSectionValues.bin`、`_sortedForHeight.bin`)或依面積/重量(`_sortedForWeight.bin`,用此版本得到近似最輕設計)。

SLS(位移控制):以虛擬載重法(virtual load methodology)在 Disp Iter 次內調整構件以滿足位移限制,先於 ULS 迭代。

LTB 長度:由程序自動計算(沿構件追蹤直到遇到 3 個以上構件相交的節點);懸臂末端 × 2;可透過 BklLenLT 手動覆寫。

截面分類:檢核 EN 1993-1-1 Class 1–4;Class 4 給出警告。

[LIT: Karamba3D v3 Manual §3.6.8, https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.8-optimize-cross-section]  
[LIT: Karamba3D v3 Manual Appendix A.4.6, https://manual.karamba3d.com/appendix/a.4-background-information/a.4.6-approach-used-for-cross-section-optimization]

### 5. LTB 檢核需要的截面性質 vs FrameCore 現有 Section

EC3 6.3.2 LTB 簡化法彈性臨界彎矩 Mcr 公式(NCCI SN003a-EN-EU / Timoshenko-Vlasov 推導):

```
Mcr = C1 × (π²EIz)/(kL)² × sqrt( Iw/Iz + (kL)²·G·It/(π²EIz) )
```

LTB 還原係數 χ_LT 計算流程:  
λ_LT = sqrt(Wy · fy / Mcr) → Φ_LT → χ_LT → M_b,Rd = χ_LT · Wy · fy / γ_M1

所需截面性質完整清單:

| 性質 | 符號 | FrameCore 現有? |
|---|---|---|
| 主軸截面模數 | Wy (Wpl,y 或 Wel,y) | 有(Zy/Zz 為塑性模數) |
| 弱軸慣性矩 | Iz | 有 |
| St. Venant 扭轉常數 | It (= J in FrameCore) | 有(J) |
| 翹曲常數 | Iw (Cw) | **缺** |
| 重心至剪力中心距 | zg | **缺** |
| 材料強度 | fy | 有 |
| 彎矩梯度係數 | C1, C2 | 需由載重形式計算(不是截面性質) |

矩形實心截面的翹曲常數 Iw 很小(≈ b³h³/144),實務上 LTB 不敏感,通常直接忽略或套簡化公式。I-形截面 Iw = Iz · (h-tf)² / 4(雙對稱),大型薄壁截面不可略。  
[LIT: EN 1993-1-1:2005 Clause 6.3.2, https://www.phd.eng.br/wp-content/uploads/2015/12/en.1993.1.1.2005.pdf]  
[LIT: Structville 2020, https://structville.com/2020/10/lateral-torsional-buckling-of-steel-beams-according-to-eurocode-3.html]  
[LIT: Dlubal RFEM Knowledge Base #001512, https://www.dlubal.com/en/support-and-learning/support/knowledge-base/001512]

AISC 360 F 章對應性質清單:Cw(= Iw)、J、Iy、ho(翼板間距)、rts(LTB 旋轉半徑)、Sx(Wpl,x 可替代)——缺漏集合與 EC3 相同,本質都需要 Iw/Cw。  
[LIT: AISC 360-22 Chapter F, Table F1-1]

### 6. 多工況 FSD:envelope D/C 標準做法

標準 stress-ratio resizing 規則多工況版本:

```
A_i^(new) = A_i^(old) × max_j( |σ_ij| / σ_allow )
```

其中 j 枚舉所有工況,|σ_ij| 是工況 j 下構件 i 的最大應力。取 max 確保在最嚴格工況下恰好達到允許應力上界。

對含壓力屈曲約束的框架,亦可對 D/C(demand/capacity 比)做同樣的 max-envelope:

```
A_i^(new) = A_i^(old) × max_j( DC_ij )
```

理論依據:多工況 FSD 定義「每根構件在至少一個工況下完全受力」;取 envelope 正是讓最嚴工況驅動截面,其餘工況自動滿足。Karamba3D 即採此路徑:所有 LCasesUtil 工況並行評估,取最大利用率。

[LIT: Mueller & Liu 2002, ASCE J. Struct. Eng. 128(6)]  
[LIT: Haftka & Gürdal 1992, §7.3]

---

## 對 FrameCore 的含義

**採用**:
- FSD 作為尺寸優化前處理或獨立工具是合理的;多工況直接取 envelope max(DC_ij) 驅動 resizing,與文獻完全一致。
- 10-bar truss 1593.2 lb 可作為 FSD 模組的標準回歸測試,我方 1593.16 lb 已通過。
- 對靜定結構可宣稱 FSD = 最優;對靜不定結構必須誠實標「啟發式,不保證全局最優」。
- 離散截面選取用 round-up(向上取整到下一號截面)是最簡單且安全的策略,符合 Karamba3D 的做法。

**避免**:
- 宣稱靜不定結構 FSD 是「最優解」——文獻明確否定。
- 不帶阻尼直接硬跳離散截面而不加任何震盪抑制機制——可能在兩個截面間無限振盪。
- 把 Karamba3D 預設迭代次數(5次)當做收斂保證——對大型靜不定結構 5 次可能不夠,需讓使用者可調。

**誠實邊界**:
- FrameCore `Section` 現只有 A/Iy/Iz/J/Asy/Asz/Zy/Zz;**缺 Iw(翹曲常數)**,無法實作完整 EC3 6.3.2 LTB 檢核或 AISC 360 完整側扭挫曲。目前只能做截面強度 D/C(彎矩/剪力/軸力)加上線性屈曲(P-Δ)。若要對標 Karamba3D 完整 EC3 合規設計,需在 `Section` 加入 Iw 欄位並提供 I 形截面閉式公式(Iw = Iz(h-tf)²/4)。矩形實心截面的 LTB 通常可忽略(Iw 極小),但 I 形截面不可略。

---

## 來源清單

- [Haftka R.T., Gürdal Z., *Elements of Structural Optimization*, Springer, 1992 — Google Books](https://books.google.com/books/about/Elements_of_Structural_Optimization.html?id=CzgIpexeh7UC)
- [NASTRAN CoFE 10-bar truss example (cites Haftka 1992 pp.238,244)](https://vtpasquale.github.io/NASTRAN_CoFE/2._Examples/b._Optimization/1._10-Bar_Truss_Sizing/)
- [Mäkinen K. et al., "Global convergence of the stress ratio method for truss sizing," *Struct. Multidiscip. Optim.* 1997](https://link.springer.com/article/10.1007/BF01742935)
- [Mueller C., Liu J., *Fully Stressed Design of Frame Structures and Multiple Load Paths*, ASCE J. Struct. Eng., 2002](https://ascelibrary.org/doi/10.1061/(ASCE)0733-9445(2002)128:6(806))
- [Burns S.A. et al., *Multiple Fully Stressed Structural Designs and the Stress Ratio Method*, CACAIE, 1995](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-8667.1995.tb00386.x)
- [Erbatur F., Al-Hussainy M.M., "Oscillation problem in element optimization with discrete sections," *Comput. Struct.* 1989](https://www.sciencedirect.com/science/article/abs/pii/004579498990240X)
- [Thanedar P.B. et al., "A pseudo-discrete rounding method for structural optimization," *Struct. Optim.* 1992](https://www.researchgate.net/publication/226284044_A_pseudo-discrete_rounding_method_for_structural_optimization)
- [Patnaik S.N., "Optimality of a Fully Stressed Design," NASA TM-1998-207411](https://ntrs.nasa.gov/api/citations/19980148007/downloads/19980148007.pdf)
- [Karamba3D v3 Manual §3.6.8 Optimize Cross Section](https://manual.karamba3d.com/3-in-depth-component-reference/3.5-algorithms/3.5.8-optimize-cross-section)
- [Karamba3D v3 Manual Appendix A.4.6 Approach Used for Cross Section Optimization](https://manual.karamba3d.com/appendix/a.4-background-information/a.4.6-approach-used-for-cross-section-optimization)
- [EN 1993-1-1:2005 Clause 6.3.2 (full standard PDF)](https://www.phd.eng.br/wp-content/uploads/2015/12/en.1993.1.1.2005.pdf)
- [Structville — Lateral-Torsional Buckling of Steel Beams per Eurocode 3](https://structville.com/2020/10/lateral-torsional-buckling-of-steel-beams-according-to-eurocode-3.html)
- [Dlubal Knowledge Base #001512 — LTB of I-section beam per EN 1993-1-1](https://www.dlubal.com/en/support-and-learning/support/knowledge-base/001512)
- [SkyCiv — Guide to Eurocode 3 Steel Design](https://skyciv.com/docs/tech-notes/other/guide-to-eurocode-3-steel-design/)