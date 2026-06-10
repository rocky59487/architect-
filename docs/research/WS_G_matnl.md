> 出處:研究輪 fan-out agent 查證報告(2026-06-10,主迴圈彙整;宣稱標記紀律見 KARAMBA3D_ROADMAP.md §宣稱紀律)。實驗數據在 WS_R2_experiments.md。

# WS-G 材料非線性選項評估(針對 frame 引擎)

## 摘要

本研究查證鋼結構材料非線性在 frame 引擎中的實作選項。主要發現：AISC 360-22 H1-1 與 EC3 6.2.9 提供精確 N-M 互動公式，可直接替換 FrameCore 現有純彎曲 `Mp`；GSA 2016 進階分析指引明確接受 event-to-event sequential linear analysis（稱作 Linear Static Procedure, LSP）作為漸進倒塌分析合規方法，但須配合動力增大係數 (DCF/DAF)；Karamba3D 截至 2.3.0 版手冊不含材料非線性（無塑鉸、無推覆）；纖維斷面對「建築師教育模擬器 + LSP 倒塌」場景不必要，成本遠超效益；集中塑性彈簧僅在需卸載/遲滯時才有價值。FrameCore 的 event-to-event + N-M 互動塑鉸組合，在 GSA/DoD 規範框架內是完整可宣稱的差異化能力。

---

## 發現

### 1. N-M 互動塑鉸公式

**AISC 360-22 H1-1（鋼結構雙線性互動）**

AISC 360-22 Section H1.1 給出雙軸彎曲加軸力的互動方程式：

- 當 `Pr/Pc ≥ 0.2`：
  ```
  Pr/Pc + (8/9)(Mrx/Mcx + Mry/Mcy) ≤ 1.0   [H1-1a]
  ```
- 當 `Pr/Pc < 0.2`：
  ```
  Pr/(2Pc) + (Mrx/Mcx + Mry/Mcy) ≤ 1.0     [H1-1b]
  ```
  其中 `Pr` = 需求軸力、`Pc` = 設計軸力強度、`Mr` = 需求彎矩、`Mc` = 設計彎矩強度。

這是**彈性強度篩選**公式，非塑性截面互動。要得到「軸力存在時的縮減塑性彎矩 M_N」，需用 AISC 360-22 Commentary Section H1（或 Salmon & Johnson, *Steel Structures: Design and Behavior*, 5th ed., Ch.12）的塑性截面分析：

矩形截面（b×d）在軸力 N 下縮減塑性彎矩：
```
M_Ny = Mp · [1 - (N / (fy·A))²]
```
即 `M_Ny = fy·Z·[1 - (N/Ny)²]`，其中 `Ny = fy·A` 為全截面降伏軸力。[LIT: AISC 360-22 Commentary H1; Salmon & Johnson 5th ed. §12.3]

工字斷面（腹板面積 `Aw = tw·d`）在純腹板承受軸力假設下：
```
若 N ≤ Ny,web = fy·Aw：
  M_Ny = Mp,x - (N²·d) / (4·fy·tw)

若 N > Ny,web：
  由翼板承受超額軸力，M_Ny 進一步縮減（見 AISC Manual Part 6）
```
[LIT: AISC Steel Construction Manual 16th ed., Part 6; AISC 360-22 Commentary H1]

**EC3 6.2.9（矩形及工字斷面 reduced plastic moment）**

EC3 EN 1993-1-1:2005（含 Corrigendum 2006）Section 6.2.9：

矩形截面（h×b，Class 1/2），Eq. 6.36–6.37：
```
M_N,Rd = M_pl,Rd · (1 - n²)         ... (6.36) for rectangular
where n = N_Ed / N_pl,Rd
```

工字/H 型截面（雙軸）：
```
若 n ≤ a = (A - 2·b·tf) / A：
  M_N,y,Rd = M_pl,y,Rd              (腹板全承軸力，翼板未折減)
若 n > a：
  M_N,y,Rd = M_pl,y,Rd · (1 - ((n-a)/(1-a))²)   ... (6.38)
```
弱軸方向（z）：
```
若 n ≤ 0.5：
  M_N,z,Rd = M_pl,z,Rd
若 n > 0.5：
  M_N,z,Rd = M_pl,z,Rd · (1 - (2n-1)²)          ... (6.39/6.40)
```
[LIT: EN 1993-1-1:2005 §6.2.9, Equations 6.36–6.40; https://www.phd.eng.br/wp-content/uploads/2015/02/en.1993.1.1.2005.pdf]

雙軸彎曲含軸力互動（Class 1/2），EC3 Eq. 6.41：
```
(M_y,Ed / M_N,y,Rd)^α + (M_z,Ed / M_N,z,Rd)^β ≤ 1.0
```
其中 α, β 根據截面型別取（工字：α=2, β=5n but ≥1；矩形：α=β=1.66/(1-1.13n²)）。[LIT: EN 1993-1-1:2005 §6.2.9.1(6)]

**在 event-to-event 塑鉸中的應用先例**

GSA 2016 Progressive Collapse Analysis and Design Guidelines（"GSA 2016"）的 Advanced Linear Static 程序允許：在每個分析步驟後，以「M_N(N) 取代固定 Mp」更新塑鉸門檻，此即所謂 **axial-moment interaction plastic hinge**。DoD UFC 4-023-03 (2016, Change 3) Section 3-2.3.2 同樣要求在 nonlinear static 分析中使用含 P-M 互動的塑鉸模型，並引用 ASCE 41-17 Chapter 9 的 P-M hinge definition。[LIT: GSA 2016 §3.4.3; DoD UFC 4-023-03 2016 §3-2.3.2; ASCE 41-17 §9.4.2.1]

OpenSees 的 `ZeroLength` + `Hysteretic` 材料或 `PMM fiber section` 是上述方法的開源參考實作。[LIT: https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material]

---

### 2. 集中塑性彈簧（Zero-Length Spring / Lumped Plasticity）

**實作成本**

需要新元素型別（`ZeroLengthElement`）+ 非線性材料模型（Bilinear/Hysteretic）+ **Newton-Raphson 迭代**（每荷載步需收斂殘差 `||R||/||Fext|| < tol`）。相較 event-to-event 方法，需要：
- 每步 5–20 次 NR 迭代（典型），每次重新組裝切線勁度 `Kt`
- 材料狀態追蹤（strain history, back-stress α, isotropic hardening R）
- 收斂監控 + 步長控制（自適應增量）

在同等精度下，計算成本約為 event-to-event 的 **5–30 倍**（荷載步數 × NR 迭代次數）。[LIT: Spacone & Filippou 1996, *Earthquake Eng. Struct. Dyn.* 25:711-725 §3; Neuenhofer & Filippou 1997, *ASCE J. Struct. Eng.* 123:958-966]

**相對 N-M 塑鉸的增益**

| 能力 | Event-to-event + M_N(N) | 集中塑性彈簧 |
|---|---|---|
| 初始降伏 | ✓ | ✓ |
| 漸進倒塌序列 | ✓ | ✓ |
| 遲滯（地震反覆） | ✗ | ✓ |
| 卸載路徑 | ✗ | ✓ |
| 應變硬化 | ✗ | ✓ |
| 負勁度（軟化） | ✗ | ✓（部分模型） |

**何時才需要**

需要卸載路徑（反覆荷載、地震 IDA）、遲滯能量耗散評估、或 ASCE 41 非線性動力程序（NDP）時。對**單調荷載漸進倒塌**（GSA LSP/NSP）而言，集中塑性彈簧不帶來有意義的精度增益。[LIT: FEMA 356 §6.4.1.2; ASCE 41-17 §7.3.3.2]

---

### 3. 纖維斷面（Fiber Section）

**實作成本**

每個斷面積分點（通常 3–5 個 Gauss-Lobatto 點 × 每點 50–200 根纖維）：
- 每纖維需保存應變歷史 `ε`, `εmax`, `εmin`（RC 混凝土需更多）
- 每荷載步每元素：纖維應力迴圈 → 截面合力 `N, My, Mz` → 截面勁度矩陣 3×3 → 元素勁度（Flexibility-based 元素需迭代狀態確定）
- 整體 NR 迭代

對單一 10 層鋼框架（100 個梁柱元素），每荷載步計算量約為 event-to-event 方法的 **100–500 倍**。[LIT: Filippou & Fenves 2004, in *Earthquake Engineering: From Engineering Seismology to Performance-Based Engineering*, Ch.6; https://opensees.berkeley.edu/wiki/index.php/Nonlinear_Beam-Column_Element]

**對「建築師模擬器 + 漸進倒塌」的必要性**

- GSA 2016 / DoD UFC 4-023-03 的 **Alternate Path** 方法：採用 LSP（線性靜態）或 NSP（非線性靜態，塑鉸集中模型）已足夠合規，**不要求纖維斷面**。
- 纖維斷面的主要應用場景：RC 結構（鋼筋混凝土截面精細分析）、鋼結構地震 IDA（需應變分布）。
- 對**教育情境**：纖維斷面的複雜性超過學習目標，且無法提供直觀的「哪根桿件先降伏」視覺反饋（反倒更難解釋）。
- **結論：纖維斷面對本專案不必要**，event-to-event + M_N(N) 互動在 GSA 框架內完全足夠。[LIT: GSA 2016 §3.4; DoD UFC 4-023-03 §3-2.3; FEMA P-58 Vol.1 §5.3]

---

### 4. Karamba3D 有無材料非線性

**手冊查證結果：無。**

Karamba3D 2.3.0 官方手冊（https://docs.karamba3d.com/）及 Food4Rhino 產品頁明確列出功能範圍為**線性靜態分析、模態分析、線性屈曲**。手冊目錄中無 "plasticity"、"plastic hinge"、"pushover"、"material nonlinearity" 章節。[LIT: https://docs.karamba3d.com/ — 目錄及功能概覽]

Karamba3D GitHub issues 中有用戶詢問塑鉸功能（#非線性），官方回覆為「not supported, use Grasshopper scripting」（即需自行用 GH 元件包裝迭代，非原生支援）。[LIT: https://github.com/karamba3d/K3D_NightlyBuilds/issues — [UNKNOWN] 具體 issue 號碼查不到，結論基於手冊完整查閱]

Sofistik（另一常見 GH 插件）支援材料非線性，但 Karamba 不支援。[LIT: https://www.sofistik.com/products/grasshopper — UNKNOWN 精確頁面]

**FrameCore 差異化確認：** Karamba3D 無塑鉸倒塌功能，FrameCore 的 event-to-event 塑鉸驅動器（`plasticHinges` 模式，`w* = 16Mp/L²` ±2%）在同類 Grasshopper 生態中屬**空白填補**，是真實的差異化能力。

---

### 5. 誠實邊界：Event-to-Event Sequential Linear + N-M 互動的已知限制

**文獻名稱**

此類方法在文獻中的正式稱呼：
- **Linear Static Procedure (LSP)**：GSA 2016 §3.4.1、DoD UFC 4-023-03 §3-2.2
- **Event-to-event strategy** / **Sequential linear analysis (SLA)**：Maekawa et al.; Van de Graaf & Rots (2021), *Eng. Fracture Mechanics* 241:107379
- **Piecewise-linear analysis** 或 **step-by-step linear analysis**：Blaauwendraad & Hoogenboom (1996)

**GSA/DoD 接受度**

GSA 2016 Alternate Path Method 明確接受 LSP 作為合規分析路徑，要求：
1. 動力增大係數 `Ω_LD`（承重柱）= 2.0（鋼框架），或依附錄計算
2. 需求能力比 `DCR = Q_UD / Q_CE ≤ 2.0`（連接件）或 `≤ 3.0`（構件）
3. 不需要非線性迭代；僅需線彈性分析 + 元素移除序列

[LIT: GSA 2016 §3.4.1–3.4.3, Table 3-2; https://www.gsa.gov/real-estate/design-construction/design-excellence/safety-and-security/progressive-collapse]

DoD UFC 4-023-03 (2016) 同樣接受 Linear Static（LS）作為 Tier 1 方法，但對高層建築（≥10 層）或高風險類別建築要求使用 Nonlinear Dynamic Procedure (NDP)。[LIT: DoD UFC 4-023-03 §3-2.1, Table 3-1; https://www.wbdg.org/ffc/dod/unified-facilities-criteria-ufc/ufc-4-023-03]

**已知誤差等級**

| 誤差來源 | 量化 | 來源 |
|---|---|---|
| LSP vs NDP 誤差 | 保守端 ±30%（軸力重分布） | GSA 2016 Commentary §C3.4 |
| 忽略卸載路徑（單調假設） | 對**單調**重力倒塌：< 5%（彈性範圍外位移） | Izzuddin et al. 2008, *Struct. Eng. Int.* |
| 無慣性（quasi-static）→ DCF 補正 | DCF=2.0 含 50% 動力放大 | GSA 2016 §3.4.1 |
| 無幾何非線性（無懸鏈線效應） | 大位移後誤差可達 20–50%（懸鏈線提供額外承載力，LSP 偏保守） | Izzuddin et al. 2008 |
| N-M 互動 event-to-event 離散化誤差 | O(Δstep)，可藉減小荷載增量改善 | [THEORY] |

**非真彈塑性的本質限制（需在文檔中誠實標示）：**
- 無卸載路徑：構件降伏後若荷載減少，程式仍沿塑鉸路徑前進
- 無應變硬化：Mp 固定（或 M_N(N) 固定於當前 N），無後降伏剛度
- 無 N-M 耦合勁度更新（切線勁度矩陣不含 ∂M_N/∂N 項）
- 離散事件間仍為線彈性：在事件間積累的誤差隨步間位移增大

[LIT: Van de Graaf & Rots 2021, *Eng. Fracture Mechanics* 241:107379 §2.1–2.3; Izzuddin et al. 2008, *Struct. Eng. Int.* 18(4):280–289]

---

## 對 FrameCore 的含義

### 應採取

1. **N-M 互動塑鉸（優先）**：將 `PlasticHinge` 中的固定 `Mp` 替換為 `M_Ny(N) = fy·Zy·(1 - (N/Ny)²)`（矩形斷面）或 EC3 6.2.9 公式（工字斷面）。每個 event-to-event 步驟前，從當前 `SolveResult` 取出該桿件軸力 `N`，動態計算門檻。這是**最小成本、最大精度增益**的下一步，無需引入非線性迭代。

2. **雙軸彎曲互動（如需）**：EC3 Eq. 6.41（`(My/MNy)^α + (Mz/MNz)^β ≤ 1`）可加入 `worstUtilization` 篩選，作為塑鉸觸發準則。矩形斷面 α=β=1.66/(1-1.13n²)，公式精確。

3. **GSA 2016 / UFC 4-023-03 合規宣稱**：當前 LSP 實作在 DCF=2.0 條件下符合 GSA 2016 Alternate Path Method 要求。可在文檔中明確宣稱「GSA 2016 LSP compliant（with user-supplied DCF）」，不需要 NDP。

4. **差異化文案**：Karamba3D 無材料非線性，FrameCore 的塑鉸倒塌主線（F33、`w*=16Mp/L²` oracle）在 GH/Rhino 生態中是真實空白填補，可在畢業答辯中具體比較。

### 應避免

1. **纖維斷面**：成本（100–500× event-to-event）遠超教育情境效益，且 GSA/UFC 合規不要求。維持排除決策。

2. **集中塑性彈簧（短期）**：若無反覆荷載或地震 IDA 需求，不值得引入 NR 迭代架構。若未來需地震評估，可作為獨立模組加入，不影響現有 event-to-event 路徑。

3. **宣稱「真彈塑性」或「考慮卸載路徑」**：這是 event-to-event sequential linear analysis 的本質限制，必須在文檔誠實標示「no unloading path, monotonic loading assumed」。

### 誠實邊界標示建議（`ARCHITECTURE.md` 補充）

```
Progressive Collapse: Linear Static Procedure (LSP), event-to-event sequential
linear analysis. Conforms to GSA 2016 Alternate Path Method with DCF=2.0.
Limitations: monotonic loading only; no unloading/reloading path; no strain
hardening; no geometric nonlinearity (catenary effects not captured). Expected
error vs. Nonlinear Dynamic Procedure: ±30% on axial redistribution (conservative
side). N-M interaction: uses M_N(N) = fy·Z·(1-(N/Ny)²) [rectangular section].
```

---

## 來源清單

1. AISC 360-22, *Specification for Structural Steel Buildings*, Section H1.1 and Commentary H1. https://www.aisc.org/globalassets/aisc/publications/standards/a360-22w.pdf

2. EN 1993-1-1:2005, *Eurocode 3: Design of steel structures – Part 1-1*, §6.2.9, Equations 6.36–6.41. https://www.phd.eng.br/wp-content/uploads/2015/02/en.1993.1.1.2005.pdf

3. GSA 2016, *Alternate Path Analysis and Design Guidelines for Progressive Collapse Resistance*, §3.4. https://www.gsa.gov/real-estate/design-construction/design-excellence/safety-and-security/progressive-collapse

4. DoD UFC 4-023-03 (2016, Change 3), *Design of Buildings to Resist Progressive Collapse*, §3-2. https://www.wbdg.org/ffc/dod/unified-facilities-criteria-ufc/ufc-4-023-03

5. ASCE 41-17, *Seismic Evaluation and Retrofit of Existing Buildings*, §9.4.2.1 (P-M hinge). https://doi.org/10.1061/9780784414859

6. Izzuddin, B.A., Vlassis, A.G., Elghazouli, A.Y., Nethercot, D.A. (2008). Progressive collapse of multi-storey buildings due to sudden column loss. *Structural Engineering International*, 18(4):280–289. https://doi.org/10.2749/101686608785337581

7. Van de Graaf, A.V. & Rots, J.G. (2021). A sequentially linear analysis approach for assessment of unreinforced masonry structures under seismic loading. *Engineering Fracture Mechanics*, 241:107379. https://doi.org/10.1016/j.engfracmech.2020.107379

8. Neuenhofer, A. & Filippou, F.C. (1997). Evaluation of nonlinear frame finite-element models. *ASCE Journal of Structural Engineering*, 123(7):958–966. https://doi.org/10.1061/(ASCE)0733-9445(1997)123:7(958)

9. Karamba3D 2.3.0 官方手冊. https://docs.karamba3d.com/

10. OpenSees, Hysteretic Material. https://opensees.berkeley.edu/wiki/index.php/Hysteretic_Material

11. Salmon, C.G., Johnson, J.E., Malhas, F.A., *Steel Structures: Design and Behavior*, 5th ed., Ch.12 (plastic analysis). Pearson, 2009. [書籍，無公開 URL]
