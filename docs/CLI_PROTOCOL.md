# FrameCore CLI 線協議(`frame_cli` — J1 文字橋)

> S6 J1 對接橋的權威線協議。`frame_cli.exe` 是 `frame::solve` 等分析的 stdin/stdout 驅動器,供
> Grasshopper / 任何外部客戶端 shell-out 驅動。單位 **N, mm, MPa**。原始碼:
> `Plugins/FrameSolver/Standalone/frame_cli.cpp`;端到端測試:`Tools/cli_roundtrip.py`(gate 第五腿)。

## 呼叫方式
把整個模型描述(每行一個 token 指令)寫到 **stdin**,以 `END` 結尾;從 **stdout** 讀回結果行。
provenance(`# frame_cli | build <sha> | …`)走 **stderr**,不污染 stdout。

```
echo "<model lines>\nEND" | frame_cli.exe
```

### Daemon 多請求模式(S6 J1.5)
驅動器對 **model BLOCK** 迴圈:每個 block 以 `END` 結束 → 求解 → 輸出 → 一行 `EOR`(end-of-response,**flush**)→
重置接下一個 block。一個保持管線開啟的客戶端可把多個模型**串流經同一行程**(讀到 `EOR` 即一塊完成)。
單發(model + END + EOF)是單 block 特例,DISP/MF 輸出逐位元不變(`EOR` 是未知 token,舊解析器忽略)。
**已驗(gate)**:同行程多 block == 各自獨立 cli 逐位元相等。`PreparedSystem`/`ReSolveSession` 跨請求重用是後續最佳化。
裸 `END`(空 block)= 純握手:只回 `VERSION` + `EOR`。

## 輸入指令(token,以空白分隔;`MAT`/`SMAT`/`SEC` 必須在引用它們的元素之前)
| 指令 | 參數 | 說明 |
|---|---|---|
| `MAT` | `E G rho [capComp capTens capShear]` | 梁材料(nu→0)。同池。**optional cap**(S6 J1b):allowable 給 SIZEOPT/D-C;省略=`make(300,300,180)`;只給 1 值→tens=comp、shear=0.6·comp |
| `SMAT` | `E nu G [capComp capTens capShear]` | 殼材料(帶 nu)。同池;cap 同上 |
| `SEC` | `A Iy Iz J cy cz Asy Asz` | 截面;`Wy=Iy/cy`、`Wz=Iz/cz` |
| `NODE` | `id x y z  fUx fUy fUz fRx fRy fRz  [pUx..pRz]` | 6 個固定旗標 0/1;optional 6 個 prescribed 位移(預設 0) |
| `MEMBER` | `id i j matIdx secIdx  refx refy refz  [active [tonly]]` | active 預設 1;**`tonly` 預設 0 = tension-only 旗標(供 `TONLY`)** |
| `SHELL` | `id n0 n1 n2 n3 matIdx t  [active]` | MITC4 平板殼 facet |
| `NLOAD` | `node Fx Fy Fz Mx My Mz` | 節點載重 |
| `UDL` | `member wx wy wz` | 桿均布(local) |
| `SPRESS` | `shellId p` | 殼橫向壓力 |
| `HINGE` | `member dof Mp` | 塑鉸(dof 4/5/10/11,signed Mp;節點側力矩由 caller 的 NLOAD 給) |
| `OPT` | `enableReleases useTimoshenko pivotTol` | 求解選項 |
| `EIGEN` | `nModes` | 附加模態頻率輸出 |
| `PDELTA` | `path` | 二階:0=凍結重用、1=K_T 參考;缺/<0=線性 |
| **`TONLY`** | `[maxIter [allowReact]]` | **S6**:tension-only 主動集 eliminator(讀 `MEMBER … tonly` 桿) |
| **`SIZEOPT`** | `Amin maxIter dcTol` | **S6**:全應力 FSD 尺寸優化(全 active 桿) |
| **`DYNC`** | `dt maxTime [rid…]` | **S6**:連續動力倒塌**摘要**;行尾 ids = initialRemovals |
| `END` | — | 結束輸入,開始求解 |

分析模式互斥:`TONLY`/`SIZEOPT`/`DYNC` 出現即設模式(後者覆蓋);皆無 → `PDELTA`(若設)→ 線性 `solve`。

## 輸出行(stdout)
| 行 | 欄位 | 說明 |
|---|---|---|
| **`VERSION`** | `sha` | **S6**:永遠是 stdout 第一行(build 握手) |
| `SINGULAR` | `0\|1` | 機構/不穩定偵測 |
| `DISP` | `nodeId ux uy uz rx ry rz` | 每節點(節點順序) |
| `RF` | `nodeId Fx Fy Fz Mx My Mz` | 反力 |
| `MF` | `id  Ni Vyi Vzi Ti Myi Mzi  Nj Vyj Vzj Tj Myj Mzj` | 每桿端內力(local;N 壓正) |
| `SF` | `id Mxx Myy Mxy Qx Qy Nxx Nyy Nxy` | 每殼合力 |
| `FREQ` | `n omega1 omega2 …` | `EIGEN` 時的模態頻率(rad/s) |
| `PDSTATUS` | `conv div iters` | `PDELTA` 時領先一行 |
| **`TONLY`** | `conv cycled iters` | **S6**;後接 `SLACK id…`(鬆弛桿),再標準 `SINGULAR/DISP/RF/MF/SF`(finalState) |
| **`SLACK`** | `id…` | **S6**:tension-only 收斂時停用的桿 |
| **`SIZEOPT`** | `conv iters singular` | **S6**;後接每桿 `AREA id A DC` + `WEIGHTVOL ΣA·L` |
| **`AREA`** | `memberId A DC` | **S6**:優化後面積 + 收斂 D/C |
| **`WEIGHTVOL`** | `value` | **S6**:材料體積 ΣA·L(mm³) |
| **`DYNC`** | `outcome nEvents nFrames Tend` | **S6**:outcome 0=Stable 1=Collapsed 2=MaxSteps 3=Invalid |
| **`DEVENT`** | `t mode nRemoved nDetached` | **S6**:每個拓撲事件摘要 |
| **`DFRAME`** | `t maxAbsU` | **S6 J1b**:逐儲存幀的時間 + 峰值\|位移\|(回放時間軸;完整 u/v 串流延後) |
| **`EOR`** | — | **S6 J1.5**:end-of-response 哨兵 + flush(每 block 一個;daemon 客戶端讀到即知該塊完成) |

數值精度 `%.12g`。未知 token 由客戶端忽略(向後相容:`Tools/opensees_compare.py` 只讀 `SINGULAR/DISP/MF`,新行不影響)。

## 範例 — 懸臂尖端載重
```
MAT 210000 80769 7850
SEC 10000 8333333.333 8333333.333 14060000 50 50 8333.333 8333.333
NODE 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0
NODE 1 2000 0 0 0 0 0 0 0 0 0 0 0 0 0 0
MEMBER 0 0 1 0 0 0 0 1
NLOAD 1 0 0 1000 0 0 0
END
```
→ `VERSION <sha>` / `SINGULAR 0` / `DISP 1 … 1.52381 …`(`δ = PL³/3EI`)。

## 範例 — tension-only X 斜撐 portal
```
…(MAT/SEC stocky/SEC brace/NODE×4/MEMBER columns+beam)…
MEMBER 3 0 3 0 1 0 1 0 1 1      ← active=1 tonly=1
MEMBER 4 1 2 0 1 0 1 0 1 1
NLOAD 2 50000 0 0 0 0 0
TONLY
END
```
→ `TONLY 1 0 2` / `SLACK 4` / 標準 state。

## 已知限制 / 未來(誠實標)
- **材料 allowable cap 目前硬編 `Capacity::make(300,300,180)`**(MAT 未帶 cap);`SIZEOPT`/D-C 用此 cap →
  尺寸結果隨之縮放。需 caller 設 cap 時,協議待加 `MAT … cap_comp cap_tens cap_shear` optional token(J1b)。
- `DYNC` 目前只出**摘要**;完整時間歷程幀(`DynCollapseFrame.u/v`)+ 碎塊交接資料(給 UE/Chaos 回放)
  屬 U 層協議,**延後 J1b**。
- **daemon 常駐**(持有 `PreparedSystem`/`ReSolveSession` 跨請求,省重複 factor)= **J1.5**;
  **C API DLL** = **J2**。皆未實作(WS_R2 §10:J1 行程開銷僅 6.7ms,~20k DOF 端到端 2.1s,J1 先夠)。
- GH `.gha` 元件 + Yak 發佈需 **Rhino 8 .NET SDK**(見 `Plugins/FrameSolver/Grasshopper/README.md`),不入引擎 gate。
