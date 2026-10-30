# performance/firstorder —— 一阶冻结核：实时顶点移位 与 稠密 IF 装配规则

本文件夹是**一阶（朴素）冻结核方法**的基准套件，兼作下游 TEMPO/GTEMPO
**MPS-IF 装配的稠密参照**。模型统一取 Toulouse（玻色 independent-boson /
费米 SIAM），所有 contour 的核条目一律直接调用 QuAPI 冻结核函数
（`bosonic_Δt/Δτ/Δm`、`fermionic_Δt/Δτ/Δm`，cell 积分闭式）。
全 K 侧二阶改善（梯形 φ、半步相位、接触对角、跳变代入）在
`../secondorder`，此处不重复。

对照基准均为严格 Trotter：实时/虚时/混合各脚本都同时输出严格一阶
（Lie）与严格二阶（Strang）cell 链的误差，用于确认"一阶冻结核确实一阶、
改善不改阶"。

---

## 1 实时顶点移位

### 1.1 为什么一阶：配对位置问题

冻结核条目本身是**精确的**：`Δt`/`Δτ` 的核闭式（胞元双积分）逐位等于严格
cell 演化 IF 的 O(V²) 系数。一阶的来源不在核，而在**装配时顶点冻结的
位置**——旧（face-value）装配把顶点对放在 cell 边界，相对严格方案滞后
偏一格，每步留下 O(V²) 残差，整体 therefore 一阶。顶点移位在**不换任何
核条目**的前提下把配对位置修正过来。

### 1.2 移位规则（方程）

两套 naive / shifted 方程（代码：`bosons/realtime.jl` 的
`frozen_retarded_A`、`vertex_shifted`；`fermions/realtime.jl` 的
`frozen_first_order`、`vertex_shifted`）。

**推迟核组合**（(1+n)−n / (1−n)+n 热抵消后的单粒子推迟核）：

- 玻色（`bosonic_Δt`）：`g[d] = −(P[d+1,1] + X[d+1,1])`，P = Δ⁺⁺，X = Δ⁺₋
- 费米（`fermionic_Δt`）：`Δf(ℓ) = X[ℓ+1,1] − P[ℓ+1,1]`
  （等价于 `fermionic_Δm` 的 −(P+X)；Δm 的 (+,−) 块 = −Δt 的 (+,−) 块）

**naive（face-value，不移位）**：

- 费米递推：`A_{m+1} = u·(A_m − Σ_{d≥1} Δf(d)·A_{m−d})`
- 玻色矩阵：`M = A0⁻¹ + KR`，`A0⁻¹ = I − u·S⁻`（开链自由飞行），
  `KR[j,k] = −(P[j,k]+X[j,k])`（k ≤ j 因果下三角），`A = M⁻¹·e₁`

**vertex-shifted（节点 lag d ↔ 条目 lag d−1，三角并入次对角）**：
记 `g₀ = Δ(0)`（等时/接触条目，naive 中弃用），`u = e^{−iϵ_d h}`：

```
M[m, m]      = 1
M[m, m−1]   −= u·(1 − g₀)          （最近邻步；m = 1 行改入源：b[1] += u(1−g₀)）
M[m, m−d]   += u·g[d−1]            （d = 2..m；m−d ≥ 1，否则 b[m] −= u·g[d−1]）
A(t_m) = (M⁻¹·b)_m
```

要点：

1. **核条目一个不换**（同一批 `Δt` 条目）；变的只有两处——
   (i) 配对改为节点 lag d 配条目 lag d−1（顶点从 cell 边界内收半格到节点）；
   (ii) 等时接触条目 g₀ 并入最近邻步，次对角成为 u(1−g₀)。
2. 矩阵形式下因果自动成立：越出初条件的 lag 落进 b（m=1 行与 m−d=0 行），
   这就是"三角形"。
3. 效果：单粒子 **V² 项严格**（V²-exact），误差常数大幅下降、阶数不变：

| 实时 T=2，N=512 | naive err | shifted err | 比值 |
|---|---|---|---|
| 玻色 | 1.345e-3 | 1.214e-4 | **~11×** |
| 费米 | 1.331e-3 | 3.490e-4 | **~3.8×** |

（数量级参考：同点严格二阶 Trotter 终点误差 ~1.2e-6，故移位后仍清晰一阶。）

### 1.3 适用范围：只有实时

顶点移位是**实时开链特有的自包含改善**。虚时与混合轮廓不走它：

- 虚时的二阶改善属 K 侧要素打包（梯形、半步、接触、跳变代入），见
  `../secondorder/*/imagtime.jl`；
- 混合轮廓把移位单独抽出与其他朴素部分组合**不自洽**（实测误差反而增大，
  见 `bosons/mixedtime.jl` 头注）——一阶→二阶只有"打包改善"一条路径。

因此其他两个 contour 的一阶文件只与严格 Trotter 1/2 阶对比，
不含移位列。

---

## 2 稠密 IF 装配规则（给定冻结核函数）

### 2.0 通用机制（四步管线）

稠密装配的输入是**某通道的完整线性系统 M·x = b**：M 把测度、动能跃迁、
核权重、跳变闭合全部吸收为一个矩阵，b 是初条件/源向量。然后：

1. **腿基**：每节点两条腿 `(a_j, ā_j) = (leg 2j−1, leg 2j)`，
   `GrassmannContourIndex(j, branch, 1, c)`，branch 取 `:τ`（虚时/混合 τ 支）
   或 `:+`（实时开链）。
2. **系数矩阵**（转置核装配，对应 IF 收缩方向 ā·C·a）：
   `C[2k, 2m−1] = −M[m,k]`。
3. **折叠**：`gif = GrassmannIF(GrassmannHybridizationMatrix(basis, C))`，
   `data` 含 2^{2n} 个系数（n = 该通道节点数；`MAXTENSORDIM = 2³⁰`
   ⟹ 单张量 n ≤ 15）。
4. **提取**（Pfaffian 余子式）：
   `(M⁻¹)[i,j] = sgn · data[全1掩码 ⊻ ((1<<(2i−1)) | (1<<(2j−2))) + 1] / data[end]`
   （掩码去掉腿 āᵢ 与 a_j）；符号 `sgn` 逐项数值标定（`target·top/data`
   必须落在整数 ±1/0 上，自校验）；**结构零**须与 M⁻¹ 的零逐项一致
   （因果三角系统的 j > i 块 minor 恒为零）。可观测量即 `x = M⁻¹·b`。

单粒子说明：Toulouse 单粒子振幅代数对玻色/费米同构（skew 双线性 +
Pfaffian 余子式），玻色统计完全包含在核条目里（推迟组合的热抵消），
腿基折叠机制与费米完全一致。下游 MPS-IF 按第 2 步的 C 矩阵逐块装配，
即可逐位复现通道解。

存储代价：单通道张量 2^{2n}；混合轮廓**按通道分段折叠**（τ 环、A、C
各一个张量），n 取各通道自己的节点数，远小于整条 L 轮廓的
2^{2Nτ+4(Nt+1)}。

### 2.1 实时（最佳装配 = 顶点移位矩阵）

按 §1.2 的 vertex-shifted 方案构造 M（g₀ 并入次对角 u(1−g₀)、节点 lag d
配条目 g[d−1]、越界项入 b），腿基 branch=:+ 开链、无回绕，
`A = M⁻¹·b`，b 只有初条件行非零。折叠与提取按 §2.0 原样执行。

> 现有 `bosons/realtime_dense.jl`、`fermions/realtime_dense.jl` 折叠的是
> face-value 矩阵（一致性检验对象是 naive 通道解）；换用最佳装配只需把
> §2.0 第 1–2 步中的 M 换成移位形式，其余规则不变。

### 2.2 虚时（周期 / 反周期热环）

**玻色**（`bosons/imagtime.jl` ↔ `imagtime_dense.jl`，核 `bosonic_Δτ`，
实数条目；无 fτ 半步——与费米约定不同，如实记录）：

```
M = (I − e^{−hϵ_d}·A_cyc) − D ，  A_cyc[mod(j,N)+1, j] = 1   （周期回绕在矩阵内，无显式源）
D[i,j] = ηⱼₖ[i−j+1] (i>j) ； ηₖⱼ[j−i+1] (i<j) ； eqt = ηⱼₖ[1]+ηₖⱼ[1] (对角)
G(τ_m) = (M⁻¹)[m+1, 1] ， m = 1..N−1 ； n_d = +Gτ[Nτ]
```

等时对角必须是**两序和** eqt：`ηⱼₖ[1]·(d†d) + ηₖⱼ[1]·(dd†)` 中不对易项
[d,d†] = 1 是 c-number，整体乘在 IF 上、在任何归一化可观测量中严格抵消
（详见 `../common/common.jl` 注释）；丢单侧是实错误。

**费米**（`fermions/imagtime.jl` ↔ `imagtime_dense.jl`，核 `fermionic_Δτ`）：

```
M = I − uτ·S⁻ − fτ·W ，  uτ = e^{−ϵ_d δτ} ， fτ = e^{−ϵ_d δτ/2}
W 同上 D 结构（朴素矩形：节点值、全权重）
反周期回绕：M[1,N] += uτ ；源 b = −uτ·e₁
Gτ = M⁻¹·b ； n_d = −Gτ[Nτ]   （KMS：n_d = −G^τ(β⁻)）
```

### 2.3 混合时间（τ 热环为输入 + 实时通道 + ζ 跨块）

框架：τ 支热环是**输入**（不被改善）→ 提取 n_d（费米另提取 γ）→
实时 A/C 通道；ζ 跨块（虚轴↔实轴耦合）进入 C 通道源。玻色与费米的
差别在 ζ 的语义。

**玻色**（`bosons/mixedtime.jl` ↔ `mixedtime_dense.jl`；核 `bosonic_Δτ` + `bosonic_Δm`）：

1. τ 环（梯形热环）：`fτ·w` 相邻节点各半权重；跳变常数（k=1 列）
   `M[m,Nτ] -= fτ·w/2`、`M[m,1] -= fτ·w/2`、`b[m] += fτ·w/2`。
   解 `Gτ = M⁻¹b`，`n_d = +Gτ[Nτ]`。
2. A 通道：`bosonic_Δm` η 块，`K[m,k] = −(P[m−k+1,1] + Xr[m−k+1,1])`
   （lag 直读、因果下三角含接触对角），`M_A = I − u·S⁻ + K`，
   `A = M_A⁻¹·e₁`。
3. C 通道：**同一个折叠张量**（同一 M_A），源
   `Svec[m+1] = Σ_j Xc[m+1, j]·Gτ[j]`（ζ 块直接收缩——玻色 ζ 语义正确），
   `C = M_A⁻¹·Svec`。
4. `G>(t) = −i·((1 + n_d)·A + C)`。

**费米**（`fermions/mixedtime.jl` ↔ `mixedtime_dense.jl`；核 `fermionic_Δm`）：

1. τ 环（朴素矩形、无 fτ、反周期）：`M = I − uτ·S⁻ − W`（W 为 ξ 块），
   `b = −uτ·e₁`；`n_d = −Gτ[Nτ]`。
2. γ 提取（朴素梯形，无消峰——一阶病态来源）：
   `Gext = [G₀; Gτ]`，`G₀ = −1 − Gτ[Nτ]`（greater 延拓），
   `γrow[k] = V_k · Σ_j δτ·(w₀Gext[j] + w₁Gext[j+1])/2 / (1 + e^{−ε_k β})`，
   `w_j = e^{−ε_k(β−τ_j)}`。
3. A 通道（半步 f、接触对角隐式）：
   `Δ(j,k) = −(Xr[j,k] + P[j,k])`（满指标），`d₀ = Δ(1,1)`，
   ```
   MA[m+1, m+1] = 1 + f·d₀/2 ；  MA[m+1, m] = −u + f·d₀/2
   MA[m+1, j]  += f·Δ(m+1, j+1)   （j ≤ m−1）
   A = MA⁻¹·e₁
   ```
   **列错位陷阱**：递推源 `s = Σ Δ(m+1,k+1)·A[k]` 的核列是 k+1、振幅列
   是 k——矩阵系数进第 j 列（不是 j+1 列），写错会漏成 O(1) 偏差。
4. C 通道：`MC = I + h·𝒦`（严格下三角），
   `𝒦(D) = Σ_k V_k² e^{−iϵ_d D}(e^{iω_k D}−1)/(iω_k)`，`ω_k = ϵ_d − ε_k`；
   源 `S_m = Σ_k −i·√α_k·J_k(t_m)·(−γrow[k])`，
   `J_k(t) = e^{−iϵ_d t}(e^{iω_k t}−1)/(iω_k)`；`C = MC⁻¹·S`。
5. `G>(t) = −i·((1 − n_d)·A + C)`。

**ζ 语义差别**：费米 ζ 跨块虽 rank-1，但直接朴素收缩存在
lesser/greater 语义错配（给出错误物理极限），故走上述 γ×𝒦 因式分解路径
（与 `../secondorder/fermions/mixedtime.jl` 的结论一致）；玻色 ζ 非退化、
直接收缩正确，无需替换。

### 2.4 已验证的一致性（本次全部实测）

稠密 IF 提取 vs 通道线性系统（机器精度即通过）：

| 文件 | 张量 | max 偏差 |
|---|---|---|
| `fermions/mixedtime_dense.jl` | τ 环 / A / C 三段 | 4e-17 / 2e-16 / 3e-17 |
| `bosons/imagtime_dense.jl` | 周期热环 ×1 | 3e-16 |
| `bosons/realtime_dense.jl` | 实时开链 ×1 | 1e-15 |
| `bosons/mixedtime_dense.jl` | τ 环 + 实时（A、C 共用） | 2e-16 / 1e-15 |
| `fermions/imagtime_dense.jl` | 反周期环 ×1 | 1e-16（此前已验证） |
| `fermions/realtime_dense.jl` | Keldysh 全腿基 | 独立一阶离散（见该文件结论） |

一阶通道解自身的误差（对照严格解，证实均一阶）：

| 通道 | naive err（末点 N） | 备注 |
|---|---|---|
| 玻色实时 | 1.35e-3 @512 | 移位后 1.21e-4，仍一阶 |
| 费米实时 | 1.33e-3 @512 | 移位后 3.49e-4，仍一阶 |
| 玻色虚时 | 2.42e-3 @256 | |
| 费米虚时 | 4.62e-4 @256 | |
| 玻色混合 | 2.60e-3 @256 | τ 热环为正确环（输入） |
| 费米混合 | 1.65e-3 @256 | γ 无消峰为一阶病态来源 |

---

## 4 BCS 浴（`performance/firstorder/bcs/`）

BCS（超导）浴版的实时快速一阶与稠密装配，模型/约定与 `../secondorder/bcs`
完全一致（复用其 `common.jl`：`bcs_toulouse`、BdG 严格参照）。Nambu 基
Φ = (d↑, c↑…; d̄↓†, c̄↓†…)，观测量 𝒜(t) = e^{−i𝓗t} 的杂质 2×2 块；推迟核
𝒩 = [cn cc; nn nc]（正常块 X−P、反常块 P−X，来自 `bcs_Δt` 四块、lag 直读），
行相位 (e^{−iϵ_d h}, e^{+iϵ_d h})。

- §1 的顶点移位规则原样适用于 2×2 矩阵条目：接触 𝒢₀ = 𝒦(0) 并入最近邻步
  P(I−𝒢₀)，节点 lag d 配条目 𝒢(d−1)，越界 lag 入 B。naive/移位均严格一阶
  （err 4.45e-3 / 3.14e-3 @N=512）；移位改善 ~1.4×（比标量费米情形小——
  Nambu 矩阵下接触项误差占比低）。
- 稠密装配（`realtime_dense.jl`）：块矩阵 M（2N 抽象节点 = N 节点 × 2 Nambu
  分量）展开到 4N 腿（每节点 band 1,2 × (a,ā)，与 `test/quapi/bcsimag.jl`
  同序）；**配对约定** C[a腿(s), ā腿(t)] = −M[s,t]（行指标↔a 腿、列指标↔ā 腿），
  提取 (M⁻¹)[s,t] 的掩码去 {a_t, ā_s}（与 §2.0 同一族约定，符号逐项标定）。
  naive/移位两套 M 均机器精度（≤3.5e-16，N ≤ 6）。
- **虚时**（`imagtime.jl` + `imagtime_dense.jl`）：Nambu Matsubara 环的三个关键结构
  （约定经 `exact_Gτ_ring` 裁判的数值扫描锁定，Δ=0 时 𝒢₂≡0 精确成立）：
  1. **每节点 4 个分量**（杂质 2×2 Nambu 块 Ĝ_{ij}）——行 1 切片不封闭：
     [Σ⊛Ĝ]_{ij} = Σ_k′ Σ̂_{i,k′}Ĝ_{k′,j} 会把行 2 分量卷进行 1 方程；
  2. **反常槽位带 −1**：𝚺̂ = [cn, −cc; −nn, nc]（IF 约定把负号吸进 Δ，
     Dyson 重排时反常自能多出 −1）；
  3. **跳变源只在 Nambu 对角分量**（b_{(i,i)} = −uτ_i，反常分量无等点反对易
     子）；行 2（空穴行）自由相位 e^{+ϵ_d δτ}。
  naive 环严格 1.00 阶、全节点采样（2.6e-2 @N=8 → 1.6e-3 @N=128）；严格 Trotter 参照用
  处置矩阵 P̃ = (I+Tᴺ)⁻¹（恒等于 M₀ = (1+e^{−β𝓗})⁻¹）：Lie 1.02 / Strang
  2.00。稠密装配 8N 腿（4 Nambu 分量 band × (a,ā)），机器精度 ≤1.1e-16。
- **混合**（`mixedtime.jl` + `mixedtime_dense.jl`）：与 secondorder 同骨架的
  naive 版——单次 `bcs_Δm` 的 ξ 块 τ 环（朴素矩形，数据 ≡ bcs_Δτ 环，逐位
  验证）→ M11/F0（跳变反解）→ naive 梯形 γ（无消峰）→ face-value march
  （无半步/无接触）→ 逐 cell 矩形 Duhamel 源（bath 相位冻结在 cell 左端点、
  −i 因子同 secondorder）。**两个必需约定**（数值扫描锁定）：
  1. 反常 march 核必须用**原始 uv/vu 加权 Δm 块**构造
     k₁₂ = Ŵ_uv − conj(Ŵ_vu)（组合块 cc/nn 的 face-value 配对不收敛——
     Δm 与 Δt 的反常推迟核是不同对象：实 Δ 下 Δt 给 sin 型虚条目、Δm 给实条目）；
  2. naive 公式差异不敏感处已对拍：γ 与 secondorder `gamma_nambu` 之差
     O(δτ) 收敛（2.9e-5@N=128）。
  naive 全链严格 1.00 阶（6.5e-2@16 → 3.7e-3@256）；稠密装配三段
  （τ 环 8Nτ 腿 + A/C 通道 4(Nt+1) 腿）机器精度 ≤2.5e-16。
- 物理路径值的 BCS 装配（含反常 ā-ā / a-a 配对的全配分函数展开）见包测试
  `test/quapi/bcsimag.jl`（GrassmannADT 乘法）。

---

## 3 文件索引

| 文件 | 内容 | 输出 |
|---|---|---|
| `bosons/realtime.jl` | naive vs 顶点移位 vs 严格 1/2 | `results/realtime.csv` |
| `bosons/imagtime.jl` | 一阶环 vs 严格 1/2 | `results/imagtime.csv` |
| `bosons/mixedtime.jl` | 朴素 E1 装配 vs 严格 1/2 | `results/mixedtime.csv` |
| `bosons/realtime_dense.jl` | 实时稠密 IF vs 线性系统 | stdout 一致性 |
| `bosons/imagtime_dense.jl` | 虚时稠密 IF vs 环方程 | stdout 一致性 |
| `bosons/mixedtime_dense.jl` | 混合稠密 IF vs 通道系统 | stdout 一致性 |
| `fermions/*` | 同名同构（费米版） | 同上 |
| `bcs/realtime.jl` | BCS Nambu 一阶（naive/移位）vs 严格 BdG + Trotter 1/2 | `bcs/results/realtime.csv` |
| `bcs/realtime_dense.jl` | BCS Nambu 稠密 IF vs 块线性系统 | stdout 一致性 |
| `bcs/imagtime.jl` | BCS Nambu Matsubara 环 vs 严格 BdG 环 + Trotter 1/2 | `bcs/results/imagtime.csv` |
| `bcs/imagtime_dense.jl` | BCS 虚时稠密 IF vs 块线性系统 | stdout 一致性 |
| `bcs/mixedtime.jl` | BCS 混合 naive（ξ 环+γ+march+Duhamel 源）vs 严格 1/2 | `bcs/results/mixedtime.csv` |
| `bcs/mixedtime_dense.jl` | BCS 混合稠密 IF vs 通道递推 | stdout 一致性 |
| `*/plot_results.py` | CSV → PNG（误差–步长双对数） | `results/*.png` |

运行方式统一为
`JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. performance/firstorder/<bosons|fermions|bcs>/<file>.jl`。
稠密文件的规模受 2^{2n} 存储/折叠代价限制，取小 N 严格对齐即可
（混合按通道分段折叠，N ≤ 10 无压力；整链稠密 N ≤ 12）。
