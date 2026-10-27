# Dense IF Assembly

本页讲三种稠密影响泛函——费米 [`QuAPI.Dense.GrassmannADT`](@ref)、玻色对角
[`QuAPI.Dense.BosonicADT`](@ref) 与玻色非对角 [`QuAPI.Dense.BosonicPT`](@ref)——如何从
[Hybridization Functions](@ref "The Hybridization Functions") 章的离散杂化函数
**组装**出来：腿/节点的排列约定、「杂化矩阵元 ⟶ 张量元」的映射公式，以及
组装引入的误差标度。记号：谱强度记作 ``\Gamma`` 或 ``\alpha``（即 ``V^{2}``，
``V`` 为系统–库耦合项的系数；谱宽记 ``D``）；误差标度以耦合系数 ``V`` 计——
``O(V^{2})`` 即「线性于谱强度」阶。

## 1. 三种稠密表示

装配都是两步：先**杂化矩阵**（把 `CorrelationMatrix` 的 Toeplitz 生成向量铺
到 contour 腿/节点上，`GrassmannHybridizationMatrix` /
`BosonicHybridizationMatrix`），再**折叠成稠密张量**（`GrassmannIF` /
`BosonicIF`）：

| 类型 | 耦合 | `basis` | `data` | 构造 |
|---|---|---|---|---|
| [`QuAPI.Dense.GrassmannADT`](@ref) | 费米（含 BCS） | `L` 条 contour 腿 | `2^L` 个 Grassmann 单项式系数 | `GrassmannIF` ∘ `GrassmannHybridizationMatrix` |
| [`QuAPI.Dense.BosonicADT`](@ref) | 玻色对角（`hyb` 对角化） | `n` 个节点 | `d^n` 个路径权重 | `BosonicIF(corr, λ)`，`λ` = 本征值 |
| [`QuAPI.Dense.BosonicPT`](@ref) | 玻色非对角（任意 `hyb`） | `n` 个节点 | `d^{2n}` 个矩阵元 | `BosonicIF(corr, hyb)`（内部 `split=true`） |

```julia
using QuAPI, QuAPI.Dense

h  = GrassmannHybridizationMatrix(corr)   # 腿基 basis + 配对矩阵 Δ
I  = GrassmannIF(h)                       # GrassmannADT
Z  = Zvalue(I)                            # Berezin 积分（全腿单项式系数）

adt = BosonicIF(bc, [0.5, -0.5])          # BosonicADT（对角耦合）
pt  = BosonicIF(bc, X)                    # BosonicPT（非对角耦合）
```

## 2. 费米：GrassmannADT

### 腿基与配对矩阵

腿的排列（`basis`）：节点按 contour 序（虚时/实时时间升序，实时每时刻
`+` 先 `−` 后，混合 τ 先 ± 后），每节点 band 内 `(a, ā)` 相邻。IF 是腿上的
Grassmann 高斯指数

```math
I[\xi_1,\dots,\xi_L] \;=\; \exp\Big(\sum_{p,q=1}^{L} C_{pq}\,\xi_p\,\xi_q\Big),
```

配对矩阵 ``C`` 的非零结构由 `_delta_terms` 决定：普通 bath 只有
``(\bar a, a)`` 块 ``C[\bar a_p, a_q] = \Delta_{pq}``（``\Delta`` 的行列序即
`contour_ordering`）；BCS 另有第二 flavor 的正常块与反常块
``(\bar a_\uparrow,\bar a_\downarrow)``、``(a_\downarrow,a_\uparrow)``。
由于 ``\xi^2=0``，指数在总阶数 ``\le L`` 后自动终止——**折叠（`_fold_grassmann`）
是有限次的单项式乘法，无任何近似**。

### 张量元、求值与边缘化

`data` 按腿 1 最快展开单项式系数，单腿占据态用位掩码索引：

```math
I \;=\; \sum_{m=0}^{2^L-1} \texttt{data}[m{+}1]\;\xi^m,
\qquad \xi^m = \prod_{k:\,m_k=1}\xi_k,
```

单项式系数带**反对易符号**（按 `basis` 升序重排计入；腿重复返回零）。
配对展开为

```math
[\xi_{i_1}\cdots\xi_{i_r}]\,I \;=\; \sum_{\text{配对方案 } M}
\operatorname{sgn}(M)\,\prod_{(\bar p,\,q)\in M} C[\bar p, q],
\qquad i_1<\cdots<i_r,
```

全腿单项式系数即 Berezin 积分（配分函数）：

```math
Z \;=\; \int\!\mathrm{d}\xi_L\cdots\mathrm{d}\xi_1\; I
\;=\; [\xi_1\cdots\xi_L]\,I \;=\; \texttt{data}[2^L],
```

即 `Zvalue`。[`QuAPI.Dense.margin`](@ref) 对未列出的腿占据态求和（Grassmann
边缘化）；两个 IF 的乘法 `A * B` 把共享腿按 Grassmann 融合、外腿并入
`basis`。

### 顶点移位后的 ``C``：τ 环

冻结核的生成向量按「顶点在 cell 边界」的滞后约定存储；稠密组装把顶点
移到节点——生成值重组为 **skew-circulant**（单 flavor，``N`` 节点）：

```math
c_0 = \eta_{kj}[2], \qquad
c_1 = \eta_{jk}[1] + \eta_{kj}[1] \;(\text{等时和}), \qquad
c_k = \eta_{jk}[k] \;(k\ge 2),
```

```math
C[\bar a_i, a_j] \;=\;
\begin{cases}
c_{\,i-j}, & i > j \\[2pt]
c_0, & i = j \\[2pt]
-\,c_{\,N+i-j}, & i < j
\end{cases}
```

（下三角为正向传播核；上三角取反周期边界。旧装配把 ``c_1`` 放在
``(i, i{-}1)`` 位——残差因此停在 ``O(V^{2})``，见第 5 节。）

BCS 把同一移位施加到四个 ``(row, col)`` 块上（正常
``(\bar a_\uparrow,a_\uparrow)``、``(\bar a_\downarrow,a_\downarrow)`` 都读
`uu`；反常 ``(\bar a_\uparrow,\bar a_\downarrow)`` 读 `ud`、
``(a_\downarrow,a_\uparrow)`` 读 `du`）。一般块移位公式（``r_i``/``s_j`` 为
块内腿位，``n`` 为节点数）：

```math
\gamma_1 = C[r_1,s_1], \quad
\gamma_k = C[r_k,s_1]\;(k\ge2), \quad
\gamma_0 = C[r_1,s_2],
```

```math
D[r_i,s_j] \;=\;
\begin{cases}
\gamma_0, & i=j \\
\gamma_{\,i-j}, & i>j \\
-\gamma_{\,n+i-j}, & i<j
\end{cases}
```

### 实时 Keldysh 块

节点 ``j_0 = 0,\dots,N``，腿位
``\mathrm{pos}(j_0,\sigma,\bar a) = 4j_0 + [\sigma=+]\cdot 0 + [\sigma=-]\cdot 2 + [\bar a]\cdot 1``。
记 ``d_{\pm} = C[\mathrm{pos}(0,\pm,\bar a),\,\mathrm{pos}(0,\pm,a)]``（Toeplitz
对角常值）。四个 ``(ā, a)`` 块的移位：

```math
\begin{aligned}
D[\bar a^{+}_{j_1}, a^{-}_{j_2}] &= C[\bar a^{+}_{j_1}, a^{-}_{j_2}],
&& \text{窗口：}\ \bar a^{+}_0\ \text{行与}\ a^{-}_N\ \text{列置零},\\[2pt]
D[\bar a^{-}_{j_1}, a^{+}_{j_2}] &= C[\bar a^{-}_{j_1}, a^{+}_{j_2}],
&& \text{窗口：}\ a^{+}_N\ \text{列与}\ \bar a^{-}_N\ \text{行置零},\\[2pt]
D[\bar a^{+}_{j_1}, a^{+}_{j_2}] &=
\begin{cases} d_{+}, & j_2 = j_1-1 \\ C[\bar a^{+}_{j_1},\, a^{+}_{j_2+1}], & \text{否则} \end{cases}
&& \text{（列移位 } a_j\!\to\! a_{j+1}\text{，动能位取本支对角）},\\[2pt]
D[\bar a^{-}_{j_1}, a^{-}_{j_2}] &=
\begin{cases} d_{-}, & j_2 = j_1+1 \\ C[\bar a^{-}_{j_1},\, a^{-}_{j_2-1}], & \text{否则} \end{cases}
&& \text{（镜像移位 } a_j\!\to\! a_{j-1}\text{）}.
\end{aligned}
```

（交叉块 ``(ā⁺,a⁻)``, ``(ā⁻,a⁺)`` 数值不变，仅把顶点窗口收紧到 cell
``1..N``；``(ā^+,a^+)``/``(ā^-,a^-)`` 的列移位把 a 算符顶点从 cell 边界搬到
节点，动能位 ``e^{-\delta t\,\epsilon_d}`` 由对角元给出。）

### 混合轮廓

``\tau`` 段与 ``\pm`` 段的顶点源映射（行源 ``\rho``、列源 ``\kappa``；无源 ⟹ 该
元为零，符号恒 ``+1``；``N_\tau`` 为虚时节点数、``N_t`` 为实时 cell 数）：

```math
\rho:\quad
\tau_j \to \begin{cases}\tau_{j-1}, & j\ge 2\\ \varnothing, & j=1\end{cases}
\qquad
+_j \to \begin{cases}\tau_{N_\tau}, & j=1\\ +_{j-1}, & j\ge2\end{cases}
\qquad
-_j \to \begin{cases}-_j, & j\le N_t\\ \varnothing, & j=N_t{+}1\end{cases}
```

```math
\kappa:\quad
\tau_k \to \tau_k
\qquad
+_k \to \begin{cases}+_k, & k\le N_t\\ \varnothing, & k=N_t{+}1\end{cases}
\qquad
-_k \to \begin{cases}\varnothing, & k=1\\ -_{k-1}, & k\ge2\end{cases}
```

```math
D[\nu, \mu] \;=\; C[\rho(\nu),\, \kappa(\mu)].
```

（数值标定：混合移位对 `exact_fermionic_Δm` 严格参照做 Richardson 外推
至 ``O(V^{2})`` 逐位一致。）

## 3. 玻色对角：BosonicADT

对角耦合 ``S`` 可对角化（本征值 ``\lambda_1,\dots,\lambda_d``）时，IF 是路径
权重的**对角嵌入**：路径 ``\nu`` 处算符取本征值 ``s_\nu\in\{\lambda_p\}``，

```math
I[s_1,\dots,s_n] \;=\;
\exp\Big(\sum_{\nu,\mu=1}^{n} s_\nu\,\Delta_{\nu\mu}\,s_\mu\Big),
```

`data` 按节点 1 最快存全部 ``d^n`` 个路径权重；配分函数是全部路径之和

```math
Z \;=\; \sum_{s_1,\dots,s_n} I[s_1,\dots,s_n] \;=\; \texttt{sum(data)}
```

（`Zvalue`）。链上相邻 ADT 的乘法收缩共享节点。**这一形式对任意
`split=false` 核是精确的等式**（第 5 节的负对照：对严格 cell 参照残差
``\sim 10^{-15}`` 平线）。

## 4. 玻色非对角：BosonicPT

``hyb`` 不可对角化（或本征基外）时，IF 是作用在系统路径上的**算符**。
代码取 ``\mathrm{op}_1 = S``、``\mathrm{op}_2 = S^{\dagger}``（`op2 = adjoint(hyb)`）。
节点 ``\nu`` 的行腿矩阵随分支取转置 / 伴随 / 原矩阵（张量指标为
``(\mathrm{out},\mathrm{in})``；+ 支元素 ``\mathrm{op}[a,m]``、− 支为
``\mathrm{op}^{\dagger}`` 的元素 ``\overline{\mathrm{op}[b,m]}``、τ 支同 + 支
但腿序相反）：

```math
X_1^{\nu} = \begin{cases}\mathrm{op}_1^{\mathsf T} & (\sigma_\nu=+)\\ \mathrm{op}_1^{\dagger} & (\sigma_\nu=-)\\ \mathrm{op}_1 & (\sigma_\nu=\tau)\end{cases},
\qquad
X_2^{\nu} = \begin{cases}\mathrm{op}_2^{\mathsf T} & (\sigma_\nu=+)\\ \mathrm{op}_2^{\dagger} & (\sigma_\nu=-)\\ \mathrm{op}_2 & (\sigma_\nu=\tau)\end{cases}.
```

``M_1 + M_2 = \Delta`` 为等时算符次序拆分（``\eta_{jk}``/``\eta_{kj}`` 侧，
第 9 节）。``H`` 按节点对 ``(\nu,\mu)`` 的分支关系组装（``j_\nu`` 为 contour
步数，``\nu`` 靠后指 ``j_\nu>j_\mu``）：

```math
\boxed{\;
\begin{aligned}
\text{同支不同节点}\;\; H &\mathrel{+}= M_1[\nu,\mu]\,X_1^{\nu}X_2^{\mu} + M_2[\mu,\nu]\,X_2^{\nu}X_1^{\mu},
&& \sigma_\nu=\sigma_\mu=+;\\
&\mathrel{+}= M_1[\nu,\mu]\,X_2^{\nu}X_1^{\mu} + M_2[\mu,\nu]\,X_1^{\nu}X_2^{\mu},
&& \sigma_\nu=\sigma_\mu=-\ \text{或}\ \tau;\\[4pt]
\text{同节点等时}\;\; H &\mathrel{+}= M_1[\nu\nu]\,X_2^{\nu}X_1^{\nu} + M_2[\nu\nu]\,X_1^{\nu}X_2^{\nu},
&& \sigma_\nu=+;\\
&\mathrel{+}= M_1[\nu\nu]\,X_1^{\nu}X_2^{\nu} + M_2[\nu\nu]\,X_2^{\nu}X_1^{\nu},
&& \sigma_\nu=-\ \text{或}\ \tau;\\[4pt]
\text{跨支}\ (+\nu,-\mu)\;\; H &\mathrel{+}= (M_1{+}M_2)[\nu,\mu]\,X_1^{\nu}X_1^{\mu}
+ (M_1{+}M_2)[\mu,\nu]\,X_2^{\nu}X_2^{\mu},
&& \text{从 } + \text{ 侧处理一次};\\[4pt]
\tau\nu\times\pm\mu\;\; H &\mathrel{+}= (M_1{+}M_2)[\nu,\mu]\,X_1^{\nu}\,Y^{\mu}
+ (M_1{+}M_2)[\mu,\nu]\,X_2^{\nu}\,Y'{}^{\mu},
&& \text{从 } \tau \text{ 侧处理一次},
\end{aligned}
\;}
```

其中混合配对的 ``\pm`` 侧算符标签与跨支相反：

```math
Y^{\mu} = \begin{cases}X_2^{\mu}, & \mu\text{ 为 } + \text{ 节点}\\ X_1^{\mu}, & \mu\text{ 为 } - \text{ 节点}\end{cases},
\qquad
Y'{}^{\mu} \text{ 取另一个}.
```

（物理来源：``b@\tau \leftrightarrow b^{\dagger}@\pm`` 的收缩把跨支的
「求和核配同算符」翻成「核留原槽位、± 侧换配对算符」。混合轮廓中
``\pm`` 节点之间的对仍用纯实时规则。）

稠密张量 = 显式矩阵指数（**不含任何 Prony/指数展开近似**）：

```math
I \;=\; \exp(H)\Big|_{\text{vec}} ,
\qquad \text{指标 } 2p{-}1 = \text{节点 } p \text{ 的 ket 路径},\ 2p = \text{bra 路径},
```

对 ``\mathrm{op}_1=\mathrm{op}_2``（厄米 ``S``）或对角 ``S``，以上退化为
``H = \sum M_1\,\mathrm{op}_1\mathrm{op}_2 + M_2\,\mathrm{op}_2\mathrm{op}_1``，
再退化为第 3 节的 ADT 对角嵌入。

## 5. 误差标度

顶点移位后，三种表示的核内容 = 严格 B-cell IF（同网格、同腿基的
exact 参照：费米为 `exact_fermionic_Δτ/Δt/Δm`，玻色为「精确 cell 演化
算符 + 单一 bath 迹」）的 **``O(V^{2})`` 系数逐位一致**（``V`` 为耦合
系数；以谱强度计即线性阶 ``O(\Gamma)``）。据此分三种情形：

1. **玻色对角 / ``\mathrm{op}_1=\mathrm{op}_2``（厄米）/ ``[\mathrm{op}_1,\mathrm{op}_2]=0``**：
   严格 IF 恰为 Gaussian ⟹ 组装**机器精确**（任意 ``\delta t``；负对照实测
   ``\sim10^{-15}`` 平线）。
2. **玻色 ``\mathrm{op}_1\neq\mathrm{op}_2``（如 ``S=\sigma^{+}``）**：严格 IF 含
   ``O(V^{4})`` 非高斯累积量（Magnus 对易子项），组装与严格 IF 的
   生成元一致到 ``O(V^{2})``，同网格残差按 ``V^{4}`` 标度。
3. **费米（冻结核）**：核闭式只保证 ``O(V^{2})`` 系数；同网格残差 =
   胞内两配对过程的 ``V^{4}`` 截断。整体 ``Z/\rho`` 另带一个
   ``O(V^{2})`` 归一化常数（「跳去又跳回」泄漏 = exact 参照的
   `prefactor`），迹归一化观测量中自动消去。

fixed-``N`` 标度实测（`performance/assembly_scaling.jl`，残差为逐 cell
局部量 ⟹ 同网格差值 = ``N\times`` 单胞截断，方案误差完全抵消）：

| 检验 | 残差比 |
|---|---|
| ``V`` 减半（固定网格；τ / RT / PT） | ``\times\,15\text{–}17``（``V^{4}``） |
| ``\delta t`` 加倍（fixed ``N``；费米 RT，``\beta`` 固定无漂移） | ``\times\,15\text{–}16``（首倍增） |
| ``\delta\tau`` 加倍（fixed ``N``；小 ``\omega`` 抑热漂移） | ``\times\,\sim14`` |
| 玻色 ``\sigma_x``（op1 = op2）负对照 | 机器精确平线 |

即**首阶同网格残差 ``\approx N\cdot V^{4}\delta t^{4}``/cell**（以谱强度
计 ``\propto\Gamma^{2}\delta t^{4}``）——塌缩参数为 ``V^{2}\delta t^{2}``
（单次胞内配对积分
``\int_{\text{cell}}\!\!\int_{\text{cell}} \Delta \propto V^{2}\delta t^{2}``）；
大 ``\delta t`` 端（``V^{2}\delta t\gtrsim0.1``）出现高阶修正使比值回落。
注意 fixed-``N`` 扫描测的是装配/核截断的 ``\delta t`` 标度；对连续极限的
全局收敛是方案误差 ``O(V^{2}\delta t^{3})/\text{cell}``（exact 核），
需固定 ``T`` 变 ``N`` 另测，不进同网格残差。

验证脚本：`performance/fermionicG4.jl`（费米 ``V`` 平台）、
`performance/assembly_scaling.jl`（fixed-``N`` ``\delta t/\delta\tau`` 扫描）、
`performance/nondiagcoupling.jl`（玻色 PT 严格 IF 对照）。
