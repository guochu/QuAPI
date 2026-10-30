# 影响泛函的 OSEE

本页说明如何从离散杂化函数出发计算**影响泛函（IF）系数态的
operator-space entanglement entropy（OSEE）**：一键入口、模式/割的约定、
费米与玻色两个分支各自的适用范围与物理含义，以及与稠密 IF 张量的
交叉验证。[Hybridization Functions](hybridization.md) 第 2 节中的 IF
指数在这里被当作一个**算符**整体向量化，其纠缠熵用于快速预估链算法
（QuAPI/GTEMPO）所需的键维。

纯数学层面（协方差快速算法的完整推导）见 GaussianStates 的文档
「理论与算法」页；数值依据是 `GTEMPO/osee/osee_report.md`（费米三方一致：
高斯协方差 ≡ 稠密 2^n SVD ≡ GrassmannMPS bond 熵，1e-14 ~ 1e-10）与
`TEMPO/osee/osee_report.md`（玻色分支的 d→∞ 极限与有限 Fock 截断的关系）。

## 1. 定义

把 [Hybridization Functions](hybridization.md) 第 2 节的 IF 指数
（QuAPI 约定，负号已吸收进 Δ）

```math
I = \exp\Big[\sum_{b_1,b_2}\sum_{i\in b_1}\sum_{j\in b_2}
\bar\psi_i^{b_1}\,\Delta^{b_1b_2}_{ij}\,\psi_j^{b_2}\Big]
```

视为算符并向量化为 doubled (ket, bra) 空间的（未归一化）纯态系数向量
`|I⟩⟩`。对 contour 节点的子集 `A`（其余节点为补集），OSEE 就是 `|I⟩⟩`
在 `A|Ā` 二分下的纠缠熵。快速算法不构造 `2^{2n}` 系数向量：IF 指数是
二次型，`|I⟩⟩` 是**高斯态**，熵由协方差矩阵的子矩阵本征值给出，
每割 ``O(n^3)``（`n` 为 contour 节点数）；态在构造时一次算出实协方差
（`covariance(s)` 读取），跨割复用。

费米与玻色的关键差别：

| | 费米 | 玻色 |
|---|---|---|
| 系数态 | Grassmann 有限和，**严格**高斯态 | 无穷维高斯算符；快速公式给 **d→∞ 高斯极限** |
| 截断 | 无（Grassmann 局域维恒 2） | 有限 Fock 截断 d 的 TEMPO PT 熵是另一端点（见 §4） |
| contour | 虚时 / 实时 / 混合（含复 Δ） | 仅实 Δ（虚时 Matsubara） |

## 2. 一键入口（GaussianStates extension）

QuAPI 对 GaussianStates 是 weak dependency：加载 GaussianStates 会触发
`QuAPIGaussianStatesExt`，给 GaussianStates v7 的两个态构造器
`FermionicGaussianState` / `BosonicGaussianState` 增加「correlation
function / 杂交矩阵 → IF 系数态」方法（态只存实协方差：费米 Majorana
`Γ = i⟨γγ⟩`、玻色 `Σ = ½⟨{R,R}⟩`）。

```julia
using QuAPI, GaussianStates

corr = fermionic_Δt(semicircular(), β = 2.0, t = 1.0, N = 4, μ = 0.2)
s = FermionicGaussianState(corr)        # IF 系数态（nmodes(s) == 2n）
n = nmodes(s) ÷ 2                       # contour 节点数
k = 2                                   # 前 k 个节点 vs 其余
entropy(s, vcat(1:k, n .+ (1:k)))       # von Neumann；Rényi 用 …; α=2
```

### 模式与割的约定

- 系数态有 `2n` 个模（doubled：模 `1..n` = ket contour 节点、
  `n+1..2n` = bra，`nmodes(s) == 2n`）；**site p 的割 = 模式对
  `[p, n+p]`**，即 `entropy(s, vcat(sites, n .+ sites))`。
- 节点序 = `contour_ordering`（`QuAPI.Dense`，见 §5 示例）：

| contour | 节点数 | `contour_ordering` |
|---|---|---|
| Δτ（虚时） | `N` | `(:τ,1), …, (:τ,N)` |
| Δt（实时） | `N+1` | `(:+,1),(:-,1),(:+,2),(:-,2), …`（时间升序，`+` 先 `−` 后） |
| Δm（混合） | `Nτ + Nt + 1` | 虚时升序接实时交替 |

「前 k 个节点」即 `contour_ordering` 的前缀；实时 contour 的一个节点
同时含 `+`、`−` 两支腿，与 GTEMPO 链的每时刻块 `[a⁺,ā⁺,a⁻,ā⁻]` 对应。
**前缀割与 GTEMPO GrassmannMPS 的 bond 割是同一约定**，两个熵逐位一致
（机器精度），因此 OSEE 可直接用作键维预估：``\chi \gtrsim e^S`` 保留
主要谱权重。

### 各 bath 类型的支持范围

- **普通费米 bath**：三种 contour 全部支持（含复 Δ；实时/混合的复 Δ
  走协方差的复数一般式）。单 flavor 时从稠密 IF 的**顶点移位后**腿矩阵
  提取 ``\Delta``（矩阵级入口 `FermionicGaussianState(h::GrassmannHybridizationMatrix)`），
  与 `GrassmannIF` 张量的 OSEE 逐位一致（§5）。
- **玻色 bath**：只支持**实 Δ**（虚时 Matsubara contour）；实时/混合
  contour 的复 Δ 会使系数态协方差复化，调用报 `ArgumentError`。
- **BCS bath**：一键入口未接（反常 ``\bar a\bar a``/``aa`` 块需要一般
  ``4N`` 腿复反对称配对，不是「每节点一对 a/ā」的布局）；组装方式见
  `GTEMPO/osee/osee_report.md` §5.6。

## 3. 两个分支的配对结构

### 费米（任意 contour）

系数态 ``|\Psi\rangle = \exp(\sum_{pq}\Delta_{pq}\,\bar a^\dagger_p a^\dagger_q)|0\rangle``
（配对块 ``\Delta[\bar a_p, a_q]``，Grassmann 幂零性保证指数是有限和），
对任意 contour（含复 Δ）恒良定义。扩展内由 SVD 闭式直接组装 Majorana
协方差（``\varepsilon_k = s_k^2/(1+s_k^2)`` 为配对模式占据，``s_k`` 为
``\Delta`` 的奇异值），不经 `purify`——系数态是 `:exp(Δ):` 型
normal-ordered Gaussian 的矢量化，仅当 ``\Delta`` 可逆时才能写成
`e^{单体温子 h}` 的 Choi 态（奇异 ``\Delta`` 无 ``h`` 可言）。

### 玻色（实 Δ）

```math
M = \exp(-\Delta), \qquad
\text{配对振幅 } M[\bar a_p, a_q]\ \text{（实对称块结构，谱} \subset (0,1)\text{）}.
```

两点值得知道：

- **为什么是 ``\exp(-\Delta)``**：玻色系数向量 ``\exp(\sum\Delta\, b^\dagger b)``
  向量化后的直接配对振幅是 ``\exp(+\Delta)``，其奇异值 > 1 使 ``I - W W^{\mathsf T}``
  不定、协方差公式数值失败。``\exp(-\Delta)`` 与之只差每个模的
  粒子–空穴反射，时间割的 Schmidt 谱相同（几何分布 ``p_n \propto t^{2n}``
  在 ``t \leftrightarrow 1/t`` 下不变），熵不变。
- **结果的含义**：这是 **d→∞ 高斯极限**的熵，不是有限 Fock 截断 d 下的
  TEMPO PT/ADT bond 熵——两者的关系见下一节。

```julia
corr_b = bosonic_Δτ(spectrum(ε -> exp(-ε), lb = 0.5, ub = 8.0), β = 2.0, N = 4)
sb = BosonicGaussianState(corr_b)                  # 纯态（辛本征值全 ½）
entropy(sb, [1, 5])                                # site 1 的割 = 模式 [1, 5]
```

另外 `BosonicGaussianState(::BosonicHybridizationMatrix)` 只接受
`split=false`；`split=true` 存的是等时算符次序拆分
`BosonicSplitDelta(M1, M2)`，不是配对振幅，传入直接报错。

## 4. 玻色 OSEE 与 Fock 截断 d 的关系

固定 δτ 时 d→∞ 极限**有限**，但收敛尺度 ``d \sim 1/\lambda_{\min}(\Delta)``
（``\lambda_{\min}`` 为 Δ 的最小正本征值；极限值量级 ``\sim \ln(1/\lambda_{\min})``）；
δτ→0 时 ``\lambda_{\min}\to 0``，极限**对数发散**——玻色 OSEE 不是
连续极限下的普适物理量。有限 d 的熵与 d→∞ 极限是**同一条曲线的两端**：
`test/osee.jl` 的「d 扫描夹逼」测试验证有限 d 的 BosonicPT 熵单调上升并
在高斯值附近穿过。

实践含义：

- 想要**与有限 d 的 TEMPO 计算严格对应**的熵：固定 d 做算符 SVD
  （`TEMPO/osee/osee_report.md` §3、§8），不要用高斯公式。
- 想要**快速扫参、键维预估**：高斯 OSEE 即可，O(模式数³)、无 Fock 截断。
  弱耦合下还可用非对角块 SVD 一阶式
  ``S_l \approx 2\sum_k s_k^2(1-\ln s_k^2)``（``s_k =`` `svdvals(Δ[1:l, l+1:N])`），
  误差约几个百分点。

## 5. 与稠密 IF 张量的交叉验证

费米侧 OSEE 与稠密 IF（`QuAPI.Dense`）的前缀割 plain SVD 熵逐位一致
（IM/RT/MX 三 contour，`atol = 1e-10`，见 `test/osee.jl`）：

```julia
using QuAPI, QuAPI.Dense, GaussianStates

corr = fermionic_Δm(semicircular(), β = 2.0, Nτ = 2, t = 0.5, Nt = 1, μ = 0.2)
T = GrassmannIF(GrassmannHybridizationMatrix(corr))     # 稠密 IF（2^2n 系数）
chain = contour_ordering(corr)
n = length(chain)
s = FermionicGaussianState(corr)

k = 1                                                   # 链前缀 k 个节点
modes = [findfirst(==(b), chain) for
         b in unique([(bi.branch, bi.j) for bi in T.basis[1:2k]])]
S_fast = entropy(s, vcat(modes, n .+ modes))
S_dense = begin                                         # 系数矩阵的 plain SVD 熵
    sv = svdvals(reshape(T.data, 2^(2k), 2^(2(n - k))))
    p = abs2.(sv); p ./= sum(p)
    -sum(x * log(x) for x in p if x > 1e-15)
end
S_fast ≈ S_dense                                        # true
```

这也反过来固定了「IF 矩阵元 ↔ 指数项」的指标约定：系数态构造取的 Δ
与 `GrassmannIF` 顶点移位后的装配完全同源（`contourindex.jl` 的
访问器），两条路径不会各说各话。

## 6. 注意事项

1. **`entropy` 直接作用在系数态上**（v7 起态只存实协方差，`entropy(s, modes)`
   对任意 `GaussianState` 有定义；纯态整态熵为 0，`ispure(s)` 可查验）。
   系数态不是密度矩阵——它本身是 doubled 空间的**纯态**（算符的矢量化），
   `entropy` 的 modes 即 doubled 物理模（site p ↔ `[p, n+p]`）。
2. **玻色须实 Δ**：实时/混合 contour 的玻色 Δ 复化后快速算法失效
   （`ArgumentError`），这不是可绕过的实现细节而是协方差的实数性边界。
3. **费米无对角项问题**：虚时费米 Δ 严格反对称（对角 ~1e-15）；玻色的
   ``\Delta_{ii}`` 只影响有限 d 精确计算的熵值（~0.4% 量级），高斯快速公式
   保留全部 Δ，无需处理。
4. **复杂度**：每割 ``O(n^3)``；系数态构造时一次算出协方差（`covariance(s)`
   读取），适合全链逐 bond 扫描（N=128 逐割秒级，见 GTEMPO 报告 §4）。
5. **Rényi**：`entropy(s, modes; α)` 支持 ``\alpha \neq 1``（两能级谱公式）。
