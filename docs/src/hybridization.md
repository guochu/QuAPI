# The Hybridization Functions

本页给出包中三个离散杂化函数 Δτ、Δt、Δm 的**数学定义**：矩阵结构、
每个元素的谱积分表达式，以及被积时间核的封闭形式。记号与
`src/{imagtime,realtime,mixedtime}` 中实现的 kernel 一一对应。

## 1. 记号与总体结构

### 输入

- 谱函数 ``A(\varepsilon)\ge 0``（`AbstractBoundedFunction`，如
  `DiscreteSpectrum`、`semicircular`、`Leggett`）；
- 逆温度 ``\beta``，化学势 ``\mu``。

### 统计因子

```math
n(\varepsilon) \equiv \frac{1}{e^{\beta(\varepsilon-\mu)}+1}, \qquad
\bar n(\varepsilon) \equiv 1-n(\varepsilon) = \frac{1}{1+e^{-\beta(\varepsilon-\mu)}}
```
（费米占据 / 空穴），以及

```math
b(\varepsilon) \equiv \frac{1}{e^{\beta(\varepsilon-\mu)}-1}, \qquad
\bar b(\varepsilon) \equiv 1+b(\varepsilon) = \frac{1}{1-e^{-\beta(\varepsilon-\mu)}}
```
（玻色占据 / 受激发射增强因子）。

### 影响泛函与符号约定

整个包把影响泛函指数中的负号吸收进杂化矩阵：

```math
I = \exp\Big[\sum_{ij}\bar\psi_i\,\Delta_{ij}\,\psi_j\Big]
```

（费米，``\bar\psi/\psi`` 为离散环境场；玻色自旋路径则 ``I=\exp[\sum_{ij}s_i\Delta_{ij}s_j]``，
``s_i=\pm 1/2``），而**不是**常见的 ``\exp[-\sum\bar\psi\Delta\psi]``。

### Toeplitz 结构

杂化矩阵是（分块）Toeplitz 的，只存两条生成向量 ``\eta_{jk}``（上三角，
``i>j``）与 ``\eta_{kj}``（下三角）：

```math
\Delta_{ij} = \begin{cases} \eta_{jk}[\,i-j+1\,], & i>j \\[2pt]
\eta_{jk}[1]+\eta_{kj}[1], & i=j \\[2pt] \eta_{kj}[\,j-i+1\,], & i<j \end{cases}
```

对角元为**两种算符排序之和**——这只对对角耦合成立（见
[`CorrelationMatrix`](@ref) 的 docstring）。

### 谱积分的一般形式

每个生成元素都是一个谱积分

```math
\eta[\cdot] = \int \mathrm{d}\varepsilon \; A(\varepsilon)\times(\text{统计因子})\times(\text{时间核}) ,
```

其中时间核是**自由传播相位在两个离散化 cell 上的精确双重积分的解析
原语**（下文各节给出封闭形式）。数值上核在 ``|\varepsilon|<10^{-6}`` 时切换到
Taylor 极限（如 ``\delta^2/2``、``\delta^2`` 等）以避免 ``0/0``；Bose 因子用
`expm1` 消除小 ``\beta\varepsilon`` 的灾难性抵消。

### 两种 μ 约定（重要）

| 函数 | 谱的处理 | 相位变量 |
|---|---|---|
| `fermionic_Δτ` / `bosonic_Δτ` | 内部 `spectrumshift`：``\bar A(\xi)=A(\xi+\mu)`` | ``\xi=\varepsilon-\mu``（相位含平移后能量） |
| `fermionic_Δt` / `bosonic_Δt` / `*_Δm` | 谱不平移 | ``\varepsilon``（绝对能量；``\mu`` 只进统计因子） |

``\mu=0`` 时两种约定一致。下面 Δτ 一节的公式写在平移变量 ``\xi`` 中
（统计因子相应取 ``n(\xi)=1/(e^{\beta\xi}+1)`` 等）；Δt、Δm 一节用未平移的
``\varepsilon``。

---

## 2. 虚时杂化函数 Δτ

网格 ``\tau_k=(k-1)\delta\tau``，``k=1,\dots,N``，``\delta\tau=\beta/N``，
cell ``I_k=[\tau_k,\tau_{k+1})``。返回 ``N\times N`` 的
`ImagCorrelationFunction`（实型 Toeplitz `CorrelationMatrix`）。

### 时间核

对分离量 ``m-1`` 个 cell（``m\ge 2``）：

```math
K_m(\xi) = \frac{2\,e^{-(m-1)\delta\tau\,\xi}\big(\cosh(\delta\tau\,\xi)-1\big)}{\xi^2}
\;=\; \int_{I_j}\!\!\int_{I_k} e^{-\xi(\tau-\tau')}\,\mathrm{d}\tau\,\mathrm{d}\tau' ,
\quad j-k=m-1
```

即 ``e^{-\xi(\tau-\tau')}`` 在两个 cell 上的精确双重积分。同一 cell 内的两
种（算符）排序各自给出三角（半 cell）积分

```math
T_-(\xi) = \frac{e^{-\delta\tau\,\xi}-(1-\delta\tau\,\xi)}{\xi^2}
= \int_0^{\delta\tau}\!\mathrm{d}\tau\!\int_0^{\tau}\! e^{-\xi(\tau-\tau')}\,\mathrm{d}\tau' ,
\qquad
T_+(\xi) = \frac{e^{+\delta\tau\,\xi}-(1+\delta\tau\,\xi)}{\xi^2} = T_-(-\xi),
```

两者之和还原整 cell 双重积分：``T_-+T_+ = K_1``。小 ``\xi`` 极限：
``K_m\to e^{-(m-1)\delta\tau\xi}\delta\tau^2``，``T_\pm\to\delta\tau^2/2``。

### 生成向量

记 ``\bar A(\xi)=A(\xi+\mu)``（`spectrumshift` 的结果）。**费米**：

```math
\begin{aligned}
\eta_{jk}[1] &= \int \mathrm{d}\xi \;\bar A(\xi)\,\bar n(\xi)\,T_-(\xi), &
\eta_{jk}[m] &= \int \mathrm{d}\xi \;\bar A(\xi)\,\bar n(\xi)\,K_m(\xi), \\[2pt]
\eta_{kj}[1] &= -\int \mathrm{d}\xi \;\bar A(\xi)\,n(\xi)\,T_+(\xi), &
\eta_{kj}[m] &= -\int \mathrm{d}\xi \;\bar A(\xi)\,n(\xi)\,K_m(-\xi),
\end{aligned}
\qquad m=2,\dots,N
```

**玻色**：

```math
\begin{aligned}
\eta_{jk}[1] &= \int \mathrm{d}\xi \;\bar A(\xi)\,\bar b(\xi)\,T_-(\xi), &
\eta_{jk}[m] &= \int \mathrm{d}\xi \;\bar A(\xi)\,\bar b(\xi)\,K_m(\xi), \\[2pt]
\eta_{kj}[1] &= +\int \mathrm{d}\xi \;\bar A(\xi)\,b(\xi)\,T_+(\xi), &
\eta_{kj}[m] &= +\int \mathrm{d}\xi \;\bar A(\xi)\,b(\xi)\,K_m(-\xi).
\end{aligned}
```

两种统计的差异集中在下三角：费米的 ``\eta_{kj}`` 整体多一个负号（反
交换统计）；对角元之和为

```math
\Delta_{ii} = \int \mathrm{d}\xi\, \bar A(\xi)\,T_-(\xi) \;\pm\; \int \mathrm{d}\xi\, \bar A(\xi)\,\sigma(\xi)\,K_1(\xi),
```

其中 ``\sigma=n``（费米，取负号）或 ``b``（玻色，取正号）——这正是费米
反交换 / 玻色交换统计的体现。

**注意**：``\mu=0`` 时玻色因子在 ``\xi=0`` 有 ``1/\xi`` 奇异性，谱支撑触及
零频（平带、Ohmic 边缘、零频 δ 峰）会使求积抛 `DomainError`。

---

## 3. 实时杂化函数 Δt

闭合 Keldysh 回路 ``t\in[0,t]``，节点 ``t_k=(k-1)\delta t``，
``k=1,\dots,N+1``（即节点 0..N），``\delta t=t/N``。返回
`RealCorrelationFunction`，含四个 ``(N+1)\times(N+1)`` Toeplitz 块
``G_{\sigma\sigma'}``（``\sigma,\sigma'=\pm``）。谱**不**平移，``\mu`` 只进
统计因子。

### 时间核

```math
R_m(\varepsilon) = \frac{2\,e^{-i(m-1)\delta t\,\varepsilon}\big(1-\cos(\delta t\,\varepsilon)\big)}{\varepsilon^2}
\;=\; \int_{I_j}\!\!\int_{I_k} e^{-i\varepsilon(t-t')}\,\mathrm{d}t\,\mathrm{d}t',
\quad j-k=m-1
```

（前向相位 ``e^{-i\varepsilon(t-t')}`` 的 cell×cell 双重积分；``R_m^*`` 对应
后向相位 ``e^{+i\varepsilon(t-t')}``），以及同一 cell 内两种排序的三角积分

```math
S_-(\varepsilon) = \frac{(1-i\,\delta t\,\varepsilon)-e^{-i\delta t\,\varepsilon}}{\varepsilon^2},
\qquad S_+(\varepsilon) = S_-(\varepsilon)^* .
```

小 ``\varepsilon`` 极限：``R_m\to e^{-i(m-1)\delta t\varepsilon}\delta t^2``，
``S_\pm\to\delta t^2/2``。

### 四个 Keldysh 块

**费米**（每格给出 ``(η_{jk}[1],\;η_{jk}[m],\;η_{kj}[1],\;η_{kj}[m])``，
``m=2,\dots,N+1``，积分号 ``\int\mathrm{d}\varepsilon\,A(\varepsilon)(\cdot)`` 略）：

| 块 | ``\eta_{jk}[1]`` | ``\eta_{jk}[m]`` | ``\eta_{kj}[1]`` | ``\eta_{kj}[m]`` |
|---|---|---|---|---|
| ``G_{++}`` | ``-\bar n\,S_-`` | ``-\bar n\,R_m`` | ``+n\,S_+`` | ``+n\,R_m^*`` |
| ``G_{+-}`` | ``+n\,S_-`` | ``+n\,R_m`` | ``+n\,S_+`` | ``+n\,R_m^*`` |
| ``G_{-+}`` | ``-\bar n\,S_-`` | ``-\bar n\,R_m`` | ``-\bar n\,S_+`` | ``-\bar n\,R_m^*`` |
| ``G_{--}`` | ``+n\,S_-`` | ``+n\,R_m`` | ``-\bar n\,S_+`` | ``-\bar n\,R_m^*`` |

**玻色**（同格式）：

| 块 | ``\eta_{jk}[1]`` | ``\eta_{jk}[m]`` | ``\eta_{kj}[1]`` | ``\eta_{kj}[m]`` |
|---|---|---|---|---|
| ``G_{++}`` | ``-\bar b\,S_-`` | ``-\bar b\,R_m`` | ``-b\,S_+`` | ``-\bar b\,R_m^*`` |
| ``G_{+-}`` | ``+b\,S_-`` | ``+b\,R_m`` | ``+b\,S_+`` | ``+b\,R_m^*`` |
| ``G_{-+}`` | ``+\bar b\,S_-`` | ``+\bar b\,R_m`` | ``+\bar b\,S_+`` | ``+\bar b\,R_m^*`` |
| ``G_{--}`` | ``-b\,S_+`` | ``-b\,R_m`` | ``-\bar b\,S_+`` | ``-\bar b\,R_m^*`` |

其中 ``G_{+-}``、``G_{-+}`` 满足 ``\eta_{kj}=\eta_{jk}^*``（Hermitian 块）。

**等时存储约定**：``G_{+-}/G_{-+}`` 的对角元按两种排序**分别**存储
（``\eta_{jk}[1]`` 与 ``\eta_{kj}[1]`` 各来自单个三角核 ``S_\mp``），与旧版
`bosonic_Δt_deprecated` 的存储约定一致；只有对角耦合时两部分才相加
（`CorrelationMatrix` 的对角装配）。

---

## 4. 混合时杂化函数 Δm（Kadanoff L 形回路）

回路由三支组成：虚时支 ``\tau\in[0,\beta)`` 取 ``N_\tau`` 个 cell
（``\delta\tau=\beta/N_\tau``），随后闭合 Keldysh 部分 ``t\in[0,t]`` 每支
``N_t+1`` 个节点（``\delta t=t/N_t``）。存储为六块：

| 存储 | 形状 | 进入的分支块 |
|---|---|---|
| ``\eta_{jk},\eta_{kj}`` | ``N_t+1`` | 实时块 ``G^{\sigma\sigma'}``（四块） |
| ``\xi_{jk},\xi_{kj}`` | ``N_\tau`` | 虚时块 ``G^{\tau\tau}`` |
| ``\zeta_{jk}`` | ``N_\tau\times(N_t+1)`` | 混合块 ``G^{\tau\sigma}`` |
| ``\zeta_{kj}`` | ``(N_t+1)\times N_\tau`` | 混合块 ``G^{\sigma\tau}`` |

### 九块的装配（`index` 的定义）

对实时指标 ``i,j``（虚时指标在 ``G^{\tau\cdot},G^{\cdot\tau}`` 中）：

```math
\begin{aligned}
G^{++}_{ij} &= \begin{cases}\eta_{jk}[i-j+1], & i>j\\ \eta_{jk}[1]+\eta_{kj}[1], & i=j\\ \eta_{kj}[j-i+1], & i<j\end{cases} &
G^{+-}_{ij} &= \begin{cases}-\eta_{kj}[i-j+1]^*, & i>j\\ -2\eta_{kj}[1], & i=j\\ -\eta_{kj}[j-i+1], & i<j\end{cases} \\[4pt]
G^{-+}_{ij} &= \begin{cases}-\eta_{jk}[i-j+1], & i>j\\ -2\eta_{jk}[1], & i=j\\ -\eta_{jk}[j-i+1]^*, & i<j\end{cases} &
G^{--}_{ij} &= \begin{cases}\eta_{kj}[i-j+1]^*, & i>j\\ \eta_{jk}[1]+\eta_{kj}[1], & i=j\\ \eta_{jk}[j-i+1]^*, & i<j\end{cases}
\end{aligned}
```

```math
G^{+\tau}_{ik} = -i\,\zeta_{kj}[i,k], \quad
G^{-\tau}_{ik} = +i\,\zeta_{kj}[i,k], \quad
G^{\tau+}_{kj} = -i\,\zeta_{jk}[k,j], \quad
G^{\tau-}_{kj} = +i\,\zeta_{jk}[k,j],
```

```math
G^{\tau\tau}_{ij} = \begin{cases}-\xi_{jk}[i-j+1], & i>j\\ -(\xi_{jk}[1]+\xi_{kj}[1]), & i=j\\ -\xi_{kj}[j-i+1], & i<j\end{cases}
```

（交叉支的相对符号与 ``\mp i`` 来自 Keldysh 回路的定向。）

### 生成向量

**实时对** ``(\eta_{jk},\eta_{kj})`` 与 Δt 的 ``G_{++}`` 块完全相同
（两种统计皆是）：

| | ``\eta_{jk}[1]`` | ``\eta_{jk}[m]`` | ``\eta_{kj}[1]`` | ``\eta_{kj}[m]`` |
|---|---|---|---|---|
| 费米 | ``-\bar n\,S_-`` | ``-\bar n\,R_m`` | ``+n\,S_+`` | ``+n\,R_m^*`` |
| 玻色 | ``-\bar b\,S_-`` | ``-\bar b\,R_m`` | ``-b\,S_+`` | ``-\bar b\,R_m^*`` |

**虚时对** ``(\xi_{jk},\xi_{kj})``（用未平移 ``\varepsilon``；与 Δτ 的生成
向量逐元素差一符号，使得 ``G^{\tau\tau}`` 恰好还原 Δτ 矩阵）：

| | ``\xi_{jk}[1]`` | ``\xi_{jk}[k]`` | ``\xi_{kj}[1]`` | ``\xi_{kj}[k]`` |
|---|---|---|---|---|
| 费米 | ``-\bar n\,T_-`` | ``-\bar n\,K_k`` | ``+n\,T_+`` | ``+n\,K_k(-\varepsilon)`` |
| 玻色 | ``-\bar b\,T_-`` | ``-\bar b\,K_k`` | ``-b\,T_+`` | ``-b\,K_k(-\varepsilon)`` |

（表中 ``K_k,T_\pm`` 即第 2 节的核，变量换成 ``\varepsilon``，统计因子含
``\mu``。）于是在 ``\mu=0`` 时有精确的边界对应
``G^{\tau\tau} = \Delta\tau`` 矩阵、``G^{\sigma\sigma'}_{\text{real}} =
\Delta t`` 的 Keldysh 块——L 形回路的实时/虚时部分分别退化为独立的
Δt / Δτ。

**混合对** ``(\zeta_{jk},\zeta_{kj})``：``\zeta_{jk}[j,k]`` 是虚时 cell
``j``（相位 ``e^{-\varepsilon\tau}``）与实时 cell ``k``（相位
``e^{+i\varepsilon t}``）的双重积分：

```math
\zeta_{jk}[j,k] = -i\!\int\mathrm{d}\varepsilon\, A(\varepsilon)\,\sigma_-(\varepsilon)\,
\frac{e^{-(j-1)\delta\tau\,\varepsilon}\,e^{i(k-1)\delta t\,\varepsilon}\,
\big(e^{-\delta\tau\,\varepsilon}-1\big)\big(e^{i\delta t\,\varepsilon}-1\big)}{\varepsilon^2},
```

```math
\zeta_{kj}[k,j] = \pm i\!\int\mathrm{d}\varepsilon\, A(\varepsilon)\,\sigma_+(\varepsilon)\,
\frac{e^{+(j-1)\delta\tau\,\varepsilon}\,e^{-i(k-1)\delta t\,\varepsilon}\,
\big(e^{\delta\tau\,\varepsilon}-1\big)\big(e^{-i\delta t\,\varepsilon}-1\big)}{\varepsilon^2},
```

其中 ``\sigma_-=\bar n,\ \sigma_+=n,\ \zeta_{kj}`` 取 ``+i``（费米）或
``\sigma_-=\bar b,\ \sigma_+=b,\ \zeta_{kj}`` 取 ``-i``（玻色）。小
``\varepsilon`` 极限下两个 ``(\ e^{\pm\cdots}-1)`` 因子各展开为一次项，
核 ``\to\delta\tau\,\delta t\,e^{\mp\cdots}``。

---

## 5. BCS 超导 bath 的杂化函数

Bogoliubov 变换已对角化的 bath，配对隙 ``\Delta``（可复数）。相干因子：

```math
E(\varepsilon)=\sqrt{\varepsilon^2+|\Delta|^2},\qquad
u(\varepsilon)=\sqrt{\tfrac{1}{2}\big(1+\tfrac{\varepsilon}{E}\big)},\qquad
v(\varepsilon)=\sqrt{\tfrac{1}{2}\big(1-\tfrac{\varepsilon}{E}\big)}\,e^{i\arg\Delta}.
```

正常谱按四种相干权重加权限色散 ``\varepsilon\mapsto E(\varepsilon)`` 复用
第 2–4 节的**费米** kernel：

```math
f \;\mapsto\; \big(u^2 f,\;\; u v\, f,\;\; v u\, f,\;\; v^2 f\big),
\qquad \text{色散 } \varepsilon\mapsto E(\varepsilon),
```

即每个块的谱积分中 ``A(\varepsilon)\to w(\varepsilon)A(\varepsilon)``（
``w\in\{u^2,uv,vu,v^2\}``），且核与统计因子中的能量全部换成
``E(\varepsilon)``。四个加权限谱的费米结果
``\Delta_{uu},\Delta_{uv},\Delta_{vu},\Delta_{vv}`` 组装成
`BCSCorrelationFunction`：

```math
\begin{aligned}
\texttt{cn} &= \Delta_{uu}-\Delta_{vv}^{\mathsf T}, &
\texttt{nc} &= \Delta_{vv}-\Delta_{uu}^{\mathsf T}, &
\texttt{cc} &= \Delta_{uv}+\Delta_{uv}^{\mathsf T}, &
\texttt{nn} &= \Delta_{vu}+\Delta_{vu}^{\mathsf T},
\end{aligned}
```

对应影响泛函指数

```math
I = \exp\Big[\bar a_\uparrow\Delta_{\uparrow\uparrow}a_\uparrow
+ \bar a_\uparrow\Delta_{\uparrow\downarrow}\bar a_\downarrow
+ a_\downarrow\Delta_{\downarrow\uparrow}a_\uparrow
+ a_\downarrow\Delta_{\downarrow\downarrow}\bar a_\downarrow\Big],
```

即正常块 ``\texttt{cn}/\texttt{nc}``（``\bar a a`` 型，字段别名
``uu/dd``）与反常块 ``\texttt{cc}/\texttt{nn}``（``\bar a\bar a``、
``aa`` 型，字段别名 ``ud/du``）。``\Delta=0`` 时
``u^2=\theta(\varepsilon)``、``v^2=\theta(-\varepsilon)``、``uv=vu=0``，
反常块消失，正常块还原普通费米结果（由测试保证）。

---

## 6. 无穷温窗（β→∞ / t→∞）

`infinite_Δτ` / `infinite_Δt` 使用与第 2、3 节**完全相同**的 kernel，取
``\beta=\infty``：此时 ``\bar n,n\to\theta(\mp\varepsilon)``（费米真空）
或 ``b\to 0,\ \bar b\to\theta(-\varepsilon)\cdot\delta(\varepsilon)^{-1}``
（玻色真空，只有 ``\varepsilon\le 0`` 的发射部分）。Toeplitz 生成向量按节点
增量生成（`compute_next!`），直到最新元素幅度小于容差 `atol`；实时 cache
通过 `branch` 按需组合出四个 Keldysh 块（共轭关系同第 3 节表格）。

---

## 7. Matsubara 频域入口

`Δiw_to_Δτ` 从频域数据直接构造虚时矩阵（推导见 `docs/quapi.pdf`）：

```math
\Delta(\tau_k) = \sum_{n=-n_{max}}^{n_{max}}
-\frac{2}{\beta\,\omega_n^2}\,\Delta(i\omega_n)\,
\big(1-\cos(\omega_n\delta\tau)\big)\,e^{-i\omega_n k\delta\tau},
\qquad \delta\tau=\beta/N,
```

要求 ``\Delta iw`` 为按 `ifrequencies(β, nmax)` 顺序的偶数个频率值，虚部
不为零时报错。
