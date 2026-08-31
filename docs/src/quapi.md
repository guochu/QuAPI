# QuAPI Algorithms

本页给出 QuAPI 路径积分算法的**参考实现**，位于 `test/quapi/*.jl`。这些
函数未导出，使用时 include 相应文件。每个函数先取 [API](api.md) 的离散杂化
矩阵，再做完整的路径求和，并与严格对角化（ED）或解析解交叉验证。

符号约定（见 [Home](index.md)）：离散路径权重为 `exp[+sᵀηs]`，IF 负号已吸收
进杂化矩阵。

## 虚时 QuAPI

### 自旋-玻色模型（增广密度张量）

```julia
quapi_adt(β, ε, Δ, ωs, λs; N) -> (Z, C, δτ)
```

模型：

```
H = ε/2 σz + Δ/2 σx + Σ_k ω_k a_k†a_k + (σz/2) Σ_k λ_k (a_k† + a_k)
```

- **参数**：`β` 逆温度；`ε`、`Δ` 自旋能隙与隧穿；`ωs`、`λs` 离散玻色
  谱频率与耦合；`N` 虚时步数（δτ = β/N）。
- **实现**：`bosonic_Δτ(DiscreteSpectrum(ωs, λs.²); β, N)` 给出 η；枚举
  2^N 个自旋路径（s_j ∈ {±1/2}，周期边界 s_N = s_0），装配秩-N、物理维 2
  的增广密度张量（ADT）；路径权重 = `exp[+sᵀηs] × Π W`（`W = e^{−δτHs}`
  为自由自旋短时传播子）。
- **返回**：配分函数 `Z`；Matsubara 关联 `C[m+1] = ⟨σz(τ_m)σz(0)⟩`，
  τ_m = m·δτ；步长 `δτ`。

### 独立玻色子极限（Δ = 0）

```julia
quapi_ibm(β, ε, ωs, λs; N) -> (Z, G, δτ)
```

- σz 守恒时唯一非平凡的 Matsubara 关联是自旋翻转（费米型）
  `G(τ_m) = ⟨d(τ_m)d†(0)⟩`，d = σ₋。由**单畴壁路径**精确给出：自旋在
  [0, τ_m) 取 +1/2、[τ_m, β) 取 −1/2，权重
  `W[1,1]^m W[2,2]^(N−m) · exp[+sᵀηs] / Z`。
- `Z` 为除以自由 bath 配分函数后的值（`Ση = βΣλ²/ω` 精确成立）；关联
  函数（比值）不受影响。
- **验证**：与 ImpurityModelBase 的 `independentbosons_Gτ` 解析解对比。

### Toulouse 共振能级模型（Grassmann 张量）

```julia
quapi_rlm_gt(β, ϵ_d, spectrum::AbstractBoundedFunction; N, μ=0.0) -> (Z, G, δτ)
quapi_rlm_gt(β, ϵ_d, ωs::AbstractVector{<:Real}, λs::AbstractVector{<:Real}; N, μ=0.0)
# 便捷方法：自动构造 DiscreteSpectrum(ωs, λs.²)
```

- **实现**：路径权重为双线性指数 `exp[Σ āᵢ (M_imp + Δ)ᵢⱼ aⱼ]`，其中
  `M_imp` 是反周期动能矩阵（`x = e^{−δτϵ_d}`，环绕元带 −x），`Δ` 来自
  `fermionic_Δτ`。表示为 2N 条腿的稠密 Grassmann 张量（系数张量存
  2^L 个 monomial 系数，避免 GrassmannTensors `hadamardproduct` 的
  O(3^L) 融合树开销）。
- **Berezin 积分读取**：`Z = (−1)^N × 全奇顶系数`；
  `G(τ_m) = (ā₀, a_m 置 even 的系数) / 顶系数`。
- **返回约定**：标准约定 `G = −⟨d(τ)d†(0)⟩`。注意 ImpurityModelBase
  的 `toulouse_Gτ` 用 raw 约定 `+⟨…⟩`，两者差一个负号。
- **验证**：自由杂质极限（λ = 0）与 `toulouse_Gτ`、耦合情形与 ED。

## 实时 QuAPI

### 自旋-玻色模型（闭合 Keldysh 回路）

```julia
quapi_rt(β, ε, Δ, ωs, λs; N, t, ρ_imp, μ=0.0) -> (σz_t, σx_t, ρ, δt)
```

- **初态**：`ρ(0) = ρ_imp ⊗ ρ_bath(β)`（热 bath，t=0 时接入相互作用）。
- **实现**：取 `bosonic_Δt` 的四个 (N+1)×(N+1) 分支块 G₊₊/G₊₋/G₋₊/G₋₋，
  对 2^N × 2^N 前后向自旋路径完整求和，权重

```
I[s⁺, s⁻] = exp[Σ_jk (s⁺_j s⁺_k G₊₊ + s⁺_j s⁻_k G₊₋
                     + s⁻_j s⁺_k G₋₊ + s⁻_j s⁻_k G₋₋)]
```

- **返回**：`⟨σz(t)⟩`、`⟨σx(t)⟩`、约化密度矩阵 `ρ`、`δt = t/N`。
- **验证**：`ed_rt_obs`（截断 Fock 空间 ED）。

### Toulouse 模型实时（Grassmann 张量）

```julia
quapi_rlm_rt_gt(β, ϵ_d, spectrum; Nt, t, μ=0.0) -> (Z, Gg, Gl, δt)
```

- 单个 Grassmann 张量，4(Nt+1) 条腿：`lp(i) = 2i+1`（ā⁺ᵢ）、
  `ap(i) = 2i+2`（a⁺ᵢ）、`lm(i) = 2(Nt+1)+2i+1`（ā⁻ᵢ）、`am(i)`（a⁻ᵢ）。
- 指数矩阵 `M` 依次包含：各节点 Berezin 测度 `−āa`；自由动能
  `x = e^{−iδtϵ_d}`（前向）/ `x̄`（后向）；回路闭合 `−ā⁻(t)a⁺(t)`；初始
  热权重 `e^{−βϵ_d} ā⁺(0)a⁻(0)`（quench）；`fermionic_Δt` 的 IF 块。
- **返回**：bath 归一化的 `Z`；**raw** 期望
  `Gg[j] = ⟨d(t_j)d†(0)⟩`、`Gl[j] = ⟨d†(0)d(t_j)⟩`（对应
  `im·G^>` / `−im·G^<`）。
- **验证**：与 `toulouse_neq_greater_lesser`、`freefermion_greater/lesser`
  对比。

## 混合时（Kadanoff L 形回路）QuAPI

### 自旋-玻色模型

```julia
quapi_mixed(β, ε, Δ, ωs, λs; Nt, Nτ, t, μ=0.0)
    -> (σz_t, σx_t, Czz, ρ, δt, δτ)
```

- **初态**：完整热态 `e^{−βH}/Z`（区别于 `quapi_rt` 的 quench 初态）。
- **回路**：三支——实时前向（+，t: 0→t）、后向（−，t→0）、虚时
  （τ: 0→β）。接点处 s⁺₀（τ=β 端）与 s⁻₀（τ=0 端）闭合，虚时传播链
  s⁻₀ → c₂ .. c_Nτ → s⁺₀。
- **实现**：`bosonic_Δm` 的九个分支块（η/ξ/ζ）拼成单个块矩阵 `G`，
  自旋路径向量 `s` 上权重 `exp(sᵀGs)`。要求 Nτ ≥ 2；自由自旋构型数
  2^(2Nt+Nτ−1)。
- **返回**：`⟨σz(t)⟩`、`⟨σx(t)⟩`、双时关联
  `Czz = ⟨σz(t)σz(0)⟩`（σz 插在接点 s⁺₀ 处）、非归一约化密度矩阵
  `ρ = e^{−iHt} e^{−βH} e^{iHt}`（对 bath 求迹）、`δt`、`δτ`。
- **验证**：`ed_mixed_obs`（ED，`Czz = tr[ρ₀ σz(t) σz]`）。

## BCS 虚时 QuAPI

### BCS bath 上的 Toulouse 模型

```julia
quapi_bcs_gt(β, ϵ_d, spectrum::AbstractBoundedFunction; N, Δ, μ=0.0) -> (Z, G, δτ)
```

模型：

```
H = Σ_σ ϵ_d d_σ†d_σ + Σ_kσ ω_k c_kσ†c_kσ + Σ_kσ V_k (d_σ†c_kσ + c_kσ†d_σ)
    + Δ Σ_k (c_k↑†c_k↓† + c_k↓c_k↑)
```

- **实现**：离散路径权重
  `exp[ā↑(M_imp+Δ↑↑)a↑ + ā↓(M_imp+Δ↓↓)a↓ + ā↑Δ↑↓ā↓ + a↓Δ↓↑a↑]`，
  Δ 块来自 `bcs_Δτ`（`corr.uu/.dd/.ud/.du`）。4N 条腿 Grassmann 张量：
  flavor ↑ 在腿 2i−1（ā）、2i（a），flavor ↓ 在腿 2N+2i−1、2N+2i。
- **Berezin 读取**：`Z = (−1)^(2N) × 顶系数`；
  `G(τ_m) = (ā↑₀, a↑_m 置 even 的系数)/顶系数`。
- **返回约定**：`G = −⟨d↑(τ_m)d↑†(0)⟩`（与 `toulouse_Gτ` 的 raw 约定差
  一负号）。
- **验证**：无解析解，参考为 `Toulouse(DiscreteBCSBath(...))` 类的严格
  BdG ED。

## BCS 实时 QuAPI（Pfaffian 精确求积）

BCS 实时路径积分用 Grassmann 高斯积分的 **Pfaffian 表示**精确求值，避免
张量收缩的组合爆炸。

### 数学工具

```julia
pfaffian(A::AbstractMatrix)                     # 反对称矩阵 Pfaffian（Parlett–Reid，O(n³)）
gauss_integral(B::AbstractMatrix, monomial::Vector{Int})
# ∫ exp(½ξᵀBξ) ξ_{i₁}···ξ_{i_{2k}} = Pf(B) · Pf( (B⁻¹)[monomial, monomial] )
```

### 变量索引

```julia
quark_single(N) / aindex_single(N)     # 单 flavor：āᵢ = base+2i+1, aᵢ = base+2i+2
                                       # base(:+) = 0, base(:-) = 2(N+1)
quark_bcs(N, σ) / aindex_bcs(N, σ)     # BCS 双 flavor：flavor ↓ 偏移 4(N+1)
```

### 反对称矩阵构造

```julia
build_B_single(N, β, ϵ_d, t, spec, μ)          # 单 flavor，4(N+1) 维
build_B_bcs(N, β, ϵ_d, t, spec, Δ, μ)          # BCS，8(N+1) 维
```

`B` 依次装配：各节点 Berezin 测度 `−āa`；前/后向跳跃 `e^{∓iδtϵ_d}`；
`fermionic_Δt` / `bcs_Δt` 的 IF 块（交叉支带 Keldysh 相对符号 −1）；
BCS 版另含反常块 `ā↑Δ↑↓ā↓`、`a↓Δ↓↑a↑`。

### 求值函数

```julia
quapi_rt_single(N, β, ϵ_d, t, spec, ρ_imp, μ=0.0) -> (ρ ./ Z, Z)
# 单 flavor 2×2 约化密度矩阵（Fock 基 |0⟩, |1⟩），边界单项式按 Keldysh
# 回路序 0⁺ → N⁺ → N⁻ → 0⁻

quapi_rt_bcs(N, β, ϵ_d, t, spec, Δ, ρ_imp, μ=0.0) -> (ρ ./ Z, Z)
# BCS 4×4 约化密度矩阵（Fock 基 |0⟩, |↑⟩, |↓⟩, |↑↓⟩，|↑↓⟩ = ā↓†ā↑†|0⟩）
```

### ED 参考

```julia
bcs_ed_imb(β, ϵ_d, ωs, Vs, Δ, nsys, ts; μ=0.0)
# ImpurityModelBase BdG 严格 ED：Toulouse(DiscreteBCSBath(ωs, Vs.², Δ; β, μ))
# + cmatrix/eigencache/separablecdm/timeevo；返回各时刻 (⟨n_↑⟩, ⟨n_↓⟩)；要求 μ=0
```

## ED 参考函数汇总

| 函数 | 位置 | 用途 |
|---|---|---|
| `ed_hamiltonian(ε, Δ, ωs, λs; d)` | `imagquapi.jl` | 截断 Fock 空间完整自旋-玻色 H |
| `ed_rt_obs(ε, Δ, ωs, λs, ρ_imp, ts; β, d=12)` | `realquapi.jl` | 实时 (σz, σx) ED |
| `ed_mixed_obs(ε, Δ, ωs, λs, ts; β, d=12)` | `mixedquapi.jl` | 混合时 (σz, σx, Czz) ED |
| `bcs_ed_imb(...)` | `bcsrealquapi.jl` | BCS BdG ED（⟨n_↑⟩,⟨n_↓⟩） |

## Green 函数约定对照

| 函数 | 约定 |
|---|---|
| `quapi_rlm_gt` / `quapi_bcs_gt` | 标准 `G = −⟨d(τ)d†(0)⟩` |
| ImpurityModelBase `toulouse_Gτ` | raw `G = +⟨d(τ)d†(0)⟩` |
| `quapi_rlm_rt_gt` 的 `Gg`/`Gl` | raw 期望 `⟨d(t)d†(0)⟩` / `⟨d†(0)d(t)⟩` |
| ImpurityModelBase `freefermion_greater/lesser` | `G^> = −im⟨d(t)d†(0)⟩` 等 |

对比时统一换算（如 `Gg ≈ im·G^>`）。
