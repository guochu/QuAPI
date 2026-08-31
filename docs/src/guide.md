# Guide

本页给出从谱函数到离散杂化矩阵的典型工作流。所有输入类型由
ImpurityModelBase 提供。

## 谱函数与 bath

谱函数（`AbstractBoundedFunction`）调用语义：`x(ε)` 在
`[lowerbound(x), upperbound(x)]` 内取 `f(ε)`，否则为 0。

常用预定义谱：

```julia
using QuAPI

# 费米（Bethe 格子半圆谱），支撑 [-t, t]
f = semicircular(t = 1.0)

# 玻色（Ohmic/sub-Ohmic），支撑 [0, 100ωc]
j = Leggett(α = 0.1, ωc = 1.0)

# 解析极限：单 δ 峰
δ = DiracDelta(ω = 0.5)

# 离散 bath（fs = λₙ² 耦合强度）
spec = DiscreteSpectrum([1.0, 2.0], [0.25, 0.09])
```

谱函数可以配上温度与化学势打包成 bath：

```julia
β, μ = 5.0, 0.3
bath = fermionicbath(f; β = β, μ = μ)     # 费米
bbath = bosonicbath(j; β = β, μ = 0.0)    # 玻色
bcs = bcsbath(f; β = β, Δ = 0.5, μ = μ)   # BCS（配对隙 Δ 可复数）
```

## 虚时（Matsubara）离散化

```julia
N = 10                              # 步数，δτ = β/N
corr = fermionic_Δτ(f; β = β, N = N, μ = μ)   # 或 bosonic_Δτ
corr = Δτ(bath; N = N)             # bath convenience，自动提取 β/μ
```

返回 [`ImagCorrelationFunction`](@ref)，包装一个 N×N 实型 Toeplitz
[`CorrelationMatrix`](@ref)。矩阵元 `corr[i, j]` 即 Δ_{ij}。

Matsubara 频域数据可以用 [`Δiw_to_Δτ`](@ref) 直接变换：

```julia
Δiw = [1.0 / (iω - 1.2) for iω in ifrequencies(β, 64)]
corr = Δiw_to_Δτ(Δiw; β = β, N = N)
```

## 实时（Keldysh）离散化

```julia
t, N = 2.0, 20                      # δt = t/N
corr = bosonic_Δt(j; β = β, N = N, t = t, μ = 0.0)
Gpp = branch(corr; b1 = :+, b2 = :+)     # 时序块 (N+1)×(N+1)
Gpm = branch(corr; b1 = :+, b2 = :-)     # lesser 型
Gmp = branch(corr; b1 = :-, b2 = :+)     # greater 型
Gmm = branch(corr; b1 = :-, b2 = :-)     # 反时序
```

四个分支也可用 `corr.G₊₊` 等字段直接访问。

## 混合时（Kadanoff L 形回路）

```julia
Nτ, Nt, t = 10, 4, 2.0             # δτ = β/Nτ, δt = t/Nt
corr = bosonic_Δm(j; β = β, Nτ = Nτ, t = t, Nt = Nt, μ = 0.0)
Gpt = branch(corr; b1 = :+, b2 = :τ)     # 混合块 (Nt+1)×Nτ
Gττ = branch(corr; b1 = :τ, b2 = :τ)     # 虚时块 Nτ×Nτ
```

## 无穷温窗

β → ∞（真空）或 t → ∞ 时 Toeplitz 生成向量可逐节点按需生成：

```julia
vac = BosonicVacuum(j)
cache = infinite_Δτ(vac; δτ = 0.01)       # 生长到相关函数衰减
compute!(cache; atol = 1e-8)              # 继续生长
n = current_size(cache)                   # 当前节点数

rcache = infinite_Δt(bath; δt = 0.05)
Gpp = branch(rcache, :+, :+)              # Keldysh 块按需组合
```

## BCS 超导 bath

Bogoliubov 变换已对角化的 bath，配对隙 `Δ` 通过相干因子
`u(ϵ) = √((1+ϵ/E)/2)`、`v(ϵ) = √((1−ϵ/E)/2)·e^{i·arg Δ}`（E = √(ϵ²+|Δ|²)）
进入四种加权限权谱 `(u², uv, vu, v²)·f`，并以色散 `ϵ ↦ E(ϵ)` 复用费米核：

```julia
corr = bcs_Δτ(f; β = β, N = N, Δ = 0.5, μ = μ)
corr.uu    # 正常块 ā↑ Δ↑↑ a↑  (= cn)
corr.ud    # 反常块 ā↑ Δ↑↓ ā↓  (= cc)
corr.du    # 反常块 a↓ Δ↓↑ a↑  (= nn)
corr.dd    # 正常块 a↓ Δ↓↓ ā↓  (= nc)
```

`Δ = 0` 时退化为普通费米结果（`ud`/`du` 块为零）。

## 注意事项

1. **Bose 零频奇异性**：`μ = 0` 且谱支撑触及 ε = 0 时 `bosonic_*` 求积抛
   `DomainError`。
2. **等时对角元**：[`CorrelationMatrix`](@ref) 对角元为两个排序之和
   `ηⱼₖ[1] + ηₖⱼ[1]`，仅对对角耦合成立。
3. **化学势**：`μ` 经 `spectrumshift` 吸收进谱轴后再积分。
4. **色散参数**：`fermionic_Δτ(f, disperse; ...)` 等方法的 `disperse`
   在积分变量上施加可选色散（默认 identity）；BCS 谱内部即以
   `ε ↦ √(ε²+|Δ|²)` 调用。
