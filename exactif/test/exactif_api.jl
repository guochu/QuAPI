# =============================================================================
# exactif/test/exactif_api.jl —— 连续谱严格核的 QuAPI 接口
#
#   exact_fermionic_Δt(f; β, N, t, μ; M, nsub)  -> RealCorrelationFunction
#   exact_fermionic_Δτ(f; β, N, μ, δτ; M)       -> ImagCorrelationFunction
#
# 输入/输出契约与包内 fermionic_Δt/Δτ 一致（f: AbstractBoundedFunction；
# 节点 0..N；minus absorbed；纯浴对象）。浴物理只经连续杂化函数 J(ε) 进入：
# 所有谱求和以 GL-M 节点做**求积**（M 是积分节点数，非浴模式——
# 收敛到连续谱由 test 的 M-收敛研究对照 ImpurityModelBase 解析解验证）。
# 条目 = Keldysh/Matsubara 轮廓上杂质模的严格两点收缩（语义锚定：
# 跨支条目 = 已验证 wick Gram 的节点化；纯浴 ⟹ 杂质自项不入 Δ；
# 虚时负 τ = 反周期延拓 −G(β+τ)，闭合逐位）。
# =============================================================================
include(joinpath(@__DIR__, "..", "brute.jl"))
include(joinpath(@__DIR__, "..", "volterra.jl"))

using LinearAlgebra: I, Diagonal, eigen, Hermitian
using QuAPI, ImpurityModelBase
using QuadGK: quadgk

"GL-M 求积节点上的谱数据（求积节点，非物理模式）。μ = 谱平移。"
function _quad_nodes(f::AbstractBoundedFunction, μ::Real, M::Int)
    lb, ub = lowerbound(f), upperbound(f)
    xs, ws = gauss_legendre(M)
    half = (ub - lb) / 2
    Es = [(lb + ub) / 2 + half * x - μ for x in xs]
    wΓ = [half * w * f((lb + ub) / 2 + half * x) for (x, w) in zip(xs, ws)]
    Vsq = sqrt.(max.(wΓ, 0.0))
    return Es, Vsq, wΓ
end

"纯浴 h₀（杂质位 1，ε=0）本征分解。"
function _bath_eig(Es, Vsq)
    M = length(Es)
    h0 = zeros(ComplexF64, M + 1, M + 1)
    for k in 1:M
        h0[1, 1+k] = Vsq[k]; h0[1+k, 1] = conj(Vsq[k]); h0[1+k, 1+k] = Es[k]
    end
    return eigen(Hermitian(h0))
end

# ------------------------------ 虚时 ------------------------------ #

"""
    exact_fermionic_Δτ(f; β, N, μ=0, δτ=β/N; M=64) -> ImagCorrelationFunction

热环收缩（模式离散版）：GL-M 节点谱 → h₀ 本征分解 →
G₀(τ>0) = Σₚ v₁ₚ² e^{−λₚτ}/(1+e^{−βλₚ})（机器精度谱闭式）；
负 τ 反周期 G(τ) = −G₀(β+τ)（闭合逐位）。节点 τᵢ = (i−1)δτ。
"""
function exact_fermionic_Δτ(f::AbstractBoundedFunction; β::Real, N::Int,
                            μ::Real=0, δτ::Real=β / N, M::Int=64)
    Es, Vsq, _ = _quad_nodes(f, μ, M)
    ev = _bath_eig(Es, Vsq)
    a2 = abs2.(ev.vectors[1, :])
    # τ>0 粒子线；τ≤0 反周期延拓（ηₖⱼ[1] = G(0⁻) = −G₀(β⁻)，语义正确）
    G0(τ) = τ > 0 ?
            sum(a2[p] * exp(-ev.values[p] * τ) / (1 + exp(-β * ev.values[p])) for p in 1:M+1) :
            -sum(a2[p] * exp(-ev.values[p] * (τ + β)) / (1 + exp(-β * ev.values[p])) for p in 1:M+1)
    ηⱼₖ = [G0((d - 1) * δτ) for d in 1:N]
    ηₖⱼ = [G0(-(d - 1) * δτ) for d in 1:N]
    return ImagCorrelationFunction(CorrelationMatrix(ηⱼₖ, ηₖⱼ))
end

# ------------------------------ 实时 ------------------------------ #

"""
    exact_fermionic_Δt(f; β, N, t, μ=0; M=64) -> (RealCorrelationFunction, blocks)

Keldysh 轮廓两点收缩，四块 (N+1)×(N+1)（一般 Matrix 存储，开链非 Toeplitz）：
  节点 (i,u) → 收缩向量 c = W†(t_辅段)·W(主段)·e₁（谱求积节点上的本征机器）；
  收缩 = Σ_l (1−n̄_l)c_v[l]*c_w[l]（v 轮廓位 ≥ w）或 −Σ n̄（w 更晚），
  l ≥ 2（纯浴：杂质位自项不入 Δ）。块内符号约定经 wick-Gram 锚定。
"""
function exact_fermionic_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real,
                            μ::Real=0, M::Int=64)
    δt = t / N
    Es, Vsq, _ = _quad_nodes(f, μ, M)
    ev = _bath_eig(Es, Vsq)
    W(θ) = θ == 0 ? Matrix{ComplexF64}(I, M+1, M+1) :
           ev.vectors * Diagonal(exp.(-im * ev.values .* θ)) * ev.vectors'
    Wd(θ) = θ == 0 ? Matrix{ComplexF64}(I, M+1, M+1) :
            ev.vectors * Diagonal(exp.(im * ev.values .* θ)) * ev.vectors'
    n̄ = vcat([0.5], [1 / (1 + exp(β * e)) for e in Es])
    n1 = N + 1
    e1 = vcat([1.0], zeros(M))
    cvec = Vector{Vector{ComplexF64}}(undef, 2 * n1)
    for i in 0:N
        cvec[i+1] = Wd(t) * W(t - i * δt) * e1          # 正向支到轮廓末端
        cvec[n1 + i + 1] = Wd(i * δt) * e1              # 反向支
    end
    pos(u, i) = u == 1 ? i : N + (N - i) + 1
    function pair(v, w)
        cv, cw = cvec[v], cvec[w]
        pv = pos(v ≤ n1 ? 1 : 2, (v - 1) % n1)
        pw = pos(w ≤ n1 ? 1 : 2, (w - 1) % n1)
        if pv == pw
            return sum((1 - n̄[l]) * conj(cv[l]) * cw[l] for l in 2:M+1)
        elseif pv > pw
            return sum(n̄[l] * cv[l] * conj(cw[l]) for l in 2:M+1)
        else
            return sum((1 - n̄[l]) * conj(cv[l]) * cw[l] for l in 2:M+1)
        end
    end
    Gpp = [pair(j + 1, i + 1) for i in 0:N, j in 0:N]
    Gpm = [pair(n1 + j + 1, i + 1) for i in 0:N, j in 0:N]
    Gmp = [pair(j + 1, n1 + i + 1) for i in 0:N, j in 0:N]
    Gmm = [pair(n1 + j + 1, n1 + i + 1) for i in 0:N, j in 0:N]
    return RealCorrelationFunction(Gpm, Gmp), (Gpp=Gpp, Gpm=Gpm, Gmp=Gmp, Gmm=Gmm)
end

"""
    exact_fidelity(f; ε_d=0, t; nsub=800) -> F(t)（Volterra，连续谱 Δ(t) 谱积分）

论文 (S10) 的保真度（含 ε_d 线项），供 U=0 Toulouse 解析参照与测试。
"""
function exact_fidelity(f::AbstractBoundedFunction; ε_d::Real=0.0, T::Real, nsub::Int=800)
    lb, ub = lowerbound(f), upperbound(f)
    Nf = max(4, round(Int, nsub * T))
    h = T / Nf
    tg = range(0, T; length=Nf + 1)
    Δ = [sum(ε -> f(ε) * exp(-im * ε * tq), lb, ub; rtol=1e-11) for tq in tg]
    F = zeros(ComplexF64, Nf + 1); F[1] = 1
    Cp = 0.0im
    for j in 0:Nf-1
        R = 0.5 * Δ[2] * F[j+1]
        for i in 0:j-1
            w = (i == 0) ? 0.5 : 1.0
            R += w * Δ[j+2-i] * F[i+1]
        end
        F[j+2] = (F[j+1] - (h / 2) * Cp - (h^2 / 2) * R - im * ε_d * h * F[j+1]) /
                 (1 + h^2 * Δ[1] / 4 + im * ε_d * h / 2)
        Cp = h * (0.5 * Δ[1] * F[j+2] + R)
    end
    Finterp(t) = begin
        x = clamp(t / h, 0, Nf)
        i = min(Nf - 1, floor(Int, x))
        F[i+1] + (x - i) * (F[i+2] - F[i+1])
    end
    return Finterp
end

# ------------------------------ 混合时间 ------------------------------ #

"""
    exact_fermionic_Δm(f; β, Nτ, t, Nt, μ=0, δτ=β/Nτ; M=64)
        -> MixedCorrelationFunction

L 形（Kadanoff-Baym）轮廓的连续谱严格核。分支序：τ 支（0→β，Nτ 节点）
最前，其后 Keldysh 部分与实时版逐位相同。九块：
  G₊₋/G₋₊：与 exact_fermionic_Δt 的交叉块一致；
  Gττ：与 exact_fermionic_Δτ 一致（环 Toeplitz）；
  G₊τ (Nt+1)×Nτ、Gτ₊ Nτ×(Nt+1)：L 轮廓 τ-交叉收缩——τ 顶点的收缩向量
  = e^{−(β−τᵢ)h₀}e₁（W(t)W†(t)=I 使 Keldysh 闭环相消）；轮廓序 τ 最早
  ⟹ ⟨C_A A_τ⟩ = Σ n̄ c_a c̄_τ、⟨A_a C_τ⟩ = Σ (1−n̄) c̄_a c_τ。
  (−,τ)/(τ,−) 的 Keldysh 相对负号由 branch 规则承担（存 G₊τ/Gτ₊）。
"""
function exact_fermionic_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real,
                            Nt::Int, μ::Real=0, δτ::Real=β / Nτ, M::Int=64)
    δt = t / Nt
    Es, Vsq, _ = _quad_nodes(f, μ, M)
    ev = _bath_eig(Es, Vsq)
    W(θ) = θ == 0 ? Matrix{ComplexF64}(I, M+1, M+1) :
           ev.vectors * Diagonal(exp.(-im * ev.values .* θ)) * ev.vectors'
    Wd(θ) = θ == 0 ? Matrix{ComplexF64}(I, M+1, M+1) :
            ev.vectors * Diagonal(exp.(im * ev.values .* θ)) * ev.vectors'
    Wim(θ) = ev.vectors * Diagonal(exp.(-ev.values .* θ)) * ev.vectors'   # 虚时演化
    n̄ = vcat([0.5], [1 / (1 + exp(β * e)) for e in Es])
    n1, e1 = Nt + 1, vcat([1.0], zeros(M))
    # 收缩向量：forward a、backward a 同 Δt
    cfwd = [Wd(t) * W(t - a * δt) * e1 for a in 0:Nt]
    cbwd = [Wd(a * δt) * e1 for a in 0:Nt]
    posf(a) = Nτ + a
    posb(a) = Nτ + Nt + (Nt - a)
    # τ-交叉块（本征基谱闭式，2 模多体 4 点逐位锚定 mx_fix12）：
    #   G₊τ[a+1,i] = Σₚ a2ₚ e^{−λₚτᵢ} e^{−iλₚt_a} fₚ   （粒子线到 C 端相位）
    #   Gτ₊[i,a+1] = −(G₊τ)ᵀ（μ=0 反对称）。
    #   注：V→0 时 G₊τ → ½（杂质自线）——系统腿结构，由调用方测度抵消
    #   （同 firstorder README §2.2 的 c-number 抵消簿记）。
    a2 = abs2.(ev.vectors[1, :])
    fk = [1 / (1 + exp(-β * ev.values[p])) for p in 1:M+1]
    G₊τ = [sum(a2[p] * exp(-ev.values[p] * ((i - 1) * δτ)) *
               exp(-im * ev.values[p] * (a * δt)) * fk[p] for p in 1:M+1)
           for a in 0:Nt, i in 1:Nτ]
    # μ=0 费米反对称 ζⱼₖ = −ζₖⱼᵀ（物理要求，docstring 语义）：
    Gτ₊ = μ == 0 ? Matrix(-transpose(G₊τ)) :
           [sum(a2[p] * exp(-ev.values[p] * ((i - 1) * δτ)) *
                exp(im * ev.values[p] * (a * δt)) * fk[p] for p in 1:M+1)
            for i in 1:Nτ, a in 0:Nt]
    # Keldysh 交叉块（与 Δt 同构造）
    function pair(v, w)
        cv, cw = v ≤ n1 ? cfwd[v] : cbwd[v - n1], w ≤ n1 ? cfwd[w] : cbwd[w - n1]
        pv = v ≤ n1 ? posf(v - 1) : posb(v - n1 - 1)
        pw = w ≤ n1 ? posf(w - 1) : posb(w - n1 - 1)
        if pv == pw
            return sum((1 - n̄[l]) * conj(cv[l]) * cw[l] for l in 2:M+1)
        elseif pv > pw
            return sum(n̄[l] * cv[l] * conj(cw[l]) for l in 2:M+1)
        else
            return sum((1 - n̄[l]) * conj(cv[l]) * cw[l] for l in 2:M+1)
        end
    end
    Gpm = [pair(n1 + j + 1, i + 1) for i in 0:Nt, j in 0:Nt]
    Gmp = [pair(j + 1, n1 + i + 1) for i in 0:Nt, j in 0:Nt]
    # ττ 块 = 环 Toeplitz（同 Δτ）
    a2 = abs2.(ev.vectors[1, :])
    G0(τ) = τ > 0 ?
            sum(a2[p] * exp(-ev.values[p] * τ) / (1 + exp(-β * ev.values[p])) for p in 1:M+1) :
            -sum(a2[p] * exp(-ev.values[p] * (τ + β)) / (1 + exp(-β * ev.values[p])) for p in 1:M+1)
    Gττ = CorrelationMatrix(ComplexF64[real(G0((d - 1) * δτ)) for d in 1:Nτ],
                            ComplexF64[real(G0(-(d - 1) * δτ)) for d in 1:Nτ])
    Gττ = CorrelationMatrix(ComplexF64.(Gττ.ηⱼₖ), ComplexF64.(Gττ.ηₖⱼ))
    return MixedCorrelationFunction(Gpm, Gmp, G₊τ, Gτ₊, Gττ)
end
