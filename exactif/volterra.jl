# =============================================================================
# exactif/volterra.jl —— 保真度 F(t) 的求解（Ng–Reichman 严格核的第一层）
#
# F(t) = [U(t)]₀₀，U(t) = exp(−it·h)，h = [[0, Vᵀ],[V̄ᵀ, diag(E)]]（ε=0 分割）。
# 满足二阶 Volterra 积分微分方程（论文式 (S10)，ε_σ = 0）：
#   F'(t) = −∫₀ᵗ Δ(t−τ) F(τ) dτ，  F(0) = 1，
#   Δ(t) = Σ_k |V_k|² e^{−iE_k t}  →  连续谱：Δ(t) = ∫ dε Γ(ε) e^{−iεt}。
#
# 两条路径：
#   1) 离散模式（Dirac 谱）：Laplace 逆变换闭式 —— F(t) = Σ_j e^{z_j t}/Q'(z_j)，
#      z_j 为 Q(z) = z + Σ_k g_k²/(z + iE_k) = 0 的根（伴随矩阵特征值），
#      机器精度（测试用）。
#   2) 连续谱：细网格乘积梯形（product trapezoid）前向推进，
#      O(h²)，配合 nsub 加密 + Richardson 可达谱积分精度。
#      M=1 解析极限（单元测试）：
#      F(t) = e^{−iEt/2}[cos(νt/2) + i(E/ν)sin(νt/2)]，ν = √(E²+4g²)。
# =============================================================================
using LinearAlgebra: eigen, exp, SymTridiagonal

"""
    fidelt_poles(Es, Vs, ts) -> Vector{ComplexF64}

离散模式的闭式 F(t)（留数法）：F(t) = Σ_j e^{z_j t}/Q'(z_j)。
Q(z) = z + Σ_k g_k²/(z+iE_k)，g_k = |V_k|（相位可重整掉，只 |V| 进入）。
根 z_j 经由结式多项式 P(z) = z·Π(z+iE_k) + Σ_k g_k² Π_{l≠k}(z+iE_l)
的伴随矩阵特征值求得。
"""
function fidelt_poles(Es::AbstractVector{<:Real}, Vs::AbstractVector{<:Number},
                      ts::AbstractVector{<:Real})
    M = length(Es)
    g2 = abs2.(Vs)
    # P(z) 的系数（升幂）：P = Π(z+iE) · [z + Σ g²/(z+iE)]
    #                    = z·Π(z+iE) + Σ_k g_k² Π_{l≠k}(z+iE)
    # 用多项式乘法构造（一元多项式 = 系数向量，升幂）
    poly_prod(factors) = reduce((p, q) -> convpoly(p, q), factors; init=[1.0 + 0im])
    function convpoly(p, q)
        r = zeros(ComplexF64, length(p) + length(q) - 1)
        for i in eachindex(p), j in eachindex(q)
            r[i+j-1] += p[i] * q[j]
        end
        return r
    end
    lin = [[im * E, 1.0 + 0im] for E in Es]   # 升幂 (iE + z)           # (z + iE)
    base = poly_prod(lin)                              # Π(z+iE)
    P = zeros(ComplexF64, M + 2)
    for i in eachindex(base)
        P[i+1] += base[i]                              # z·Π(z+iE)
    end
    for k in 1:M
        others = poly_prod(lin[1:end .!= k])
        for i in eachindex(others)
            P[i] += g2[k] * others[i]                  # g_k² Π_{l≠k}
        end
    end
    # 结式根：伴随矩阵
    n = length(P) - 1
    comp = zeros(ComplexF64, n, n)
    for i in 1:n-1
        comp[i+1, i] = 1
    end
    for i in 1:n
        comp[i, n] = -P[i] / P[end]
    end
    z = eigen(comp).values
    Qp(zj) = 1 - sum(g2[k] / (zj + im * Es[k])^2 for k in 1:M)
    return [sum(exp(zj * t) / Qp(zj) for zj in z) for t in ts]
end

"""
    fidelt_continuous(Γ, lb, ub, ts; nsub) -> Vector{ComplexF64}

连续谱 Γ(ε)（支撑 [lb,ub]）的 F(t)：细网格乘积梯形求解 Volterra。
内点数 = nsub·(最大 t)/h。Δ(t) 在细网格上一次算好（GL 谱积分）。
"""
function fidelt_continuous(Γ::Function, lb::Real, ub::Real,
                           ts::AbstractVector{<:Real}; nsub::Int=200, npts_spec::Int=400)
    Tmax = maximum(ts)
    Tmax > 0 || return [1.0 + 0im for _ in ts]
    # 细网格（含端点）：h = Tmax/(nsub·Tmax/h...) —— 直接取 nsub 点/单位……
    # 采用：每单位时间 nsub 个子步
    Nf = max(2, round(Int, nsub * Tmax))
    h = Tmax / Nf
    tg = range(0, Tmax; length=Nf + 1)
    # Δ(t) 在细网格（GL 谱积分）
    xs, ws = gauss_legendre(npts_spec)
    half = (ub - lb) / 2
    Es = [(lb + ub) / 2 + half * x for x in xs]
    wΓ = [half * w * Γ(ω) for (ω, w) in zip(Es, ws)]
    Δ = [sum(wΓ[k] * exp(-im * Es[k] * t) for k in eachindex(Es)) for t in tg]
    # 乘积梯形（Volterra 二阶的正确形式）：F' = −C，C(t) = ∫₀ᵗ Δ(t−τ)F(τ)dτ
    # F_{j+1} = F_j − (h/2)(C_j + C_{j+1})，C_{j+1} 的被积函数梯形含隐式 ½Δ₀F_{j+1}
    F = zeros(ComplexF64, Nf + 1)
    F[1] = 1
    Cprev = 0.0 + 0.0im                       # C(t_0) = 0
    for j in 0:Nf-1
        R = 0.5 * Δ[j+2] * F[1]               # τ=0 端点（j=0 时即唯一外端点）
        for i in 1:j
            R += Δ[j+2-i] * F[i+1]            # 内点（含 i=j，权 1）
        end
        F[j+2] = (F[j+1] - (h / 2) * Cprev - (h^2 / 2) * R) / (1 + h^2 * Δ[1] / 4)
        Cprev = h * (0.5 * Δ[1] * F[j+2] + R)
    end
    return [F[round(Int, t / h) + 1] for t in ts]
end

# Gauss–Legendre 节点（Golub–Welsch）
function gauss_legendre(n::Int)
    n == 1 && return [0.0], [2.0]
    β = [0.5 / sqrt(1 - 1 / (2 * k)^2) for k in 1:n-1]
    F = eigen(SymTridiagonal(zeros(n), β))
    return F.values, 2 .* (F.vectors[1, :] .^ 2)
end
