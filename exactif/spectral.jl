# =============================================================================
# exactif/spectral.jl —— 谱积分层：φ, D, E, K 与小矩阵 A、B（严格核的代数核心）
#
# 输入：离散模式（Es, |Vs|²）或连续谱 Γ(ε)（支撑 [lb,ub]，GL 求积）。
#
# 量（n = 1..N）：
#   φ_n(ε) = ∫₀^{nΔt} e^{−iε(nΔt−τ)} F(τ) dτ        （W 的行/列向量生成函数）
#   D[m,n] = Σ_k |V_k|² φ_m(E_k) conj(φ_n(E_k))      （连续：∫Γ φ_m φ̄_n）
#   E[m,n] = Σ_k |V_k|² f(E_k) φ_m(E_k) conj(φ_n(E_k))（连续：∫Γ f φ_m φ̄_n）
#   K：Krylov 递推矩阵 C 的行向量组：C[m,1] = −F(mΔt)，C[m,m+1] = 1，
#      K[i,:] = (Cⁱ)₁,· —— {aₘ} 基下行向量 r⁽ⁱ⁾ = a·U_Bⁱ 的展开系数。
#   A = K D K†，B = K E K† —— 论文式 (S31)/(S33) 的 Gram 矩阵
#      （⟨v⁽ⁱ⁾|v⁽ʲ⁾⟩ 与 ⟨v⁽ⁱ⁾|(1+e^{βh})⁻¹|v⁽ʲ⁾⟩）。
#   Z_N/Z = det(I − B)（(S34)）。
#
# 验证（runtests.jl）：离散模式下 A、B、Z_N/Z 与 brute.jl 的显式矩阵
# 计算（r⁽ⁱ⁾·v⁽ʲ⁾、r⁽ⁱ⁾·f·v⁽ʲ⁾）逐位机器精度一致。
# =============================================================================
using LinearAlgebra: I, det, eigen

# ---------------------- F(t) 与 φ_n（离散模式：闭式） ---------------------- #

"""
    spectral_discrete(Es, Vs, β, N, Δt) -> NamedTuple

离散模式的严格核代数层（全部闭式、机器精度）：
F(nΔt) 与 φ_n(E_k) 用留数闭式；D、E、K、A、B、ZNoZ 直接组装。
"""
function spectral_discrete(Es, Vs, β, N, Δt)
    M = length(Es)
    ts = [n * Δt for n in 1:N]
    Fs = fidelt_poles(Es, Vs, ts)                 # F(nΔt)，闭式
    # 留数参数（供 φ 闭式）
    g2 = abs2.(Vs)
    P, z, Qp = _pole_data(Es, g2)
    # φ_n(E_k) = Σ_j [e^{z_j nΔt} − e^{−iE_k nΔt}]/(z_j + iE_k)/Q'(z_j)
    φ = zeros(ComplexF64, N, M)
    for n in 1:N, k in 1:M
        s = zero(ComplexF64)
        for j in eachindex(z)
            s += (exp(z[j] * n * Δt) - exp(-im * Es[k] * n * Δt)) /
                 (z[j] + im * Es[k]) / Qp(z[j])
        end
        φ[n, k] = s
    end
    f = [1 / (1 + exp(β * Es[k])) for k in 1:M]
    g1 = [abs2(Vs[k]) for k in 1:M]
    D = [sum(g1[k] * φ[m, k] * conj(φ[n, k]) for k in 1:M) for m in 1:N, n in 1:N]
    E = [sum(g1[k] * f[k] * φ[m, k] * conj(φ[n, k]) for k in 1:M) for m in 1:N, n in 1:N]
    K = _krylov_K(Fs, N)
    A = K * D * K'
    B = K * E * K'
    return (; Fs, φ, D, E, K, A, B, ZNoZ=det(I - B), f)
end

function _pole_data(Es, g2)
    M = length(Es)
    convpoly(p, q) = begin
        r = zeros(ComplexF64, length(p) + length(q) - 1)
        for i in eachindex(p), j in eachindex(q)
            r[i+j-1] += p[i] * q[j]
        end
        r
    end
    lin = [[im * E, 1.0 + 0im] for E in Es]   # 升幂 (iE + z)
    base = reduce((p, q) -> convpoly(p, q), lin; init=[1.0 + 0im])
    P = zeros(ComplexF64, M + 2)
    for i in eachindex(base)
        P[i+1] += base[i]
    end
    for k in 1:M
        others = reduce((p, q) -> convpoly(p, q), lin[1:end .!= k]; init=[1.0 + 0im])
        for i in eachindex(others)
            P[i] += g2[k] * others[i]
        end
    end
    n = length(P) - 1
    comp = zeros(ComplexF64, n, n)
    for i in 1:n-1
        comp[i+1, i] = 1
    end
    for i in 1:n
        comp[i, n] = -P[i] / P[end]
    end
    z = (comp isa Matrix{ComplexF64} ? eigen(comp) : eigen(comp)).values
    Qp(zj) = 1 - sum(g2[k] / (zj + im * Es[k])^2 for k in 1:M)
    return P, z, Qp
end

# ---------------------- Krylov 递推矩阵 K ---------------------- #

"""
    _krylov_K(Fs, N) -> N×N

K 的第 i 行 = e₁ᵀCⁱ⁻¹（i = 1..N），C[m,1] = −F(mΔt)，C[m,m+1] = 1。
行向量组 {r⁽ⁱ⁾ = a·U_B^{i−1}} 在基 {aₘ} 下的展开：
  r⁽ⁱ⁾ = Σ_n K[i,n]·a_n，aₘ,k = −iV_k φ_m(E_k)。
"""
function _krylov_K(Fs::AbstractVector, N::Int)
    C = zeros(ComplexF64, N, N)
    for m in 1:N
        C[m, 1] = -Fs[m]
        m < N && (C[m, m+1] = 1)
    end
    K = zeros(ComplexF64, N, N)
    row = zeros(ComplexF64, 1, N)      # 1×N 行向量
    row[1] = 1
    K[1, :] .= vec(row)
    for i in 2:N
        row = row * C
        K[i, :] .= vec(row)
    end
    return K
end

# ---------------------- 连续谱版本 ---------------------- #

"""
    spectral_continuous(Γ, lb, ub, β, N, Δt; nsub, npts_spec, npts_phi)

连续谱 Γ(ε)：F 用细网格 Volterra（volterra.jl），φ_n 在同一细网格上
梯形积分，D、E 用 GL 谱积分。收敛性由 nsub/npts_* 控制（runtests 检验）。
"""
function spectral_continuous(Γ::Function, lb::Real, ub::Real, β, N, Δt;
                             nsub::Int=200, npts_spec::Int=400)
    Ts = [n * Δt for n in 1:N]
    Tmax = maximum(Ts)
    # 细网格 F
    Nf = max(2, round(Int, nsub * Tmax))
    h = Tmax / Nf
    tg = range(0, Tmax; length=Nf + 1)
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
    Fs = [F[round(Int, T / h) + 1] for T in Ts]
    # φ_n(ε_q) 在细网格上梯形：∫₀^{nΔt} e^{−iε(T−τ)}F(τ)dτ
    φ = zeros(ComplexF64, N, length(Es))
    for n in 1:N
        np = round(Int, Ts[n] / h)                # 细网格点数
        for q in eachindex(Es)
            ε = Es[q]
            s = zero(ComplexF64)
            for i in 0:np
                w = (i == 0 || i == np) ? 0.5 : 1.0
                s += w * exp(-im * ε * (Ts[n] - tg[i+1])) * F[i+1]
            end
            φ[n, q] = h * s
        end
    end
    f = [1 / (1 + exp(β * ε)) for ε in Es]
    D = [sum(wΓ[q] * φ[m, q] * conj(φ[n, q]) for q in eachindex(Es))
         for m in 1:N, n in 1:N]
    E = [sum(wΓ[q] * f[q] * φ[m, q] * conj(φ[n, q]) for q in eachindex(Es))
         for m in 1:N, n in 1:N]
    K = _krylov_K(Fs, N)
    A = K * D * K'
    B = K * E * K'
    return (; Fs, φ, D, E, K, A, B, ZNoZ=det(I - B), f, Es, wΓ)
end
