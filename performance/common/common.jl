# =============================================================================
# benchmark/bosons/common.jl
#
# 无相互作用玻色子 Toulouse 杂质模型的 QuAPI 离散化精度基准——公共部分。
#
# 模型：H = ϵ_d d†d + Σ_k ω_k b†_k b_k + Σ_k V_k (d†b_k + b_k†d)
#       （ImpurityModelBase 的 Toulouse 模型；全线可严格对角化）
#
# 严格解（ImpurityModelBase）：
#   · 虚时    toulouse_Gτ(Toulouse(bath; ϵ_d))                → freebosons_Gτ
#   · 实/混合 toulouse_greater_lesser(Toulouse(bath; ϵ_d))     → freebosons_greater_lesser
#
# 演化器（高效：单粒子（Gaussian）表象，无 Fock 截断，O(L³) / 步）：
#   把 H 拆成  A = 杂质自由部分 (ϵ_d d†d)  与  B = 浴自由部分 + 耦合 (Σω b†b + H_int)
#   1 阶：  e^{-τA} e^{-τB}
#   2 阶：  e^{-τA/2} e^{-τB} e^{-τA/2}                （Strang）
#   4 阶：  Suzuki S₄(τ) = S₂(θτ)² S₂(γτ) S₂(θτ)²
#   每个因子都是 (L+1)×(L+1) 单粒子矩阵，链 = 矩阵连乘；观测量由单粒子
#   Green 函数公式给出（见下）。
#
# QuAPI 杂化矩阵：bosonic_Δτ / bosonic_Δt / bosonic_Δm。本文件用独立的
# “模式求和”重建同一离散核并逐元对照，确认 bosonic_Δ* 给出的正是严格的
# 双 cell 积分（即“离散 hybridization matrix”本身无误差）。
# =============================================================================
using QuAPI, ImpurityModelBase, LinearAlgebra, Printf

# ------------------------------------------------------------------ #
# 模型构造
# ------------------------------------------------------------------ #
"""
    boson_toulouse(; β, ϵd, Γ, wlo, whi, Lmodes, μ)

在 [wlo, whi] 上取平坦谱 J(ω) = Γ（gapped，避开 ω=0 的 1/ε 奇点），离散成
`Lmodes` 个浴模。返回 NamedTuple，含离散浴、严格模型、单粒子矩阵等。
"""
function boson_toulouse(; β=4.0, ϵd=1.0, Γ=0.15, wlo=0.5, whi=3.0, Lmodes=48, μ=0.0)
    J = bounded(w -> Γ, wlo, whi)
    freqs = collect(range(wlo, whi; length=Lmodes + 1))
    bath = discretebath(Boson, freqs, J; β=β, μ=μ)
    model = Toulouse(bath; ϵ_d=ϵd)
    ws, αs = frequencies(bath), spectrumvalues(bath)
    L = length(ws)
    hS = zeros(L + 1, L + 1); hS[1, 1] = ϵd
    hb = zeros(L + 1, L + 1)                 # 浴自由 + 耦合（单粒子）
    hi = zeros(L + 1, L + 1)
    for k in 1:L
        hb[1 + k, 1 + k] = ws[k]
        hi[1, 1 + k] = sqrt(αs[k])
        hi[1 + k, 1] = sqrt(αs[k])
    end
    (; β=β, ϵd=ϵd, Γ=Γ, μ=μ, bath, model, ws, αs, L, hS, hb, hi,
       hfull=hS + hb + hi, Id=Matrix{Float64}(I, L + 1, L + 1))
end

"严格虚时 Green 函数 G(τ)（ImpurityModelBase，解析 Matsubara 求和 / ED）。"
exact_Gτ(mdl) = toulouse_Gτ(mdl.model)

"""
    exact_greater_lesser(mdl) -> (G>, G<)

严格实时 greater/lesser（热平衡）。约定与 freebosons_greater_lesser 一致：
G>(t) = -i⟨d(t)d†(0)⟩,  G<(t) = +i⟨d†(0)d(t)⟩。
"""
exact_greater_lesser(mdl) = toulouse_greater_lesser(mdl.model)

# ------------------------------------------------------------------ #
# Trotter 单粒子步
# ------------------------------------------------------------------ #
const SUZUKI_θ = 1 / (4 - 4^(1 / 3))
const SUZUKI_γ = 1 - 4 * SUZUKI_θ

"虚时单粒子步矩阵 e^{-h·(…)}(order)。"
function imag_step(mdl, order::Int, h::Real)
    eA(t) = exp(-t * mdl.hS)
    eB(t) = exp(-t * (mdl.hb + mdl.hi))
    s2(t) = eA(t / 2) * eB(t) * eA(t / 2)
    order == 1 && return eA(h) * eB(h)
    order == 2 && return s2(h)
    order == 4 && return (s2(SUZUKI_θ * h)^2) * s2(SUZUKI_γ * h) * (s2(SUZUKI_θ * h)^2)
    throw(ArgumentError("order 必须为 1/2/4"))
end

"实时单粒子步矩阵 e^{-i h·(…)}(order)。"
function real_step(mdl, order::Int, h::Real)
    eA(t) = exp(-im * t * mdl.hS)
    eB(t) = exp(-im * t * (mdl.hb + mdl.hi))
    s2(t) = eA(t / 2) * eB(t) * eA(t / 2)
    order == 1 && return eA(h) * eB(h)
    order == 2 && return s2(h)
    order == 4 && return (s2(SUZUKI_θ * h)^2) * s2(SUZUKI_γ * h) * (s2(SUZUKI_θ * h)^2)
    throw(ArgumentError("order 必须为 1/2/4"))
end

# ------------------------------------------------------------------ #
# 观测量（单粒子 Gaussian 公式）
#   Gτ(τ)  = [st^{Nτ-m} (I - U_β)^{-1}]_{11},  U_β = st^{Nτ}
#   G>(t)  = -i [R(t) (I + γᵀ)]_{11},          γ = (U_β^{-1} - I)^{-1}
#   G<(t)  = +i [R(t) γᵀ]_{11}
# ------------------------------------------------------------------ #
"虚时链：返回 (Gτ 在 τ=m·h, m=ms) 以及 h。"
function chain_Gτ(mdl, order::Int, N::Int; ms=1:N-1)
    h = mdl.β / N
    st = imag_step(mdl, order, h)
    U = st^N
    S = inv(mdl.Id - U)                 # = I + γ
    return [ (st^m * S)[1, 1] for m in ms ], ms, h
end

"虚时链给出的 γ = (U_β^{-1}-I)^{-1}（用于实时/混合轮廓的平衡初态）。"
function chain_gamma(mdl, order::Int, N::Int)
    h = mdl.β / N
    U = imag_step(mdl, order, h)^N
    inv(inv(U) - mdl.Id)
end

"严格热平衡单粒子密度矩阵 γ_ij = ⟨d_i†d_j⟩（由 H 的单粒子矩阵对角化）。"
function exact_gamma(mdl)
    ev = eigen(Hermitian(mdl.hfull))
    nk = [1 / (exp(mdl.β * max(λ, 1e-14)) - 1) for λ in ev.values]
    ev.vectors * Diagonal(nk) * ev.vectors'
end

"""
    chain_greater_lesser(mdl, order, Nt, T; γ, ms)

实时链（总时长 T，Nt 步）。返回 (G>, G<) 在 t=m·T/Nt，m=ms。
"""
function chain_greater_lesser(mdl, order::Int, Nt::Int, T::Real;
                              γ=chain_gamma(mdl, 2, 256), ms=1:Nt-1)
    h = T / Nt
    st = real_step(mdl, order, h)
    idim = size(mdl.Id, 1)
    γt = transpose(γ)
    gt = Vector{ComplexF64}(undef, length(ms))
    lt = Vector{ComplexF64}(undef, length(ms))
    for (idx, m) in enumerate(ms)
        R = st^m
        gt[idx] = -im * (R * (mdl.Id + γt))[1, 1]
        lt[idx] = im * (R * γt)[1, 1]
    end
    return gt, lt, ms .* h
end

# ------------------------------------------------------------------ #
# QuAPI 杂化矩阵：模式求和重建（用于逐元验证 bosonic_Δτ/Δt/Δm）
#   约定：cell i = [τ_i, τ_{i+1}]（左端节点 τ_i=(i-1)h），核取双 cell 积分。
# ------------------------------------------------------------------ #
nb(w, β) = 1 / (exp(β * w) - 1)

"虚时 ηⱼₖ / ηₖⱼ（与 bosonic_Δτ 同序）"
function delta_tau_modesum(ws, αs, β, N)
    h = β / N
    ηjk = zeros(N); ηkj = zeros(N)
    for (w, α) in zip(ws, αs)
        n = nb(w, β)
        f2 = (2 * sinh(w * h / 2) / w)^2          # |∫_cell e^{-wτ}|²
        # 不等时：ηjk[d] 滞后 (d-1)h；ηkj[d] 取周期回绕（等价 n·e^{-w(β-lag)}）
        for d in 2:N
            lag = (d - 1) * h
            ηjk[d] += α * (1 + n) * exp(-w * lag) * f2
            # 反向核：e^{-(β-lag)w}（周期回绕，等价 n·e^{+lag·w} 的另一种写法）
            ηkj[d] += α * (1 + n) * exp(-w * (β - lag)) * f2
        end
        # 等时：两种序各取半胞三角积分
        ηjk[1] += α * (1 + n) * (exp(-w * h) - 1 + w * h) / w^2
        ηkj[1] += α * n * (exp(w * h) - 1 - w * h) / w^2
    end
    return ηjk, ηkj
end

"""
    verify_delta_tau(mdl, N) -> (maxdiff, ηjk_q, ηkj_q)

对照 bosonic_Δτ 与模式求和重建值，返回最大偏差。
"""
function verify_delta_tau(mdl, N)
    corr = bosonic_Δτ(mdl.bath.spectrum, β=mdl.β, N=N, μ=mdl.μ)
    ηjk, ηkj = delta_tau_modesum(mdl.ws, mdl.αs, mdl.β, N)
    d = max(maximum(abs.(ηjk .- corr.data.ηⱼₖ)), maximum(abs.(ηkj .- corr.data.ηₖⱼ)))
    return d, corr.data.ηⱼₖ, corr.data.ηₖⱼ
end

"实时四个 Keldysh 块（与 bosonic_Δt 同序、同符号）的模式求和重建。"
function delta_t_modesum(ws, αs, β, N, T)
    h = T / N
    n = N + 1
    blocks = Dict{String,Tuple{Vector{ComplexF64},Vector{ComplexF64}}}()
    for key in ("++", "+-", "-+", "--")
        blocks[key] = (zeros(ComplexF64, n), zeros(ComplexF64, n))
    end
    for (w, α) in zip(ws, αs)
        b = nb(w, β); bb = b + 1
        # 不等时核：e^{-iωΔk h}·2(1-cos ωh)/ω²（= |∫_cell e^{iωτ}|² 带相位）
        c2 = 2 * (1 - cos(w * h)) / w^2
        # 等时核
        gm = (1 - im * w * h - exp(-im * w * h)) / w^2    # _gⱼⱼ_r
        gp = (1 + im * w * h - exp(im * w * h)) / w^2     # _gₖₖ_r
        fjj = bb * gm
        fkk = b * gp
        for d in 1:N
            lag = d * h
            ker = exp(-im * w * lag) * c2
            fjk = bb * ker          # _bosonic_fⱼₖ_r
            fkj = b * ker           # _bosonic_fₖⱼ_r
            # G₊₊
            blocks["++"][1][1 + d] += -α * fjk
            blocks["++"][2][1 + d] += -α * (b * exp(-im * w * (-d * h)) * c2)
            # G₊₋
            blocks["+-"][1][1 + d] += α * fkj
            blocks["-+"][1][1 + d] += α * fjk
            # G₋₋
            blocks["--"][1][1 + d] += -α * fkj
            blocks["--"][2][1 + d] += -α * (bb * exp(-im * w * (-d * h)) * c2)
        end
        # 等时项
        blocks["++"][1][1] += -α * fjj
        blocks["++"][2][1] += -α * fkk
        blocks["+-"][1][1] += α * b * gm
        blocks["+-"][2][1] += α * fkk
        blocks["-+"][1][1] += α * fjj
        blocks["-+"][2][1] += α * bb * gp
        blocks["--"][1][1] += -α * conj(fkk)
        blocks["--"][2][1] += -α * conj(fjj)
    end
    # 共轭对称的交叉项（bosonic_Δt 中 ηₖⱼ[i+1] = ηⱼₖ[i+1]'）
    for key in ("+-", "-+")
        jk, kj = blocks[key]
        for d in 1:N
            kj[1 + d] = conj(jk[1 + d])
        end
    end
    return blocks
end

"对照 bosonic_Δt 的四个块，返回最大偏差。"
function verify_delta_t(mdl, N, T)
    corr = bosonic_Δt(mdl.bath.spectrum, β=mdl.β, N=N, t=T, μ=mdl.μ)
    ref = delta_t_modesum(mdl.ws, mdl.αs, mdl.β, N, T)
    d = 0.0
    for (key, (b1, b2)) in (("++", (:+, :+)), ("+-", (:+, :-)), ("-+", (:-, :+)), ("--", (:-, :-)))
        q = branch(corr; b1=b1, b2=b2)
        jk, kj = ref[key]
        d = max(d, maximum(abs.(jk .- q.ηⱼₖ)), maximum(abs.(kj .- q.ηₖⱼ)))
    end
    return d
end

"""
    verify_delta_m(mdl, Nτ, Nt, T) -> (dτ, d±)

混合轮廓 Δm 的一致性：其 τ 块应等于 Δτ，其 ± 块应等于 Δt（同一网格）。
"""
function verify_delta_m(mdl, Nτ, Nt, T)
    m = bosonic_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=Nτ, t=T, Nt=Nt, μ=mdl.μ)
    τq = bosonic_Δτ(mdl.bath.spectrum, β=mdl.β, N=Nτ, μ=mdl.μ)
    rq = bosonic_Δt(mdl.bath.spectrum, β=mdl.β, N=Nt, t=T, μ=mdl.μ)
    matof(x) = [ComplexF64(x[i, j]) for i in 1:size(x, 1), j in 1:size(x, 1)]
    dτ = maximum(abs.(matof(branch(m; b1=:τ, b2=:τ)) .- matof(τq.data)))
    dpm = maximum(maximum(abs.(matof(branch(m; b1=b1, b2=b2)) .-
                              matof(branch(rq; b1=b1, b2=b2))))
                  for (b1, b2) in ((:+, :+), (:+, :-), (:-, :+), (:-, :-)))
    return dτ, dpm
end

# ------------------------------------------------------------------ #
# 冻结顶点（vertex-frozen）QuAPI 路径积分——玻色 Gaussian 精确求逆
#   G = (G0^{-1} - Δ̄)^{-1}，G0 为节点自由杂化 Green 函数
#
# 等时项（[d, d†] = 1，非对角杂化）的正确处理（推导见下）：
#   算符表示：ηⱼₖ[1]·(d†d) + ηₖⱼ[1]·(d d†)          ← 两种算符次序，必须都保留
#            = (ηⱼₖ[1] + ηₖⱼ[1])·n̂ + ηₖⱼ[1]·[d, d†]
#            = (ηⱼₖ[1] + ηₖⱼ[1])·n̂ + ηₖⱼ[1] .
#   末项 [d,d†] = 1 是 c-number：它整体乘在 IF 上，在任何归一化可观测量
#   （G = 分子/分母，cell 数相同）中严格抵消。因此等时项对路径的系数恰为
#   ηⱼₖ[1] + ηₖⱼ[1] —— 这不是“忽略不对易”，不对易关系正是那个被丢掉的常数。
#   丢单侧（只用 ηⱼₖ[1] 或 ηₖⱼ[1]）会把 n̂ 的系数改掉 ηₖⱼ[1]，是实打实的错误。
#
#   若耦合算符的不对易子不是 c-number（如 σ⁺: [σ⁺,σ⁻] = σ_z），则该项是算符，
#   必须在 PT 里按 M1/M2 分开保留（BosonicHybridizationMatrix split 模式）。
#   对角/厄米耦合下两种次序的算符积相同，直接相加即可。
#
# 该规则把 cell 双积分核配节点处冻结的路径值，是 1 阶规则（与等时项无关）。
# ------------------------------------------------------------------ #
function frozen_dyson_Gτ(mdl, N; cellrule::Int=1)
    h = mdl.β / N
    A = zeros(N, N)
    for j in 1:N
        A[mod(j, N) + 1, j] = 1.0
    end
    G0 = inv(Matrix{Float64}(I, N, N) - exp(-h * mdl.ϵd) * A)
    corr = bosonic_Δτ(mdl.bath.spectrum, β=mdl.β, N=N, μ=mdl.μ)
    eqt = corr.data.ηⱼₖ[1] + corr.data.ηₖⱼ[1]     # 两种次序的系数之和（常数项已约去）
    D = zeros(N, N)
    for i in 1:N, j in 1:N
        if i > j
            D[i, j] = corr.data.ηⱼₖ[i - j + 1]
        elseif i < j
            D[i, j] = corr.data.ηₖⱼ[j - i + 1]
        else
            D[i, j] = eqt
        end
    end
    T = cellrule_T(N, cellrule)
    G = inv(inv(G0) - transpose(T) * D * T)
    return [G[m + 1, 1] for m in 1:N-1], h
end

"""
    cellrule_T(N, p)

cell 内系统路径的重构矩阵：`ẑ_i = Σ_j T[i,j] z_j`（cell i = [τ_i, τ_{i+1}]）。
p=1 节点值（矩形），p=2 中点（线性），p=4 三次插值后取 cell 平均。
用于检验“把杂交项的重构阶数提高”能否改善冻结顶点规则的收敛阶。
"""
function cellrule_T(N, p)
    T = zeros(N, N)
    if p == 1
        for i in 1:N
            T[i, i] = 1.0
        end
    elseif p == 2
        for i in 1:N
            T[i, i] = 0.5
            T[i, mod(i, N) + 1] = 0.5
        end
    elseif p == 4
        w = (-1 / 24, 13 / 24, 13 / 24, -1 / 24)
        for i in 1:N, (k, off) in zip(1:4, (-1, 0, 1, 2))
            T[i, mod(i - 1 + off, N) + 1] = w[k]
        end
    else
        throw(ArgumentError("cellrule 必须为 1/2/4"))
    end
    return T
end

"""
    equal_time_identity_check(ηjk1, ηkj1; nmax=30)

数值验证等时项算符恒等式（截断 Fock 空间）：
    exp(ηjk1·n̂) exp(ηkj1·(n̂+1)) = e^{ηkj1}·exp((ηjk1+ηkj1)·n̂) .
即“两种次序分开保留”与“系数相加”只差一个整体常数 e^{ηkj1}，归一化量不变。
返回 (max|Osplit − e^{ηkj1}Osum|, 归一化期望之比 / e^{ηkj1} − 1)。
"""
function equal_time_identity_check(ηjk1, ηkj1; nmax=30)
    n = Matrix(Diagonal(0.0:nmax))
    Iden = Matrix{Float64}(I, nmax + 1, nmax + 1)
    Osplit = exp(ηjk1 .* n) * exp(ηkj1 .* (n + Iden))
    Osum = exp((ηjk1 + ηkj1) .* n)
    d = maximum(abs.(Osplit .- exp(ηkj1) .* Osum))
    w = [exp(-0.7 * k) for k in 0:nmax]
    ratio = sum(w .* diag(Osplit)) / sum(w .* diag(Osum))
    return d, ratio / exp(ηkj1) - 1
end

# ------------------------------------------------------------------ #
# 实时（实轴、直积初态）的冻结核求解器
#
# 直积初态 ρ₀ = ρ_S(n_sys) ⊗ ρ_bath(β) 下，精确关系为
#     G>(t) = -i(1+n_sys) A(t)，  G<(t) = +i n_sys A(t)，
# 其中 A(t) = [e^{-iHt}]₁₁ 是杂质回归振幅（H 为杂质+浴的单粒子矩阵）。
# A(t) 满足推迟单粒子 Dyson 方程，冻结核离散化为
#     (A0^{-1} - K) A = 1,   A0^{-1} = I - e^{-iϵh}S⁻,
#     K[j,k] = -Σ_α α ∫_{cell j}∫_{cell k} e^{-iω(t-t')} (t>t'，对角用三角形)，
# 即用 bosonic_Δt 的 cell 双积分核（≈ -Δ^{++}）驱动。复杂度 O(N³)（三角求解）。
# ------------------------------------------------------------------ #
function gauss01(n)
    # Golub–Welsch：Jacobi 矩阵非对角元 β_k = 1/√(4k²−1)（此前误写 0.5/√(4−1/k²)，
    # 节点全错、权重碰巧正确——已修正并对 ∫u¹⁰、∫cos 3u 自检）
    β = [1 / sqrt(4 * k^2 - 1) for k in 1:n-1]
    ev = eigen(SymTridiagonal(zeros(n), β))
    return (ev.values .+ 1) ./ 2, ev.vectors[1, :] .^ 2
end
const _GX, _GW = gauss01(6)

"cell 对 ∫∫ e^{-iω(t-t')}；j==k 时只取 t>t' 三角形（Duffy）。"
function _cellint_rt(w, dj, h)
    s = 0.0 + 0.0im
    if dj == 0
        for (x, wx) in zip(_GX, _GW), (y, wy) in zip(_GX, _GW)
            s += wx * wy * h^2 * x * exp(-im * w * h * x * (1 - y))
        end
    else
        for (x, wx) in zip(_GX, _GW), (y, wy) in zip(_GX, _GW)
            s += wx * wy * h^2 * exp(-im * w * (dj * h + h * x - h * y))
        end
    end
    s
end

"精确回归振幅 A(t) = [e^{-iHt}]₁₁（单粒子对角化）。"
function exact_return_amplitude(mdl)
    ev = eigen(Hermitian(mdl.hfull))
    return t -> sum(ev.vectors[1, k]^2 * exp(-im * ev.values[k] * t)
                    for k in 1:length(ev.values))
end

"""
    frozen_retarded_A(mdl, N, T) -> (A_t, h)

实轴冻结核求解器：**直接由 `bosonic_Δt` 的块驱动**，解推迟单粒子 Dyson 方程

    (A0⁻¹ + K^R) A = 1,    A0⁻¹ = I − e^{-iϵ_d h}S⁻,
    K^R[j,k] = −(Δ⁺⁺ + Δ⁺₋)[j,k]   (t_j ≥ t_k 的因果下三角，含等时),

其中 Δ⁺⁺/Δ⁺₋ 为 `bosonic_Δt` 的 CorrelationMatrix 块（直接 getindex，含等时
ηⱼₖ[1]+ηₖⱼ[1] 约定）。该组合是**精确恒等式**：lag τ>0 时 Δ⁺⁺ = −(1+n)R、
Δ⁺₋ = +nR（R = Σα e^{-iωτ}∫∫cell），−(Δ⁺⁺+Δ⁺₋) = R，热因子 (1+n)−n 严格
相消——正是 retarded 核与态无关的定义性质；等时对角 (−(1+n)tri−n·tri′) +
(n·tri+n·tri′) = −tri 同样严格成立（三角形 t>t′ 约定）。返回 A(t_j)，j=0..N。
"""
function frozen_retarded_A(mdl, N, T)
    h = T / N
    corr = bosonic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+)     # Δ⁺⁺
    X = branch(corr; b1=:+, b2=:-)     # Δ⁺₋
    A0inv = Matrix{ComplexF64}(I, N + 1, N + 1)
    u = exp(-im * mdl.ϵd * h)
    for j in 2:N+1
        A0inv[j, j - 1] = -u
    end
    KR = zeros(ComplexF64, N + 1, N + 1)
    for j in 1:N+1, k in 1:j           # 因果：只填下三角（含对角）
        KR[j, k] = -(P[j, k] + X[j, k])
    end
    A = inv(A0inv + KR)
    return [A[j, 1] for j in 1:N+1], h
end

"""
    chain_return_amplitude(mdl, order, N, T) -> (A_t, h)

严格 cell 的实时链给出的 A(t)（用于与冻结核对照）。
"""
function chain_return_amplitude(mdl, order::Int, N::Int, T::Real)
    h = T / N
    st = real_step(mdl, order, h)
    return [ (st^m)[1, 1] for m in 0:N ], h
end

# ------------------------------------------------------------------ #
# 误差 / 收敛阶
# ------------------------------------------------------------------ #
"对 (err at Ns) 拟合逐级比值与等效阶 log2(err_i/err_{i+1})。"
function order_table(io, name, Ns, errs)
    println(io, name)
    prev = 0.0
    for (N, e) in zip(Ns, errs)
        if prev == 0.0
            @printf(io, "   N=%4d   err=%.4e\n", N, e)
        else
            @printf(io, "   N=%4d   err=%.4e   ratio=%.2f   order=%.2f\n",
                    N, e, prev / e, log2(prev / e))
        end
        prev = e
    end
end
