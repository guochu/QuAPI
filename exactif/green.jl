# =============================================================================
# exactif/green.jl —— 严格离散 IF 的接点-环路（Wick）表示
#
# 表示（本文件用多体矩阵实现并与 brute.jl 的 X-链张量互检；
# 单粒子 Wick 归约见 wick_if —— 收缩元全在 (F,A,B) 代数内）：
#
#   I[{s}] = Tr_{d,B}[ 𝒥ₛₜₐᵣₜ ⋅ (bra 链: Ŵ† 与接点) ⋯ 𝒥_turn ⋅ (ket 链: Ŵ 与接点) ⋯ ρ_B ]
#
# 接点字典（杂质模费米算符）：
#   ket 内点 m（时刻 mΔt）：𝒥ₘ = |αₘ⁺⟩⟨βₘ⁺| ∈ {1, 1−n̂, d̂†, d̂}
#   bra 内点 m：            𝒥'ₘ = |βₘ⁻⟩⟨αₘ⁻|
#   转折（t=NΔt）：         |β_N⁻⟩⟨β_N⁺|
#   起点（t=0，过 ρ_B）：    |α₀⁺⟩⟨α₀⁻|
# 全部为单粒子层级算符 ⟹ 整迹 = 高斯期望 ⟹ Wick 行列式。
# V=0 检验：环路退化为 δ-链 ✓（runtests L7）。
# =============================================================================
using LinearAlgebra: exp, tr, I, diagm

"""杂质模算符 |a⟩⟨b| 在 (d+bath) 联合基的矩阵（模 1 = 杂质 = 位 0）。"""
function dblock(a::Int, b::Int, M::Int)
    dim = 1 << (M + 1)
    O = zeros(ComplexF64, dim, dim)
    for jb in 0:(1 << M)-1          # 基：态 = jb*2 + imp（杂质 = 低位，同 mb_generator）
        O[2jb + a + 1, 2jb + b + 1] = 1
    end
    return O
end

"""
    loop_path_values(Es, Vs, β, N, Δt) -> Vector{ComplexF64}

接点-环路形式的过程张量（多体矩阵实现，与 if_full_tensor 互检用）。
顶点位同 (S39)：α_n⁺=4n+1, α_n⁻=4n+2, β_{n+1}⁺=4n+3, β_{n+1}⁻=4n+4。
"""
function loop_path_values(Es, Vs, β, N, Δt)
    M = length(Es)
    h = single_particle_h(Es, Vs)
    Wmb = exp(-im * Δt * mb_generator(h))
    nb = [1 / (1 + exp(β * Es[k])) for k in 1:M]
    ρB = diagm(0 => [prod((1 - nb[k], nb[k])[(j >> (k - 1)) & 1 + 1] for k in 1:M)
                     for j in 0:(1 << M)-1])
    # 联合基 ρ = 𝟙_d ⊗ ρ_B（未归一；I[∅] 归一由字典吸收）
    dimJ = 1 << (M + 1)
    ρ = kron(ρB, Matrix{ComplexF64}(I, 2, 2))
    nvert = 4N
    val(S, p) = (S >> (p - 1)) & 1
    out = Vector{ComplexF64}(undef, 1 << nvert)
    for S in 0:(1 << nvert)-1
        # 环路乘积 L（左=轮廓晚，右=轮廓早；Tr[L·ρ]，ρ 之后最早作用 W₁）：
        # 晚→早序列：𝒥_start, W†₁,𝒥'₁, W†₂,…,𝒥'_{N−1}, W†_N, 𝒥_turn,
        #            W_N, 𝒥_{N−1}, …, 𝒥₁, W₁
        L = Matrix{ComplexF64}(I, dimJ, dimJ)
        # 晚→早逐项右乘：𝒥_start, [W†·𝒥'_m]_{m=1..N−1}, W†(胞N), 𝒥_turn,
        #                [W·𝒥_{m}]_{m=N−1..1}, W(胞1)
        L = L * dblock(val(S, 1), val(S, 2), M)              # 𝒥_start=|α₀⁺⟩⟨α₀⁻|
        for m in 1:N-1
            L = L * Wmb' * dblock(val(S, 4m), val(S, 4m + 2), M)   # 𝒥'_m=|β_m⁻⟩⟨α_m⁻|
        end
        L = L * Wmb'
        L = L * dblock(val(S, 4N), val(S, 4N - 1), M)         # 𝒥_turn=|β_N⁻⟩⟨β_N⁺|
        for m in N-1:-1:1
            L = L * Wmb * dblock(val(S, 4m + 1), val(S, 4m - 1), M) # 𝒥_m=|α_m⁺⟩⟨β_m⁺|
        end
        L = L * Wmb
        out[S+1] = tr(L * ρ)
    end
    return out
end

# =============================================================================
# Wick 单粒子归约：I[{s}] = Σ_接点展开 ± Pf(收缩矩阵)
#
# 接点字典（费米算符，杂质模）：
#   |1⟩⟨1| = n̂ = d̂†d̂    |0⟩⟨0| = 1 − n̂    |1⟩⟨0| = d̂†    |0⟩⟨1| = d̂
# 高斯态 = (𝟙_d⊗ρ_B)/Tr：⟨f̂_i f̂_j†⟩ = δ_ij(1−n̄_j)（杂质位 n̄ = 1/2），
# 反常收缩 = 0（数守恒）⟹ 只有 (湮灭,产生) 配对非零。
#
# 轮廓（晚→早）：𝒥_start, [W†·𝒥'_m]_{m=N−1..1}, W†(胞N), 𝒥_turn,
#                [W·𝒥_m]_{m=N−1..1}, W(胞1)；共 2N 个元胞 + 2N+2 个接点。
# 缝隙 g（0..2N，早→晚编号）的单粒子演化 P_g = 从该缝隙到最晚端的元胞积。
# =============================================================================

"""
    wick_path_values(Es, Vs, β, N, Δt) -> Vector{ComplexF64}

Wick 归约的过程张量（单粒子代数，无多体矩阵）。与 loop/X-链机器精度互检。
"""
function wick_path_values(Es, Vs, β, N, Δt)
    M = length(Es)
    h = single_particle_h(Es, Vs)
    W1 = exp(-im * Δt * h)
    # 轮廓元胞（早→晚）：E₁..E_N = ket W，E_{N+1}..E_{2N} = bra W†
    E = vcat([W1 for _ in 1:N], [W1' for _ in 1:N])
    n̄ = [0.5; [1 / (1 + exp(β * Es[k])) for k in 1:M]]
    js = junction_list(N)                     # 晚→早接点序列
    # 每接点的"晚侧累计" T_j：从该接点向晚方向穿越的全部元胞积（晚左）
    # 接点 j（晚→早第 j 个）的晚侧元胞 = 轮廓早→晚的最后 (j−1) 个：E_{2N−j+2..2N}
    # T_j = E_{2N}···E_{2N−j+2}（j=1 ⟹ I）
    # 各接点的晚侧元胞数（js 晚→早 = [start, 𝒥'_{N−1..1}, turn, 𝒥_{N−1..1}]）：
    #   start: 0；𝒥'_m: m；turn: N；𝒥_m: N + (N−m)
    counts = Int[0]
    for m in 1:N-1; push!(counts, m); end
    push!(counts, N)
    for m in N-1:-1:1; push!(counts, N + (N - m)); end
    # T_j = 晚侧元胞积（最晚居左：E[2N]·E[2N−1]⋯）
    Ts = Vector{Matrix{ComplexF64}}(undef, length(js))
    for (j, cnt) in enumerate(counts)
        acc = Matrix{ComplexF64}(I, M + 1, M + 1)
        for cidx in 0:cnt-1
            acc = acc * E[2N - cidx]
        end
        Ts[j] = acc
    end
    val(S, p) = (S >> (p - 1)) & 1
    Trρ = 2.0
    out = Vector{ComplexF64}(undef, 1 << (4N))
    for S in 0:(1 << (4N)) - 1
        terms = [(1.0 + 0.0im, Tuple{Int,Int}[])]   # (系数, [(接点序 j, +1=C†/−1=A)])
        for (j, (pa, pb)) in enumerate(js)
            va, vb = val(S, pa), val(S, pb)
            nt = eltype(terms)[]
            if va == 1 && vb == 1
                for (c, ops) in terms; push!(nt, (c, vcat(ops, [(j, 1), (j, -1)]))); end
            elseif va == 0 && vb == 0
                for (c, ops) in terms
                    push!(nt, (c, ops))
                    push!(nt, (-c, vcat(ops, [(j, 1), (j, -1)])))
                end
            elseif va == 1 && vb == 0
                for (c, ops) in terms; push!(nt, (c, vcat(ops, [(j, 1)]))); end
            else
                for (c, ops) in terms; push!(nt, (c, vcat(ops, [(j, -1)]))); end
            end
            terms = nt
        end
        total = 0.0 + 0.0im
        for (c, ops) in terms
            isempty(ops) && (total += c; continue)
            k = length(ops)
            # 演化后向量：A(湮灭) = conj(T_j)·e₁；C(产生) = T_j 的第 1 行
            # 产生算符：T 的第 1 列；湮灭算符：conj(T 的第 1 列)
            vecs = [ops[i][2] == 1 ? collect(Ts[ops[i][1]][:, 1]) :
                                   collect(conj(Ts[ops[i][1]][:, 1])) for i in 1:k]
            MM = zeros(ComplexF64, k, k)
            for i in 1:k, j2 in i+1:k
                ti, tj = ops[i][2], ops[j2][2]
                v = if ti == 1 && tj == -1          # ⟨C A⟩ = Σ v_l u_l n̄_l
                    sum(vecs[i][l] * vecs[j2][l] * n̄[l] for l in 1:M+1)
                elseif ti == -1 && tj == 1          # ⟨A C⟩ = Σ u_l v_l (1−n̄_l)
                    sum(vecs[i][l] * vecs[j2][l] * (1 - n̄[l]) for l in 1:M+1)
                else
                    0.0im
                end
                MM[i, j2] = v
                MM[j2, i] = -v
            end
            total += c * pfaffian(MM)
        end
        out[S+1] = Trρ * total
    end
    return out
end

"""缝隙表：每个接点（按轮廓晚→早序列）→ (顶点位 ket 侧, bra 侧 or 末端)。"""
junction_list(N) = begin
    # 晚→早：𝒥_start(α₀⁺←α₀⁻), [𝒥'_m=β_m⁻←α_m⁻]_{m=1..N−1} 按晚→早即 m=N−1..1,
    #        𝒥_turn(β_N⁻←β_N⁺), [𝒥_m=α_m⁺←β_m⁺]_{m=N−1..1}
    js = Tuple{Int,Int}[]   # (bra/ket 侧顶点位 a, b)  接点=|a⟩⟨b|
    push!(js, (1, 2))                    # 𝒥_start = |α₀⁺⟩⟨α₀⁻|
    for m in 1:N-1                       # 晚→早：𝒥'_1, 𝒥'_2, …（𝒥'_m 在 braCell_m 与 braCell_{m+1} 之间）
        push!(js, (4m, 4m + 2))          # 𝒥'_m = |β_m⁻⟩⟨α_m⁻|
    end
    push!(js, (4N, 4N - 1))              # 𝒥_turn = |β_N⁻⟩⟨β_N⁺|
    for m in N-1:-1:1
        push!(js, (4m + 1, 4m - 1))      # 𝒥_m = |α_m⁺⟩⟨β_m⁺|
    end
    js
end

# =============================================================================
# 谱闭式 Wick（无 W 矩阵）：连续热库严格离散 IF 的实现
#
# 观察：Pfaffian 只依赖 Gram 矩阵
#   Gf[i,j] = ½·w_i[1]conj(w_j[1]) + Σ_k f_k·w_i[k+1]conj(w_j[k+1])
#   G1[i,j] = ½·w_i[1]conj(w_j[1]) + Σ_k (1−f_k)·w_i[k+1]conj(w_j[k+1])
# 其中 w_j = T_j e₁（T_j = 接点晚侧元胞积）。由于 W = U(Δt) 是 h 的函数，
# W^c = U(cΔt)，全部 w_j 分量有 F/φ/ψ 闭式（V 实规范，ψ 为 Duhamel 项）：
#   (a) c ≤ N：w = [conj(F_c); +iV_k conj(φ_c(E_k))]
#   (b) c = N+s：w = U(NΔt)†U(sΔt)e₀，分量经 (S9) 第四式展开（含 ψ_N）
# 连续热库：Σ_k → ∫dε Γ(ε)（谱求积，与 spectral.jl 同构）。
# =============================================================================

"""ψ_t(E_k,E_l) = ∫₀ᵗ (e^{−iE_kτ}−e^{−iE_lτ})/(E_k−E_l) F(t−τ) dτ（留数闭式）。"""
function psi_poles(Es, g2, z, Qp, t)
    M = length(Es)
    ψ = zeros(ComplexF64, M, M)
    for k in 1:M, l in 1:M
        s = 0.0 + 0.0im
        for j in eachindex(z)
            zj = z[j]
            if k != l
                inner = ((1 - exp(-(zj + im * Es[k]) * t)) / (zj + im * Es[k]) -
                         (1 - exp(-(zj + im * Es[l]) * t)) / (zj + im * Es[l])) / (Es[k] - Es[l])
            else
                # 对角项（数值锁定：带额外因子 −i）
                a = zj + im * Es[k]
                inner = -im * (1 - (1 + a * t) * exp(-a * t)) / a^2
            end
            s += exp(zj * t) / Qp(zj) * inner
        end
        ψ[k, l] = s
    end
    return ψ
end

"""
    wick_path_values_spec(Es, Vs, β, N, Δt; Vreal=true)

谱闭式 Wick 求值（不构造 W 矩阵；Vreal=true 时要求 Vs 实非负——
重相位规范总能做到）。与 wick_path_values 机器精度对照（runtests L8）。
"""
function wick_path_values_spec(Es, Vs, β, N, Δt; Vreal::Bool=true)
    M = length(Es)
    Vreal || error("复 V 请先做模式重相位规范（物理只依赖 |V_k|²）")
    # F/φ/ψ 闭式数据
    ts = [n * Δt for n in 1:N]
    Fs = fidelt_poles(Es, Vs, ts)
    φ = zeros(ComplexF64, N, M)
    g2 = abs2.(Vs)
    P, z, Qp = _pole_data_pub(Es, g2)
    for n in 1:N, k in 1:M
        s = 0.0im
        for j in eachindex(z)
            s += (exp(z[j] * n * Δt) - exp(-im * Es[k] * n * Δt)) /
                 (z[j] + im * Es[k]) / Qp(z[j])
        end
        φ[n, k] = s
    end
    ψN = psi_poles(Es, g2, z, Qp, N * Δt)
    f = [1 / (1 + exp(β * Es[k])) for k in 1:M]
    # 接点表（与 wick_path_values 一致）
    js = junction_list(N)
    counts = Int[0]
    for m in 1:N-1; push!(counts, m); end
    push!(counts, N)
    for m in N-1:-1:1; push!(counts, N + (N - m)); end
    nJ = length(js)
    # w_j 的闭式分量
    w0 = zeros(ComplexF64, nJ)              # 杂质位分量
    wV = zeros(ComplexF64, nJ, M)           # 浴分量（已含 V_k 因子）
    for (j, c) in enumerate(counts)
        if c <= N
            w0[j] = c == 0 ? 1.0 : conj(Fs[c])
            for k in 1:M
                wV[j, k] = c == 0 ? 0.0 : im * Vs[k] * conj(φ[c, k])
            end
        else
            s = c - N
            w0[j] = conj(Fs[N]) * (s == 0 ? 1.0 : Fs[s]) +
                    (s == 0 ? 0.0 : sum(g2[k] * conj(φ[N, k]) * φ[s, k] for k in 1:M))
            Fs_ = s == 0 ? 1.0 : Fs[s]
            φs = s == 0 ? zeros(M) : φ[s, :]
            for k in 1:M
                t1 = im * Vs[k] * conj(φ[N, k]) * Fs_
                t2 = -im * Vs[k] * exp(im * Es[k] * N * Δt) * φs[k]
                t3 = sum(Vs[k] * abs2(Vs[l]) * φs[l] * conj(ψN[l, k]) for l in 1:M)
                wV[j, k] = t1 + t2 + t3
            end
        end
    end
    # Gram 矩阵
    Gf = zeros(ComplexF64, nJ, nJ)
    G1 = zeros(ComplexF64, nJ, nJ)
    for i in 1:nJ, j2 in 1:nJ
        Gf[i, j2] = 0.5 * w0[i] * conj(w0[j2]) + sum(f[k] * wV[i, k] * conj(wV[j2, k]) for k in 1:M)
        G1[i, j2] = 0.5 * w0[i] * conj(w0[j2]) + sum((1 - f[k]) * wV[i, k] * conj(wV[j2, k]) for k in 1:M)
    end
    return wick_core_gram(js, Gf, G1, N)
end

"""Pfaffian 求值内核：由 Gram 矩阵 (Gf, G1) 与路径位串输出张量。"""
function wick_core_gram(js, Gf, G1, N)
    val(S, p) = (S >> (p - 1)) & 1
    out = Vector{ComplexF64}(undef, 1 << (4N))
    for S in 0:(1 << (4N)) - 1
        terms = [(1.0 + 0.0im, Tuple{Int,Int}[])]
        for (j, (pa, pb)) in enumerate(js)
            va, vb = val(S, pa), val(S, pb)
            nt = eltype(terms)[]
            if va == 1 && vb == 1
                for (c, ops) in terms; push!(nt, (c, vcat(ops, [(j, 1), (j, -1)]))); end
            elseif va == 0 && vb == 0
                for (c, ops) in terms
                    push!(nt, (c, ops))
                    push!(nt, (-c, vcat(ops, [(j, 1), (j, -1)])))
                end
            elseif va == 1 && vb == 0
                for (c, ops) in terms; push!(nt, (c, vcat(ops, [(j, 1)]))); end
            else
                for (c, ops) in terms; push!(nt, (c, vcat(ops, [(j, -1)]))); end
            end
            terms = nt
        end
        total = 0.0im
        for (c, ops) in terms
            isempty(ops) && (total += c; continue)
            k = length(ops)
            MM = zeros(ComplexF64, k, k)
            for i in 1:k, j2 in i+1:k
                (i1, t1), (i2, t2) = ops[i], ops[j2]
                v = if t1 == 1 && t2 == -1
                    Gf[i1, i2]
                elseif t1 == -1 && t2 == 1
                    G1[i2, i1]        # ⟨A_i C_j⟩ = conj(G1[i,j]) = G1[j,i]（Hermitian）
                else
                    0.0im
                end
                MM[i, j2] = v; MM[j2, i] = -v
            end
            total += c * pfaffian(MM)
        end
        out[S+1] = 2 * total
    end
    return out
end

# 公开极点数据（供 ψ/φ 使用）
function _pole_data_pub(Es, g2)
    M = length(Es)
    conv(p, q) = begin
        r = zeros(ComplexF64, length(p) + length(q) - 1)
        for i in eachindex(p), j2 in eachindex(q); r[i+j2-1] += p[i] * q[j2]; end
        r
    end
    lin = [[im * E, 1.0 + 0im] for E in Es]
    base = reduce(conv, lin; init=[1.0 + 0im])
    P = zeros(ComplexF64, M + 2)
    for i in eachindex(base); P[i+1] += base[i]; end
    for k in 1:M
        others = reduce(conv, lin[1:end .!= k]; init=[1.0 + 0im])
        for i in eachindex(others); P[i] += g2[k] * others[i]; end
    end
    n = M + 1
    comp = zeros(ComplexF64, n, n)
    for i in 1:n-1; comp[i+1, i] = 1; end
    for i in 1:n; comp[i, n] = -P[i] / P[n+1]; end
    z = eigen(comp).values
    Qp(zj) = 1 - sum(g2[k] / (zj + im * Es[k])^2 for k in 1:M)
    return P, z, Qp
end

"""
    wick_path_values_cont(Γ, lb, ub, β, N, Δt; nsub, npts_spec)

连续热库严格离散 IF：F(t) 用 Volterra 细网格，φ_n(ε)、ψ_N(ε,ε')
用时间求积，全部 Gram 为 ε-节点谱求积。V 实规范自动成立（Γ=|V|²ρ）。
"""
function wick_path_values_cont(Γ::Function, lb::Real, ub::Real, β, N, Δt;
                               nsub::Int=400, npts_spec::Int=400)
    M = npts_spec
    xs, ws = gauss_legendre(M)
    half = (ub - lb) / 2
    Es = [(lb + ub) / 2 + half * x for x in xs]
    wΓ = [half * w * Γ(ε) for (ε, w) in zip(Es, ws)]
    # F 细网格（Volterra）
    Tmax = N * Δt
    Nf = max(2, round(Int, nsub * Tmax))
    h = Tmax / Nf
    tg = range(0, Tmax; length=Nf + 1)
    Δ = [sum(wΓ[q] * exp(-im * Es[q] * t) for q in 1:M) for t in tg]
    F = zeros(ComplexF64, Nf + 1)
    F[1] = 1
    Cp = 0.0im
    for j in 0:Nf-1
        R = 0.5 * Δ[2] * F[j+1]                  # i = j 项（权 ½），与 volterra.jl 一致
        for i in 0:j-1
            w = (i == 0) ? 0.5 : 1.0
            R += w * Δ[j+2-i] * F[i+1]
        end
        F[j+2] = (F[j+1] - (h / 2) * Cp - (h^2 / 2) * R) / (1 + h^2 * Δ[1] / 4)
        Cp = h * (0.5 * Δ[1] * F[j+2] + R)
    end
    Fs = [F[round(Int, n * Δt / h) + 1] for n in 1:N]
    # φ_n(ε_q) = ∫₀^{nΔt} e^{−iε(T−τ)}F(τ)dτ（梯形）
    φ = zeros(ComplexF64, N, M)
    for n in 1:N, q in 1:M
        np = round(Int, n * Δt / h)
        T = n * Δt
        s = 0.0im
        for i in 0:np
            w = (i == 0 || i == np) ? 0.5 : 1.0
            s += w * exp(-im * Es[q] * (T - tg[i+1])) * F[i+1]
        end
        φ[n, q] = h * s
    end
    # ψ_N(ε_q, ε_l)（时间求积；对角带 −i 因子——数值锁定）
    ψN = zeros(ComplexF64, M, M)
    tN = N * Δt
    for q in 1:M, l in 1:M
        s = 0.0im
        for i in 0:Nf
            τ = tg[i+1]
            w = (i == 0 || i == Nf) ? 0.5 : 1.0
            g = q == l ? -im * τ * exp(-im * Es[q] * τ) :
                (exp(-im * Es[q] * τ) - exp(-im * Es[l] * τ)) / (Es[q] - Es[l])
            s += w * g * F[Nf - i + 1]
        end
        ψN[q, l] = h * s
    end
    f = [1 / (1 + exp(β * ε)) for ε in Es]
    js = junction_list(N)
    counts = Int[0]
    for m in 1:N-1; push!(counts, m); end
    push!(counts, N)
    for m in N-1:-1:1; push!(counts, N + (N - m)); end
    nJ = length(js)
    w0 = zeros(ComplexF64, nJ)
    wV = zeros(ComplexF64, nJ, M)
    Vsq = sqrt.(wΓ)
    for (j, c) in enumerate(counts)
        if c <= N
            w0[j] = c == 0 ? 1.0 : conj(Fs[c])
            for k in 1:M
                wV[j, k] = c == 0 ? 0.0 : im * Vsq[k] * conj(φ[c, k])
            end
        else
            s = c - N
            w0[j] = conj(Fs[N]) * (s == 0 ? 1.0 : Fs[s]) +
                    (s == 0 ? 0.0 : sum(wΓ[k] * conj(φ[N, k]) * φ[s, k] for k in 1:M))
            Fs_ = s == 0 ? 1.0 : Fs[s]
            φs = s == 0 ? zeros(M) : φ[s, :]
            for k in 1:M
                t1 = im * Vsq[k] * conj(φ[N, k]) * Fs_
                t2 = -im * Vsq[k] * exp(im * Es[k] * tN) * φs[k]
                t3 = sum(Vsq[k] * wΓ[l] * φs[l] * conj(ψN[l, k]) for l in 1:M)
                wV[j, k] = t1 + t2 + t3
            end
        end
    end
    Gf = zeros(ComplexF64, nJ, nJ)
    G1 = zeros(ComplexF64, nJ, nJ)
    for i in 1:nJ, j2 in 1:nJ
        Gf[i, j2] = 0.5 * w0[i] * conj(w0[j2]) +
                    sum(f[k] * wV[i, k] * conj(wV[j2, k]) for k in 1:M)
        G1[i, j2] = 0.5 * w0[i] * conj(w0[j2]) +
                    sum((1 - f[k]) * wV[i, k] * conj(wV[j2, k]) for k in 1:M)
    end
    return wick_core_gram(js, Gf, G1, N)
end
