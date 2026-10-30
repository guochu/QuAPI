# =============================================================================
# firstorder/bcs/imagtime.jl —— BCS 浴虚时：一阶冻结核 Nambu Matsubara 环
# vs 严格 BdG 环 vs 严格 Trotter 1/2。输出 results/imagtime.csv
#
# 模型与约定与 secondorder/bcs 完全一致（common.jl 复用）。观测对象：杂质
# Nambu 2×2 Matsubara 块 Ĝ(τ_m) = −[e^{−𝓗τ}M₀]_{imp,imp}（M₀ = ⟨ΦΦ†⟩ =
# (1+e^{−β𝓗})⁻¹），行 1 = [greater 𝒢₁, 反常 𝒢₂]（exact_Gτ_ring 同款）。
#
# 环方程（4 分量/节点：(i,j) = (1,1),(1,2),(2,1),(2,2)，行 i = 杂质 Nambu 行）：
#   自由相位（行 i）：uτ_i = e^{∓ϵ_d δτ}，fτ_i = e^{∓ϵ_d δτ/2}（行 2 = 空穴行）；
#   核 𝚺̂ 的 lag 权重（naive 矩形，同正常费米环）：k<m 用 ηⱼₖ lag、k>m 用
#   ηₖⱼ lag、对角 ηⱼₖ[1]+ηₖⱼ[1]；Nambu 槽位
#       𝚺̂ = [cn  −cc ; −nn  nc]        （反常槽位带 −1：IF 约定把负号吸进
#       Δ，Dyson 重排时反常自能多出 −1——符号由 exact_Gτ_ring 裁判扫描锁定）
#   卷积按行混排：[𝚺⊛Ĝ]_{ij} = Σ_k′ 𝚺̂_{i,k′} Ĝ_{k′,j}（行 2 分量参与行 1
#   方程——这是 2 分量切片不封闭的根源）；
#   跳变源只在 Nambu 对角分量：b[(i,i)] = −uτ_i（{Φ_i,Φ_i†}=1；反常
#   (1,2)/(2,1) 无等点反对易子，源为零）。
#   提取：Ĝ_11 = 𝒢₁、Ĝ_12 = 𝒢₂（行 1）。
#
#   Δ=0 退化：𝒢₂ ≡ 0 精确成立（已验证）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bcs/imagtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "bcs", "common.jl"))

"""
    bcs_tau_ring_naive(mdl, N) -> (M, b)

BCS Nambu Matsubara 一阶冻结核环（naive 矩形）。未知堆叠 s = 4(m−1)+idx(i,j)。
"""
function bcs_tau_ring_naive(mdl, N)
    δτ = mdl.β / N
    uτ = [exp(-mdl.ϵd * δτ), exp(mdl.ϵd * δτ)]
    fτ = [exp(-mdl.ϵd * δτ / 2), exp(mdl.ϵd * δτ / 2)]
    corr = bcs_Δτ(mdl.bath.spectrum; β=mdl.β, N=N, μ=mdl.μ, Δ=mdl.Δ)
    ηj(X) = X.data.ηⱼₖ
    ηk(X) = X.data.ηₖⱼ
    NJ(l) = [ηj(corr.uu)[l]  -ηj(corr.ud)[l];
             -ηj(corr.du)[l]  ηj(corr.nc)[l]]
    NK(l) = [ηk(corr.uu)[l]  -ηk(corr.ud)[l];
             -ηk(corr.du)[l]  ηk(corr.nc)[l]]
    Nd0 = NJ(1) + NK(1)
    M = zeros(ComplexF64, 4N, 4N)
    b = zeros(ComplexF64, 4N)
    idx(i, j) = 2(i - 1) + j
    for m in 1:N, i in 1:2, j in 1:2
        s = 4 * (m - 1) + idx(i, j)
        M[s, s] += 1.0
        if m ≥ 2
            M[s, 4*(m-2)+idx(i, j)] -= uτ[i]
        else
            M[s, 4*(N-1)+idx(i, j)] += uτ[i]           # 反周期回绕（per Nambu 行）
            i == j && (b[s] = -uτ[i])                  # 跳变源只在 Nambu 对角
        end
        for k in 1:N
            W = k < m ? NJ(m-k+1) : k > m ? NK(k-m+1) : Nd0
            for ip in 1:2                              # [Σ⊛Ĝ]_{ij} = Σ_ip Σ̂_{i,ip} Ĝ_{ip,j}
                M[s, 4*(k-1)+idx(ip, j)] -= fτ[i] * W[i, ip]
            end
        end
    end
    return M, b
end

"""
    strict_chain_bcs_tau(mdl, N, order) -> (Ĝ[1..N], h)

严格 Trotter 虚时链：T = e^{−h hA}e^{−h hB}（Lie）/ eA·eB·eA（Strang），
处置矩阵 P̃ = (I+Tᴺ)⁻¹，Ĝ(τ_m) = −[Tᵐ P̃]_{imp,imp}。
"""
function strict_chain_bcs_tau(mdl, N, order::Int)
    h = mdl.β / N
    𝓗 = bcs_bdg(mdl)
    n = mdl.L + 1
    hA = zeros(ComplexF64, 2n, 2n)
    hA[1, 1] = mdl.ϵd
    hA[n+1, n+1] = -mdl.ϵd
    hB = 𝓗 - hA
    eA = exp(-h / 2 * hA)
    eB = exp(-h * hB)
    T = order == 1 ? exp(-h * hA) * eB : eA * eB * eA
    P = inv(Matrix{ComplexF64}(I, 2n, 2n) + T^N)
    idx = [1, n+1]
    out = Vector{Matrix{ComplexF64}}(undef, N)
    Tm = Matrix{ComplexF64}(I, 2n, 2n)
    for m in 1:N
        Tm = Tm * T
        out[m] = -(Tm * P)[idx, idx]
    end
    return out, h
end

"严格 BdG 环：Ĝ(τ_m) = −[e^{−𝓗τ_m}M₀]_{imp,imp}（2×2 块，机器精度参照）。"
function exact_tau_block(mdl, N)
    V, λ, Vinv = bcs_eig(mdl)
    n = mdl.L + 1
    M0 = bcs_M0(mdl)
    h = mdl.β / N
    idx = [1, n+1]
    return [-(V * (Diagonal(exp.(-λ .* (m * h))) * (Vinv * M0)))[idx, idx] for m in 1:N], h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Δ=0.5, Ns=[8, 16, 32, 64, 128], fracs=[0.25, 0.5, 0.75])
    mdl = bcs_toulouse(; Lmodes=Lmodes, Δ=Δ)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "BCS 浴虚时：一阶冻结核 Nambu Matsubara 环 vs 严格 BdG 环 vs 严格 Trotter")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，Δ=%.2f\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64}[]

    println(io, "\n[0] Δ=0 退化：𝒢₂ ≡ 0（反常分量）")
    mdl0 = bcs_toulouse(; Lmodes=Lmodes, Δ=0.0)
    for N in (4, 8)
        M0, b0 = bcs_tau_ring_naive(mdl0, N)
        G0 = M0 \ b0
        @printf(io, "  N=%2d: max|Ĝ12| = %.1e\n", N,
                maximum(abs.(G0[4*(m-1)+2] for m in 1:N)))
    end

    println(io, "\n[1] 一阶冻结核 Nambu 环 vs 严格 BdG 环（杂质 2×2 块全元、全节点）")
    es = Float64[]
    for N in Ns
        M, b = bcs_tau_ring_naive(mdl, N)
        G = M \ b
        ex, h = exact_tau_block(mdl, N)
        e = 0.0
        for m in 1:N
            e = max(e, maximum(abs.([G[4*(m-1)+2*(i-1)+j] - ex[m][i, j]
                                     for i in 1:2, j in 1:2])))
        end
        push!(es, e)
        push!(rows, (N, h, e, NaN, NaN))
    end
    order_table(io, "    一阶冻结核 Nambu 环（naive 矩形）", Ns, es)

    println(io, "\n[2] 严格 Trotter 参照（处置矩阵 P̃ = (I+Tᴺ)⁻¹）")
    for order in (1, 2)
        es = Float64[]
        for (i, N) in enumerate(Ns)
            Gs, h = strict_chain_bcs_tau(mdl, N, order)
            ex, _ = exact_tau_block(mdl, N)
            e = maximum(maximum(abs.([Gs[m][i2, j2] - ex[m][i2, j2]
                                      for i2 in 1:2, j2 in 1:2])) for m in ms_of(N))
            push!(es, e)
            rows[i] = (rows[i][1], rows[i][2], rows[i][3],
                       order == 1 ? e : rows[i][4], order == 2 ? e : rows[i][5])
        end
        order_table(io, "    strict order=$order", Ns, es)
    end

    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "imagtime.csv"), "w") do fh
        println(fh, "N,dtau,err_frozen1,err_strict1,err_strict2")
        for (N, h, a, c, d) in rows
            @printf(fh, "%d,%.6e,%.6e,%.6e,%.6e\n", N, h, a, c, d)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "imagtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
