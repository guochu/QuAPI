# =============================================================================
# firstorder/bcs/realtime.jl —— BCS 浴实时：一阶冻结核（Nambu 矩阵）naive vs
# 顶点移位 vs 严格 1/2 阶 Trotter。对照：严格 BdG e^{−i𝓗t} 杂质 2×2 块。
# 输出 results/realtime.csv
#
# 模型与约定与 secondorder/bcs 完全一致（common.jl 复用）：Nambu 基
# Φ = (d↑, c↑…; d̄↓†, c̄↓†…)，杂质 𝒜(t) = e^{−i𝓗t} 的 2×2 块；推迟核
# 𝒩 = [cn cc; nn nc]（正常块 X−P、反常块 P−X，lag 直读），
# 行相位（u, uI）= (e^{−iϵ_d h}, e^{+iϵ_d h})。
#
#   · 一阶（naive，face-value，无半步/无接触）：𝒜_{m+1} = P(𝒜_m − Σ_{d≥1} 𝒦(d)𝒜_{m−d})；
#   · 顶点移位（仍一阶，常数↓）：接触 𝒢₀ 并入最近邻步 P(I−𝒢₀)，节点 lag d 配
#     条目 𝒢(d−1)——与 firstorder/fermions/realtime.jl 同款规则的 Nambu 矩阵化，
#     写成块线性系统 M·𝒜 = B（可直接稠密折叠）；
#   · Δ=0 退化为两个解耦正常通道（已验证）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bcs/realtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "bcs", "common.jl"))

"BCS 推迟 Nambu 核条目（bcs_Δt 四块 → 𝒩 = [cn cc; nn nc]，lag 直读）。"
function bcs_kernel_lag(mdl, N, T)
    corr = bcs_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ, Δ=mdl.Δ)
    P1, X1 = branch(corr.uu, :+, :+), branch(corr.uu, :+, :-)
    P2, X2 = branch(corr.dd, :+, :+), branch(corr.dd, :+, :-)
    P3, X3 = branch(corr.ud, :+, :+), branch(corr.ud, :+, :-)
    P4, X4 = branch(corr.du, :+, :+), branch(corr.du, :+, :-)
    return d -> [-(X1[d+1, 1] + P1[d+1, 1])     P3[d+1, 1] + X3[d+1, 1];
                 P4[d+1, 1] + X4[d+1, 1]     -(X2[d+1, 1] + P2[d+1, 1])]
end

"杂质行相位（无半步）：P = diag(e^{−iϵ_d h}, e^{+iϵ_d h})。"
row_phase(mdl, h) = Diagonal([exp(-im * mdl.ϵd * h), exp(im * mdl.ϵd * h)])

"""
    bcs_naive_system(mdl, N, T) -> (M, B)

一阶 naive 的块线性系统（未知 𝒜_1..𝒜_N，每个 2×2 按 2 行堆叠；
B = 初条件源 e₁⊗I₂）：𝒜_n = P(𝒜_{n−1} − Σ_{d=1}^{n−2} 𝒦(d)𝒜_{n−1−d})。
"""
function bcs_naive_system(mdl, N, T)
    h = T / N
    P = row_phase(mdl, h)
    K = bcs_kernel_lag(mdl, N, T)
    M = zeros(ComplexF64, 2N, 2N)
    B = zeros(ComplexF64, 2N, 2)
    M[1:2, 1:2] = Matrix{ComplexF64}(I, 2, 2)
    B[1:2, :] = Matrix{ComplexF64}(I, 2, 2)
    for n in 2:N
        M[2n-1:2n, 2n-1:2n] = Matrix{ComplexF64}(I, 2, 2)
        M[2n-1:2n, 2n-3:2n-2] = -P
        for d in 1:n-2
            M[2n-1:2n, 2(n-1-d)-1:2(n-1-d)] += P * K(d)
        end
    end
    return M, B
end

"""
    bcs_shifted_system(mdl, N, T) -> (M, B)

顶点移位的块线性系统：𝒜_m = P(I−𝒢₀)𝒜_{m−1} − P·Σ_{d≥2} 𝒢(d−1)𝒜_{m−d}，
𝒜_0 = I。接触 𝒢₀ = 𝒦(0) 并入最近邻步；越界 lag 落进 B。
"""
function bcs_shifted_system(mdl, N, T)
    h = T / N
    P = row_phase(mdl, h)
    K = bcs_kernel_lag(mdl, N, T)
    G0 = K(0)
    M = zeros(ComplexF64, 2N, 2N)
    B = zeros(ComplexF64, 2N, 2)
    M[1:2, 1:2] = Matrix{ComplexF64}(I, 2, 2)
    B[1:2, :] = P * (I - G0)
    for m in 2:N
        M[2m-1:2m, 2m-1:2m] = Matrix{ComplexF64}(I, 2, 2)
        M[2m-1:2m, 2m-3:2m-2] = -P * (I - G0)
        for d in 2:m-1
            M[2m-1:2m, 2(m-d)-1:2(m-d)] += P * K(d - 1)
        end
        B[2m-1:2m, :] = P * K(m - 1)          # d = m → 𝒜_0 = I
    end
    return M, B
end

"""
    strict_chain_bcs(mdl, N, T, order) -> (𝒜[1..N], h)

严格 cell 链参照：BdG 𝓗 = hA + hB（hA = 杂质 Nambu 块 diag(ϵ_d,−ϵ_d)），
order=1 Lie（eA·eB），order=2 Strang（eA·eB·eA）。
"""
function strict_chain_bcs(mdl, N, T, order::Int)
    h = T / N
    𝓗 = bcs_bdg(mdl)
    n = mdl.L + 1
    hA = zeros(ComplexF64, 2n, 2n)
    hA[1, 1] = mdl.ϵd
    hA[n+1, n+1] = -mdl.ϵd
    hB = 𝓗 - hA
    eA = exp(-im * h / 2 * hA)
    eB = exp(-im * h * hB)
    st = order == 1 ? exp(-im * h * hA) * eB : eA * eB * eA
    idx = [1, n+1]
    return [(st^m)[idx, idx] for m in 1:N], h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Δ=0.5, T=2.0, Ns=[32, 64, 128, 256, 512],
              fracs=[0.25, 0.5, 0.75])
    mdl = bcs_toulouse(; Lmodes=Lmodes, Δ=Δ)
    ex = exact_retarded_Nambu(mdl)          # −i𝒢^R_𝓝 = e^{−i𝓗t} 杂质块
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "BCS 浴实时：一阶冻结核（Nambu 矩阵）naive vs 顶点移位 vs 严格参照")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，Δ=%.2f，T=%.1f\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64,Float64}[]
    err_of(𝒜s, N) = maximum(maximum(abs.([𝒜s[m][p, q] - (ex(m * T / N) ./ (-im))[p, q]
                                          for m in ms_of(N), p in 1:2, q in 1:2])))
    println(io, "\n[1] 一阶冻结核 vs 精确 e^{−i𝓗t} 杂质块")
    es1 = Float64[]; es2 = Float64[]
    for N in Ns
        M, B = bcs_naive_system(mdl, N, T)
        st = M \ B
        naive = [st[2m-1:2m, :] for m in 1:N]
        M2, B2 = bcs_shifted_system(mdl, N, T)
        st2 = M2 \ B2
        shifted = [st2[2m-1:2m, :] for m in 1:N]
        push!(es1, err_of(naive, N))
        push!(es2, err_of(shifted, N))
        push!(rows, (N, T / N, es1[end], es2[end], NaN, NaN))
    end
    order_table(io, "    一阶冻结核（naive，face-value）", Ns, es1)
    order_table(io, "    顶点移位改善（仍一阶，常数↓）", Ns, es2)
    println(io, "\n[2] 严格 Trotter 参照（BdG 分裂，杂质 Nambu 块）")
    for order in (1, 2)
        es = Float64[]
        for (i, N) in enumerate(Ns)
            𝒜s, _ = strict_chain_bcs(mdl, N, T, order)
            e = err_of(𝒜s, N)
            push!(es, e)
            rows[i] = (rows[i][1], rows[i][2], rows[i][3], rows[i][4],
                       order == 1 ? e : rows[i][5], order == 2 ? e : rows[i][6])
        end
        order_table(io, "    strict order=$order", Ns, es)
    end
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "realtime.csv"), "w") do fh
        println(fh, "N,dt,err_frozen1,err_shifted,err_strict1,err_strict2")
        for (N, h, a, b, c, d) in rows
            @printf(fh, "%d,%.6e,%.6e,%.6e,%.6e,%.6e\n", N, h, a, b, c, d)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "realtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
