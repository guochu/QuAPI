# =============================================================================
# firstorder/bosons/mixedtime.jl —— 玻色混合轮廓：一阶冻结核 vs 严格 1/2 阶
# Trotter。对照：严格解（全平衡 G>）。输出 results/mixedtime.csv
#
#   · 一阶（E1 冻结核装配，face-value）：τ 支用正确环解作热背景（n_d、φτ
#     ——环是热背景输入），实时通道 = QuAPI 一阶装配 A = (A0⁻¹+K)⁻¹·源
#     （A0⁻¹ = I−uS⁻，K = face-value 因果下三角含接触对角；无半步 f、
#     无梯形），ζ 跨块朴素收缩（行 +1、G 节点值权重、无 ×f）；
#   · 全 K 侧改善（梯形 φ、半步相位 f、接触对角、跳变完整代入、跨块终止
#     cell + ×f）= 直接二阶，见 secondorder/bosons/mixedtime.jl（firstorder
#     不重复）。注：仅顶点移位单独抽出与其他朴素部分组合不自洽（实测误差
#     反而增大），故混合轮廓的一阶→二阶只有"打包改善"一条路径。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bosons/mixedtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "common", "common.jl"))

"""
    frozen_mixed_naive(mdl, Nτ, Nt, T) -> G>[t_0..t_Nt]

玻色混合一阶冻结核（朴素 E1 装配）。`G[m+1] ≈ G>(t_m)`。
"""
function frozen_mixed_naive(mdl, Nτ, Nt, T)
    h = T / Nt
    δτ = mdl.β / Nτ
    u = exp(-im * mdl.ϵd * h)
    fτ = exp(-mdl.ϵd * δτ / 2)
    dm = bosonic_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=Nτ, t=T, Nt=Nt, μ=mdl.μ)
    P = branch(dm; b1=:+, b2=:+)
    Xr = branch(dm; b1=:+, b2=:-)
    Xc = branch(dm; b1=:+, b2=:τ)
    d0 = -(P[1, 1] + Xr[1, 1])                  # 接触对角（冻结核本体条目）
    # ① τ 支热背景：正确环解（bosonic_Δτ 驱动——热背景是输入而非被改善对象）
    corr = bosonic_Δτ(mdl.bath.spectrum; β=mdl.β, N=Nτ, μ=mdl.μ)
    ηjk, ηkj = corr.data.ηⱼₖ, corr.data.ηₖⱼ
    d0τ = ηjk[1] + ηkj[1]
    M = zeros(Float64, Nτ, Nτ)
    b = zeros(Float64, Nτ)
    for m in 1:Nτ
        M[m, m] += 1.0
        if m ≥ 2
            M[m, m-1] -= exp(-mdl.ϵd * δτ)
        else
            b[m] += exp(-mdl.ϵd * δτ)
            M[m, Nτ] -= exp(-mdl.ϵd * δτ)
        end
        for k in 1:Nτ
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0τ
            if k == 1
                M[m, Nτ] -= fτ * w * 0.5
                M[m, 1] -= fτ * w * 0.5
                b[m] += fτ * w * 0.5
            else
                M[m, k-1] -= fτ * w * 0.5
                M[m, k] -= fτ * w * 0.5
            end
        end
    end
    Gτ = M \ b
    n_d = Gτ[Nτ]
    # ② 实时支 A/C：E1 冻结核装配 A = (A0⁻¹+K)⁻¹·源（隐式、稳定）
    #    A0⁻¹ = I − u·S⁻（full-step 自由飞行）；K = face-value 因果下三角
    #    K[m,k] = −(P+X)[m−k+1]（lag m−k 直读、含对角接触 [1,1]）
    A0inv = Matrix{ComplexF64}(I, Nt + 1, Nt + 1)
    for j in 2:Nt+1; A0inv[j, j-1] = -u; end
    K = zeros(ComplexF64, Nt + 1, Nt + 1)
    for m in 1:Nt+1, k in 1:m
        K[m, k] = -(P[m-k+1, 1] + Xr[m-k+1, 1])
    end
    A = (A0inv + K) \ [1.0; zeros(ComplexF64, Nt)]          # A 通道（初值 a₀=1）
    # C 通道：同骨架 + ζ 跨块源（朴素：行 +1、G 节点值权重、无 ×f）
    Svec = zeros(ComplexF64, Nt + 1)
    for m in 1:Nt
        Svec[m+1] = sum(Xc[m+1, j] * Gτ[j] for j in 1:Nτ)
    end
    C = (A0inv + K) \ Svec
    return (-im) .* ((1 + n_d) .* A .+ C)
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[16, 32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    GL = exact_greater_lesser(mdl)
    γ = exact_gamma(mdl)
    println(io, "=" ^ 74)
    println(io, "玻色混合：一阶冻结核 vs 严格 1/2 阶 Trotter（全平衡 G>(t)）")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f（Nτ=Nt 同步）\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "  全 K 侧二阶改善见 secondorder/bosons/mixedtime.jl")
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64}[]
    println(io, "\n[1] 一阶冻结核 vs 精确 G>(t)")
    es = Float64[]
    for N in Ns
        G = frozen_mixed_naive(mdl, N, N, T)
        ms = ms_of(N)
        e = maximum(abs.([G[m+1] - GL(m * T / N)[1] for m in ms]))
        push!(es, e)
        push!(rows, (N, T / N, e, NaN, NaN))
    end
    order_table(io, "    一阶冻结核（朴素 E1 装配）", Ns, es)
    println(io, "\n[2] 严格 Trotter 参照（γ 取 exact_gamma）")
    for order in (1, 2)
        es = Float64[]
        for (i, N) in enumerate(Ns)
            h = T / N
            gt, _, _ = chain_greater_lesser(mdl, order, N, T; γ=γ, ms=ms_of(N))
            e = maximum(abs.(gt .- [GL(m * T / N)[1] for m in ms_of(N)]))
            push!(es, e)
            rows[i] = (rows[i][1], rows[i][2], rows[i][3],
                       order == 1 ? e : rows[i][4], order == 2 ? e : rows[i][5])
        end
        order_table(io, "    strict order=$order", Ns, es)
    end
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "mixedtime.csv"), "w") do f
        println(f, "N,dt,err_frozen1,err_strict1,err_strict2")
        for (N, h, a, c, d) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e,%.6e\n", N, h, a, c, d)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "mixedtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
