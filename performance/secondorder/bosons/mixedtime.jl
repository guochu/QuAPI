# =============================================================================
# secondorder/bosons/mixedtime.jl —— U1GTEMPO2 "冻结核不变 + K 侧二阶"的混合
# 时间复现（L 形轮廓：τ 环热制备 + 实时演化，玻色 Toulouse 平衡 G>(t)），
# 并输出画图数据 results/mixedtime.csv（Nτ=Nt 同步，dt=dτ）
#
# 精确单粒子分解：G>(t)/(−i) = (1+n_d)·A(t) + cross(t)，三个源通道共享同一
# 实时隐式梯形系统（realtime.jl 骨架）：
#   · 通道 1（retarded）：A_m —— 结点初始振幅，bosonic_Δm 的 ++/+- 块推迟组合，
#     梯形 φ + 半步相位 f + 接触对角（同实时三要素）；
#   · 通道 2（占据）：n_d = G^τ(β⁻)（ξ 块环端点，玻色 KMS）；
#   · 通道 3（初关联）：rhs_m = f·Σ_j ζ⁺ᵗ[m,j]·φ^τ_j —— (+:τ) 跨块 × τ 环解
#     梯形值（含双侧跳变），×吸收半步相位 f，块行取终止 cell（块为起始 cell
#     约定，行偏移 −1）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/bosons/mixedtime.jl
# =============================================================================
include(joinpath(@__DIR__, "imagtime.jl"))    # 传递引入 common.jl（τ 支改由 Δm ξ 块驱动，不再调 weyl_k_imagtime）

"""
    weyl_k_mixed(mdl, Nτ, Nt, T) -> G>[t_0..t_Nt]

混合轮廓 K 侧二阶。`G[m+1] ≈ G>(t_m)`。
"""
function weyl_k_mixed(mdl, Nτ, Nt, T)
    h = T / Nt
    δτ = mdl.β / Nτ
    u = exp(-im * mdl.ϵd * h)          # 实时支相位
    f = exp(-im * mdl.ϵd * h / 2)
    uτ = exp(-mdl.ϵd * δτ)             # τ 支相位（实衰减，虚时无 i）
    fτ = exp(-mdl.ϵd * δτ / 2)
    mix = bosonic_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=Nτ, t=T, Nt=Nt, μ=mdl.μ)
    P = branch(mix; b1=:+, b2=:+)
    Xr = branch(mix; b1=:+, b2=:-)
    Xc = branch(mix; b1=:+, b2=:τ)              # (Nt+1)×Nτ 跨块（起始 cell 约定）
    d0 = -(P[1, 1] + Xr[1, 1])
    # ① τ 支环（ξ 块 = bosonic_Δτ 同数据，逐位验证；四要素：梯形 φ、半步相位
    #    fτ、接触对角、跳变 G₀ = G_N+1 完整代入——自由项 b += u 与 φ₁ 左半 b += fτ·w/2）
    ξ = branch(mix; b1=:τ, b2=:τ)                # branch(:τ,:τ) = CM(−ξⱼₖ,−ξₖⱼ) = Δτ 条目
    ηjk, ηkj = ξ.ηⱼₖ, ξ.ηₖⱼ
    d0τ = ηjk[1] + ηkj[1]
    M = zeros(Float64, Nτ, Nτ)
    b = zeros(Float64, Nτ)
    for m in 1:Nτ
        M[m, m] += 1.0
        if m ≥ 2
            M[m, m-1] -= uτ
        else
            b[m] += uτ                          # −u·G₀ = −u·(G_N+1)：常数 +u
            M[m, Nτ] -= uτ
        end
        for k in 1:Nτ
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0τ
            if k == 1                           # φ₁ 左半 = (G_N+1)/2
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
    Gτ0 = Gτ[Nτ] + 1
    φτ = [((j == 1 ? Gτ0 : Gτ[j-1]) + Gτ[j]) / 2 for j in 1:Nτ]
    # 实时双通道
    A = zeros(ComplexF64, Nt + 1); A[1] = 1.0
    C = zeros(ComplexF64, Nt + 1)
    for m in 1:Nt
        s1 = 0.0 + 0.0im
        s3 = 0.0 + 0.0im
        for k in 1:m-1
            Δ = -(P[m+1, k+1] + Xr[m+1, k+1])
            s1 += Δ * (A[k] + A[k+1]) / 2
            s3 += Δ * (C[k] + C[k+1]) / 2
        end
        rhs = (1 ≤ m ≤ Nt + 1) ? f * sum(Xc[m, j] * φτ[j] for j in 1:Nτ) : 0.0im
        A[m+1] = (u*A[m] - f*s1 - f*d0*A[m]/2) / (1 + f*d0/2)
        C[m+1] = (u*C[m] - f*s3 - f*d0*C[m]/2 + rhs) / (1 + f*d0/2)
    end
    return (-im) .* ((1 + n_d) .* A .+ C)
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[16, 32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    GL = exact_greater_lesser(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    gexact(N, m) = GL(m * T / N)[1]
    println(io, "=" ^ 74)
    println(io, "混合 Weyl-K 二阶（冻结核 + K 侧）vs 严格 cell 链：玻色 Toulouse G>(t)")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f（Nτ=Nt 同步）\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64}[]
    println(io, "\n[1] Weyl-K 二阶 vs 精确 G>(t)")
    es = Float64[]
    for N in Ns
        G = weyl_k_mixed(mdl, N, N, T)
        e = maximum(abs.([G[m+1] - gexact(N, m) for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, T / N, e, NaN, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照（γ 取 exact_gamma）")
    for order in (1, 2)
        es = Float64[]
        for (i, N) in enumerate(Ns)
            gt, _, _ = chain_greater_lesser(mdl, order, N, T; γ=exact_gamma(mdl), ms=ms_of(N))
            e = maximum(abs.(gt .- [gexact(N, m) for m in ms_of(N)]))
            push!(es, e)
            rows[i] = (rows[i][1], rows[i][2], rows[i][3],
                       order == 1 ? e : rows[i][4], order == 2 ? e : rows[i][5])
        end
        order_table(io, "    strict order=$order", Ns, es)
    end
    # CSV 输出（画图用）
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "mixedtime.csv"), "w") do f
        println(f, "N,dt,err_weylk,err_strict1,err_strict2")
        for (N, h, e1, e2, e3) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e,%.6e\n", N, h, e1, e2, e3)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "mixedtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
