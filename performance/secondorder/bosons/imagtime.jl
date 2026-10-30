# =============================================================================
# secondorder/bosons/imagtime.jl —— U1GTEMPO2 "冻结核不变 + K 侧二阶"的虚时复现
# （玻色 Toulouse，Matsubara 环 G(τ)），并输出画图数据 results/imagtime.csv
#
# 方案（环形隐式梯形，冻结核 = bosonic_Δτ 条目原样、face-value cell 配对）：
#   G_m − u·G_{m−1} − f·Σ_k D[m,k]·φ_k = 0,  G₀ := G_N + 1（双侧跳变，c-number）
#   u = e^{−ϵ_d h}；f = e^{−ϵ_d h/2}；φ_k = (G_{k−1}+G_k)/2（环梯形，φ₁ 左半取
#   G_N+1）；D = 前向 ηⱼₖ / 环绕 ηₖⱼ / 对角等时和 d₀（接触全量对角，隐式）。
# 跳变 = [d,d†]=1（与耦合无关），完整代入自由项与 φ₁ 是环上二阶的决定性要素。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/bosons/imagtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "common", "common.jl"))

"""
    weyl_k_imagtime(mdl, N) -> (G[1..N], h)

虚时 K 侧二阶（Weyl-K 环）。`G[m] ≈ G(τ_m)`，`G[N] = G(β⁻)`。
"""
function weyl_k_imagtime(mdl, N)
    h = mdl.β / N
    u = exp(-mdl.ϵd * h)
    f = exp(-mdl.ϵd * h / 2)
    corr = bosonic_Δτ(mdl.bath.spectrum; β=mdl.β, N=N, μ=mdl.μ)
    ηjk, ηkj = corr.data.ηⱼₖ, corr.data.ηₖⱼ
    d0 = ηjk[1] + ηkj[1]                        # 等时（两序三角之和）
    M = zeros(Float64, N, N)
    b = zeros(Float64, N)
    for m in 1:N
        M[m, m] += 1.0
        if m ≥ 2
            M[m, m-1] -= u
        else
            b[m] += u                           # −u·G₀ = −u·(G_N+1)：常数 +u
            M[m, N] -= u
        end
        for k in 1:N
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0
            if k == 1                           # φ₁ 左半 = (G_N+1)/2
                M[m, N] -= f * w * 0.5
                M[m, 1] -= f * w * 0.5
                b[m] += f * w * 0.5
            else
                M[m, k-1] -= f * w * 0.5
                M[m, k] -= f * w * 0.5
            end
        end
    end
    return M \ b, h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Ns=[16, 32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Gex = exact_Gτ(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "虚时 Weyl-K 二阶（冻结核 + K 侧）vs 严格 cell 链：玻色 Toulouse G(τ)")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64}[]
    println(io, "\n[1] Weyl-K 二阶 vs 精确 G(τ)（τ=β/4,β/2,3β/4）")
    es = Float64[]
    for N in Ns
        G, h = weyl_k_imagtime(mdl, N)
        e = maximum(abs.([G[m] for m in ms_of(N)] .- [Gex(m * h) for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, h, e, NaN, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照")
    for order in (1, 2)
        es = Float64[]
        for (i, N) in enumerate(Ns)
            vals, ms, h = chain_Gτ(mdl, order, N; ms=ms_of(N))
            e = maximum(abs.(vals .- [Gex(mod(m * h, mdl.β)) for m in ms]))
            push!(es, e)
            rows[i] = (rows[i][1], rows[i][2], rows[i][3],
                       order == 1 ? e : rows[i][4], order == 2 ? e : rows[i][5])
        end
        order_table(io, "    strict order=$order", Ns, es)
    end
    # CSV 输出（画图用）
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "imagtime.csv"), "w") do f
        println(f, "N,dtau,err_weylk,err_strict1,err_strict2")
        for (N, h, e1, e2, e3) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e,%.6e\n", N, h, e1, e2, e3)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "imagtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
