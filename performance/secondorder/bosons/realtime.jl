# =============================================================================
# secondorder/bosons/realtime.jl —— U1GTEMPO2 "冻结核不变 + K 侧二阶"的实时复现
# （玻色 Toulouse，A(t) = [e^{-iHt}]₁₁），并输出画图数据 results/realtime.csv
#
# 方案（隐式梯形，冻结核 = bosonic_Δt 条目原样、face-value cell 配对）：
#   A_j − u·A_{j−1} + f·Σ_{k<j} Δ[j,k]·φ_k + f·d₀·φ_j = 0,  A₀ = 1
#   u = e^{−iϵ_d h}；f = e^{−iϵ_d h/2}（半步相位）；φ_k = (A_{k−1}+A_k)/2（Weyl ½
#   源平均）；Δ = −(Δ⁺⁺+Δ⁺⁻)（推迟 cell 对积分）；d₀ = 等时三角（接触全量对角）。
# 机制：与严格 Strang 链相差 O(h²)，从而继承其二阶。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/bosons/realtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "common", "common.jl"))

"""
    weyl_k_realtime(mdl, N, T) -> (A[1..N], h)

实时 K 侧二阶（Weyl-K 隐式梯形）。`A[m] ≈ A(t_m)`。
"""
function weyl_k_realtime(mdl, N, T)
    h = T / N
    u = exp(-im * mdl.ϵd * h)
    f = exp(-im * mdl.ϵd * h / 2)
    corr = bosonic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+)
    X = branch(corr; b1=:+, b2=:-)
    d0 = -(P[1, 1] + X[1, 1])                   # 等时接触三角
    A = zeros(ComplexF64, N + 1)
    A[1] = 1.0
    for m in 1:N
        s = 0.0 + 0.0im
        for k in 1:m-1
            s += -(P[m+1, k+1] + X[m+1, k+1]) * (A[k] + A[k+1]) / 2
        end
        A[m+1] = (u * A[m] - f * s - f * d0 * A[m] / 2) / (1 + f * d0 / 2)
    end
    return A[2:end], h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[16, 32, 64, 128, 256, 512],
             fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Aex = exact_return_amplitude(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "实时 Weyl-K 二阶（冻结核 + K 侧）vs 严格 cell 链：玻色 Toulouse A(t)")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64}[]
    println(io, "\n[1] Weyl-K 二阶 vs 精确 A(t)")
    es = Float64[]
    for N in Ns
        A, h = weyl_k_realtime(mdl, N, T)
        e = maximum(abs.([A[m] for m in ms_of(N)] .- [Aex(m * h) for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, h, e, NaN, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照")
    for order in (1, 2)
        es = Float64[]
        for (i, N) in enumerate(Ns)
            st = real_step(mdl, order, T / N)
            e = maximum(abs.([(st^m)[1, 1] for m in ms_of(N)] .-
                              [Aex(m * T / N) for m in ms_of(N)]))
            push!(es, e)
            rows[i] = (rows[i][1], rows[i][2], rows[i][3],
                       order == 1 ? e : rows[i][4], order == 2 ? e : rows[i][5])
        end
        order_table(io, "    strict order=$order", Ns, es)
    end
    # CSV 输出（画图用）
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "realtime.csv"), "w") do f
        println(f, "N,dt,err_weylk,err_strict1,err_strict2")
        for (N, h, e1, e2, e3) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e,%.6e\n", N, h, e1, e2, e3)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "realtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
