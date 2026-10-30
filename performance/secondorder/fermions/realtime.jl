# =============================================================================
# secondorder/fermions/realtime.jl —— 费米 Toulouse 实时 K 侧二阶（Weyl-K）
# 输出画图数据 results/realtime.csv
#
# 与玻色侧（../bosons/realtime.jl）同骨架，唯核组合不同：
#   玻色推迟核 = −(Δ⁺⁺+Δ⁺₋)（热因子 (1+n)−n 相消）
#   费米推迟核 = −(Δ⁺⁺−Δ⁺₋) = Δ⁺₋−Δ⁺⁺（热因子 (1−n)+n 相消；Pauli 反映射：
#               Δ⁺⁺ = −(1−n_F)R、Δ⁺₋ = +n_F R，单模验尸实测）
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/fermions/realtime.jl
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

"""
    weyl_k_realtime(mdl, N, T) -> (A[1..N], h)

费米实时 K 侧二阶（Weyl-K 隐式梯形，核 = X−P 推迟组合）。`A[m] ≈ A(t_m)`。
"""
function weyl_k_realtime(mdl, N, T)
    h = T / N
    u = exp(-im * mdl.ϵd * h)
    f = exp(-im * mdl.ϵd * h / 2)
    corr = fermionic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+)
    X = branch(corr; b1=:+, b2=:-)
    Δ(j, k) = -(X[j, k] + P[j, k])              # 费米推迟组合（统一约定；≡ 旧 X−P）
    d0 = Δ(1, 1)
    A = zeros(ComplexF64, N + 1)
    A[1] = 1.0
    for m in 1:N
        s = 0.0 + 0.0im
        for k in 1:m-1
            s += Δ(m+1, k+1) * (A[k] + A[k+1]) / 2
        end
        A[m+1] = (u*A[m] - f*s - f*d0*A[m]/2) / (1 + f*d0/2)
    end
    return A[2:end], h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[16, 32, 64, 128, 256, 512],
             fracs=[0.25, 0.5, 0.75])
    mdl = fermi_toulouse(; Lmodes=Lmodes)
    Aex = exact_return_amplitude(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "费米实时 Weyl-K 二阶（冻结核 + K 侧）vs 严格 cell 链：A(t)")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，T=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64}[]
    println(io, "\n[1] Weyl-K 二阶 vs 精确 A(t)")
    es = Float64[]
    for N in Ns
        A, h = weyl_k_realtime(mdl, N, T)
        e = maximum(abs.([A[m] for m in ms_of(N)] .- [Aex(m * h) for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, h, e, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照（费米扇区 Lie≡Strang 简并同玻色）")
    es = Float64[]
    for (i, N) in enumerate(Ns)
        st = real_step(mdl, 2, T / N)
        e = maximum(abs.([(st^m)[1, 1] for m in ms_of(N)] .-
                          [Aex(m * T / N) for m in ms_of(N)]))
        push!(es, e)
        rows[i] = (rows[i][1], rows[i][2], rows[i][3], e)
    end
    order_table(io, "    strict order=2", Ns, es)
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "realtime.csv"), "w") do f
        println(f, "N,dt,err_weylk,err_strict2")
        for (N, h, e1, e2) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e\n", N, h, e1, e2)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "realtime.csv"))")
end

# real_step 需在费米模型上工作（common.jl 的 boson 版本同构）
function real_step(mdl, order::Int, h::Real)
    eA(t) = exp(-im * t * mdl.hS)
    eB(t) = exp(-im * t * (mdl.hb + mdl.hi))
    s2(t) = eA(t / 2) * eB(t) * eA(t / 2)
    order == 1 && return eA(h) * eB(h)
    order == 2 && return s2(h)
    throw(ArgumentError("order 必须为 1/2"))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
