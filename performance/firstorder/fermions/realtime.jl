# =============================================================================
# firstorder/fermions/realtime.jl —— 费米实时：一阶冻结核装配 vs 顶点移位改善
# 对照：严格解 A(t) 与严格二阶 Trotter。输出 results/realtime.csv
#
#   · 一阶（face-value 配对，不移位）：K[j,k] = (Δ⁺₋−Δ⁺⁺)[j,k]（费米推迟
#     组合，lag j−k 直读）+ 杂质相位每步 u=e^{−iϵ_d h}；
#   · 顶点移位改善：节点 lag d ↔ 条目 lag d−1 + 三角并入次对角（同批条目，
#     误差常数下降，阶数不变）。
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "fermions", "common.jl"))

function frozen_first_order(mdl, N, T)
    h = T / N
    u = exp(-im * mdl.ϵd * h)
    corr = fermionic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+); X = branch(corr; b1=:+, b2=:-)
    Δf(ℓ) = -(X[ℓ+1, 1] + P[ℓ+1, 1])
    A = zeros(ComplexF64, N + 1); A[1] = 1.0
    for m in 1:N
        s = 0.0 + 0.0im
        for d in 1:m-1
            s += Δf(d) * A[m-d]            # face-value：lag d ↔ 条目 d（不移位）
        end
        A[m+1] = u * (A[m] - s)
    end
    return A[2:end], h
end

function vertex_shifted(mdl, N, T)
    h = T / N
    u = exp(-im * mdl.ϵd * h)
    corr = fermionic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+); X = branch(corr; b1=:+, b2=:-)
    Δf(ℓ) = -(X[ℓ+1, 1] + P[ℓ+1, 1])
    g0 = Δf(0)
    g = [Δf(d) for d in 1:N-1]
    M = Matrix{ComplexF64}(I, N, N); b = zeros(ComplexF64, N)
    for m in 1:N
        if m ≥ 2; M[m, m-1] -= u * (1 - g0); else; b[m] += u * (1 - g0); end
        for d in 2:m
            if m - d ≥ 1; M[m, m-d] += u * g[d-1]; else; b[m] -= u * g[d-1]; end
        end
    end
    return M \ b, h
end

function main(; Lmodes=48, T=2.0, Ns=[32, 64, 128, 256, 512], fracs=[0.25, 0.5, 0.75])
    mdl = fermi_toulouse(; Lmodes=Lmodes)
    Aex = exact_return_amplitude(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "费米实时：一阶冻结核 vs 顶点移位改善 vs 严格参照")
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64,Float64}[]
    es1 = Float64[]; es2 = Float64[]; es3 = Float64[]; es4 = Float64[]
    for N in Ns
        h = T / N
        e1 = maximum(abs.([frozen_first_order(mdl, N, T)[1][m] for m in ms_of(N)] .-
                           [Aex(m * h) for m in ms_of(N)]))
        e2 = maximum(abs.([vertex_shifted(mdl, N, T)[1][m] for m in ms_of(N)] .-
                           [Aex(m * h) for m in ms_of(N)]))
        eA = exp(-im * h / 2 * mdl.hS)
        eB = exp(-im * h * (mdl.hb + mdl.hi))
        st1 = eA^2 * eB                    # order1: eA(1)·eB = e^{-iεh}·eB
        st1 = exp(-im * h * mdl.hS) * eB
        st2 = eA * eB * eA
        e3 = maximum(abs.([(st1^N)[1, 1] - Aex(T)]))
        e4 = maximum(abs.([(st2^N)[1, 1] - Aex(T)]))
        push!(es1, e1); push!(es2, e2); push!(es3, e3); push!(es4, e4)
        push!(rows, (N, h, e1, e2, e3, e4))
    end
    order_table(io, "    一阶冻结核（face-value）", Ns, es1)
    order_table(io, "    顶点移位改善（仍一阶，常数↓）", Ns, es2)
    order_table(io, "    严格一阶 Trotter（终点，参照）", Ns, es3)
    order_table(io, "    严格二阶 Trotter（终点，参照）", Ns, es4)
    out = joinpath(@__DIR__, "results"); mkpath(out)
    open(joinpath(out, "realtime.csv"), "w") do f
        println(f, "N,dt,err_frozen1,err_shifted,err_strict1,err_strict2")
        for (N, h, a, b, c, d) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e,%.6e,%.6e\n", N, h, a, b, c, d)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "realtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
