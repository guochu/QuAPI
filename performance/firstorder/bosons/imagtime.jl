# =============================================================================
# firstorder/bosons/imagtime.jl —— 玻色虚时：一阶冻结核环 vs 严格 Trotter 1/2 阶
# 对照：严格解 G(τ)。输出 results/imagtime.csv
#
#   · 一阶（frozen_dyson_Gτ 同款）：bosonic_Δτ 条目 + 等时两序和，旧装配；
#   · 严格一阶/二阶 Trotter 热环参照（chain_Gτ）。
# 注：顶点移位是实时轮廓特有的自包含改善；虚时的二阶改善属 K 侧要素打包
# （梯形+半步+接触+跳变代入），见 secondorder/bosons/imagtime.jl。
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "common", "trotter_frozen.jl"))

function main(; Lmodes=48, Ns=[32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Gex = exact_Gτ(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "玻色虚时：一阶冻结核环 vs 严格 1/2 阶 Trotter 参照")
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64}[]
    es1 = Float64[]; es3 = Float64[]; es4 = Float64[]
    for N in Ns
        h = mdl.β / N
        v1, _ = frozen_dyson_Gτ(mdl, N)
        ex1 = [Gex(mod(m * h, mdl.β)) for m in 1:N-1]
        e1 = maximum(abs.(v1 .- ex1))
        vals, ms2, _ = chain_Gτ(mdl, 1, N; ms=ms_of(N))
        e3 = maximum(abs.(vals .- [Gex(mod(m * h, mdl.β)) for m in ms2]))
        vals4, ms4, _ = chain_Gτ(mdl, 2, N; ms=ms_of(N))
        e4 = maximum(abs.(vals4 .- [Gex(mod(m * h, mdl.β)) for m in ms4]))
        push!(es1, e1); push!(es3, e3); push!(es4, e4)
        push!(rows, (N, h, e1, e3, e4))
    end
    order_table(io, "    一阶冻结核环", Ns, es1)
    order_table(io, "    严格一阶 Trotter 环（参照）", Ns, es3)
    order_table(io, "    严格二阶 Trotter 环（参照）", Ns, es4)
    out = joinpath(@__DIR__, "results"); mkpath(out)
    open(joinpath(out, "imagtime.csv"), "w") do f
        println(f, "N,dtau,err_frozen1,err_strict1,err_strict2")
        for (N, h, a, c, d) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e,%.6e\n", N, h, a, c, d)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "imagtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
