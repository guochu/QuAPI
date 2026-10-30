# =============================================================================
# firstorder/fermions/imagtime.jl —— 费米虚时：一阶冻结核环 vs 严格 Trotter 1/2 阶
# 对照：严格解 G(τ)（greater 约定谱求和）。输出 results/imagtime.csv
#
#   · 一阶（朴素矩形）：反周期双符号闭合、节点值权重、无梯形、无跳变常数
#     代入 φ₁；
#   · 严格参照：精确 cell 热环 (I+T^N)⁻¹（费米分母，对离散链精确）。
# 注：顶点移位是实时轮廓特有的自包含改善；虚时的二阶改善属 K 侧要素打包
# （梯形+半步+接触+跳变常数代入），见 secondorder/fermions/imagtime.jl。
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "fermions", "common.jl"))

"""
    frozen_first_order_ring(mdl, N) -> (G[1..N], h)

费米虚时一阶冻结核环（朴素矩形）。`G[m] ≈ G(τ_m)`（greater 约定）。
"""
function frozen_first_order_ring(mdl, N)
    h = mdl.β / N
    u = exp(-mdl.ϵd * h)
    f = exp(-mdl.ϵd * h / 2)
    corr = fermionic_Δτ(mdl.bath.spectrum; β=mdl.β, N=N, μ=mdl.μ)
    ηjk, ηkj = corr.data.ηⱼₖ, corr.data.ηₖⱼ
    d0 = ηjk[1] + ηkj[1]
    M = zeros(ComplexF64, N, N)
    b = zeros(ComplexF64, N)
    for m in 1:N
        M[m, m] += 1.0
        if m ≥ 2; M[m, m-1] -= u
        else; b[m] -= u; M[m, N] += u
        end
        for k in 1:N
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0
            M[m, k] -= f * w                     # 矩形：节点值、无跳变常数代入
        end
    end
    return M \ b, h
end

"严格参照（greater 约定谱求和）与精确 cell 热环。"
function exact_refs(mdl, N)
    ev = eigen(Hermitian(mdl.hfull))
    h = mdl.β / N
    ex = [-sum(ev.vectors[1, p]^2 * exp(-ev.values[p] * m * h) /
               (1 + exp(-mdl.β * ev.values[p])) for p in 1:length(ev.values)) for m in 1:N]
    T = exp(-h * mdl.hfull)
    P = inv(mdl.Id + T^N)
    strict = [-(T^m * P)[1, 1] for m in 1:N]     # 精确 cell 热环（greater 约定）
    return ex, strict, h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Ns=[16, 32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = fermi_toulouse(; Lmodes=Lmodes)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "费米虚时：一阶冻结核环 vs 严格 1/2 阶 Trotter 参照：G(τ)")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64}[]
    println(io, "\n[1] 一阶冻结核环 vs 精确 G(τ)")
    es = Float64[]
    for N in Ns
        G, h = frozen_first_order_ring(mdl, N)
        ex, _ = exact_refs(mdl, N)
        push!(es, maximum(abs.([G[m] - ex[m] for m in ms_of(N)])))
        push!(rows, (N, h, es[end], NaN))
    end
    order_table(io, "    一阶冻结核环（朴素矩形）", Ns, es)
    println(io, "\n[2] 严格 cell 热环参照")
    es = Float64[]
    for (i, N) in enumerate(Ns)
        _, h = frozen_first_order_ring(mdl, N)
        ex, _ = exact_refs(mdl, N)
        T = exp(-h * mdl.hfull)
        P = inv(mdl.Id + T^N)
        vals = [-(T^m * P)[1, 1] for m in 1:N]
        e = maximum(abs.([vals[m] - ex[m] for m in ms_of(N)]))
        push!(es, e)
        rows[i] = (rows[i][1], rows[i][2], rows[i][3], e)
    end
    order_table(io, "    strict order=2（精确 cell 热环）", Ns, es)
    out = joinpath(@__DIR__, "results"); mkpath(out)
    open(joinpath(out, "imagtime.csv"), "w") do f
        println(f, "N,dtau,err_frozen1,err_strict2")
        for (N, h, a, b) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e\n", N, h, a, b)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "imagtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
