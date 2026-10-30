# =============================================================================
# secondorder/fermions/imagtime.jl —— 费米 Toulouse 虚时 K 侧二阶（Weyl-K 环）
# 输出画图数据 results/imagtime.csv
#
# 与玻色侧（../bosons/imagtime.jl）同骨架，两处费米特异性：
#   · 跳变：反周期 + δ 源 ⟹ G₀ = −1 − G_N（玻色为 G₀ = G_N + 1）；符号翻转
#     同时出现在右端常数与 G_N 列系数（φ₁ 左半 = (−1−G_N)/2）；
#   · 约定：解为 Matsubara greater 约定 G(τ) = −⟨d(τ)d†(0)⟩（V=0 极限
#     −(1−n̄)e^{−ϵτ} 逐点机器精度验证）；严格参照取谱求和同约定。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/fermions/imagtime.jl
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

"""
    weyl_k_imagtime(mdl, N) -> (G[1..N], h)

费米虚时 K 侧二阶（Weyl-K 环，greater 约定）。`G[m] ≈ G(τ_m)`，`G[N] = G(β⁻)`。
"""
function weyl_k_imagtime(mdl, N)
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
        if m ≥ 2
            M[m, m-1] -= u
        else
            b[m] -= u                           # −u·G₀ = −u·(−1−G_N) = +u + u·G_N
            M[m, N] += u
        end
        for k in 1:N
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0
            if k == 1                           # φ₁ = (−1−G_N+G₁)/2
                M[m, N] += f * w * 0.5
                M[m, 1] -= f * w * 0.5
                b[m] -= f * w * 0.5
            else
                M[m, k-1] -= f * w * 0.5
                M[m, k] -= f * w * 0.5
            end
        end
    end
    return M \ b, h
end

"严格参照（greater 约定谱求和）：G(τ_m) = −Σ|U₁ₚ|²e^{−λ_pτ}/(1+e^{−βλ_p})。"
function exact_Gτ(mdl, N)
    ev = eigen(Hermitian(mdl.hfull))
    h = mdl.β / N
    return [(-1)^0 * -sum(ev.vectors[1, p]^2 * exp(-ev.values[p] * m * h) /
                          (1 + exp(-mdl.β * ev.values[p])) for p in 1:length(ev.values))
            for m in 1:N], h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Ns=[16, 32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = fermi_toulouse(; Lmodes=Lmodes)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "费米虚时 Weyl-K 二阶（冻结核 + K 侧）vs 严格 cell 链：G(τ)")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64}[]
    println(io, "\n[1] Weyl-K 二阶 vs 精确 G(τ)（greater 约定，τ=β/4,β/2,3β/4）")
    es = Float64[]
    for N in Ns
        G, h = weyl_k_imagtime(mdl, N)
        ex, _ = exact_Gτ(mdl, N)
        e = maximum(abs.([G[m] - ex[m] for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, h, e, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照（热环公式，(I+T^N)⁻¹ 费米分母）")
    es = Float64[]
    for (i, N) in enumerate(Ns)
        h = mdl.β / N
        T = exp(-h * mdl.hfull)
        P = inv(mdl.Id + T^N)
        vals = [-(T^m * P)[1, 1] for m in 1:N]   # greater 约定（与参照同号）
        ex, _ = exact_Gτ(mdl, N)
        e = maximum(abs.([vals[m] - ex[m] for m in ms_of(N)]))
        push!(es, e)
        rows[i] = (rows[i][1], rows[i][2], rows[i][3], e)
    end
    order_table(io, "    strict order=2（精确 cell 热环）", Ns, es)
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "imagtime.csv"), "w") do f
        println(f, "N,dtau,err_weylk,err_strict2")
        for (N, h, e1, e2) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e\n", N, h, e1, e2)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "imagtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
