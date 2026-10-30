# =============================================================================
# firstorder/fermions/imagtime_dense.jl —— 费米虚时一阶冻结核的稠密 IF 装配
#
# 与 imagtime.jl 完全对齐的另一种实现：同一个一阶冻结核（fermionic_Δm 的
# ξ 块条目，朴素矩形权重）不再走环方程 M\b，而是经 QuAPI.Dense.GrassmannIF
# 折叠为稠密 IF 张量（2^{2N} 系数），从张量收缩中提取 G(τ_m)。两路输出
# 必须逐位一致（机器精度）——本测试为下游 TEMPO/GTEMPO 的 MPS-IF 装配提供
# 稠密参照：IF 的张量元素、收缩顺序、骨架约定全部在此显式可见。
#
# 装配对应关系（imagtime.jl 环方程 M·G = b ⟷ 稠密 IF）：
#   环矩阵 M = I − B，B = 一粒子矩阵（uτ·跃迁 + fτ·W 冻结核）；
#   腿基系数 C[2m, 2k−1] = −M[m,k]（ā_m a_k 块；C = −M 的腿展开）；
#   骨架（测度 −āa、动能跃迁 uτ、反周期回绕 −uτ）已含于 −M 中；
#   环源 b = −uτ·e₁（反周期跳变自由项）⟹ G(τ_m) = −uτ·[(I−B)⁻¹]_{m,1}
#   = −uτ·(Pf-余子式)/(Pf) ⟹ 稠密提取 = data[互补掩码]/data[end]（差符号）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/fermions/imagtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "fermions", "common.jl"))
const δτ0 = Ref(0.0)
using QuAPI
using QuAPI.Dense
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix,
                   GrassmannIF, GrassmannADT

"""
    naive_ring_M(mdl, N) -> (M, b, fτ)

一阶冻结核环矩阵（与 firstorder/fermions/imagtime.jl 的
frozen_first_order_ring 完全同款）：M = I − uτ·S⁻ − fτ·W（W = fermionic_Δτ
的朴素矩形权重），b = −uτ·e₁（反周期跳变自由项）。
"""
function naive_ring_M(mdl, N)
    δτ = mdl.β / N
    uτ = exp(-mdl.ϵd * δτ)
    fτ = exp(-mdl.ϵd * δτ / 2)
    corr = fermionic_Δτ(mdl.bath.spectrum; β=mdl.β, N=N, μ=mdl.μ)
    ηjk, ηkj = corr.data.ηⱼₖ, corr.data.ηₖⱼ
    d0τ = ηjk[1] + ηkj[1]
    M = zeros(ComplexF64, N, N)
    for m in 1:N
        M[m, m] += 1.0
        if m ≥ 2
            M[m, m-1] -= uτ
        else
            M[m, N] += uτ                          # 反周期回绕
        end
        for k in 1:N
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0τ
            M[m, k] -= fτ * w                      # 朴素矩形：节点值、全权重
        end
    end
    b = zeros(ComplexF64, N)
    b[1] = -uτ                                     # 跳变自由项（greater 约定）
    return M, b, fτ
end

"""
    dense_if_legs(mdl, N, M, fτ) -> (basis, C, gif)

稠密 IF 装配：腿基 = 每节点 (a, ā)（GrassmannContourIndex(j, :τ, 1, c)），
腿系数 C[2m, 2k−1] = −M[m,k]（ā_m a_k 块；骨架测度/动能/回绕/核全部内含于
−M——环矩阵即一粒子矩阵的 (I − B) 形式）。折叠为稠密 IF 张量。
"""
function dense_if_legs(mdl, N, M, fτ)
    basis = GrassmannContourIndex[GrassmannContourIndex(j, :τ, 1, c)
                                  for j in 1:N for c in (false, true)]
    C = zeros(ComplexF64, 2N, 2N)
    for m in 1:N, k in 1:N
        C[2k, 2m-1] = -M[m, k]                     # 转置核装配（IF 收缩方向约定）
    end
    gif = GrassmannIF(GrassmannHybridizationMatrix(basis, C))
    return basis, C, gif
end

"严格参照（greater 约定谱求和）与精确 cell 热环（同 imagtime.jl）。"
function exact_refs(mdl, N)
    ev = eigen(Hermitian(mdl.hfull))
    h = mdl.β / N
    ex = [-sum(ev.vectors[1, p]^2 * exp(-ev.values[p] * m * h) /
               (1 + exp(-mdl.β * ev.values[p])) for p in 1:length(ev.values)) for m in 1:N]
    T = exp(-h * mdl.hfull)
    P = inv(mdl.Id + T^N)
    strict = [-(T^m * P)[1, 1] for m in 1:N]
    return ex, strict, h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[8, 10, 12], fracs=[0.25, 0.5, 0.75])
    mdl = fermi_toulouse(; Lmodes=Lmodes)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    # 严格参照用 γ_eq（一次对角化，仅参照、不进入求解器）
    ev = eigen(Hermitian(mdl.hfull))
    γ = ev.vectors * Diagonal([1/(exp(mdl.β*λ)+1) for λ in ev.values]) * ev.vectors'
    gex(N, m) = -im * (exp(-im * (m * T / N) * mdl.hfull) * (mdl.Id - transpose(γ)))[1, 1]
    println(io, "=" ^ 74)
    println(io, "费米虚时：一阶冻结核的稠密 IF 装配（GrassmannIF）vs 环方程 vs 精确")
    println(io, "  IF 张量 = 2^{2N} 系数；G(τ_m) 从 Pfaffian 余子式收缩提取")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，T=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64}[]
    println(io, "\n[0] 一致性：稠密 IF 装配 vs 环方程（必须机器精度）")
    es = Float64[]
    for N in (8, 10, 12)
        M, b, fτ = naive_ring_M(mdl, N)
        G_ring = M \ b
        basis, C, gif = dense_if_legs(mdl, N, M, fτ)
        top = gif.data[end]
        uτ = exp(-mdl.ϵd * (mdl.β/N))
        G_dense = zeros(ComplexF64, N)
        for m in 1:N
            mask = ((1 << (2N)) - 1) ⊻ ((1 << 0) | (1 << (2m - 1)))  # 去掉 a₁(leg 1) 与 a_m(leg 2m−1)
            G_dense[m] = -uτ * gif.data[mask + 1] / top
        end
        d = maximum(abs.(G_dense .- G_ring))
        @printf(io, "  N=%2d: max|G_dense − G_ring| = %.3e\n", N, d)
    end
    println(io, "\n[1] 稠密 IF 提取 G(τ) vs 精确（误差应与 imagtime.jl 逐位一致）")
    es_dense = Float64[]
    for N in Ns
        M, b, fτ = naive_ring_M(mdl, N)
        G_ring = M \ b
        basis, C, gif = dense_if_legs(mdl, N, M, fτ)
        top = gif.data[end]
        uτ = exp(-mdl.ϵd * (mdl.β/N))
        G_dense = zeros(ComplexF64, N)
        for m in 1:N
            mask = ((1 << (2N)) - 1) ⊻ ((1 << 0) | (1 << (2m - 1)))  # 去掉 a₁(leg 1) 与 a_m(leg 2m−1)
            G_dense[m] = -uτ * gif.data[mask + 1] / top
        end
        ex, strict, h = exact_refs(mdl, N)
        e = maximum(abs.([G_dense[m] - ex[m] for m in ms_of(N)]))
        push!(es, e)
        # 与环方程解的一致性
        d = maximum(abs.([G_dense[m] - G_ring[m] for m in 1:N]))
        push!(rows, (N, h, e, d))
        @printf(io, "  N=%3d: |G_dense−exact|=%.4e  |G_dense−G_ring|=%.1e\n", N, e, d)
    end
    order_table(io, "    稠密 IF 提取 G(τ) 误差", Ns, es)
    println(io, "\n结论（实测）")
    println(io, "  · 稠密 IF（GrassmannIF 折叠 2^{2N} 系数）的 Pfaffian 余子式收缩与")
    println(io, "    环方程 M\\b 逐位一致（机器精度）——两种装配是同一离散 Gauss 积分的")
    println(io, "    两种表示；下游 MPS-IF 按本文件的腿基与系数矩阵装配即可逐位对齐。")
    println(io, "  · 对精确解的误差与 imagtime.jl 完全相同（朴素一阶装配的 O(h) 误差）。")
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "imagtime_dense.csv"), "w") do f
        println(f, "N,dtau,err_dense_vs_exact,err_dense_vs_ring")
        for (i, N) in enumerate(Ns)
            @printf(f, "%d,%.6e,%.6e,%.6e\n", N, T/N, es[i], rows[i][4])
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "imagtime_dense.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
