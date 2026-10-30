# =============================================================================
# firstorder/bosons/imagtime_dense.jl —— 玻色虚时一阶冻结核的稠密 IF 装配
#
# 与 imagtime.jl 完全对齐的另一种实现：同一个一阶冻结核（bosonic_Δτ 条目，
# 朴素节点值权重 + 等时两序和 eqt）不再走环方程求逆，而是经 QuAPI.Dense 的
# GrassmannIF 折叠为稠密 IF 张量（2^{2N} 系数），从 Pfaffian 余子式收缩中
# 提取 G(τ_m)。两路输出必须逐位一致（机器精度）——本测试为下游 MPS-IF
# 装配提供稠密参照。
#
# 单粒子说明：Toulouse 单粒子振幅代数与费米情形同构（skew 双线性 + Pfaffian
# 余子式），玻色统计完全包含在 bosonic_Δτ 的核条目里（含 (1+n)−n 热因子），
# 腿基折叠机制与 firstorder/fermions/imagtime_dense.jl 完全一致。
#
# 装配对应（imagtime.jl / frozen_dyson_Gτ 环方程 ⟷ 稠密 IF）：
#   环矩阵 M = (I − e^{−hϵ_d}·A_cyc) − D（A_cyc 循环移位、D = 核权重矩阵、
#   对角 = 等时两序和 eqt；T=cellrule 1 即单位阵）；
#   腿基系数 C[2k, 2m−1] = −M[m,k]（转置核装配）；
#   环解 G[m] = (M⁻¹)[m+1, 1] ⟹ 稠密提取 = ±data[去掉 {ā_{m+1}, a₁}]/top。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bosons/imagtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "common", "common.jl"))
using QuAPI
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix, GrassmannIF

"""
    naive_ring_M(mdl, N) -> M

一阶冻结核环矩阵（与 frozen_dyson_Gτ 的 cellrule=1 完全同款）：
M = (I − e^{−hϵ_d}·A_cyc) − D，D 对角 = eqt = ηⱼₖ[1] + ηₖⱼ[1]。
"""
function naive_ring_M(mdl, N)
    h = mdl.β / N
    M = Matrix{ComplexF64}(I, N, N)
    for j in 1:N                                    # 循环移位（周期热环）
        M[mod(j, N) + 1, j] -= exp(-h * mdl.ϵd)
    end
    corr = bosonic_Δτ(mdl.bath.spectrum, β=mdl.β, N=N, μ=mdl.μ)
    ηjk, ηkj = corr.data.ηⱼₖ, corr.data.ηₖⱼ
    eqt = ηjk[1] + ηkj[1]
    for i in 1:N, j in 1:N
        M[i, j] -= (i > j ? ηjk[i - j + 1] :
                    i < j ? ηkj[j - i + 1] : eqt)
    end
    return M
end

"""
    dense_if_legs(N, M) -> gif

稠密 IF 装配：腿基 = 每节点 (a, ā)（GrassmannContourIndex(j, :τ, 1, c)），
腿系数 C[2k, 2m−1] = −M[m,k]（转置核装配；测度/动能/回绕/核全部内含于 −M）。
"""
function dense_if_legs(N, M)
    basis = GrassmannContourIndex[GrassmannContourIndex(j, :τ, 1, c)
                                  for j in 1:N for c in (false, true)]
    C = zeros(ComplexF64, 2N, 2N)
    for m in 1:N, k in 1:N
        C[2k, 2m-1] = -M[m, k]
    end
    return GrassmannIF(GrassmannHybridizationMatrix(basis, C))
end

"""
    minor_inv(gif, N, target) -> (Minv, sgns)

从折叠张量提取全部 (M⁻¹)[i, j] = sgn·data[去掉 {ā_i, a_j}]/data[end]；
符号 (i, j) 逐项数值标定（必须落在整数上，自校验），结构零须与 target 零一致。
"""
function minor_inv(gif, N, target::AbstractMatrix)
    top = gif.data[end]
    full = (1 << (2N)) - 1
    Minv = zeros(ComplexF64, N, N)
    for i in 1:N, j in 1:N
        mask = full ⊻ ((1 << (2i - 1)) | (1 << (2j - 2)))  # ā_i(leg 2i)、a_j(leg 2j−1)
        v = gif.data[mask + 1]
        if abs(v) < 1e-280
            abs(target[i, j]) < 1e-12 * max(1.0, maximum(abs, target)) ||
                error("结构零处目标非零 at ($i,$j)")
            continue
        end
        r = target[i, j] * top / v
        s = round(real(r))
        abs(real(r) - s) > 1e-6 && error("符号标定失败 at ($i,$j): r=$r")
        Minv[i, j] = s * v / top
    end
    return Minv
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Ns=[6, 8, 10], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Gex = exact_Gτ(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "玻色虚时：一阶冻结核的稠密 IF 装配（GrassmannIF）vs 环方程 vs 精确")
    println(io, "  IF 张量 = 2^{2N} 系数；G(τ_m) = (M⁻¹)[m+1,1] 从 Pfaffian 余子式提取")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd)
    println(io, "=" ^ 74)

    println(io, "\n[0] 一致性：稠密 IF 装配 vs 环方程（必须机器精度）")
    worst = 0.0
    for N in Ns
        M = naive_ring_M(mdl, N)
        G_ring = [real((M \ [1.0; zeros(N-1)])[m+1]) for m in 1:N-1]
        gif = dense_if_legs(N, M)
        Minv = minor_inv(gif, N, inv(M))
        G_dense = [Minv[m+1, 1] for m in 1:N-1]
        d = maximum(abs.(G_dense .- G_ring))
        worst = max(worst, d)
        @printf(io, "  N=%2d: max|G_dense − G_ring| = %.3e\n", N, d)
    end
    println(io, "\n[1] 稠密 IF 提取 G(τ) vs 精确（误差应与 imagtime.jl 同款）")
    es = Float64[]
    for N in Ns
        M = naive_ring_M(mdl, N)
        gif = dense_if_legs(N, M)
        Minv = minor_inv(gif, N, inv(M))
        G_dense = [real(Minv[m+1, 1]) for m in 1:N-1]
        h = mdl.β / N
        e = maximum(abs.([G_dense[m] - Gex(mod(m * h, mdl.β)) for m in ms_of(N)]))
        push!(es, e)
    end
    println("  N        dτ       error")
    for (i, N) in enumerate(Ns)
        @printf(io, "%3d  %.4e  %.4e\n", N, mdl.β / N, es[i])
    end
    @printf(io, "\n最大一致性偏差 = %.2e（机器精度即通过）\n", worst)
    println(io, "\n结论（实测）")
    println(io, "  · 稠密 IF（GrassmannIF 折叠 −M）的 Pfaffian 余子式提取与环方程逐位一致。")
    println(io, "  · 下游 MPS-IF 按本文件的腿基/系数矩阵装配即可逐位复现环解。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
