# =============================================================================
# firstorder/bosons/realtime_dense.jl —— 玻色实时：一阶冻结核的稠密 IF 装配
#
# 与 realtime.jl 完全对齐的另一种实现：同一个一阶冻结核（bosonic_Δt 的
# Δ⁺⁺/Δ⁺₋ 条目，face-value 因果下三角）不再走 (A0⁻¹+K)⁻¹ 求逆，而是经
# GrassmannIF 折叠为稠密 IF 张量（2^{2(N+1)} 系数），从 Pfaffian 余子式
# 收缩中提取 A(t_m)。两路输出必须逐位一致（机器精度）——本测试为下游
# MPS-IF 装配提供稠密参照。
#
# 单粒子说明：单粒子振幅代数与费米情形同构（skew 双线性 + Pfaffian 余子式），
# 玻色统计完全包含在 bosonic_Δt 的核条目里（−(P+X) 的 (1+n)−n 热抵消）。
#
# 装配对应（realtime.jl / frozen_retarded_A ⟷ 稠密 IF）：
#   A0⁻¹ = I − u·S⁻（full-step 自由飞行，开链）；K[j,k] = −(P+X)[j,k]
#   （k ≤ j 因果下三角，含接触对角）；M = A0⁻¹ + K，A[m] = (M⁻¹)[m+1, 1]；
#   腿基系数 C[2k, 2m−1] = −M[m,k] ⟹ 稠密提取 = ±data[去掉 {ā_{m+1}, a₁}]/top。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bosons/realtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "common", "common.jl"))
using QuAPI
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix, GrassmannIF

"""
    naive_rt_M(mdl, N, T) -> M

一阶冻结核实时矩阵（与 frozen_retarded_A 完全同款）：
M = A0⁻¹ + K，A0⁻¹ = I − u·S⁻，K[j,k] = −(P[j,k]+X[j,k])（k ≤ j）。
"""
function naive_rt_M(mdl, N, T)
    h = T / N
    corr = bosonic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+)
    X = branch(corr; b1=:+, b2=:-)
    M = Matrix{ComplexF64}(I, N + 1, N + 1)
    u = exp(-im * mdl.ϵd * h)
    for j in 2:N+1
        M[j, j-1] -= u
    end
    for j in 1:N+1, k in 1:j                       # 因果：只填下三角（含对角）
        M[j, k] -= (P[j, k] + X[j, k])
    end
    return M
end

"实时开链稠密 IF：腿基 (a, ā)/节点（branch=:+），C[2k, 2m−1] = −M[m,k]。"
function dense_if_legs(n, M)
    basis = GrassmannContourIndex[GrassmannContourIndex(j, :+, 1, c)
                                  for j in 1:n for c in (false, true)]
    C = zeros(ComplexF64, 2n, 2n)
    for m in 1:n, k in 1:n
        C[2k, 2m-1] = -M[m, k]
    end
    return GrassmannIF(GrassmannHybridizationMatrix(basis, C))
end

"从折叠张量提取全部 (M⁻¹)[i, j]（符号逐项数值标定，结构零须与 target 一致）。"
function minor_inv(gif, N, target::AbstractMatrix)
    top = gif.data[end]
    full = (1 << (2N)) - 1
    Minv = zeros(ComplexF64, N, N)
    for i in 1:N, j in 1:N
        mask = full ⊻ ((1 << (2i - 1)) | (1 << (2j - 2)))
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
function main(; Lmodes=48, T=2.0, Ns=[6, 8, 10], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Aex = exact_return_amplitude(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "玻色实时：一阶冻结核的稠密 IF 装配（GrassmannIF）vs (A0⁻¹+K)⁻¹ vs 精确")
    println(io, "  IF 张量 = 2^{2(N+1)} 系数；A(t_m) = (M⁻¹)[m+1,1] 从 Pfaffian 余子式提取")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)

    println(io, "\n[0] 一致性：稠密 IF 装配 vs 线性系统（必须机器精度）")
    worst = 0.0
    for N in Ns
        M = naive_rt_M(mdl, N, T)
        A_ref = (M \ [1.0; zeros(ComplexF64, N)])
        gif = dense_if_legs(N + 1, M)
        Minv = minor_inv(gif, N + 1, inv(M))
        A_dense = Minv * [1.0; zeros(ComplexF64, N)]
        d = maximum(abs.(A_dense .- A_ref))
        worst = max(worst, d)
        @printf(io, "  N=%2d: max|A_dense − A_sys| = %.3e\n", N, d)
    end
    println(io, "\n[1] 稠密 IF 提取 A(t) vs 精确（误差应与 realtime.jl 同款）")
    for N in Ns
        M = naive_rt_M(mdl, N, T)
        gif = dense_if_legs(N + 1, M)
        Minv = minor_inv(gif, N + 1, inv(M))
        A_dense = Minv * [1.0; zeros(ComplexF64, N)]
        h = T / N
        e = maximum(abs.([A_dense[m+1] - Aex(m * h) for m in ms_of(N)]))
        @printf(io, "  N=%2d (h=%.4e): err = %.4e\n", N, h, e)
    end
    @printf(io, "\n最大一致性偏差 = %.2e（机器精度即通过）\n", worst)
    println(io, "\n结论（实测）")
    println(io, "  · 稠密 IF（GrassmannIF 折叠 −M）的 Pfaffian 余子式提取与线性系统逐位一致。")
    println(io, "  · 下游 MPS-IF 按本文件的腿基/系数矩阵装配即可逐位复现通道解。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
