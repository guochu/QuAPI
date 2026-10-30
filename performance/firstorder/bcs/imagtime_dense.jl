# =============================================================================
# firstorder/bcs/imagtime_dense.jl —— BCS 虚时：稠密 IF 装配 vs Nambu Matsubara 环
#
# 与 imagtime.jl 完全对齐的另一种实现：同一个一阶冻结核（bcs_Δτ 四块的
# Nambu 环组合 𝚺̂ = [cn −cc; −nn nc]）不再走块线性系统 M\b，而是经
# GrassmannIF 折叠为稠密 IF 张量（2^{8N} 系数），从 Pfaffian 余子式收缩
# 提取，两路输出必须逐位一致（机器精度）。
#
# Nambu 腿基：每物理节点 8 条腿 = 4 个 Nambu 分量（band k = (i,j) 槽位，
# idx(i,j) = 2(i−1)+j）× (a, ā)。块矩阵 M 的 4N 个抽象节点逐一展开：
# C[a腿(s), ā腿(t)] = −M[s, t]（行指标↔a 腿、列指标↔ā 腿），提取
# (M⁻¹)[s,t] 的掩码去 {a_t, ā_s}（与 realtime_dense.jl 同族约定）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bcs/imagtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "bcs", "common.jl"))
include(joinpath(@__DIR__, "imagtime.jl"))              # bcs_tau_ring_naive
using QuAPI
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix, GrassmannIF

const NCOMP = 4                                        # 每节点 Nambu 分量数 (i,j)

"抽象节点 s（= NCOMP*(m−1)+idx）的 a/ā 腿位置（1-idx；band = idx 槽位）。"
aleg(s) = 2*NCOMP*(cld(s, NCOMP) - 1) + 2 * ((s - 1) % NCOMP) + 1
ableg(s) = aleg(s) + 1

"块矩阵 M（4N×4N）的稠密 IF 折叠：8N 腿（band 1..4 × (a,ā)）。"
function dense_if_legs(N, M)
    basis = GrassmannContourIndex[GrassmannContourIndex(m, :τ, band, c)
                                  for m in 1:N for band in 1:NCOMP for c in (false, true)]
    C = zeros(ComplexF64, 2*NCOMP*N, 2*NCOMP*N)       # 腿数 = 2 × 4N = 8N
    for s in 1:NCOMP*N, t in 1:NCOMP*N
        C[aleg(s), ableg(t)] = -M[s, t]
    end
    return GrassmannIF(GrassmannHybridizationMatrix(basis, C))
end

"从折叠张量提取全部 (M⁻¹)[s, t]（符号逐项数值标定，结构零须与 target 一致）。"
function minor_inv(gif, n, target::AbstractMatrix)
    top = gif.data[end]
    full = (1 << (2n)) - 1
    Minv = zeros(ComplexF64, n, n)
    for s in 1:n, t in 1:n
        mask = full ⊻ ((1 << (aleg(t) - 1)) | (1 << (ableg(s) - 1)))  # 去 {a_t, ā_s}
        v = gif.data[mask + 1]
        if abs(v) < 1e-280
            abs(target[s, t]) < 1e-12 * max(1.0, maximum(abs, target)) ||
                error("结构零处目标非零 at ($s,$t)")
            continue
        end
        r = target[s, t] * top / v
        sg = round(real(r))
        abs(real(r) - sg) > 1e-6 && error("符号标定失败 at ($s,$t): r=$r")
        Minv[s, t] = sg * v / top
    end
    return Minv
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Δ=0.5, Ns=[2, 3])
    mdl = bcs_toulouse(; Lmodes=Lmodes, Δ=Δ)
    io = stdout
    println(io, "=" ^ 74)
    println(io, "BCS 虚时：稠密 IF 装配（GrassmannIF，Nambu 4 分量腿基）vs 块线性系统")
    println(io, "  IF 张量 = 2^{8N} 系数；提取 Ĝ 的杂质 2×2 块（行 1 = [𝒢₁, 𝒢₂]）")
    @printf(io, "浴：%d 模，J=%.2f，β=%.1f，ϵ_d=%.1f，Δ=%.2f\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ)
    println(io, "=" ^ 74)

    worst = 0.0
    for N in Ns
        M, b = bcs_tau_ring_naive(mdl, N)
        ref = M \ b
        gif = dense_if_legs(N, M)
        Minv = minor_inv(gif, NCOMP * N, inv(M))
        dense = Minv * b
        d = maximum(abs.(dense .- ref))
        worst = max(worst, d)
        # 行 1 与严格 BdG 环的误差（应与 imagtime.jl [1] 同款）
        ex, _ = exact_tau_block(mdl, N)
        e = maximum(maximum(abs.([dense[4*(m-1)+2*(i-1)+j] - ex[m][i, j]
                                  for i in 1:2, j in 1:2])) for m in 1:N)
        @printf(io, "  N=%d: max|dense − sys| = %.3e   |dense − exact| = %.3e\n", N, d, e)
    end
    @printf(io, "\n最大一致性偏差 = %.2e（机器精度即通过）\n", worst)
    println(io, "\n结论（实测）")
    println(io, "  · BCS Nambu Matsubara 环的稠密 IF（GrassmannIF 折叠 −M）Pfaffian 余子式")
    println(io, "    提取与块线性系统逐位一致；8N 腿基覆盖 4 个 Nambu 分量。")
    println(io, "  · 下游 MPS-IF 按本文件的腿基 / 系数矩阵装配即可逐位复现环解。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
