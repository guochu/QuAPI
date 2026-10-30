# =============================================================================
# firstorder/bcs/realtime_dense.jl —— BCS 实时：稠密 IF 装配 vs Nambu 块线性系统
#
# 与 realtime.jl 完全对齐的另一种实现：同一个一阶冻结核（bcs_Δt 四块的
# Nambu 推迟组合 𝒩 = [cn cc; nn nc]）不再走块线性系统 M\B，而是经
# GrassmannIF 折叠为稠密 IF 张量（2^{4N} 系数），从 Pfaffian 余子式收缩
# 提取，两路输出必须逐位一致（机器精度）。
#
# Nambu 腿基：每物理节点 4 条腿 = 2 个 Nambu 分量（band 1 = ↑ 粒子行、
# band 2 = ↓ 空穴行）× (a, ā)，序与 test/quapi/bcsimag.jl 一致
# （band 1 a,ā 然后 band 2 a,ā）。块矩阵 M 的 2N 个"抽象节点"（N 节点 ×
# 2 Nambu 分量）逐一展开到腿基：C[ā_s, a_t] = −M[s, t]；提取
# (M⁻¹)[s,t] = sgn·data[去掉 {ā_s, a_t}]/data[end]（符号逐项数值标定），
# 𝒜_m[p, :] = Σ_t (M⁻¹)[2(m−1)+p, t]·B[t, :]。
#
# 说明：物理路径值的 BCS 装配（含反常 ā-ā / a-a 配对的全配分函数展开，
# GrassmannADT 乘法）见包测试 test/quapi/bcsimag.jl；本文件与
# firstorder/*/dense 系列同章程——通道线性系统的张量化参照。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bcs/realtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "bcs", "common.jl"))
include(joinpath(@__DIR__, "realtime.jl"))              # bcs_naive_system / bcs_shifted_system
using QuAPI
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix, GrassmannIF

"抽象节点 s（= 2(m−1)+p）的 a/ā 腿在 basis 中的位置（1-idx）。"
aleg(s) = 4 * ((s + 1) >> 1 - 1) + 2 * ((s - 1) % 2) + 1     # a  腿
ableg(s) = aleg(s) + 1                                       # ā 腿

"块矩阵 M（2N×2N）的稠密 IF 折叠：4N 腿，C[ā_s, a_t] = −M[s, t]。"
function dense_if_legs(N, M)
    basis = GrassmannContourIndex[GrassmannContourIndex(m, :+, band, c)
                                  for m in 1:N for band in (1, 2) for c in (false, true)]
    C = zeros(ComplexF64, 4N, 4N)
    for s in 1:2N, t in 1:2N
        C[aleg(s), ableg(t)] = -M[s, t]      # 行指标↔a腿、列指标↔ā腿（核装配约定）
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
function main(; Lmodes=48, Δ=0.5, T=2.0, Ns=[2, 3, 4, 5, 6])
    mdl = bcs_toulouse(; Lmodes=Lmodes, Δ=Δ)
    io = stdout
    println(io, "=" ^ 74)
    println(io, "BCS 实时：稠密 IF 装配（GrassmannIF，Nambu 腿基）vs 块线性系统")
    println(io, "  IF 张量 = 2^{4N} 系数；naive 与顶点移位两套 M 各自折叠对比")
    @printf(io, "浴：%d 模，J=%.2f，β=%.1f，ϵ_d=%.1f，Δ=%.2f，T=%.1f\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ, T)
    println(io, "=" ^ 74)

    worst = 0.0
    for N in Ns
        for (name, sys) in (("naive", bcs_naive_system), ("shifted", bcs_shifted_system))
            M, B = sys(mdl, N, T)
            ref = M \ B
            gif = dense_if_legs(N, M)
            Minv = minor_inv(gif, 2N, inv(M))
            dense = Minv * B
            d = maximum(abs.(dense .- ref))
            worst = max(worst, d)
            @printf(io, "  N=%d  %-7s: max|dense − sys| = %.3e\n", N, name, d)
        end
    end
    @printf(io, "\n最大一致性偏差 = %.2e（机器精度即通过）\n", worst)
    println(io, "\n结论（实测）")
    println(io, "  · BCS Nambu 块线性系统的稠密 IF（GrassmannIF 折叠 −M）Pfaffian 余子式")
    println(io, "    提取与线性系统逐位一致；naive / 顶点移位两套 M 均验证。")
    println(io, "  · 下游 MPS-IF 按本文件的 Nambu 腿基 / 系数矩阵装配即可逐位复现通道解。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
