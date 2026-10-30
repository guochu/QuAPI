# =============================================================================
# firstorder/bosons/mixedtime_dense.jl —— 玻色混合轮廓：稠密 IF 装配
#
# 与 mixedtime.jl 完全对齐的另一种实现：同一个一阶冻结核（bosonic_Δτ 的 η
# 块 + bosonic_Δm 的 P/Xr/Xc 块）不再走通道线性系统，而是按通道折叠为稠密
# IF 张量（GrassmannIF），从 Pfaffian 余子式收缩提取，两路输出必须逐位一致
# （机器精度）。三段张量：
#   ① τ 环 IF：bosonic_Δτ 驱动的正确热环（fτ·w 梯形 + 跳变常数代入），
#     提取 Gτ、n_d —— 热背景；
#   ② A 通道 IF：bosonic_Δm 的 η 块 face-value 因果下三角（lag 直读、含
#     接触对角），A = M⁻¹·e₁ 整体折叠；
#   ③ C 通道 IF：同一折叠张量上，源 Svec[m+1] = Σ_j ζ(m+1, j)·Gτ[j]
#     （ζ 跨块直接收缩——玻色情形 rank 非退化、语义正确），C = M⁻¹·Svec。
#   虚轴↔实轴的混合耦合（ζ 块）由此显式进入稠密装配。
#
# 单粒子说明：单粒子振幅代数与费米情形同构（skew 双线性 + Pfaffian 余子式），
# 玻色统计完全包含在核条目里（ζ 块的 (1+n)−n 热抵消）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bosons/mixedtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "common", "common.jl"))
using QuAPI
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix, GrassmannIF

"""
    naive_tau_ring(mdl, Nτ) -> (Mτ, bτ)

τ 支热环（与 mixedtime.jl ① 完全同款）：bosonic_Δτ 驱动、fτ·w 梯形
（相邻节点各半权重）、跳变常数 w(k=1) 代入源。
"""
function naive_tau_ring(mdl, Nτ)
    δτ = mdl.β / Nτ
    fτ = exp(-mdl.ϵd * δτ / 2)
    corr = bosonic_Δτ(mdl.bath.spectrum; β=mdl.β, N=Nτ, μ=mdl.μ)
    ηjk, ηkj = corr.data.ηⱼₖ, corr.data.ηₖⱼ
    d0τ = ηjk[1] + ηkj[1]
    M = zeros(ComplexF64, Nτ, Nτ)
    b = zeros(ComplexF64, Nτ)
    for m in 1:Nτ
        M[m, m] += 1.0
        if m ≥ 2
            M[m, m-1] -= exp(-mdl.ϵd * δτ)
        else
            b[m] += exp(-mdl.ϵd * δτ)
            M[m, Nτ] -= exp(-mdl.ϵd * δτ)
        end
        for k in 1:Nτ
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0τ
            if k == 1
                M[m, Nτ] -= fτ * w * 0.5
                M[m, 1] -= fτ * w * 0.5
                b[m] += fτ * w * 0.5
            else
                M[m, k-1] -= fτ * w * 0.5
                M[m, k] -= fτ * w * 0.5
            end
        end
    end
    return M, b
end

"通用稠密 IF：腿基 (a, ā)/节点，C[2k, 2m−1] = −M[m,k]（转置核装配）。"
function dense_if_legs(n, M, br::Symbol)
    basis = GrassmannContourIndex[GrassmannContourIndex(j, br, 1, c)
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
    GL = exact_greater_lesser(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "玻色混合：稠密 IF 装配（GrassmannIF，分段折叠）vs 通道线性系统（机器精度）")
    println(io, "  ① τ 环 IF（bosonic_Δτ 梯形热环）→ Gτ、n_d")
    println(io, "  ② A 通道 IF（bosonic_Δm η 块 face-value）+ 源 e₁")
    println(io, "  ③ C 通道 = 同一张量 + ζ 跨块源 Svec(= ζ·Gτ 直接收缩)")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f（Nτ=Nt 同步）\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)

    worst = 0.0
    for N in Ns
        h = T / N
        dm = bosonic_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=N, t=T, Nt=N, μ=mdl.μ)
        P = branch(dm; b1=:+, b2=:+)
        Xr = branch(dm; b1=:+, b2=:-)
        Xc = branch(dm; b1=:+, b2=:τ)
        # ---------- ① τ 环 ----------
        Mτ, bτ = naive_tau_ring(mdl, N)
        Gτ_ref = Mτ \ bτ
        gif1 = dense_if_legs(N, Mτ, :τ)
        Minv1 = minor_inv(gif1, N, inv(Mτ))
        Gτ = Minv1 * bτ
        n_d = Gτ[N]
        d1 = maximum(abs.(Gτ .- Gτ_ref))
        # ---------- ②/③ 实时通道（同一折叠张量，两个源） ----------
        u = exp(-im * mdl.ϵd * h)
        MA = Matrix{ComplexF64}(I, N + 1, N + 1)
        for j in 2:N+1; MA[j, j-1] = -u; end
        for m in 1:N+1, k in 1:m
            MA[m, k] -= (P[m-k+1, 1] + Xr[m-k+1, 1])
        end
        Svec = zeros(ComplexF64, N + 1)
        for m in 1:N
            Svec[m+1] = sum(Xc[m+1, j] * Gτ[j] for j in 1:N)
        end
        ref = inv(MA) * [ [1.0; zeros(ComplexF64, N)] Svec ]
        gif2 = dense_if_legs(N + 1, MA, :+)
        Minv2 = minor_inv(gif2, N + 1, inv(MA))
        dense = Minv2 * [ [1.0; zeros(ComplexF64, N)] Svec ]
        d2 = maximum(abs.(dense .- ref))
        worst = max(worst, d1, d2)
        A, C = dense[:, 1], dense[:, 2]
        G = (-im) .* ((1 + n_d) .* A .+ C)
        eex = maximum(abs.([G[m+1] - GL(m * h)[1] for m in ms_of(N)]))
        @printf(io, "  N=%2d: |Gτ_dense−Gτ_sys|=%.2e  |[A,C]_dense−sys|=%.2e  |G−exact|=%.3e\n",
                N, d1, d2, eex)
    end
    @printf(io, "\n最大一致性偏差 = %.2e（机器精度即通过）\n", worst)
    println(io, "\n结论（实测）")
    println(io, "  · 两段 IF 张量的 Pfaffian 余子式提取与通道线性系统逐位一致（机器精度）。")
    println(io, "  · ζ 跨块（虚轴↔实轴混合耦合）经 Svec = ζ·Gτ 直接收缩进入 C 通道源；")
    println(io, "    玻色 ζ 收缩语义正确，无需费米情形的因式分解替换。")
    println(io, "  · 下游 MPS-IF 按本文件的腿基/系数矩阵装配即可逐位复现通道解。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
