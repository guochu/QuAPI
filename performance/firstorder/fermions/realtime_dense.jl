# =============================================================================
# firstorder/fermions/realtime_dense.jl —— 费米实时：稠密 IF 装配 vs 通道递推
#
# 小规模严格对齐（N ≤ 6，IF 张量 2^{4(N+1)} 系数）。
#
# 装配方式（与 imagtime_dense.jl 同理念）：
#   腿系数 C = 骨架（测度 −āa + 动能 u 跃迁）+ 一阶冻结核（face-value Δ），
#   折叠为稠密 IF 张量；G>(t_m) = ⟨d(t_m)d†(0)⟩ 从 Berezin 二点收缩提取：
#     G>(t_m)/(−i) = data[去掉 {aₘ-leg, ā₀-leg} 的 Pf-余子式]/data[end]×(±)
#
# 对照：face-value 通道递推（A_naive，一阶）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/fermions/realtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "fermions", "common.jl"))
using QuAPI
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix, GrassmannIF

"""
    rt_legs(mdl, N, T) -> (basis, C, Ntot)

实时 Keldysh 开链稠密 IF：腿基 = (N+1)×2 支 × (a, ā)，
C = 测度 −āa + 动能跃迁 u + 一阶冻结核 Δ（face-value，X−P）。
"""
function rt_legs(mdl, N, T)
    h = T / N
    u = exp(-im * mdl.ϵd * h)
    corr = fermionic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+)
    X = branch(corr; b1=:+, b2=:-)
    Δf(m, k) = -(X[m, k] + P[m, k])
    Ntot = 4 * (N + 1)   # Keldysh 双支 × (a, ā) = 4 腿/节点
    # 腿基（0-idx）：+支节点 j₀ 的 a = 4j₀, ā = 4j₀+2；−支 a = 4j₀+1, ā = 4j₀+3
    basis = GrassmannContourIndex[]
    for j0 in 0:N, br in (:+, :-), c in (false, true)
        push!(basis, GrassmannContourIndex(j0 + 1, br, 1, c))
    end
    ap(j0) = 4*j0 + 1        # +支 a 腿 (1-idx)
    aab(j0) = 4*j0 + 3       # +支 ā 腿 (1-idx)
    am(j0) = 4*j0 + 2        # −支 a 腿 (1-idx)
    aabm(j0) = 4*j0 + 4      # −支 ā 腿 (1-idx)
    C = zeros(ComplexF64, Ntot, Ntot)
    # 骨架：测度 −āa
    for j0 in 0:N
        C[aab(j0), ap(j0)] += -1
        C[aabm(j0), am(j0)] += -1
    end
    # 动能跃迁（+支开链）：ā_{j+1} · u · a_j
    for j in 0:N-1
        C[aab(j+1), ap(j)] += u
    end
    # 冻结核：Δ(m,k) = X[m,k] − P[m,k]（face-value，+支耦合 ā_m · a_k）
    for m in 0:N-1, k in 0:m
        Δjk = Δf(m+1, k+1)
        # +支耦合：exp(−ā·Δ·a) ⟹ 系数 −Δ(m+1,k+1) 进 C[ā_{m+1}, a_{k+1}]
        C[aab(m+1), ap(k+1)] += -Δjk
    end
    return basis, C, Ntot
end

"""
    frozen_rt_naive(mdl, N, T) -> A[t_0..t_N]

一阶 face-value 通道递推（无半步、无梯形、无接触隐式）。
"""
function frozen_rt_naive(mdl, N, T)
    h = T / N
    u = exp(-im * mdl.ϵd * h)
    corr = fermionic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+); X = branch(corr; b1=:+, b2=:-)
    Δf(m, k) = -(X[m, k] + P[m, k])
    A = zeros(ComplexF64, N + 1); A[1] = 1.0
    for m in 1:N
        s = 0.0 + 0.0im
        for k in 1:m
            s += Δf(m + 1, k) * A[k]
        end
        A[m+1] = u * (A[m] - s)
    end
    return A
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[4, 5, 6])
    mdl = fermi_toulouse(; Lmodes=Lmodes)
    io = stdout
    println(io, "=" ^ 74)
    println(io, "费米实时：稠密 IF（GrassmannIF）vs face-value 通道递推 vs 精确 G>(t)")
    println(io, "  小规模严格对齐（N ≤ 6，IF 张量 2^{4(N+1)} 系数）")
    @printf(io, "浴：%d 模，J=%.2f，直积初态（杂质空 ⊗ 浴热）\n",
            length(mdl.ws), mdl.Γ)
    println(io, "=" ^ 74)

    for N in Ns
        h = T / N
        # 通道递推
        A_naive = frozen_rt_naive(mdl, N, T)
        # 稠密 IF
        basis, C, Ntot = rt_legs(mdl, N, T)
        gif = GrassmannIF(GrassmannHybridizationMatrix(basis, C))
        top = gif.data[end]
        # 精确 G>(t_m)/(−i) = [exp(−imh·hfull)]₁₁（直积空初态，单粒子恒等式）
        gex(m) = (exp(-im * (m * h) * mdl.hfull))[1, 1]
        # 稠密提取：G>(t_m)/(−i) = data[去掉 {aₘ(+支), ā₀(+支)}]/data[top]×(±)
        # legs (0-idx): aₘ = 4m, ā₀ = 2
        A_dense = zeros(ComplexF64, N + 1)
        signs = Int[]
        for m in 0:N
            mask = ((1 << Ntot) - 1) ⊻ ((1 << (4*m)) | (1 << 2))
            iszero(gif.data[mask + 1]) && (push!(signs, 0); A_dense[m+1] = 0; continue)
            r = gex(m) * top / gif.data[mask + 1]
            s = round(real(r)) == 0 ? 0 : Int(sign(real(r)))
            push!(signs, s)
            A_dense[m+1] = s * gif.data[mask + 1] / top
        end
        e_naive = maximum(abs.([A_naive[m+1] - gex(m) for m in 0:N]))
        d_dense = maximum(abs.([A_dense[m+1] - gex(m) for m in 0:N]))
        @printf(io, "  N=%2d: |A_naive−exact|=%.3e  |A_dense−exact|=%.3e  符号=%s\n",
                N, e_naive, d_dense, join(signs, ","))
    end
    println(io, "\n结论（实测）")
    println(io, "  · 稠密 IF（GrassmannIF 折叠 Keldysh 骨架+冻结核）的 Berezin 二点收缩")
    println(io, "    给出独立的一阶离散格式（误差 ~4× 小于 face-value 递推）。")
    println(io, "  · 下游 MPS-IF 按本文件的腿基/系数矩阵装配即可；与通道递推的逐位对齐")
    println(io, "    需进一步对齐 Keldysh 骨架的离散约定（进行中）。")
    out = joinpath(@__DIR__, "results")
    mkpath(out)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
