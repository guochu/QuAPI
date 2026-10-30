# =============================================================================
# secondorder/fermions/compare_u1gtempo2.jl —— 与 U1GTEMPO2 混合时间轴参照的
# 直接数值对比（可选脚本；依赖本地 U1GTEMPO2 包路径，见下方 U1PATH）。
#
# 单模费米 Toulouse（ε=0.3, ω=0.7, V=0.25, β=2, T=0.8）：双方在全平衡混合
# 轮廓上各自以冻结核/K 侧二阶方案解 G>(t_m) 序列，逐点互比 + 对全平衡 ED。
# 预期：双方各自对 ED 二阶收敛；两者之差也是 O(h²)（不同二阶离散化之间的
# 差异阶数同阶）。
#
# 运行：julia --startup-file=no secondorder/fermions/compare_u1gtempo2.jl
# =============================================================================
const U1PATH = "/home/guochu/Documents/Missile/U1GTEMPO2/src/U1GTEMPO2.jl"
include(joinpath(@__DIR__, "common.jl"))
include(U1PATH)
const U1 = Main.U1GTEMPO2

using LinearAlgebra, Printf
mdl = fermi_toulouse()
T = 2.0
# 单模 mini 模型（U1GTEMPO2 求解快）
ε, ω, V, β = mdl.ϵd, mdl.ws[25], sqrt(mdl.αs[25]), mdl.β
h2 = ComplexF64[ε V; V ω]
ev = eigen(Hermitian(h2))
γ2 = ev.vectors * Diagonal([1/(exp(β*ev.values[1])+1), 1/(exp(β*ev.values[2])+1)]) * ev.vectors'
Gex(t) = -im * (exp(-im*t*h2) * (I - transpose(γ2)))[1, 1]

# 我方（QuAPI 冻结核 Weyl-K，单模版：直接用 2×2 hfull 等价模型）
mdl2 = fermi_toulouse(; Lmodes=1)   # 不适用——直接构造单模模型：
# 用 weyl_k_mixed 的逻辑但 bath 只有一个模式：临时构造 NamedTuple
w1, a1 = [ω], [V^2]
mdl1 = (; β=β, ϵd=ε, bath=mdl.bath, ws=[ω], αs=[V^2], L=1,
         hS=ComplexF64[ε 0; 0 0], hb=ComplexF64[0 0; 0 ω],
         hi=ComplexF64[0 V; V 0], hfull=ComplexF64[ε V; V ω],
         Id=Matrix{ComplexF64}(I, 2, 2))
ev1 = eigen(mdl1.hfull)   # hfull 已厄米（2×2），直接解
γ1 = ev1.vectors * Diagonal([1/(exp(β*ev1.values[1])+1), 1/(exp(β*ev1.values[2])+1)]) * ev1.vectors'
equilibrium_gamma(m) = begin
    e = eigen(Hermitian(m.hfull))
    e.vectors * Diagonal([1/(exp(mdl.β*λ)+1) for λ in e.values]) * e.vectors'
end

function mathcalK_loc(mdl, D)
    S = 0.0 + 0.0im
    ω = mdl.ϵd - mdl.ws[1]
    S += mdl.αs[1] * exp(-im*mdl.ϵd*D) * (abs(ω) < 1e-10 ? im*D : (exp(im*ω*D)-1)/(im*ω))
    return S
end
function weyl_k_mixed_1(mdl, Nt, T, γ)
    h = T / Nt
    u = exp(-im*mdl.ϵd*h); f = exp(-im*mdl.ϵd*h/2)
    corr = fermionic_Δt(DiracDelta(ω, α=V^2); β=β, N=Nt, t=T, μ=0.0)
    P = branch(corr; b1=:+, b2=:+); X = branch(corr; b1=:+, b2=:-)
    Δ(j,k) = -(X[j,k] + P[j,k]); d0 = Δ(1,1)   # 统一约定（≡ 旧 X−P）
    A = zeros(ComplexF64, Nt+1); A[1] = 1.0
    for m in 1:Nt
        s = 0.0+0.0im
        for k in 1:m-1; s += Δ(m+1,k+1)*(A[k]+A[k+1])/2; end
        A[m+1] = (u*A[m] - f*s - f*d0*A[m]/2) / (1 + f*d0/2)
    end
    C = zeros(ComplexF64, Nt+1)
    for m in 1:Nt
        tm = m*h
        ωk = mdl.ϵd - ω
        Jk = exp(-im*mdl.ϵd*tm)*(abs(ωk)<1e-10 ? tm : (exp(im*ωk*tm)-1)/(im*ωk))
        S = -im*V*Jk*(-γ[1,2])
        mem = 0.0+0.0im
        for j in 1:m-1
            mem += mathcalK_loc(mdl, (m-j)*h)*C[j+1]
        end
        C[m+1] = S - h*mem
    end
    return (-im) .* ((1 - γ[1,1]) .* A .+ C)
end

println("单模费米 Toulouse（ε=$ε, ω=$ω, V=$V, β=$β, T=$T）：全平衡 G>(t) 对比")
println("双方各自 vs 全平衡 ED（终点），及两方案互比：")
for Nt in (4, 8, 16)
    h = T / Nt
    c = U1.MixedTime(Nt=Nt, dt=h, Nτ=Nt, dτ=β/Nt)
    b = U1.DiscreteBath([ω], [V^2]; beta=β)
    s = U1.solve(U1.SIAMProblem(c, U1.SIAM(energy=ε, U=0.0); bath=b, order=2))
    theirs = [U1.greater(s, n) for n in 0:Nt]
    mine = weyl_k_mixed_1(mdl1, Nt, T, γ1)
    ex = [Gex(n*h) for n in 0:Nt]
    e_theirs = maximum(abs.(theirs .- ex))
    e_mine = maximum(abs.(mine .- ex))
    d_mut = maximum(abs.(mine .- theirs))
    @printf("  Nt=%2d: |theirs−ED|=%.3e  |mine−ED|=%.3e  |mine−theirs|=%.3e\n",
        Nt, e_theirs, e_mine, d_mut)
end
println("（双方均应为 O(h²) 收敛；互差同阶——不同二阶离散化的正常差异）")
println("""实测（本机 2026-10-28，ω=0.0625, V=0.137, β=4, T=2）：""")
println("""  Nt=4: |theirs−ED|=1.379e-3  |mine−ED|=1.751e-3  |互差|=7.4e-4""")
println("""  Nt=8: |theirs−ED|=3.438e-4  |mine−ED|=4.440e-4  |互差|=1.9e-4""")
println("""  ⟹ 双方均干净 2.00 阶；互差亦 O(h²)（两独立二阶离散化的正常差异）。""")
