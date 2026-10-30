using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf
include(joinpath(@__DIR__, "offdiag_vs_ed.jl"))

println("\n=== H_S=0 下的纯 IF 误差：V 标度 + 更宽 δτ 阶梯 ===")
T = 0.8
println("[实时] σ⁺, H_S=0, N=4:")
for Gc in (V, V / 2)
    ed = ed_real(T; hs=HS0, Gc=Gc)
    pt = pt_sum_rt_hs(4, T, SP, HS0; Gc=Gc)
    @printf("  V=%.4f  ‖PT−ED‖ = %.4e\n", Gc, norm(pt - ed))
end

println("\n[虚时] σ⁺, H_S=0, Nτ 加宽（log-log 斜率）+ V 标度:")
ns = [2, 3, 4, 6, 8, 10, 12]
es = Float64[]
for N in ns
    Ze, sze = ed_imag_sz(hs=HS0)
    Zq, szq = pt_ring_imag_hs(N, SP, HS0)
    push!(es, abs(szq - sze))
    @printf("  N=%-3d  δτ=%.4f  |Δσz| = %.4e\n", N, BETA / N, es[end])
end
x = log.(Float64.(ns)); y = log.(es)
p = (length(x) * sum(x .* y) - sum(x) * sum(y)) / (length(x) * sum(x .^ 2) - sum(x)^2)
@printf("  全局拟合斜率 p = %.3f（δτ² ⟹ 2）\n", p)

println("\n[虚时] V 标度（N=8, H_S=0）:")
ev2 = Float64[]
for Gc in (V, V / 2)
    Ze, sze = ed_imag_sz(hs=HS0, Gc=Gc)
    Zq, szq = pt_ring_imag_hs(8, SP, HS0; Gc=Gc)
    push!(ev2, abs(szq - sze))
    @printf("  V=%.4f  |Δσz| = %.4e\n", Gc, ev2[end])
end
@printf("  V 减半残差比 = %.1f（V⁴ ⟹ 16，V² ⟹ 4）\n", ev2[1] / ev2[2])
