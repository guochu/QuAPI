# 探针：提取 PT 隐含的单粒子配对 Σ_eff = A0⁻¹ − A_PT⁻¹，与 Δt 核对照
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf
include(joinpath(@__DIR__, "offdiag_toulouse.jl"))

h = 0.3
Nmax = 3
# A_PT(t_j)：逐前缀跑 PT（N'=j 个 cell，T'=j·h）
Apt = ones(ComplexF64, Nmax + 1)
for j in 1:Nmax
    rhoG = zeros(ComplexF64, DF, DF); rhoG[2, 1] = 1.0
    Apt[j + 1] = tr(pt_sum_rt(j, j * h, adag(); rho0=rhoG) * adown())
end
Apt_mat = zeros(ComplexF64, Nmax + 1, Nmax + 1)
for j in 0:Nmax, k in 0:j
    Apt_mat[j + 1, k + 1] = Apt[j - k + 1]    # 平移不变：A(t_j,t_k)=A((j-k)h)
end
# A0⁻¹（节点自由传播子之逆）
u = exp(-im * EPSD * h)
A0inv = Matrix{ComplexF64}(I, Nmax + 1, Nmax + 1)
for j in 1:Nmax; A0inv[j + 1, j] = -u; end
Σeff = A0inv - inv(Apt_mat)
println("Σ_eff = A0⁻¹ − A_PT⁻¹（下三角）：")
for j in 1:Nmax+1
    print("  row ", j - 1, ": ")
    for k in 1:j
        print(round(Σeff[j, k], digits=6), "  ")
    end
    println()
end
println()
println("对照：Δt 闭式核 −(Δ⁺⁺+Δ⁺₋)（lag 含等时）：")
KR = -KR_closed(Nmax, Nmax * h)   # 取负号与 realtime_frozen 的 K^R=−(Δ⁺⁺+Δ⁺₋) 对齐
for j in 1:Nmax+1
    print("  row ", j - 1, ": ")
    for k in 1:j
        print(round(KR[j, k], digits=6), "  ")
    end
    println()
end
