# HO 杂质（op1 = a† ≠ op2 = a，[a†,a] = −1 c 数换位子）—— PT 严格性检验
#
# 结论（2026-10 数值）：
# · 校准 L=2 hyb=σx（op1=op2）：‖链−PT‖ ~1e-10（bath 截断 nb 极限）✓ 机器一致；
#   hyb=a†(L=2≡σ⁺ JC)：G⁴ 标度（比值 14.6≈16）✓ 复现 nondiagcoupling 附录。
# · HO 大 L（6→12）：残差收敛到结构平台（L=12: 4.7e-4@G=0.4 / 3.1e-5@G=0.2，
#   比值 15.4≈16）——**非截断伪影**：当前包的 PT 对 op1≠op2 一律只有 O(G²) 精确，
#   a†/a 的 c 数换位子情形也不例外。
# · 数学注记：a†/a 情形 IF 恰为双线性（cell = 位移 D(γ(τ)a) × 算符相位 e^{ic(τ)a†a}，
#   c 数封闭 ⟹ 无 σ⁺ 型原则障碍）——G⁴ 残差是实现缺口（疑缺 op1≠op2 分支公式的
#   高阶双线性内容，如对角算符相位的完整闭式；该相位在 op1=op2 时退化为 c 数，
#   故从未暴露）。玻色侧可修的真实缺口。
include(joinpath(@__DIR__, "nondiagcoupling.jl"))   # guard 保护，不跑 main
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf

function pt_sum_L(Tarr, L, HS, ρ0, lsys, lbath; ltail = 0.0)
    M = length(lbath)
    Afull = [exp(-im * lsys[k] * HS) for k in 1:M]
    Utail = exp(-im * ltail * HS)
    ρ = zeros(ComplexF64, L, L)
    inner = L^(2M + 1)
    for cfgP in 0:(inner-1)
        vP = [((cfgP ÷ L^(k-1)) % L) + 1 for k in 1:2M+1]
        amp = one(ComplexF64)
        for k in 1:M; amp *= Afull[k][vP[2k], vP[2k-1]]; end
        for cfgM in 0:(inner-1)
            vM = [((cfgM ÷ L^(k-1)) % L) + 1 for k in 1:2M+1]
            amp2 = one(ComplexF64)
            for k in 1:M; amp2 *= Afull[k][vM[2k], vM[2k-1]]; end
            for a in 1:L, b in 1:L
                f = Utail[a, vP[2M+1]] * conj(Utail[b, vM[2M+1]])
                iszero(f) && continue
                inds = Int[]
                for k in 1:M
                    push!(inds, vP[2k], k < M ? vP[2k+1] : vP[2M+1])
                    push!(inds, vM[2k], k < M ? vM[2k+1] : vM[2M+1])
                end
                ρ[a, b] += ρ0[vP[1], vM[1]] * amp * conj(amp2) * f * Tarr[inds...]
            end
        end
    end
    return ρ
end

# mode = :spinx（L=2, hyb=σx, 耦合 G σx⊗(b+b†)）| :adag（hyb=a†, 耦合 G(a†⊗b+a⊗b†)）
function ho_pt_check(; mode = :adag, nmax = 3, Ω = 0.8, ω = 1.0, G = 0.4, β = 1.0,
                     N = 2, T = 0.9, nb = 8)
    L = nmax + 1
    a = Matrix(diagm(1 => sqrt.(1.0:nmax)))
    ad = a'
    HS = Matrix(Ω * (ad * a))
    ρ0 = exp(-HS + 0.5 * (a + ad)); ρ0 ./= tr(ρ0)
    bth = Matrix(diagm(1 => sqrt.(1.0:(nb-1))))
    hb = ω * Diagonal(0.0:(nb-1))
    if mode === :spinx
        hyb = [0.0 1.0; 1.0 0.0]
        coup = G * (kron(hyb, bth) + kron(hyb, bth'))
        HS = [.23 .31; .31 -.17]              # 一般非对角 H_S（同 nondiagcoupling 精神）
    else
        hyb = ad
        coup = G * (kron(ad, bth) + kron(a, bth'))
    end
    A(τ) = kron(exp(-im * τ * HS), Matrix(I, nb, nb))
    B(τ) = exp(-im * τ * (kron(Matrix(I, L, L), Matrix(hb)) + coup))
    rb = exp.(-β * ω * (0.0:(nb-1))); rb ./= sum(rb)
    dt = T / N
    U = A(dt)
    for j in 1:N
        U = B(dt) * U
        U = A(j == N ? 0.0 : dt) * U        # 末段系统演化长度 0（同 chain_specs order=1）
    end
    R = U * kron(ρ0, Diagonal(rb)) * U'
    R4 = reshape(R, nb, L, nb, L)
    ρch = [sum(R4[α, i, α, j] for α in 1:nb) for i in 1:L, j in 1:L]
    corr = bosonic_Δt(DiracDelta(ω, α = G^2), β = β, N = N, t = T, μ = 0.0)
    hmat = BosonicHybridizationMatrix(corr; split = true)
    keep = trues(2 * (N + 1)); keep[1] = keep[2] = false
    idx = findall(keep)
    h2 = BosonicHybridizationMatrix(pt_basis(N),
                                    QuAPI.Dense.BosonicSplitDelta(hmat.Δ.M1[idx, idx],
                                                                  hmat.Δ.M2[idx, idx]))
    P = BosonicIF(h2, hyb)
    Tarr = reshape(P.data, ntuple(_ -> L, 2 * length(P.basis)))
    ρpt = pt_sum_L(Tarr, L, HS, ρ0, fill(dt, N), fill(dt, N))
    return norm(ρch .- ρpt), norm(ρch)
end

println("== 校准 L=2（nb=20）==")
println("-- hyb = σx（op1=op2，应 ~1e-13）--")
for (N, G) in ((2, 0.4), (2, 0.2), (3, 0.4))
    d, nρ = ho_pt_check(; mode = :spinx, nmax = 1, N = N, G = G, nb = 20)
    @printf("  N=%d G=%.2f:  ‖链−PT‖ = %.3e  (‖ρ‖=%.3f)\n", N, G, d, nρ)
end
println("-- hyb = a† (L=2 ≡ σ⁺ JC，应 G⁴ 标度 比值≈16) --")
for G in (0.4, 0.2)
    d, nρ = ho_pt_check(; mode = :adag, nmax = 1, N = 2, G = G, nb = 20)
    @printf("  N=2 G=%.2f:  ‖链−PT‖ = %.3e  (‖ρ‖=%.3f)\n", G, d, nρ)
end
println("== HO：L 扫描（残差随 L 收敛到平台 ⟹ 结构性；随 L 消失 ⟹ 截断伪影）==")
for nmax in (5, 7, 9, 11)
    for G in (0.4, 0.2)
        d, nρ = ho_pt_check(; mode = :adag, nmax = nmax, N = 1, G = G, nb = 20)
        @printf("  L=%d N=1 G=%.2f:  ‖链−PT‖ = %.3e  (‖ρ‖=%.3f)\n", nmax + 1, G, d, nρ)
    end
end
