# =============================================================================
# benchmark/bosons/offdiag_tempo.jl
#
# 原计划：用 TEMPO 包的 PT-MPO 做 σ⁺ 非对角耦合的长时间 QuAPI 离散化标度
# 测试（offdiag_order.jl 的大 N 推广：实时固定 T 扫 δt、虚时固定 β 扫 δτ，
# 检验 V⁴δt² 是否在多步数下保持）。
#
# 实际结果：该测试目前被一个 **TEMPO NonDiagonalHyb 实时组装的约定问题阻塞**，
# 本文件保留为诊断报告 + 阻塞证据：
#
# [0] G=0 对照：PT 复现自由演化到机器精度 ⟹ 门/边界/观测管线正确。
# [1] σx 对照（NonAdditiveHyb + PartialIF）：|PT − Lie| = 1e-12 ⟹ 复现
#     second_order_check_nondiag.jl 的恒等式，管线与源文件等价。
# [2] σ⁺（NonDiagonalHyb）：与 σ⁺⊗b+σ⁻⊗b† 的 Lie 链之差 = 0.2926，
#     对 IF 构造路线（TDVP2IF+ExactConstructor / XTRGIF+Exact / XTRGIF+Prony）
#     和 TDVP 流步长 δ 完全不变 ⟹ 不是流/截断误差，是组装内容的系统差。
#     与通道交换链（σ⁺⊗b†+σ⁻⊗b）差 V²δt¹（随 N 以 ratio 2.00 下降），
#     即 PT 收敛到一个既非原链也非交换链的连续极限；与原链之差随 N 饱和。
# [3] 佐证：QuAPI 的同一 IF 张量（BosonicIF split 核，crosscheck 逐位一致）
#     用 offdiag_order.jl 的门约定收缩后，与 σ⁺-Lie 链只差 V⁴δt² ——
#     张量内容正确，问题出在 TEMPO 侧 NonDiagonalHyb 的组装/配对。
# [4] 阻塞：长时间 V⁴δt² 标度测试需先在 TEMPO 包内解决该约定。
#     参照：InfiniteTEMPO 的 freebosons 实时测试（独立实现）以 ~6e-4
#     通过 Toulouse 严格解，可作为正确行为的对照实现。
#
# 运行（用 TEMPO 的环境）：
#   JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no \
#       --project=/home/guochu/Documents/Missile/TEMPO \
#       benchmark/bosons/offdiag_tempo.jl
# =============================================================================
using TEMPO, ImpurityModelBase, LinearAlgebra, Printf

BLAS.set_num_threads(4)

# ------------------------------ 模型（与 offdiag_order.jl 完全一致） ------------------------------ #
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]
const HS = .37X + .19Y + .23Z
const SP = ComplexF64[0 1; 0 0]
const SM = ComplexF64[0 0; 1 0]
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5Z) / 2
const OMEGA, V, BETA = 1.3, 0.35, 1.7
const NB = 28

coup(Gc, swap=false) = Gc * (swap ? kron(SM, diagm(1 => sqrt.(1.0:(NB-1)))) +
                                       kron(SP, diagm(-1 => sqrt.(1.0:(NB-1)))) :
                                    kron(SP, diagm(1 => sqrt.(1.0:(NB-1)))) +
                                       kron(SM, diagm(-1 => sqrt.(1.0:(NB-1)))))
therm() = (p = exp.(-BETA * OMEGA * (0.0:(NB-1))); p ./ sum(p))
ptrace2(R) = [sum(reshape(R, NB, 2, NB, 2)[i,a,i,b] for i in 1:NB) for a in 1:2, b in 1:2]

"严格 Lie 链（实时）：U = (A·B)^n，B 精确含耦合 cell。"
function strict_real(n, T; Gc=V, swap=false)
    h = T / n
    hb = OMEGA * Diagonal(0.0:(NB-1))
    A = kron(exp(-im * h * HS), Matrix(I, NB, NB))
    B = exp(-im * h * (kron(Matrix(I, 2, 2), Matrix(hb)) + coup(Gc, swap)))
    U = (A * B)^n
    ptrace2(U * kron(RHO0, Diagonal(therm())) * U')
end

"ED 精确（对角化全空间）。"
function ed_real(T; nb=NB)
    hb = OMEGA * Diagonal(0.0:(nb-1))
    H = kron(HS, Matrix(I, nb, nb)) + kron(Matrix(I, 2, 2), Matrix(hb)) + coup(V)
    ev, Vv = eigen(Hermitian(H))
    U = Vv * Diagonal(exp.(-im * T * ev)) * Vv'
    p = exp.(-BETA * OMEGA * (0.0:(nb-1))); p ./= sum(p)
    [sum(reshape(U * kron(RHO0, Diagonal(p)) * U', nb, 2, nb, 2)[i,a,i,b] for i in 1:nb) for a in 1:2, b in 1:2]
end

# ------------------------------ TEMPO PT ------------------------------ #
_trc(D) = truncdimcutoff(D=D, ϵ=1.0e-13, add_back=0)

"route: :partial（仅 NonAdditiveHyb）| :tdvp2（TDVP2IF+ExactConstructor）| :xtrg（XTRGIF+ExactConstructor）"
function pt_real(n, T, hybstyle; Gc=V, D=256, route=:tdvp2, δ=0.05)
    lat = PTLattice(N=n, δt=T / n, contour=:real)
    bath = bosonicbath(DiracDelta(OMEGA, α=Gc^2), β=BETA)
    corr = correlationfunction(bath, lat)
    mi = route === :partial ? hybriddynamics(lat, corr, hybstyle, trunc=_trc(D)) :
         route === :tdvp2 ? hybriddynamics(lat, corr, hybstyle,
              TDVP2IF(alg_expan=ExactConstructor(), trunc=_trc(D), δ=δ)) :
         hybriddynamics(lat, corr, hybstyle,
              XTRGIF(k=max(ceil(Int, log2(n)), 3), fast=true,
                     alg_mult=DMRG1(trunc=_trc(D)), alg_expan=ExactConstructor()))
    k = sysdynamics(lat, ImpurityHamiltonian(HS), trunc=_trc(D))
    m = mult(k, mi, trunc=_trc(D))
    initialstate!(m, lat, RHO0)
    return rdm(lat, m), bonddim(mi)
end

fit_order(dts, errs) = begin
    x = log.(dts); y = log.(errs)
    (length(x) * sum(x .* y) - sum(x) * sum(y)) / (length(x) * sum(x .^ 2) - sum(x)^2)
end

# ------------------------------- 主流程 ------------------------------- #
function main()
    t0 = time()
    io = stdout
    n0, T = 8, 1.6
    println(io, "=" ^ 78)
    println(io, "σ⁺ 非对角耦合 × TEMPO PT：长时间标度测试的诊断（详见文件头）")
    @printf(io, "H_S 一般非对角；浴：单模 ω=%.2f, V=%.2f, β=%.2f, T=%.1f\n", OMEGA, V, BETA, T)
    println(io, "=" ^ 78)

    # ---------- [0] G=0 对照 ----------
    println(io, "\n[0] G=0（N=$n0）：PT vs 自由 Lie 链")
    r, _ = pt_real(n0, T, NonDiagonalHyb(SP); Gc=0.0)
    @printf(io, "    |PT − Lie| = %.2e（应机器精度）\n", norm(r - strict_real(n0, T; Gc=0.0)))

    # ---------- [1] σx 对照（管线验证） ----------
    println(io, "\n[1] σx（NonAdditiveHyb，PartialIF）：恒等式复现")
    for route in (:partial, :tdvp2, :xtrg)
        route === :partial || continue
        rx, _ = pt_real(n0, T, NonAdditiveHyb(X); route=route)
        hb = OMEGA * Diagonal(0.0:(NB-1)); b = diagm(1 => sqrt.(1.0:(NB-1)))
        h = T / n0
        A = kron(exp(-im*h*HS), Matrix(I, NB, NB))
        B = exp(-im*h*(kron(Matrix(I,2,2), Matrix(hb)) + V*kron(X, b+b')))
        U = (A*B)^n0
        rl = ptrace2(U * kron(RHO0, Diagonal(therm())) * U')
        @printf(io, "    PartialIF   |PT − Lie| = %.2e（secondorder 恒等式 ✓）\n", norm(rx - rl))
    end
    for (tag, route) in (("TDVP2+Exact ", :tdvp2), ("XTRGIF+Exact ", :xtrg))
        rx, _ = pt_real(n0, T, NonAdditiveHyb(X); route=route)
        hb = OMEGA * Diagonal(0.0:(NB-1)); b = diagm(1 => sqrt.(1.0:(NB-1)))
        h = T / n0
        A = kron(exp(-im*h*HS), Matrix(I, NB, NB))
        B = exp(-im*h*(kron(Matrix(I,2,2), Matrix(hb)) + V*kron(X, b+b')))
        U = (A*B)^n0
        rl = ptrace2(U * kron(RHO0, Diagonal(therm())) * U')
        @printf(io, "    %s|PT − Lie| = %.2e（该路线自身精度水平）\n", tag, norm(rx - rl))
    end

    # ---------- [2] σ⁺ 诊断：路线/δ 不变性 ----------
    println(io, "\n[2] σ⁺（NonDiagonalHyb）vs σ⁺-Lie 链（σ⁺⊗b+σ⁻⊗b†）：路线与 δ 不变性")
    for (tag, kw) in (("TDVP2 δ=0.05 ", (route=:tdvp2,)),
                      ("TDVP2 δ=0.01 ", (route=:tdvp2, δ=0.01)),
                      ("XTRGIF+Exact ", (route=:xtrg,)))
        r, bond = pt_real(n0, T, NonDiagonalHyb(SP); kw...)
        @printf(io, "    %s |Δ| = %.4e  (bond=%d)\n", tag, norm(r - strict_real(n0, T)), bond)
    end

    # ---------- [3] σ⁺ vs 交换链：δt¹ 收敛（收敛到不同连续极限） ----------
    println(io, "\n[3] σ⁺ vs 交换链（σ⁺⊗b†+σ⁻⊗b）与 vs ED：随 N 的行为")
    local prev = 0.0
    for n in (8, 16, 32)
        r, _ = pt_real(n, T, NonDiagonalHyb(SP))
        eSwap = norm(r - strict_real(n, T; swap=true))
        eED = norm(r - ed_real(T))
        @printf(io, "    N=%-3d  |PT−swap|=%.4e%s   |PT−ED|=%.4e（饱和）\n", n, eSwap,
                prev > 0 ? @sprintf("  ratio=%.2f", prev / eSwap) : "", eED)
        prev = eSwap
    end

    println(io, "\n[4] 结论：长时间 V⁴δt² 标度测试被阻塞 —— NonDiagonalHyb 实时组装收敛到")
    println(io, "    既非 σ⁺b+σ⁻b† 也非其通道交换的连续极限；IF 张量本身正确（QuAPI 侧")
    println(io, "    V⁴δt²）。需先在 TEMPO 包内排查 NonDiagonalHyb 的组装/配对（对照：")
    println(io, "    InfiniteTEMPO freebosons 实时测试以 ~6e-4 通过）。")
    println(io, "\n总耗时 ", round(time() - t0, digits=1), " s")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
