# =============================================================================
# benchmark/bosons/offdiag_vs_ed.jl
#
# 补上 offdiag_order.jl 缺失的一环：把「基于 BosonicIF 的路径求和装配」直接
# 与 ED（全空间对角化，2*NB = 56 维）对比，判断 V⁴δt² 是 IF 的真实截断误差，
# 还是我自己装配（腿序 / ρ₀ 绑定 / A 因子位置）的假象。
#
# 关键判据 —— 令 H_S = 0：
#   A = exp(-i h H_S) = I  ⟹  链 U = (B·A)^N = B^N = exp(-i T H) 与 ED 逐位相同。
#   于是 ‖PT − ED‖ 中**不含任何分裂误差**，只剩装配 + IF 内容。
#   · 装配正确 ⟹ ‖PT − ED‖ 应 ~ V⁴δt⁴（单胞）/ 随 N 下降；
#   · 装配错误 ⟹ ‖PT − ED‖ 饱和在 O(1) 或 O(V²)（正是 TEMPO NonDiagonalHyb
#     那条路径的症状：0.29 对路线/步长完全不变）。
#   H_S = 0 的代价：σ⁺ 只由耦合驱动，仍能激发 |↓⟩→|↑⟩，物理不退化。
#
# 两部分：
#   [A] H_S = 0：‖PT − ED‖（纯装配/IF），实时 + 虚时（虚时同时校准 τ 环腿序）
#   [B] H_S ≠ 0（原模型）：三角分解 ‖PT−ED‖ / ‖chain−ED‖ / ‖PT−chain‖，
#       验证两条路线都收敛到 ED，且 PT−chain 是最小量（阶数分离）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. \
#           benchmark/bosons/offdiag_vs_ed.jl
# =============================================================================
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf
include(joinpath(@__DIR__, "offdiag_order.jl"))

const HS0 = zeros(ComplexF64, 2, 2)

# ------------------------------ ED 精确解（56 维） ------------------------------ #
_Hsys(hs) = kron(hs, Matrix(I, NB, NB)) + kron(Matrix(I, 2, 2), Matrix(OMEGA * Diagonal(0.0:(NB - 1))))

function ed_real(T; hs=HS, op=SP, op2=SM, Gc=V)
    H = Hermitian(_Hsys(hs) + coup_mat(op, op2, Gc))
    ev, Vv = eigen(H)
    U = Vv * Diagonal(exp.(-im * T * ev)) * Vv'
    ptrace2(U * kron(RHO0, Diagonal(thermal_bath())) * U')
end

function ed_imag_sz(; hs=HS, op=SP, op2=SM, Gc=V)
    H = Hermitian(_Hsys(hs) + coup_mat(op, op2, Gc))
    ev, Vv = eigen(H)
    w = exp.(-BETA * ev)
    O = kron(SZ, Matrix(I, NB, NB))
    sz = sum(real(dot(Vv[:, i], O * Vv[:, i])) * w[i] for i in eachindex(ev))
    return real(sum(w)), real(sz / sum(w))
end

# ------------------------------ 链（可指定 H_S） ------------------------------ #
function chain_real(N, T; hs=HS, op=SP, op2=SM, Gc=V)
    h = T / N
    hb = OMEGA * Diagonal(0.0:(NB - 1))
    A = kron(exp(-im * h * hs), Matrix(I, NB, NB))
    B = exp(-im * h * (kron(Matrix(I, 2, 2), Matrix(hb)) + coup_mat(op, op2, Gc)))
    U = (B * A)^N
    ptrace2(U * kron(RHO0, Diagonal(thermal_bath())) * U')
end

function chain_imag(N; hs=HS, op=SP, op2=SM, Gc=V)
    h = BETA / N
    hb = OMEGA * Diagonal(0.0:(NB - 1))
    A = kron(exp(-h * hs), Matrix(I, NB, NB))
    B = exp(-h * (kron(Matrix(I, 2, 2), Matrix(hb)) + coup_mat(op, op2, Gc)))
    U = (B * A)^N
    O = kron(SZ, Matrix(I, NB, NB))
    return real(tr(U)), real(tr(O * U) / tr(U))
end

# ------------------------------ PT 求和（可指定 H_S） ------------------------------ #
"pt_sum_rt 的可变 H_S 版本；其余与 offdiag_order.jl 逐字一致。"
function pt_sum_rt_hs(N, T, hyb, hs; Gc=V)
    corr = bosonic_Δt(DiracDelta(OMEGA, α=Gc^2); β=BETA, N=N, t=T, μ=0.0)
    hmat = BosonicHybridizationMatrix(corr; split=true)
    keep = trues(2 * (N + 1)); keep[1] = keep[2] = false
    idx = findall(keep)
    basis = [BosonicContourIndex(k + 1, b) for k in 1:N for b in (:+, :-)]
    h2 = BosonicHybridizationMatrix(basis,
        QuAPI.Dense.BosonicSplitDelta(hmat.Δ.M1[idx, idx], hmat.Δ.M2[idx, idx]))
    P = BosonicIF(h2, hyb)
    M = N
    Tarr = reshape(P.data, ntuple(_ -> 2, 4M))
    Af = exp(-im * (T / N) * hs)
    A2 = [Af[i, j] for i in 1:2, j in 1:2]
    ρ = zeros(ComplexF64, 2, 2)
    inner = 2^(2M + 1)
    @inbounds for cfgP in 0:(inner-1)
        amp = one(ComplexF64)
        for k in 1:M
            amp *= A2[((cfgP >> (2k - 1)) & 1) + 1, ((cfgP >> (2k - 2)) & 1) + 1]
        end
        iszero(amp) && continue
        aP = ((cfgP >> (2M)) & 1) + 1
        for cfgM in 0:(inner-1)
            r00 = RHO0[((cfgP >> 0) & 1) + 1, ((cfgM >> 0) & 1) + 1]
            iszero(r00) && continue
            amp2 = one(ComplexF64)
            for k in 1:M
                amp2 *= A2[((cfgM >> (2k - 1)) & 1) + 1, ((cfgM >> (2k - 2)) & 1) + 1]
            end
            lin = 1; sh = 0
            for k in 1:M
                lin += ((cfgP >> (2k - 1)) & 1) << sh; sh += 1
                lin += (k < M ? ((cfgP >> (2k)) & 1) : ((cfgP >> (2M)) & 1)) << sh; sh += 1
                lin += ((cfgM >> (2k - 1)) & 1) << sh; sh += 1
                lin += (k < M ? ((cfgM >> (2k)) & 1) : ((cfgM >> (2M)) & 1)) << sh; sh += 1
            end
            bM = ((cfgM >> (2M)) & 1) + 1
            ρ[aP, bM] += r00 * amp * conj(amp2) * Tarr[lin]
        end
    end
    return ρ
end

"pt_ring_imag 的可变 H_S 版本；其余与 offdiag_order.jl 逐字一致。"
function pt_ring_imag_hs(N, hyb, hs; Gc=V)
    corr = bosonic_Δτ(DiracDelta(OMEGA, α=Gc^2); β=BETA, N=N, μ=0.0)
    P = BosonicIF(BosonicHybridizationMatrix(corr; split=true), hyb)
    Tarr = reshape(P.data, ntuple(_ -> 2, 2N))
    A = exp(-(BETA / N) * hs)
    AO = SZ * A
    Zq = 0.0 + 0.0im
    Nq = 0.0 + 0.0im
    for x in 0:(2^N-1), y in 0:(2^N-1)
        amp = one(ComplexF64); ampo = one(ComplexF64)
        for p in 1:N
            xp = ((x >> (p-1)) & 1) + 1
            yp = ((y >> (p-1)) & 1) + 1
            amp *= A[yp, xp]
            ampo *= (p == 1 ? AO[yp, xp] : A[yp, xp])
        end
        (iszero(amp) && iszero(ampo)) && continue
        inds = Int[]
        for p in 1:N
            push!(inds, ((x >> mod(p, N)) & 1) + 1, ((y >> (p-1)) & 1) + 1)
        end
        t = Tarr[inds...]
        Zq += amp * t
        Nq += ampo * t
    end
    return real(Zq), real(Nq / Zq)
end

fit_order(ns, es) = begin
    x = log.(ns[2:end] ./ ns[1:end-1])
    (length(x) * sum(x .* log.(es[2:end] ./ es[1:end-1])) -
     sum(x) * sum(log.(es[2:end] ./ es[1:end-1]))) / (length(x) * sum(x .^ 2) - sum(x)^2)
end

# ------------------------------- 主流程 ------------------------------- #
function main()
    io = stdout
    println(io, "=" ^ 76)
    println(io, "BosonicIF 路径求和装配 vs ED（offdiag_order.jl 缺失的独立校验）")
    @printf(io, "单模浴 ω=%.2f, V=%.2f, β=%.2f, NB=%d（全空间 %d 维，ED 精确）\n",
            OMEGA, V, BETA, NB, 2NB)
    println(io, "=" ^ 76)

    # ================= [A] H_S = 0：纯装配/IF 检验（无分裂误差） =================
    println(io, "\n[A1] 实时 σ⁺, H_S=0：链 ≡ ED（B^N = e^{-iTH}），‖PT−ED‖ 只剩装配+IF")
    println(io, "     N    ‖chain−ED‖(应~1e-15)   ‖PT−ED‖        ‖PT−chain‖")
    T = 0.8
    ea = Float64[]; Nb = Int[]
    for N in (2, 3, 4, 5)
        ed = ed_real(T; hs=HS0)
        ch = chain_real(N, T; hs=HS0)
        pt = pt_sum_rt_hs(N, T, SP, HS0)
        push!(ea, norm(pt - ed)); push!(Nb, N)
        @printf(io, "     %-4d %.3e            %.4e        %.4e\n",
                N, norm(ch - ed), norm(pt - ed), norm(pt - ch))
    end
    println(io, "     判读：chain−ED ~ 1e-15 确认无分裂误差；PT−ED 若随 N 下降（而非饱和）⟹ 装配正确。")
    if length(ea) >= 3 && all(ea .> 0)
        @printf(io, "     ‖PT−ED‖ 末端两点比值 = %.2f（δt² ⟹ 4，δt⁴ ⟹ 16，饱和 ⟹ 1）\n",
                ea[end-1] / ea[end])
    end

    println(io, "\n[A2] 虚时 σ⁺, H_S=0：链 ≡ ED，同时独立校准 τ 环腿序")
    println(io, "     N    |chain−ED|(Z)  |Δσz|(chain)   |Δσz|(PT−ED)")
    eb = Float64[]
    for N in (2, 3, 4, 6, 8)
        Ze, sze = ed_imag_sz(hs=HS0)
        Zc, szc = chain_imag(N; hs=HS0)
        Zq, szq = pt_ring_imag_hs(N, SP, HS0)
        push!(eb, abs(szq - sze))
        @printf(io, "     %-4d %.3e       %.3e       %.4e\n", N, abs(Zc - Ze) / abs(Ze), abs(szc - sze), abs(szq - sze))
    end
    if length(eb) >= 3 && all(eb .> 0)
        @printf(io, "     |Δσz|(PT−ED) 末端两点比值 = %.2f（δτ² ⟹ 4）\n", eb[end-1] / eb[end])
    end

    # ================= [B] H_S ≠ 0：三角分解 =================    println(io, "\n[B1] 实时 σ⁺, H_S≠0, T=$(T)：三条边")
    println(io, "     N    ‖chain−ED‖     ‖PT−ED‖       ‖PT−chain‖")
    ec = Float64[]; ep = Float64[]; epc = Float64[]
    for N in (2, 3, 4, 5)
        ed = ed_real(T)
        ch = chain_real(N, T)
        pt = pt_sum_rt(N, T, SP)
        push!(ec, norm(ch - ed)); push!(ep, norm(pt - ed)); push!(epc, norm(pt - ch))
        @printf(io, "     %-4d %.4e       %.4e      %.4e\n", N, ec[end], ep[end], epc[end])
    end
    println(io, "     判读：若 ‖PT−ED‖ ≈ ‖chain−ED‖ ≫ ‖PT−chain‖，则两条路线同步逼近 ED，")
    println(io, "           共享的 O(δt) Lie 误差主导，PT−chain 才是 IF 截断（V⁴δt²）。")
    println(io, "           若 ‖PT−ED‖ 明显大于 ‖chain−ED‖ 且不下降 ⟹ 装配有偏。")

    println(io, "\n[B2] 虚时 σ⁺, H_S≠0：三条边（σz）")
    println(io, "     N    |σz_chain−ED|  |σz_PT−ED|    |σz_PT−chain|")
    for N in (2, 3, 4, 6, 8)
        Ze, sze = ed_imag_sz()
        Zc, szc = chain_imag(N)
        Zq, szq = pt_ring_imag(N, SP)
        @printf(io, "     %-4d %.4e        %.4e      %.4e\n", N,
                abs(szc - sze), abs(szq - sze), abs(szq - szc))
    end

    # ================= [C] σx 对照（装配的已知正确情形） =================
    println(io, "\n[C] σx, H_S=0：‖PT−ED‖（NG≡0，应只剩数值噪声）")
    for N in (2, 4)
        ed = ed_real(T; hs=HS0, op=X, op2=X)
        pt = pt_sum_rt_hs(N, T, X, HS0)
        @printf(io, "     N=%d  ‖PT−ED‖ = %.3e\n", N, norm(pt - ed))
    end

    println(io, "\n结论区（按输出判读）：")
    println(io, "  · [A] 若 ‖PT−ED‖ 在 H_S=0 时随 N 下降 ⟹ 装配正确，")
    println(io, "    offdiag_order 的 V⁴δt² 是真实 IF 截断，不是装配假象。")
    println(io, "  · [A] 若 ‖PT−ED‖ 饱和 ⟹ 装配有问题（腿序/ρ₀/因子），需先修装配。")
    println(io, "  · [B] 给出 PT 与 chain 相对 ED 的绝对位置，判断‘共同分裂误差’假设是否成立。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
