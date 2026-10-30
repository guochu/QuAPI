# =============================================================================
# benchmark/bosons/offdiag_toulouse.jl
#
# Toulouse 模型（杂质=玻色模，耦合 H_int = V(d†⊗b + d⊗b†)，非对角共轭对）下
# QuAPI PT（BosonicIF + split 核）与冻结核单粒子解法的关系判定。
#
# 与 offdiag_order.jl 的 σ⁺（自旋）模型区别：Toulouse 的 [d†,d] = 1 是 c-number，
# 双重嵌套对易子 [d†,[d†,d]] = 0 —— 非高斯内容可能与自旋情形（[σ⁺,σ⁻]=σ_z）
# 结构不同。dF=2 截断会把 [d†,d] 截成 σ_z（即退化为 σ⁺ 模型），故用 dF=3。
#
# 步骤：
#  [1] PT vs 严格全 Fock 链（同一 dF 截断、同一 Lie 分裂）⟹ 阶数判定
#  [2] G⁺(t) 三方对比：PT 提取 / 冻结单粒子解 (A0⁻¹+K^R)A=1 / 严格链
#  [3] （视 [1][2] 结果）提取 PT 隐含的单粒子配对 Σ_eff = A0⁻¹ − A_PT⁻¹，
#      与 Δt 核对照，找出冻结核高效解法与 PT 完全对上所需的顶点移位。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. benchmark/bosons/offdiag_toulouse.jl
# =============================================================================
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf

const OMEGA, V, BETA = 1.3, 0.35, 1.7
const NB = 28
const DF = 3                      # 杂质 Fock 截断
const EPSD = 0.6                  # ϵ_d

# 杂质算符（dF 截断）
adag(dF=DF) = begin
    m = zeros(ComplexF64, dF, dF)
    for n in 0:dF-2; m[n+2, n+1] = sqrt(n + 1); end
    m
end
adown(dF=DF) = adjoint(adag(dF))
n_op(dF=DF) = Diagonal(0.0:(dF-1))
HS = EPSD .* Matrix(n_op())       # 对角 H_S = ϵ_d n̂

coup_full(Gc) = Gc * (kron(adag(), diagm(1 => sqrt.(1.0:(NB-1)))) +
                      kron(adown(), diagm(-1 => sqrt.(1.0:(NB-1)))))

# ------------------------------ 严格全 Fock 链（实时） ------------------------------ #
# U = (B·A)^N（A 先、B 后，与 PT 收缩的节点次序一致），初态 ρ_S = |0⟩⟨0| ⊗ 浴热态
function ref_chain(N, T; Gc=V)
    h = T / N
    hb = OMEGA * Diagonal(0.0:(NB - 1))
    A = kron(exp(-im * h * HS), Matrix(I, NB, NB))
    B = exp(-im * h * (kron(Matrix(I, DF, DF), Matrix(hb)) + coup_full(Gc)))
    U = (B * A)^N
    rb = exp.(-BETA * OMEGA * (0.0:(NB - 1))); rb ./= sum(rb)
    ρ0 = zeros(ComplexF64, DF * NB, DF * NB)
    ρ0[1:NB, 1:NB] = Diagonal(rb)                  # 杂质 |0⟩ ⊗ 浴
    R = U * ρ0 * U'
    return [sum(reshape(R, NB, DF, NB, DF)[i,a,i,b] for i in 1:NB) for a in 1:DF, b in 1:DF]
end

# ------------------------------ QuAPI PT（实时，路径求和，一般 dF） ------------------------------ #
function pt_tensor_rt(N, T, hyb; Gc=V)
    corr = bosonic_Δt(DiracDelta(OMEGA, α=Gc^2); β=BETA, N=N, t=T, μ=0.0)
    hmat = BosonicHybridizationMatrix(corr; split=true)
    keep = trues(2 * (N + 1)); keep[1] = keep[2] = false
    idx = findall(keep)
    basis = [BosonicContourIndex(k + 1, b) for k in 1:N for b in (:+, :-)]
    h2 = BosonicHybridizationMatrix(basis,
        QuAPI.Dense.BosonicSplitDelta(hmat.Δ.M1[idx, idx], hmat.Δ.M2[idx, idx]))
    return BosonicIF(h2, hyb)
end

# 路径求和：v = [n₁,m₁,…,n_N,m_N,o]（dF 进制）。ρ₀ 因子支持非对角初态（ket≠bra），
# 用于 G⁺ 提取（ket 起 |1⟩=d†|0⟩，bra 起 |0⟩）。
function pt_sum_rt(N, T, hyb; Gc=V, rho0=Matrix{ComplexF64}(Diagonal([1.0, 0.0, 0.0])))
    P = pt_tensor_rt(N, T, hyb; Gc=Gc)
    M = N
    Tarr = reshape(P.data, ntuple(_ -> DF, 4M))   # 2N 节点 × (ket,bra) = 4N 指标
    Af = exp(-im * (T / N) * HS)
    ρ = zeros(ComplexF64, DF, DF)
    inner = DF^(2M + 1)
    digs(v, k) = (v ÷ DF^(k - 1)) % DF + 1     # 第 k 位（1 起）
    for cfgP in 0:(inner-1)
        amp = one(ComplexF64)
        for k in 1:M
            amp *= Af[digs(cfgP, 2k), digs(cfgP, 2k - 1)]
        end
        iszero(amp) && continue
        aP = digs(cfgP, 2M + 1)
        for cfgM in 0:(inner-1)
            r00 = rho0[digs(cfgP, 1), digs(cfgM, 1)]
            iszero(r00) && continue
            amp2 = one(ComplexF64)
            for k in 1:M
                amp2 *= Af[digs(cfgM, 2k), digs(cfgM, 2k - 1)]
            end
            # 线性指标：inds = [P₁k,P₁b,M₁k,M₁b,P₂k,…]，首个指标最低位（Julia 列主序）
            lin = 0; sh = 0
            for k in 1:M
                lin += (digs(cfgP, 2k) - 1) * DF^sh; sh += 1
                lin += (digs(cfgP, k < M ? 2k + 1 : 2M + 1) - 1) * DF^sh; sh += 1
                lin += (digs(cfgM, 2k) - 1) * DF^sh; sh += 1
                lin += (digs(cfgM, k < M ? 2k + 1 : 2M + 1) - 1) * DF^sh; sh += 1
            end
            ρ[aP, digs(cfgM, 2M + 1)] += r00 * amp * conj(amp2) * Tarr[lin + 1]
        end
    end
    return ρ
end

# ------------------------------ 冻结单粒子解法（realtime_frozen.jl 的核） ------------------------------ #
"K^R = −(Δ⁺⁺+Δ⁺₋)（因果下三角，含等时），闭式逐模式。"
function KR_closed(N, T; Gc=V)
    h = T / N
    b = 1 / (exp(BETA * OMEGA) - 1)
    KR = zeros(ComplexF64, N + 1, N + 1)
    for j in 1:N+1, k in 1:j
        lag = j - k
        s = lag == 0 ? (1 - im * OMEGA * h - exp(-im * OMEGA * h)) / OMEGA^2 :
                       exp(-im * OMEGA * lag * h) * 2 * (1 - cos(OMEGA * h)) / OMEGA^2
        KR[j, k] = Gc^2 * s
    end
    return KR
end

frozen_A(N, T; Gc=V) = begin
    h = T / N
    A0inv = Matrix{ComplexF64}(I, N + 1, N + 1)
    u = exp(-im * EPSD * h)
    for j in 2:N+1; A0inv[j, j - 1] = -u; end
    [ (inv(A0inv + KR_closed(N, T; Gc=Gc)))[j, 1] for j in 1:N+1 ]
end

# ------------------------------ 严格单粒子链 ------------------------------ #
strict_A(N, T; Gc=V) = begin
    h = T / N
    H2 = [EPSD Gc; Gc OMEGA]          # 单粒子 2×2（单模浴）
    U = (exp(-im * h * [0 Gc; Gc OMEGA]) * exp(-im * h * [EPSD 0; 0 0]))^N
    [ (U^1)[1,1] for _ in 1:1 ]  # 占位，实际在 main 里直接算
end

function main()
    io = stdout
    println(io, "=" ^ 74)
    println(io, "Toulouse（d† 耦合，dF=$DF）：QuAPI PT vs 严格链 vs 冻结单粒子解")
    @printf(io, "浴：单模 ω=%.2f, V=%.2f, β=%.2f, ϵ_d=%.2f\n", OMEGA, V, BETA, EPSD)
    println(io, "=" ^ 74)

    # ---------- [1] PT vs 严格全 Fock 链：n(t) 误差的阶 ----------
    T = 0.8
    println(io, "\n[1] PT(BosonicIF) vs 严格链（同一 dF 截断、同一 Lie 分裂），n(t) 误差")
    Ns = [2, 3]
    ds = Float64[]
    for N in Ns
        r = real(tr(n_op() * ref_chain(N, T)))
        q = real(tr(n_op() * pt_sum_rt(N, T, adag())))
        push!(ds, abs(q - r))
        @printf(io, "    N=%d  δt=%.4f  |n_PT − n_strict| = %.4e\n", N, T / N, ds[end])
    end
    length(ds) > 1 && @printf(io, "    order = %.2f\n",
        log(ds[1] / ds[2]) / log(Ns[2] / Ns[1]))
    # V 标度
    dv = Float64[]
    for Gc in (V, V / 2)
        r = real(tr(n_op() * ref_chain(2, T; Gc=Gc)))
        q = real(tr(n_op() * pt_sum_rt(2, T, adag(); Gc=Gc)))
        push!(dv, abs(q - r))
    end
    @printf(io, "    V 减半残差比 = %.2f（V⁴⟹16；V²⟹4；机器精度⟹噪声）\n", dv[1] / dv[2])

    # ---------- [2] G⁺ 三方：PT 提取 vs 冻结 vs 严格单粒子 ----------
    println(io, "\n[2] G⁺(T)/(-i) 三方对比（N=2, T=0.8）：A_PT / A_frozen / A_strict")
    N, T2 = 2, 0.8
    rhoG = zeros(ComplexF64, DF, DF); rhoG[2, 1] = 1.0     # ⟨d(t)d†⟩：ket 起 |1⟩，bra 起 |0⟩
    Gp = pt_sum_rt(N, T2, adag(); rho0=rhoG)
    A_pt = tr(Gp * adown())      # 实测归一：G⁺ 提取在此约定下直接给出 A（校准于严格链）
    h = T2 / N
    H2 = ComplexF64[EPSD V; V OMEGA]
    Usp = (exp(-im * h * ComplexF64[0 V; V OMEGA]) * exp(-im * h * ComplexF64[EPSD 0; 0 0]))^N
    A_strict_sp = Usp[1, 1]
    A_fr = frozen_A(N, T2)[end]
    @printf(io, "    A_PT      = %s\n", string(round(A_pt, digits=8)))
    @printf(io, "    A_strict  = %s\n", string(round(A_strict_sp, digits=8)))
    @printf(io, "    A_frozen  = %s\n", string(round(A_fr, digits=8)))
    @printf(io, "    |A_PT − A_strict| = %.2e    |A_frozen − A_strict| = %.2e\n",
            abs(A_pt - A_strict_sp), abs(A_fr - A_strict_sp))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
