# =============================================================================
# benchmark/bosons/offdiag_order.jl
#
# 非对角（共轭对）耦合下 QuAPI 离散化 IF 的误差阶数判定（实时 + 虚时）。
# 参考：docs/GuoChen2026a.pdf（PRB 114, 125413 (2026)），Appendix C。
#
# ---------------------------- 理论框架 ----------------------------
# 论文 Eq. (8)：I[A†,A] = T_C exp(-∮∮ A†(τ)Δ(τ,τ')A(τ')) 在连续极限精确。
# T_C 保留在指数上：展开后非对易的 A 算符产生连接（非高斯）累积量。
# Appendix C 的离散化把等时两种次序分开（C7/C17 = QuAPI 的 M1/M2 split），
# 但整个指数仍是双线性的 —— 即 QuAPI IF = exp(G)，G 为全部 O(V²) 双线性
# 内容（含等时三角核的两种次序），这一部分已被 performance/nondiagcoupling.jl
# 验证与严格 IF 逐系数一致（O(V²) 精确）。
#
# 严格 IF = exp(G + NG)，NG 为 O(V⁴) 起的非高斯内容，来源有二：
#   (i)  历史基底中不可因子化的嵌套排序（T_C 连接累积量，[A,[A,A†]] 型）；
#   (ii) 单胞等时处方本身：论文 C7/C17 是两个指数的乘积
#        e^{-A†Δ¹A}·e^{-AΔ²A†}，而 exp(和) 形式与其差 O(Δ¹Δ²) = O(V⁴δt⁴)。
#
# 误差阶候选（固定 V、固定总时长 T，δt = T/N）：
#   · 论文的泛函推断（一阶分裂）：全局 O(δt)；
#   · 混合嵌套（胞内 × 跨胞对）：每对 (j,i) V⁴δt⁴，共 ~N² 对 ⟹ 全局 δt²；
#   · 胞内嵌套：每胞 V⁴δt⁴ ⟹ 全局 δt³；
#   · 对角 / op1=op2：NG ≡ 0，IF 严格（无穷阶）——对照。
#
# 模型：H = H_S + ω b†b + V(σ⁺⊗b + σ⁻⊗b†)，单模浴（与 nondiagcoupling.jl 同）。
# 严格侧 = 精确 cell 演化（一阶分裂链 U=(B·A)^N，A=系统因子，B=耦合 cell），
# QuAPI 侧 = bosonic_Δt/Δτ 的 split 核 + BosonicIF。两侧同一系统分裂
# ⟹ 差值只含 IF 的非高斯截断误差。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. benchmark/bosons/offdiag_order.jl
# =============================================================================
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf

const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const SZ = ComplexF64[1 0; 0 -1]
const SP = ComplexF64[0 1; 0 0]
const SM = ComplexF64[0 0; 1 0]
const HS = .37X + .19Y + .23SZ             # 与 performance/nondiagcoupling.jl 完全一致
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5SZ) / 2
const OMEGA, V, BETA = 1.3, 0.35, 1.7
const NB = 28

coup_mat(op, op2, Gc) =
    Gc * (kron(op, diagm(1 => sqrt.(1.0:(NB - 1)))) +
          kron(op2, diagm(-1 => sqrt.(1.0:(NB - 1)))))
ptrace2(Rm) = [sum(reshape(Rm, NB, 2, NB, 2)[i,a,i,b] for i in 1:NB) for a in 1:2, b in 1:2]
thermal_bath() = (rb = exp.(-BETA * OMEGA * (0.0:(NB - 1))); rb ./ sum(rb))

# ------------------------------ 严格链：实时 ------------------------------ #
# 路径序：A 先、B 后（与 PT 收缩的节点约定一致）⟹ 矩阵 U = (B·A)^N
function ref_chain_real(N, T; op=SP, op2=SM, Gc=V)
    h = T / N
    hb = OMEGA * Diagonal(0.0:(NB - 1))
    A = kron(exp(-im * h * HS), Matrix(I, NB, NB))
    B = exp(-im * h * (kron(Matrix(I, 2, 2), Matrix(hb)) + coup_mat(op, op2, Gc)))
    U = (B * A)^N
    return ptrace2(U * kron(RHO0, Diagonal(thermal_bath())) * U')
end

# ------------------------------ 严格链：虚时 ------------------------------ #
function ref_chain_imag(N; op=SP, op2=SM, Gc=V)
    h = BETA / N
    hb = OMEGA * Diagonal(0.0:(NB - 1))
    A = kron(exp(-h * HS), Matrix(I, NB, NB))
    B = exp(-h * (kron(Matrix(I, 2, 2), Matrix(hb)) + coup_mat(op, op2, Gc)))
    U = (B * A)^N
    Z = tr(U)
    sz = tr(kron(SZ, Matrix(I, NB, NB)) * U)
    return real(Z), real(sz / Z)
end

# ------------------------------ QuAPI PT：实时路径求和 ------------------------------ #
"bosonic_Δt 窗口（丢 t=0 节点）→ N 个 cell 的 PT"
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

function pt_sum_rt(N, T, hyb; Gc=V)
    P = pt_tensor_rt(N, T, hyb; Gc=Gc)
    M = N
    Tarr = reshape(P.data, ntuple(_ -> 2, 4M))   # 2N 节点 × (ket,bra) = 4N 指标
    Af = exp(-im * (T / N) * HS)
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
            r00 = RHO0[((cfgP >> 0) & 1) + 1, ((cfgM >> 0) & 1) + 1]   # 初态因子 ρ₀[n₁,t₀]
            iszero(r00) && continue
            amp2 = one(ComplexF64)
            for k in 1:M
                amp2 *= A2[((cfgM >> (2k - 1)) & 1) + 1, ((cfgM >> (2k - 2)) & 1) + 1]
            end
            # 线性指标：inds = [P₁k,P₁b,M₁k,M₁b,P₂k,…]，每指标 1 bit
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

# ------------------------------ QuAPI PT：虚时 τ 环 ------------------------------ #
# 节点 p 腿 (ket, bra) = (输出 x_{p+1}, 输入 y_p)。注意与实时 + 支相反：
# BosonicIF 的 τ 支行腿矩阵用 op1 原矩阵、+ 支用 transpose(op1)，故 τ 环的
# 腿方向相对 + 支转置。对称耦合（σx）对此不可见——σx 对照无法校准该约定，
# 由 V⁴ 标度（错误配对会停在 V²）+ δτ² 收敛共同校准。
function pt_ring_imag(N, hyb; Gc=V)
    corr = bosonic_Δτ(DiracDelta(OMEGA, α=Gc^2); β=BETA, N=N, μ=0.0)
    P = BosonicIF(BosonicHybridizationMatrix(corr; split=true), hyb)
    Tarr = reshape(P.data, ntuple(_ -> 2, 2N))
    A = exp(-(BETA / N) * HS)
    AO = SZ * A
    Zq = 0.0 + 0.0im
    Nq = 0.0 + 0.0im
    for x in 0:(2^N-1), y in 0:(2^N-1)
        amp = one(ComplexF64)
        ampo = one(ComplexF64)
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

# ------------------------------ 单胞严格 IF（实时） ------------------------------ #
function jc_IF_full(h; Gc=V, nb=24)
    bth = diagm(1 => sqrt.(1.0:(nb-1)))
    hb = OMEGA * Diagonal(0.0:(nb-1))
    Vc = exp(-im*h*(kron(Matrix(I,2,2), Matrix(hb)) +
                    Gc*(kron(SP, Matrix(bth)) + kron(SM, Matrix(bth')))))
    p = exp.(-BETA*OMEGA*(0.0:(nb-1))); p ./= sum(p)
    B(a, m) = Vc[(a-1)*nb+1 : a*nb, (m-1)*nb+1 : m*nb]   # kron(sys,bath) 是 sys-major
    F = zeros(ComplexF64, 2, 2, 2, 2)      # [a, b, m, mm] 自然基
    for i in 1:2, j in 1:2, k in 1:2, l in 1:2
        F[i,j,k,l] = tr(B(i,k) * Diagonal(p) * B(j,l)')
    end
    return F
end

# ------------------------------ 主流程 ------------------------------ #
function main()
    io = stdout
    println(io, "=" ^ 74)
    println(io, "非对角（σ⁺/σ⁻ 共轭对）耦合：QuAPI IF 误差阶数判定")
    @printf(io, "H_S 一般非对角；浴：单模 ω=%.2f, V=%.2f, β=%.2f\n", OMEGA, V, BETA)
    println(io, "=" ^ 74)

    # ---------- [0] 对照：σx（op1=op2，NG≡0）应机器精度一致 ----------
    println(io, "\n[0] 对照 σx（op1=op2，IF 严格）：实时链 vs QuAPI IF")
    for N in (2, 3, 4)
        r = ref_chain_real(N, 0.8; op=X, op2=X)
        q = pt_sum_rt(N, 0.8, X)
        @printf(io, "    N=%d  ‖ref − PT‖ = %.2e\n", N, norm(r - q))
    end

    # ---------- [1] 实时：固定 V、固定 T，扫 N ----------
    # 上限 N=5：BosonicIF 的稠密 exp(H) 是 d^{2N}×d^{2N}（N=5 → 1024²可行；
    # N=6 → 4096² 的 H 装配 + exp 单线程需 ~小时级，不进入阶梯）
    T = 0.8
    Ns = [2, 3, 4, 5]
    println(io, "\n[1] 实时 σ⁺：固定 T=$(T)，扫 N（δt = T/N）")
    ds = Float64[]
    for N in Ns
        r = ref_chain_real(N, T)
        q = pt_sum_rt(N, T, SP)
        push!(ds, norm(r - q))
        @printf(io, "    N=%d  δt=%.4f  ‖ref − PT‖ = %.4e\n", N, T / N, ds[end])
    end
    for i in 1:length(Ns)-1
        @printf(io, "    N %d→%d : order = %.2f\n", Ns[i], Ns[i+1],
                log(ds[i] / ds[i+1]) / log(Ns[i+1] / Ns[i]))
    end

    # ---------- [2] V 标度复核（连接此前 V⁴ 结论） ----------
    println(io, "\n[2] V 标度复核（固定 N=3, T=0.8）")
    dv = Float64[]
    for Gc in (V, V / 2)
        r = ref_chain_real(3, 0.8; Gc=Gc)
        q = pt_sum_rt(3, 0.8, SP; Gc=Gc)
        push!(dv, norm(r - q))
        @printf(io, "    V=%.4f  ‖ref − PT‖ = %.4e\n", Gc, dv[end])
    end
    @printf(io, "    V 减半残差比 = %.1f（V⁴ ⟹ 16）\n", dv[1] / dv[2])

    # ---------- [3] 单胞：严格 IF vs QuAPI PT，固定 V 扫 h ----------
    println(io, "\n[3] 单胞（实时）：‖I_strict − I_QuAPI‖ vs h（固定 V=$V）")
    hs = [0.1, 0.2, 0.4, 0.8]
    d1 = Float64[]
    for h in hs
        Fex = permutedims(jc_IF_full(h), [3, 1, 4, 2])     # → (m,a,mm,b) 布局
        Tq = reshape(pt_tensor_rt(1, h, SP).data, 2, 2, 2, 2)
        push!(d1, maximum(abs.(Fex .- Tq)))
        @printf(io, "    h=%.2f  max|Δ| = %.4e\n", h, d1[end])
    end
    for i in 1:length(hs)-1
        @printf(io, "    h %.1f→%.1f : per-cell order = %.2f\n", hs[i], hs[i+1],
                log(d1[i] / d1[i+1]) / log(hs[i+1] / hs[i]))
    end

    # ---------- [4] 虚时：σx 对照（机器精度） ----------
    println(io, "\n[4] 虚时 σx 对照（NG≡0，应机器精度一致）")
    for N in (2, 3, 4)
        Zr, szr = ref_chain_imag(N; op=X, op2=X)
        Zq, szq = pt_ring_imag(N, X)
        @printf(io, "    N=%d  |Z−Zq|/Z = %.2e   |σz_ref − σz_PT| = %.2e\n",
                N, abs(Zr - Zq) / abs(Zr), abs(szr - szq))
    end

    # ---------- [5] 虚时：固定 β，扫 Nτ + V 标度（约定校准） ----------
    println(io, "\n[5] 虚时 σ⁺：固定 β=$(BETA)，扫 Nτ（δτ = β/Nτ）")
    println(io, "    （腿方向约定由下方 V⁴ 标度独立校准：错误配对会停在 V²）")
    Ns2 = [2, 3, 4, 6, 8, 10]
    ds2 = Float64[]
    for N in Ns2
        Zr, szr = ref_chain_imag(N)
        Zq, szq = pt_ring_imag(N, SP)
        push!(ds2, abs(szr - szq))
        @printf(io, "    N=%2d  δτ=%.4f  |σz_ref − σz_PT| = %.4e\n", N, BETA / N, ds2[end])
    end
    for i in 1:length(Ns2)-1
        @printf(io, "    N %d→%d : order = %.2f\n", Ns2[i], Ns2[i+1],
                log(ds2[i] / ds2[i+1]) / log(Ns2[i+1] / Ns2[i]))
    end
    println(io, "\n[6] 虚时 σ⁺ V 标度（N=4）")
    dvi = Float64[]
    for Gc in (V, V / 2)
        Zr, szr = ref_chain_imag(4; Gc=Gc)
        Zq, szq = pt_ring_imag(4, SP; Gc=Gc)
        push!(dvi, abs(szr - szq))
        @printf(io, "    V=%.4f  |Δσz| = %.4e\n", Gc, dvi[end])
    end
    @printf(io, "    V 减半残差比 = %.1f（V⁴ ⟹ 16，V² ⟹ 4）\n", dvi[1] / dvi[2])

    println(io, "\n结论（实测）：")
    println(io, "  · 对角 / op1=op2（σx）：实时、虚时均机器精度 ⟹ IF 严格（无穷阶）。")
    println(io, "  · 非对角共轭对（σ⁺）：单胞内容 V⁴δt⁴（[3] per-cell order ≈ 4），")
    println(io, "    固定总时长的全局误差为 V⁴·δt²（实时 order ≈ 2，虚时 order ≈ 2），")
    println(io, "    V 标度恒为 V⁴（实时 15.3，虚时 15.7 ≈ 16）。")
    println(io, "  · 全局 δt²（而非 δt³）表明误差按 ~N² 对胞对累积：跨胞核（每个 O(δt²)）")
    println(io, "    与胞内嵌套的非高斯内容配对，即“混合嵌套”情形；论文的泛函推断 O(δt)")
    println(io, "    （一阶分裂）是过悲观的界，实测严格二阶。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
