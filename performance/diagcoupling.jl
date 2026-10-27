# =============================================================================
# QuAPI 玻色 IF 离散化验证 —— 对角耦合（ADT）：1/2/4 阶 Trotter 链 vs IF 路径求和
#
# 目的：证明 QuAPI 的玻色影响泛函离散化本身不引入额外误差——
#   "系统侧按 1/2/4 阶 Trotter 展开的全空间演化 RDM" 与
#   "基于相应 QuAPI 离散化 IF 的路径求和 RDM" 逐位一致（机器精度 ~1e-15），
#   即 IF 精确对应参考链在每个 cell 内把耦合算符冻结为该 cell 路径值的离散化。
#
# 模型（单模热浴，谱为 DiracDelta，可全空间精确求解作参考链）：
#   H = H_S + Σ ω b†b + S ⊗ G(b+b')，  S = diag(λ₁,λ₂)（对角耦合）
#   谱密度 A(ω') = G² δ(ω'-ω)，β 有限。
#
# 参考链（全空间矩阵乘，无任何近似）：
#   U = A(lsys[1]) · B(l₁) · A(lsys[2]) · B(l₂) ⋯ B(l_M) · A(ltail)
#   A(τ) = e^{-iτH_S}⊗1，B(τ) = expm(-iτ(1⊗h_b + S⊗G(b+b')))
#   （cell 内耦合冻结：对角 S 时 B 自动谱分解为位移振子，cell kernel 为精确二重积分）
#   · 1 阶  (B(h)A(h))^N          : lsys = [h×N; 0]
#   · 2 阶  (u·B(h)·u)^N, u=A(h/2): lsys = [h/2, h×(N-1), h/2]
#   · 4 阶  Suzuki S₄ = S₂(θh)²S₂(γh)S₂(θh)²，θ=1/(4-4^{1/3})，γ=1-4θ<0：
#     每 cell 浴段 [θh,θh,γh,θh,θh]，系统因子（含跨界合并：中间 cell 尾因子为 θh）
#
# IF 侧（显式双重路径枚举，非 MPS）：
#   ρ[a,b] = Σ ρ₀[s₀,t₀]·amp⁺·conj(amp⁻)·exp(v⁺'D₊₊v⁺+v⁻'D₋₋v⁻+v⁺'D₊₋v⁻+v⁻'D₋₊v⁺)·tail
#   · 等距 cell（1/2 阶）：四块直接取自 QuAPI bosonic_Δt 的 [2:N+1,2:N+1] Toeplitz 窗口
#     （(N+1) 节点矩阵的任意连续 N×N 窗口等价；t=0 节点不参与）
#   · 非等距 cell（4 阶 Suzuki，含负长度 γ 段）：手工闭式核（S∓ 三角核 + K 中心核），
#     脚本内先对等距情形断言其与 bosonic_Δt 四块逐位一致
#
# 运行：julia --startup-file=no --project=.. performance/diagcoupling.jl
# =============================================================================
using QuAPI, ImpurityModelBase, LinearAlgebra, Printf

# ------------------------------ 模型参数 ------------------------------
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]
const HS = .37X + .19Y + .23Z                 # 杂质哈密顿量
const SV = [-0.7, 0.9]                        # 对角耦合算符 S 的本征值 λ
const SMAT = ComplexF64[SV[1] 0; 0 SV[2]]
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5Z) / 2
const OMEGA, G, BETA = 1.3, 0.35, 1.7         # 浴模频率 / 耦合 / 温度
const NB = 28                                  # 参考链玻色 Fock 截断

bathmats(nb) = begin
    b = diagm(1 => sqrt.(1.0:(nb-1)))
    hb = OMEGA * Diagonal(0.0:(nb-1))
    bb = G * (b + b')
    p = exp.(-BETA * OMEGA * (0.0:(nb-1))); p ./= sum(p)
    Matrix(hb), Matrix(bb), Diagonal(p)
end
ptrace(Rm, nb) = [sum(reshape(Rm, nb, 2, nb, 2)[i,a,i,b] for i in 1:nb) for a in 1:2, b in 1:2]

# ------------------------------ 参考链（全空间） ------------------------------
function ref_chain(lsys, lbath; nb=NB, Gc=G)
    bth = diagm(1 => sqrt.(1.0:(nb-1)))
    hb = OMEGA * Diagonal(0.0:(nb-1))
    bb = Gc * (bth + bth')
    rb = exp.(-BETA * OMEGA * (0.0:(nb-1))); rb ./= sum(rb)
    A(τ) = kron(exp(-im*τ*HS), Matrix(I, nb, nb))
    B(τ) = exp(-im*τ*(kron(I(2), Matrix(hb)) + kron(SMAT, bb)))
    M = length(lbath)
    U = A(lsys[1])
    for j in 1:M
        U = B(lbath[j]) * U
        U = A(j < M ? lsys[j+1] : 0.0) * U
    end
    U = A(lsys[end]) * U
    return ptrace(U * kron(RHO0, Diagonal(rb)) * U', nb)
end

# ------------------------------ 手工核（任意 cell 长度，可为负） ------------------------------
boseb(ω) = 1 / (exp(BETA*ω) - 1)
fcell(ω, l) = 2sin(ω*l/2) / ω
Kk(ω, cp, cq, lp, lq) = exp(-im*ω*(cp-cq)) * fcell(ω,lp) * conj(fcell(ω,lq))
Sm(ω, l) = (1 - im*l*ω - exp(-im*l*ω)) / ω^2
Sp(ω, l) = (1 + im*l*ω - exp(im*l*ω)) / ω^2
function manual_delta(lb; Gc=G)
    M = length(lb)
    b = boseb(OMEGA); bb = b + 1
    edges = vcat(0.0, cumsum(lb)); cen = [(edges[j]+edges[j+1])/2 for j in 1:M]
    out = Dict{String,Matrix{ComplexF64}}()
    for key in ("++", "+-", "-+", "--")
        Mx = zeros(ComplexF64, M, M)
        for p in 1:M, q in 1:M
            if p == q
                s₋, s₊ = Sm(OMEGA, lb[p]), Sp(OMEGA, lb[p])
                key == "++" && (Mx[p,p] = -Gc^2*(bb*s₋ + b*s₊))
                key == "+-" && (Mx[p,p] = +Gc^2*(b*s₋ + b*s₊))
                key == "-+" && (Mx[p,p] = +Gc^2*(bb*s₋ + bb*s₊))
                key == "--" && (Mx[p,p] = -Gc^2*(b*s₋ + bb*s₊))
            else
                k = Kk(OMEGA, cen[p], cen[q], lb[p], lb[q])
                key == "++" && (Mx[p,q] = (p > q ? -Gc^2*bb : -Gc^2*b) * k)
                key == "+-" && (Mx[p,q] = +Gc^2*b * k)
                key == "-+" && (Mx[p,q] = +Gc^2*bb * k)
                key == "--" && (Mx[p,q] = (p > q ? -Gc^2*b : -Gc^2*bb) * k)
            end
        end
        out[key] = Mx
    end
    return out
end

# 自检 1：手工核（等距）≡ bosonic_Δt 四块（任意 N×N Toeplitz 窗口）
let N = 4, T = 1.1
    corr = bosonic_Δt(DiracDelta(OMEGA, α=G^2), β=BETA, N=N, t=T, μ=0.0)
    Dm = manual_delta(fill(T/N, N))
    for (key, (b1, b2)) in (("++", (:+, :+)), ("+-", (:+, :-)), ("-+", (:-, :+)), ("--", (:-, :-)))
        Bq = Matrix(branch(corr, b1=b1, b2=b2))[2:N+1, 2:N+1]
        d = norm(Bq - Dm[key])
        @printf("手工核自检 %s: |manual - bosonic_Δt 窗口| = %.3e\n", key, d)
        d < 1e-13 || error("手工核与 QuAPI 核不一致: $key")
    end
end

# ------------------------------ IF 路径求和（显式枚举，非 MPS） ------------------------------
function ifsum(lsys, lbath, D)
    M = length(lbath)
    Afull = [exp(-im*lsys[j]*HS) for j in 1:M]
    Utail = exp(-im*lsys[end]*HS)
    ρ = zeros(ComplexF64, 2, 2)
    for cfgP in 0:(2^(M+1)-1)
        s0 = (cfgP & 1) + 1
        σP = [((cfgP >> j) & 1) + 1 for j in 1:M]
        v⁺ = [SV[σP[j]] for j in 1:M]
        amp = one(ComplexF64); prev = s0
        for j in 1:M; amp *= Afull[j][σP[j], prev]; prev = σP[j]; end
        for cfgM in 0:(2^(M+1)-1)
            t0 = (cfgM & 1) + 1
            σM = [((cfgM >> j) & 1) + 1 for j in 1:M]
            v⁻ = [SV[σM[j]] for j in 1:M]
            amp2 = one(ComplexF64); prev = t0
            for j in 1:M; amp2 *= Afull[j][σM[j], prev]; prev = σM[j]; end
            e = exp(v⁺'D["++"]*v⁺ + v⁻'D["--"]*v⁻ + v⁺'D["+-"]*v⁻ + v⁻'D["-+"]*v⁺)
            for a in 1:2, b in 1:2
                ρ[a,b] += RHO0[s0,t0] * amp * conj(amp2) * e * Utail[a,σP[M]] * conj(Utail[b,σM[M]])
            end
        end
    end
    return ρ
end

# QuAPI 等距核的 N×N 窗口四块
quapi_blocks(N, T) = begin
    corr = bosonic_Δt(DiracDelta(OMEGA, α=G^2), β=BETA, N=N, t=T, μ=0.0)
    Dict("++" => Matrix(branch(corr, b1=:+, b2=:+))[2:N+1, 2:N+1],
         "+-" => Matrix(branch(corr, b1=:+, b2=:-))[2:N+1, 2:N+1],
         "-+" => Matrix(branch(corr, b1=:-, b2=:+))[2:N+1, 2:N+1],
         "--" => Matrix(branch(corr, b1=:-, b2=:-))[2:N+1, 2:N+1])
end

# ------------------------------ 各阶链的 (lsys, lbath) ------------------------------
# Suzuki 4 阶：N 个外层 cell，每个 5 段浴；相邻 cell 的 A(θh/2)·A(θh/2) 合并为 A(θh)
# lsys[k+1] = 浴段 k 与 k+1 之间的系统因子；lsys[1] 首 lead、lsys[M+1] 末尾因子
const θ = 1 / (4 - 4^(1/3))
const γ = 1 - 4θ
function suzuki_segments(N, h)
    seg = [θ*h, θ*h, γ*h, θ*h, θ*h]
    lbath = repeat(seg, N)
    interior = [θ*h, (θ+γ)*h/2, (θ+γ)*h/2, θ*h]   # cell 内 s₁→s₂ … s₄→s₅ 的 4 个因子
    lsys = Vector{Float64}(undef, 5N + 1)
    lsys[1] = θ*h/2                                 # 首 lead
    for c in 1:N
        base = 5(c-1)
        c > 1 && (lsys[base+1] = θ*h)               # 跨界合并（前尾 θh/2 + 本首 θh/2）
        lsys[base+2 : base+5] = interior
    end
    lsys[5N+1] = θ*h/2                              # 仅末 cell 尾因子保持 θh/2
    return lsys, lbath
end

chain_specs(order, N, T) = begin
    h = T / N
    if order == 1
        (vcat(fill(h, N), 0.0), fill(h, N), quapi_blocks(N, T))
    elseif order == 2
        (vcat(h/2, fill(h, N-1), h/2), fill(h, N), quapi_blocks(N, T))
    else
        lsys, lbath = suzuki_segments(N, h)
        (lsys, lbath, manual_delta(lbath))
    end
end

# ------------------------------ 主流程 ------------------------------
function main()
    println("=" ^ 78)
    println("对角耦合（ADT）：Trotter 链 vs QuAPI IF 路径求和 —— 机器精度一致性验证")
    println("H_S 有一般非对角元; S = diag($(SV)); 浴: 单模 ω=$OMEGA, G=$G, β=$BETA")
    println("=" ^ 78)
    cases = [(1, 6, 0.9), (1, 5, 1.3), (2, 6, 0.9), (2, 5, 1.3), (4, 1, 0.8), (4, 2, 0.9)]
    # 噪声底对照：G=0（完全解耦，e 权重≡1，纯酉路径和）下同一最大 M 链的残差
    lsys0, lbath0, _ = chain_specs(4, 2, 0.9)
    noise = norm(ref_chain(lsys0, lbath0; Gc=0.0) - ifsum(lsys0, lbath0, manual_delta(lbath0; Gc=0.0)))
    @printf("双精度路径求和噪声底（G=0 对照, M=%d）: %.2e\n\n", length(lbath0), noise)
    tol(M) = max(1e-12, 10*noise)
    ok = true
    for (order, N, T) in cases
        lsys, lbath, D = chain_specs(order, N, T)
        M = length(lbath)
        r = ref_chain(lsys, lbath)
        ρ = ifsum(lsys, lbath, D)
        d = norm(r - ρ)
        ok &= d < tol(M)
        @printf("%d 阶  N=%d  M=%2d  T=%.2f:  ‖链 − IF‖_F = %.3e   (‖ρ‖=%.3f)\n",
                order, N, M, T, d, norm(r))
    end
    println()
    println(ok ? "结论：全部通过 —— QuAPI 的玻色 IF 离散化与对应阶 Trotter 链逐位一致（小 M 时 ≈1e-15；" *
                 "M=10 时残差为双精度路径求和的消消噪声，已由 G=0 对照标定），无额外离散误差。" :
                 "结论：存在超差项，见上表。")
    return ok
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
