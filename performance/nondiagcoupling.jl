# =============================================================================
# QuAPI 玻色 IF 离散化验证 —— 非对角耦合（PT）：1/2/4 阶 Trotter 链 vs Dense IF 收缩
#
# 目的：证明 QuAPI 的玻色影响泛函离散化对非对角耦合算符同样不引入额外误差——
#   耦合算符取 S = σx（非对角、厄米，op1 = op2 = S，即 NonAdditiveHyb 情形；
#   TEMPO 包自带的对 PT 路线做机器精度阶数验证时也采用这类算符），
#   "系统侧按 1/2/4 阶 Trotter 展开的全空间演化 RDM" 与
#   "基于相应 QuAPI 离散化 IF（QuAPI.Dense 的 BosonicIF，显式路径收缩，非 MPS）
#    算出的 RDM" 逐位一致（机器精度 ~1e-15）。
#
# 模型：H = H_S + 1⊗h_b + G·σx ⊗ (b+b')，单模热浴（DiracDelta 谱），β 有限。
#
# IF 侧（QuAPI.Dense）：
#   BosonicIF(BosonicHybridizationMatrix(corr; split=true), σx) = exp(H)，
#   H 按分支组装（+ 支行腿矩阵 opᵀ、− 支 op'，跨支求和核配同算符，详见 models.jl
#   docstring；对 op1 = op2 退化为 Σ_pq [M1·op_p·op_q + M2·op_q·op_p]），
#   节点按 contour 序 (+,1),(-,1),(+,2),(-,2),… 交错排列，每节点 (ket 腿, bra 腿)。
#   · 等距 cell（1/2 阶）：M1/M2 取 bosonic_Δt 的 [2:N+1] Toeplitz 窗口（丢 t=0 节点）
#   · 非等距 cell（4 阶 Suzuki，含负长度 γ 段）：手工 M1/M2（等时 S∓ 三角核按
#     ηⱼₖ/ηₖⱼ 拆分 + K 中心核），脚本内先断言等距情形与 QuAPI 的 M1/M2 逐位一致
#
# 收缩（显式枚举前后向路径，每支 [n₁,m₁,n₂,m₂,…] 共 2M+1 位）：
#   ρ₀[n₁,ñ₁] → A₁[m₁,n₁] → node₁[m₁|n₂] → A₂[m₂,n₂] → … → node_M[m_M|o] → A_tail[a,o]
#   节点 (+,k) 的腿 = (m_k, n_{k+1})（门输入 m、输出 n），(-,k) 同理（conj 链）
#
# 链结构（与对角耦合脚本相同）：
#   · 1 阶 (B(h)A(h))^N              : lsys = h×N, ltail = 0
#   · 2 阶 (u·B(h)·u)^N, u = A(h/2)  : lsys = [h/2, h×(N-1)], ltail = h/2
#   · 4 阶 Suzuki（跨界 A(θh/2)²→A(θh) 合并）
#
# 附录（σ⁺ 非厄米耦合，op1 = σ⁺ ≠ op2 = σ⁻，[op1,op2] ≠ 0）：
#   以严格 IF（精确 cell 耦合演化算符 + 偏迹环境，双实现互证）为参照，验证
#   BosonicIF 的分支修正：O(G²) 系数一致（提取噪声级），全链残差 O(G⁴)
#   （非高斯累积量，G 减半残差 ×16）——对 [op1,op2]≠0 严格 IF 无 exp(双线性)
#   闭式，Gaussian 因子化误差不因修复消失，但被压到 O(G⁴)（详见文末打印）。
#
# 运行：julia --startup-file=no --project=.. performance/nondiagcoupling.jl
# =============================================================================
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf

# ------------------------------ 模型参数 ------------------------------
const X = ComplexF64[0 1; 1 0]
const Y = ComplexF64[0 -im; im 0]
const Z = ComplexF64[1 0; 0 -1]
const HS = .37X + .19Y + .23Z                 # 杂质哈密顿量（含非对角元）
const SP = ComplexF64[0 1; 0 0]                # σ⁺（附录用）
const SM = ComplexF64[0 0; 1 0]                # σ⁻
const RHO0 = (Matrix{ComplexF64}(I, 2, 2) + .3X + .4Y + .5Z) / 2
const OMEGA, G, BETA = 1.3, 0.35, 1.7
const NB = 28

ptrace(Rm, nb) = [sum(reshape(Rm, nb, 2, nb, 2)[i,a,i,b] for i in 1:nb) for a in 1:2, b in 1:2]

# ------------------------------ 参考链（全空间 expm） ------------------------------
# op = σx：非对角厄米耦合；op = σ⁺：JC 型 H_int = G(σ⁺⊗b + σ⁻⊗b†)（附录）
function ref_chain(op, lsys, lbath; nb=NB, Gc=G, op2=op)
    bth = diagm(1 => sqrt.(1.0:(nb-1)))
    hb = OMEGA * Diagonal(0.0:(nb-1))
    coup = Gc * (kron(op, Matrix(bth)) + kron(op2, Matrix(bth')))
    rb = exp.(-BETA * OMEGA * (0.0:(nb-1))); rb ./= sum(rb)
    A(τ) = kron(exp(-im*τ*HS), Matrix(I, nb, nb))
    B(τ) = exp(-im*τ*(kron(I(2), Matrix(hb)) + coup))
    M = length(lbath)
    U = A(lsys[1])
    for j in 1:M
        U = B(lbath[j]) * U
        U = A(j < M ? lsys[j+1] : 0.0) * U
    end
    U = A(lsys[end]) * U
    return ptrace(U * kron(RHO0, Diagonal(rb)) * U', nb)
end

# ------------------------------ 链结构（各阶） ------------------------------
const θ = 1 / (4 - 4^(1/3))
const γ = 1 - 4θ
function suzuki_segments(N, h)
    lbath = repeat([θ*h, θ*h, γ*h, θ*h, θ*h], N)
    interior = [θ*h, (θ+γ)*h/2, (θ+γ)*h/2, θ*h]
    lsys = Vector{Float64}(undef, 5N + 1)
    lsys[1] = θ*h/2
    for c in 1:N
        base = 5(c-1)
        c > 1 && (lsys[base+1] = θ*h)               # 跨界合并
        lsys[base+2 : base+5] = interior
    end
    lsys[5N+1] = θ*h/2
    return lsys, lbath
end
chain_specs(order, N, T) = begin
    h = T / N
    if order == 1
        (vcat(fill(h, N), 0.0), fill(h, N))
    elseif order == 2
        (vcat(h/2, fill(h, N-1), h/2), fill(h, N))
    else
        suzuki_segments(N, h)
    end
end

# ------------------------------ 手工 M1/M2（任意 cell 长度，可为负） ------------------------------
boseb(ω) = 1 / (exp(BETA*ω) - 1)
fcell(ω, l) = 2sin(ω*l/2) / ω
Kk(ω, cp, cq, lp, lq) = exp(-im*ω*(cp-cq)) * fcell(ω,lp) * conj(fcell(ω,lq))
Sm(ω, l) = (1 - im*l*ω - exp(-im*l*ω)) / ω^2
Sp(ω, l) = (1 + im*l*ω - exp(im*l*ω)) / ω^2
# 节点序：(+,1),(-,1),(+,2),(-,2),…（与 contour_ordering 一致）；p = 2(k-1)+(b==:+ ? 1 : 2)
function manual_M12(lb; Gc=G)
    M = length(lb)
    b = boseb(OMEGA); bbar = b + 1
    edges = vcat(0.0, cumsum(lb)); cen = [(edges[j]+edges[j+1])/2 for j in 1:M]
    n = 2M
    M1 = zeros(ComplexF64, n, n); M2 = zeros(ComplexF64, n, n)
    for p in 1:n, q in 1:n
        kp, bp = (p+1)÷2, isodd(p) ? :+ : :-
        kq, bq = (q+1)÷2, isodd(q) ? :+ : :-
        if kp == kq
            s₋, s₊ = Sm(OMEGA, lb[kp]), Sp(OMEGA, lb[kp])
            if (bp, bq) == (:+, :+)
                M1[p,q] = -Gc^2*bbar*s₋;  M2[p,q] = -Gc^2*b*s₊
            elseif (bp, bq) == (:+, :-)
                M1[p,q] = +Gc^2*b*s₋;     M2[p,q] = +Gc^2*b*s₊
            elseif (bp, bq) == (:-, :+)
                M1[p,q] = +Gc^2*bbar*s₋;  M2[p,q] = +Gc^2*bbar*s₊
            else
                M1[p,q] = -Gc^2*b*s₋;     M2[p,q] = -Gc^2*bbar*s₊
            end
        else
            k = Kk(OMEGA, cen[kp], cen[kq], lb[kp], lb[kq])
            later = kp > kq
            if (bp, bq) == (:+, :+)
                later ? (M1[p,q] = -Gc^2*bbar*k) : (M2[p,q] = -Gc^2*b*k)
            elseif (bp, bq) == (:+, :-)
                later ? (M1[p,q] = +Gc^2*b*k) : (M2[p,q] = +Gc^2*b*k)
            elseif (bp, bq) == (:-, :+)
                later ? (M1[p,q] = +Gc^2*bbar*k) : (M2[p,q] = +Gc^2*bbar*k)
            else
                later ? (M1[p,q] = -Gc^2*b*k) : (M2[p,q] = -Gc^2*bbar*k)
            end
        end
    end
    return M1, M2
end

# 自检：手工 M1/M2（等距）≡ QuAPI 窗口 M1/M2
let N = 3, T = 0.9
    corr = bosonic_Δt(DiracDelta(OMEGA, α=G^2), β=BETA, N=N, t=T, μ=0.0)
    hmat = BosonicHybridizationMatrix(corr; split=true)
    keep = trues(2*(N+1)); keep[1] = keep[2] = false
    idx = findall(keep)
    M1q, M2q = hmat.Δ.M1[idx, idx], hmat.Δ.M2[idx, idx]
    M1m, M2m = manual_M12(fill(T/N, N))
    d = max(norm(M1q - M1m), norm(M2q - M2m))
    @printf("手工 M1/M2 自检（等距 vs QuAPI 窗口）: max‖Δ‖ = %.3e\n", d)
    d < 1e-13 || error("手工 M1/M2 与 QuAPI 不一致")
end

# ------------------------------ PT 张量与收缩 ------------------------------
pt_basis(M) = [BosonicContourIndex(k+1, b) for k in 1:M for b in (:+, :-)]

function pt_tensor(M1, M2, hyb)
    M = size(M1, 1) ÷ 2
    h2 = BosonicHybridizationMatrix(pt_basis(M),
          QuAPI.Dense.BosonicSplitDelta(M1, M2))
    return BosonicIF(h2, hyb)
end
# 等距 cell：直接取 bosonic_Δt 的窗口 M1/M2
function pt_tensor_uniform(N, T, hyb; Gc=G)
    corr = bosonic_Δt(DiracDelta(OMEGA, α=Gc^2), β=BETA, N=N, t=T, μ=0.0)
    hmat = BosonicHybridizationMatrix(corr; split=true)
    keep = trues(2*(N+1)); keep[1] = keep[2] = false
    idx = findall(keep)
    return pt_tensor(hmat.Δ.M1[idx, idx], hmat.Δ.M2[idx, idx], hyb)
end

# 收缩：枚举 [n₁,m₁,…,n_M,m_M,o]（每支 2M+1 位）
#   链: ρ₀[n₁,ñ₁] → A₁[m₁,n₁] → node₁[m₁|n₂] → … → node_M[m_M|o] → A_tail[a,o]
#   节点(+,k) 腿 = (m_k, n_{k+1})；节点(-,k) 腿 = (m̃_k, ñ_{k+1})
pt_sum(P, lsys, lbath; ltail) =
    pt_sum_arr(reshape(P.data, ntuple(_ -> 2, 2*length(P.basis))), lsys, lbath; ltail)
function pt_sum_arr(Tarr, lsys, lbath; ltail)
    M = length(lbath)
    Afull = [exp(-im*lsys[k]*HS) for k in 1:M]
    Utail = exp(-im*ltail*HS)
    ρ = zeros(ComplexF64, 2, 2)
    inner = 2^(2M+1)
    for cfgP in 0:(inner-1)
        vP = [((cfgP >> (k-1)) & 1) + 1 for k in 1:2M+1]   # [n₁,m₁,…,n_M,m_M,o]
        amp = one(ComplexF64)
        for k in 1:M; amp *= Afull[k][vP[2k], vP[2k-1]]; end
        for cfgM in 0:(inner-1)
            vM = [((cfgM >> (k-1)) & 1) + 1 for k in 1:2M+1]
            amp2 = one(ComplexF64)
            for k in 1:M; amp2 *= Afull[k][vM[2k], vM[2k-1]]; end
            for a in 1:2, b in 1:2
                f = Utail[a, vP[2M+1]] * conj(Utail[b, vM[2M+1]])
                iszero(f) && continue
                inds = Int[]
                for k in 1:M
                    push!(inds, vP[2k], k < M ? vP[2k+1] : vP[2M+1])   # node(+,k)
                    push!(inds, vM[2k], k < M ? vM[2k+1] : vM[2M+1])   # node(-,k)
                end
                ρ[a,b] += RHO0[vP[1],vM[1]] * amp * conj(amp2) * f * Tarr[inds...]
            end
        end
    end
    return ρ
end

# ================= 附录工具：σ⁺（op1 ≠ op2）的严格 IF 与修正 Gaussian 张量 =================
# 严格 cell 块：V = expm(-ih(h_b + G(σ⁺⊗b + σ⁻⊗b†)))，按 (系统出 a, bath α, 系统入 m, bath β) 取块
function jc_vblocks(h; Gc=G, nb=24)
    bth = diagm(1 => sqrt.(1.0:(nb-1)))
    hb = OMEGA * Diagonal(0.0:(nb-1))
    V = exp(-im*h*(kron(Matrix(I,2,2), Matrix(hb)) +
                   Gc*(kron(SP, Matrix(bth)) + kron(SM, Matrix(bth')))))
    permutedims(reshape(V, nb, 2, nb, 2), (2, 1, 4, 3))     # Vb[a,α,m,β]
end
# 严格 N-cell IF 张量的 O(G²) 系数（自然基 (a₁..,b₁..,m₁..,mm₁..)）：
#   F = Tr[B⁺_{a_N m_N}···B⁺_{a₁ m₁} ρ_b (B⁻_{b₁ mm₁})†···(B⁻_{b_N mm_N})†]
function jc_exact_c2(N, h; Gc=0.02, nb=24)
    Vb = jc_vblocks(h; Gc, nb)
    p = exp.(-BETA*OMEGA*(0.0:(nb-1))); p ./= sum(p)
    ρb = Diagonal(p)
    B(a, m) = Vb[a, :, m, :]
    out = zeros(ComplexF64, ntuple(_ -> 2, 4N))
    for idx in CartesianIndices(out)
        I = Tuple(idx)
        a = I[1:N]; b_ = I[N+1:2N]; m = I[2N+1:3N]; mm = I[3N+1:4N]
        id0 = all(a[q] == m[q] && b_[q] == mm[q] for q in 1:N) ? 1.0 : 0.0   # G=0 恒等张量
        Mx = B(a[N], m[N])
        for q in N-1:-1:1; Mx = Mx * B(a[q], m[q]); end
        Mx = Mx * ρb
        for q in 1:N; Mx = Mx * B(b_[q], mm[q])'; end
        out[idx] = (tr(Mx) - id0) / Gc^2
    end
    return out
end
# 严格单 cell IF 超算符（全空间）与 JC 激发数分块解析（双实现交叉验证）
function jc_IF_full(h; Gc=G, nb=24)
    Vb = jc_vblocks(h; Gc, nb)
    p = exp.(-BETA*OMEGA*(0.0:(nb-1))); p ./= sum(p)
    ρb = Diagonal(p)
    F = zeros(ComplexF64, 2, 2, 2, 2)
    for i in 1:2, j in 1:2, k in 1:2, l in 1:2
        F[i,j,k,l] = tr(Vb[i,:,k,:] * ρb * Vb[j,:,l,:]')
    end
    F
end
function jc_IF_sector(h; Gc=G, nb=24)
    # 扇区 n≥1: 块基 (|↑,n-1⟩,|↓,n⟩)，H_n = [(n-1)ω G√n; G√n nω]；扇区 0: |↓,0⟩
    R = Vector{Matrix{ComplexF64}}(undef, nb)
    R[1] = reshape([one(ComplexF64)], 1, 1)
    for n in 1:nb-1
        R[n+1] = exp(-im*h*ComplexF64[(n-1)*OMEGA Gc*sqrt(n); Gc*sqrt(n) n*OMEGA])
    end
    pn = exp.(-BETA*OMEGA*(0:nb-1)); pn ./= sum(pn)
    uu(m,n) = (m==n)   ? R[n+2][1,1] : 0.0im
    ud(m,n) = (m==n-1) ? R[n+1][1,2] : 0.0im
    du(m,n) = (m==n+1) ? R[n+2][2,1] : 0.0im
    dd(m,n) = n==0 ? (m==0 ? R[1][1,1] : 0.0im) : (m==n ? R[n+1][2,2] : 0.0im)
    el(s,m,si,n) = s==1 ? (si==1 ? uu(m,n) : ud(m,n)) : (si==1 ? du(m,n) : dd(m,n))
    F = zeros(ComplexF64, 2, 2, 2, 2)
    N = nb - 2
    for i in 1:2, j in 1:2, k in 1:2, l in 1:2
        s = 0.0 + 0.0im
        for n in 0:N, m in 0:N
            a = el(i, m, k, n); c = el(j, m, l, n)
            (iszero(a) || iszero(c)) && continue
            s += pn[n+1] * a * conj(c)
        end
        F[i,j,k,l] = s
    end
    F
end
# 修正 Gaussian 生成元（严格到 O(G²)）：
# 元素矩阵 E[r,c] = O[c,r]（门的 (入 r, 出 c) 腿元素），节点序 (+,1),(-,1),(+,2),…
#   同 cell:   +b|f|²·Eσ⁺₊Eσ⁺₋ + (1+b)|f|²·Eσ⁻₊Eσ⁻₋
#              -(1+b)S₋·(σ⁺σ⁻)₊D₋ - bS₊·(σ⁻σ⁺)₊D₋ -(1+b)S₊·D₊(σ⁺σ⁻)₋ - bS₋·D₊(σ⁻σ⁺)₋
#   跨支:      +b·e^{-iω(c_p-c_q)}|f|²·Eσ⁺ₚEσ⁺_q + (1+b)·e^{-iω(c_q-c_p)}|f|²·Eσ⁻ₚEσ⁻_q
#   同支(p晚): -(1+b)e^{-iωΔ}|f|²·Eσ⁺ₚEσ⁻_q - b·conj·Eσ⁻ₚEσ⁺_q（bra-bra 共轭）
const Eσp = SP'; const Eσm = SM'; const EpmE = (SP*SM)'; const EmpE = (SM*SP)'
const Dd = Matrix{Float64}(I, 2, 2)
function _addpair!(H, ν, μ, w, E1, E2)
    for I in CartesianIndices(H)
        ok = true
        for η in 1:ndims(H) ÷ 2
            (η ≠ ν && η ≠ μ && I[2η-1] ≠ I[2η]) && (ok = false; break)
        end
        ok || continue
        H[I] += w * E1[I[2ν-1], I[2ν]] * E2[I[2μ-1], I[2μ]]
end
end
function corrected_pt_gen(N, T, Gc)
    h = T / N
    bth, bb = boseb(OMEGA), 1 + boseb(OMEGA)
    f2 = abs2((1 - exp(-im*OMEGA*h)) / (im*OMEGA))
    S₋, S₊ = Sm(OMEGA, h), Sp(OMEGA, h)
    nn = 2N
    H = zeros(ComplexF64, ntuple(_ -> 2, 2nn))
    nodeP(p) = 2p - 1; nodeM(p) = 2p
    for p in 1:N, q in 1:N
        if p == q
            _addpair!(H, nodeP(p), nodeM(p), Gc^2*bth*f2,  Eσp, Eσp)
            _addpair!(H, nodeP(p), nodeM(p), Gc^2*bb*f2,   Eσm, Eσm)
            _addpair!(H, nodeP(p), nodeM(p), -Gc^2*bb*S₋,  EpmE, Dd)
            _addpair!(H, nodeP(p), nodeM(p), -Gc^2*bth*S₊, EmpE, Dd)
            _addpair!(H, nodeM(p), nodeP(p), -Gc^2*bb*S₊,  EpmE, Dd)
            _addpair!(H, nodeM(p), nodeP(p), -Gc^2*bth*S₋, EmpE, Dd)
        else
            Kpq = exp(-im*OMEGA*(p-q)*h) * abs2(fcell(OMEGA, h))   # e^{-iω(c_p-c_q)}|f|²（σ⁺σ⁺ 用）
            _addpair!(H, nodeP(p), nodeM(q), Gc^2*bth*Kpq,  Eσp, Eσp)
            _addpair!(H, nodeP(p), nodeM(q), Gc^2*bb*conj(Kpq), Eσm, Eσm)
            p > q || continue
            K = exp(-im*OMEGA*(p-q)*h) * abs2(fcell(OMEGA, h))
            _addpair!(H, nodeP(p), nodeP(q), -Gc^2*bb*K,   Eσp, Eσm)
            _addpair!(H, nodeP(p), nodeP(q), -Gc^2*bth*K', Eσm, Eσp)
            _addpair!(H, nodeM(p), nodeM(q), -Gc^2*bb*K',  Eσp, Eσm)
            _addpair!(H, nodeM(p), nodeM(q), -Gc^2*bth*K,  Eσm, Eσp)
        end
    end
    return H
end
function corrected_pt_tensor(N, T, Gc)
    H = corrected_pt_gen(N, T, Gc)
    nn = 2N
    perm = vcat([2ν-1 for ν in 1:nn], [2ν for ν in 1:nn])
    Hm = reshape(permutedims(H, perm), 2^nn, 2^nn)
    Tm = reshape(exp(Hm), ntuple(_ -> 2, 2nn))
    return reshape(permutedims(Tm, sortperm(perm)), ntuple(_ -> 2, 2nn))
end
# 自然基 (a,b,m,mm) → 张量布局 (m,a,mm,b)（每 cell 交错）的置换
natural2tensor_perm(N) = begin
    perm = Int[]
    for p in 1:N; append!(perm, [2N+p, p, 3N+p, N+p]); end
    perm
end

# ------------------------------ 主流程 ------------------------------
function main()
    println("=" ^ 78)
    println("非对角耦合（PT, S = σx, op1 = op2）：Trotter 链 vs QuAPI.Dense IF —— 一致性验证")
    println("H_S 有一般非对角元; 浴: 单模 ω=$OMEGA, G=$G, β=$BETA")
    println("=" ^ 78)
    cases = [(1, 3, 0.9), (1, 4, 1.2), (2, 3, 0.9), (2, 4, 1.2), (4, 1, 0.8)]
    ok = true
    for (order, N, T) in cases
        lsys, lbath = chain_specs(order, N, T)
        M = length(lbath)
        P = M == N ? pt_tensor_uniform(N, T, X) :
            pt_tensor(manual_M12(lbath)..., X)
        r = ref_chain(X, lsys, lbath)
        ρ = pt_sum(P, lsys, lbath; ltail=lsys[end])
        d = norm(r - ρ)
        tol = M ≤ 5 ? 1e-12 : 1e-10
        ok &= d < tol
        @printf("%d 阶  N=%d  M=%2d  T=%.2f:  ‖链 − PT‖_F = %.3e   (‖ρ‖=%.3f)  [%s]\n",
                order, N, M, T, d, norm(r), M == N ? "QuAPI 窗口张量" : "手工非等距张量")
    end
    println()
    println(ok ? "σx 主检验全部通过：QuAPI.Dense 的 PT IF 与对应阶 Trotter 链逐位一致（≈1e-15），无额外离散误差。" :
                 "σx 主检验存在超差项，见上表。")

    # ---------------- 附录：σ⁺（op1 ≠ op2）—— 严格 IF 对照与 PT 修复验证 ----------------
    println()
    println("-" ^ 78)
    println("附录：非厄米耦合 S = σ⁺（op1 = σ⁺ ≠ op2 = σ⁻）—— 严格 IF 与 PT 修复验证")
    println("-" ^ 78)
    h = 0.3
    # (1) 严格单 cell IF：expm 耦合演化算符 + 偏迹环境，并与 JC 激发数分块解析交叉验证
    dsec = maximum(abs.(jc_IF_full(h) .- jc_IF_sector(h)))
    @printf("(1) 严格单 cell IF 双实现（全空间 expm+偏迹 vs JC 分块解析）: max|Δ| = %.2e\n", dsec)
    dsec < 1e-12 || error("严格 IF 双实现不一致")
    # (2) BosonicIF(σ⁺) 的 O(G²) 系数 vs 严格提取：修复后应为提取噪声级
    #     （修复前旧表达式 H=Σ M1·op1_p·op2_q + M2·op2_p·op1_q 在此为 9.98e-02 结构错配）
    Gsm = 0.02
    c2ex = permutedims(jc_exact_c2(1, h; Gc=Gsm), natural2tensor_perm(1))
    Tq = reshape(pt_tensor_uniform(1, h, SP; Gc=Gsm).data, 2,2,2,2)
    cq = [(Tq[m,a,mm,b] - ((m==a && mm==b) ? 1.0 : 0.0im)) / Gsm^2
          for m in 1:2, a in 1:2, mm in 1:2, b in 1:2]
    dq = maximum(abs.(c2ex .- cq))
    @printf("(2) BosonicIF(σ⁺) 的 O(G²) 系数 vs 严格: max|Δ| = %.2e（提取噪声级；修复前 9.98e-02）\n", dq)
    dq < 1e-4 || error("BosonicIF(σ⁺) 的 O(G²) 系数与严格 IF 不符")
    # (3) 手工修正生成元（独立实现）== 严格 O(G²) 提取
    for N in (1, 2)
        c2exN = permutedims(jc_exact_c2(N, h; Gc=Gsm), natural2tensor_perm(N))
        genN = corrected_pt_gen(N, N*h, 1.0)
        @printf("(3) N=%d: 手工修正生成元 vs 严格 O(G²) 提取 max|Δ| = %.2e\n",
                N, maximum(abs.(c2exN .- genN)))
    end
    # (4) 全链对照（σ⁺ 耦合参照 σ⁺⊗b + σ⁻⊗b†）：残差应按 O(G⁴) 标度（G 减半 ×16）
    for (N, T) in ((2, 0.6), (3, 0.9))
        lsys, lbath = chain_specs(1, N, T)
        dG = Float64[]
        for Gc in (G, G/2)
            r = ref_chain(SP, lsys, lbath; Gc, op2=SM)
            ρq = pt_sum(pt_tensor_uniform(N, T, SP; Gc), lsys, lbath; ltail=0.0)
            push!(dG, norm(r - ρq))
            @printf("(4) N=%d G=%.3f:  ‖ref − PT‖ = %.2e\n", N, Gc, dG[end])
        end
        @printf("    G 减半残差比 = %.1f（O(G⁴) ⟹ ≈16）\n", dG[1] / dG[2])
        10 < dG[1] / dG[2] || error("全链残差未按 O(G⁴) 标度")
    end
    println("""
    结论：
    · 严格 IF 总可以通过「精确算出 cell 耦合演化算符 + 偏迹环境」得到（全空间 expm 与 JC 激发数
      分块解析双实现，互证到 1e-15）。对 [op1,op2]≠0 它不含任何 exp(双线性) 闭式。
    · BosonicIF 的 PT 构建已按分支修正：+ 支行腿矩阵用 opᵀ、− 支用 op'，同支块 − 分支 M1↔M2
      互换，跨支块用求和核 (M1+M2) 配同算符配对——跨支对应为 op1×op1*（b 权重）与 op2×op2*
      （1+b 权重），即 bra 支 (B⁻)† 中的 a† 恰与 σ⁺ 元素相配。修复前旧式把跨支权重放到
      op1×op2 错误通道（O(G²) 结构错配 9.98e-02，产生不存在的相干混合通道、缺失布居转移通道）。
    · 修复后 BosonicIF(σ⁺) 与严格 IF 一致到 O(G²)（(2) 提取噪声级），全链残差按 G⁴ 标度
      （(4) G 减半残差×16）——即 O(G⁴) 的非高斯累积量（Magnus 对易子项），低于一/二阶
      Trotter 误差，不破坏收敛阶。
    · 对角耦合与 op1=op2（σx / a+a′）时修正公式退化为原公式（彼时严格，主检验 ≈1e-15）。""")
    return ok
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
