# =============================================================================
# QuAPI 费米 IF 严格性验证 —— 单模费米 bath：严格准绝热链 IF（Fock→相干态）vs QuAPI.Dense
#
# 问题：QuAPI.Dense 的费米 IF（fermionic_Δt 离散化 → GrassmannIF，I = exp(Σ āΔa)）
#   是否是"严格 IF"——是否与严格准绝热链的 IF 在相干态（Grassmann）表象下逐位一致
#   （玻色对角耦合已验证逐位一致，见 diagcoupling.jl）。
#
# 方法（单模 bath，谱 A(ε) = g²δ(ε−ε₀)）：
#   模型  H = ε_d d†d + ε₀ c†c + g(d†c + c†d)，bath 热态 ρ_b = diag(1−n, n)
#   链    U = (B(h)A(h))^N，A = e^{−ihε_d d†d}，
#         B = expm[−ih(ε₀c†c + g(d†c+c†d))]   ←—— 杂质-环境耦合项的演化算符，严格计算
#   1) 严格 IF（Fock 表象）：杂质 Fock 态路径插入 + 环境自由度严格求迹
#        I_Fock[s⁺_{0..N}; r⁻_{0..N}]
#          = Tr_b[B_{s⁺_N s⁺_{N-1}}⋯B_{s⁺_1 s⁺_0} ρ_b (B_{r⁻_1 r⁻_0})†⋯(B_{r⁻_N r⁻_{N-1}})†]
#      （B_{s's} = ⟨s'|B|s⟩ 为 2×2 bath 矩阵；用于自检 A 的 ρ(T) 收缩）
#   2) 相干态表象转换（外代数 bitmask 实现，逐 cell 相干矩阵元乘积求迹）：
#        ⟨η̄⁺_k|B|η⁺_{k-1}⟩ = Σ_{s',s} (−ā⁺_k)^{s'} B[s',s] (−a⁺_{k-1})^s（2×2 多项式矩阵）
#        I_coh = Tr_b[Π_k⟨η̄⁺_k|B|η⁺_{k-1}⟩ ρ_b Π_k⟨η̄⁻_{k-1}|B†|η⁻_k⟩]
#                / Π_k(1+ā⁺_k a⁺_{k-1})(1+ā⁻_{k-1} a⁻_k)
#      （分母为 g=0 自由重叠 ⟹ I(g=0)≡1，与 QuAPI IF 的归一化一致；注意须用完整相干
#        矩阵元——Fock 路径版只在 Berezin 积分后等价，作为多项式不含部分占据重叠项）
#   3) 与 QuAPI.Dense 对比：GrassmannIF 全系数（链侧腿子集 2^{4N} 个单项式）逐位比较
#
# 腿约定（0 基节点 j=0..N，时间 j·h；两分支每节点 (a, ā)）：
#   链侧自然腿（B 因子相干矩阵元产生）：ā⁺(j) j=1..N（cell j 输出）、a⁺(j) j=0..N−1
#   （cell j+1 输入）；ā⁻(j) j=0..N−1、a⁻(j) j=1..N（后向时间对应）。
#   QuAPI 的 Δ 腿对把节点 p 映到节点右侧的 cell（block[p+1,p+1] = 同 cell 双重积分）。
#   ⟹ 两种对齐都做对比：
#     (A) 直接对齐（包内 quapi_rlm_rt_gt 的使用约定）
#     (B) ā 腿左移一格对齐（链侧 ā(node k) ↔ QuAPI 行 k：同 cell 系数 ↔ QuAPI 对角元）
#
# 自检（须全部机器精度通过，否则报错）：
#   A. Fock 路径收缩 ρ(T) == 直接 4 维严格演化（B、JW 记号、路径求和全对）
#   B. g=0 ⇒ I_coh ≡ 1（重叠除法正确）
#   D. g=0 ⇒ Grassmann 全积分（测度×动能×初态×末态投影）== 自由杂质解析演化
#      （费米相干态密度矩阵核有 ⟨−η̄|⋯|η⟩ 型负号约定：通道符号 (−1)^{a(1−b)} 由 D 校准）
#   E. 小 g 二次核 == 块迹闭式（与多项式机器完全独立的矩阵代数路径）
#   C.（信息项）I_exact × 离散动能骨架 的 ρ(T)：组合本身是 O(δt²) 离散化，非机器精度
#
# 结论预告（脚本实测证实；T(ε)=∫₀^h ds∫₀^s ds' e^{−iε(s−s')} 三角核，
#         R(m)=∫∫_{cell 对, 中心分离 m·h} e^{−iε(t−t')} 矩形核；均取排序单项式符号约定）：
#   · 同分支二次核 O(g²) 主阶逐项一致：跨 cell 双方 = +n̄g²R(m)（吻合至 O(g⁴)）；
#     同 cell 链 = g²[n̄T(ε₀) − n·δt²/2] = tr((B00−B11)ρ_b)，QuAPI = g²[n̄T − nT*(ε₀)]
#     （n 部分仅差 g²n[T*−δt²/2] = O(g²δt³)）。
#   · 常数项（主导差异）：链 IF = tr(B00ᴺρB00†ᴺ) = 1 − O(N·g²δt²)——cell 内"跳去又跳回"
#     虚过程对节点相干插入不可见（B00 块非幺正的概率泄漏）；QuAPI 恒为 1。
#   · 混合分支（G₊₋/G₋₋ 型）与同节点系数带相位/符号差；高次单项式 O(g⁴) 非 Gaussian。
#   ⟹ 非逐位严格（区别于玻色对角耦合 diagcoupling.jl 的逐位一致）：费米杂质算符无法
#     在路径扇内冻结，IF 带内在 O(g²δt²)/cell 离散误差；ρ(T) 层面保持一阶收敛。
#
# 运行：julia --startup-file=no --project=../test performance/fermionicIF.jl
#   （根目录无 Manifest；test 环境含 QuAPI/ImpurityModelBase 的 path 依赖）
# =============================================================================
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf

# ------------------------------ 模型参数 ------------------------------
const EPSD = 0.6                                  # 杂质能级 ε_d
const EPS0 = 0.9                                  # 单 bath 模能量 ε₀（谱 g²δ(ε−ε₀)）
const GC   = 0.35                                 # 杂质-环境耦合 g
const BETA = 1.2                                  # 逆温度（μ=0）
const RHO0 = ComplexF64[0.62 0.2-0.15im; 0.2+0.15im 0.38]   # 一般混合初态

nocc()  = 1 / (exp(BETA * EPS0) + 1)              # bath 模占据 n
nbar()  = 1 - nocc()

# ------------------------------ 双模费米代数（4 维 Fock） ------------------------------
# 基底 |n_d n_c⟩ = (d†)^{n_d}(c†)^{n_c}|0⟩，线性指标 = 2n_d + n_c + 1（杂质为慢指标）。
# d† 无 JW 串；c† 带 (−1)^{n_d}（d† 在 c† 左侧的次序约定）。
function cellB(h; g = GC)
    ddag = zeros(ComplexF64, 4, 4); ddag[3, 1] = 1; ddag[4, 2] = 1
    cdag = zeros(ComplexF64, 4, 4); cdag[2, 1] = 1; cdag[4, 3] = -1
    dd, c = ddag', cdag'
    return exp(-im * h * (EPS0 * (cdag * c) + g * (ddag * c + cdag * dd)))
end
cellA(h) = kron(exp(-im * h * EPSD * ComplexF64[0 0; 0 1]), Matrix(I, 2, 2))
rhob() = Diagonal([1 - nocc(), nocc()])

# 求迹 bath（快指标）：R[(n_d,n_c),(n_d',n_c')]，线性指标 = 2n_d + n_c + 1
ptrace4(R) = [sum(reshape(R, 2, 2, 2, 2)[i, a, i, b] for i in 1:2) for a in 1:2, b in 1:2]

# ------------------------------ 外代数（Grassmann）工具 ------------------------------
# 多项式 = Vector{ComplexF64}，长 2^L；data[m+1] 为单项式 ∏_{k∈m}ξ_k（腿位升序）的系数。
# 腿位（1 基，与 QuAPI basis 顺序一致）：pos(节点 j0, 分支, ā?)
pos(j0::Int, br::Symbol, ā::Bool) = 4 * j0 + (br === :+ ? 0 : 2) + (ā ? 2 : 1)

# 乘积序腿序列 → (mask, 重排符号)；腿序列无重复（链侧结构保证）
function mono_from_seq(seq::Vector{Int})
    s = 0
    for i in 1:length(seq), j in (i + 1):length(seq)
        seq[i] > seq[j] && (s += 1)
    end
    mask = 0
    for p in seq; mask |= 1 << (p - 1); end
    return mask, isodd(s) ? -1 : 1
end

# ξ^A·ξ^B（均升序、不相交）的合并符号：(−1)^{Σ_{b∈B}#{a∈A: a>b}}
function merge_sign(A::Int, B::Int)
    s = 0
    while B > 0
        lb = B & -B; B -= lb
        s += count_ones(A & ~(2 * lb - 1))
    end
    return isodd(s) ? -1 : 1
end

# 稀疏因子：[(腿序列, 系数)...] → [(mask, 系数·重排符号)...]
function sf(terms::Vector{Tuple{Vector{Int},ComplexF64}})
    out = Tuple{Int,ComplexF64}[]
    for (seq, c) in terms
        if isempty(seq)
            push!(out, (0, c))
        else
            m, s = mono_from_seq(seq)
            push!(out, (m, c * s))
        end
    end
    return out
end

# 稠密 × 稀疏（out = P·S）
function gds(P::Vector{ComplexF64}, S::Vector{Tuple{Int,ComplexF64}})
    out = zeros(ComplexF64, length(P))
    for (m, c) in S
        iszero(c) && continue
        for m2 in 0:length(P)-1
            v = P[m2+1]
            iszero(v) && continue
            (m & m2) == 0 || continue
            out[(m | m2) + 1] += c * v * merge_sign(m2, m)
        end
    end
    return out
end

kval(data::Vector{ComplexF64}, ps::Int...) =
    data[foldl((a, p) -> a | (1 << (p - 1)), ps; init = 0) + 1]

# ------------------------------ 严格链：Fock 表象 IF 与 ρ(T) ------------------------------
# I_Fook[sp, rp]：sp/rp 的第 j 位（j=0..N）= 节点 j 的杂质占据 s_j / r_j
function fock_IF(N::Int, T::Real; g = GC)
    h = T / N
    B = cellB(h; g = g)
    Blk = [B[2s1+1:2s1+2, 2s0+1:2s0+2] for s1 in 0:1, s0 in 0:1]   # ⟨s'|B|s⟩
    rb = rhob()
    n1 = N + 1
    IF = zeros(ComplexF64, 2^n1, 2^n1)
    for sp in 0:2^n1-1, rp in 0:2^n1-1
        X = Matrix{ComplexF64}(I, 2, 2)
        for k in N:-1:1
            X = X * Blk[((sp >> k) & 1) + 1, ((sp >> (k - 1)) & 1) + 1]
        end
        Y = Matrix{ComplexF64}(I, 2, 2)
        for k in 1:N
            Y = Y * Blk[((rp >> k) & 1) + 1, ((rp >> (k - 1)) & 1) + 1]'
        end
        IF[sp+1, rp+1] = tr(X * rb * Y)
    end
    return IF
end

function exact_rho(N::Int, T::Real; g = GC)
    h = T / N
    B = cellB(h; g = g); A4 = cellA(h)
    U = Matrix{ComplexF64}(I, 4, 4)
    for _ in 1:N; U = B * A4 * U; end
    return ptrace4(U * kron(RHO0, Matrix(rhob())) * U')
end

# 自检 A 用的 Fock 路径收缩（含 A 因子与初态）
function fock_rho(N::Int, T::Real, IF; g = GC)
    h = T / N; x = exp(-im * h * EPSD)
    n1 = N + 1
    ρ = zeros(ComplexF64, 2, 2)
    for a in 0:1, b in 0:1
        acc = 0.0 + 0.0im
        for sp in 0:2^n1-1, rp in 0:2^n1-1
            ((sp >> N) & 1) == a && ((rp >> N) & 1) == b || continue
            amp = RHO0[(sp & 1) + 1, (rp & 1) + 1]
            amp *= x^count_ones(sp & (2^N - 1))           # Π_k e^{−ihε_d s_{k-1}}
            amp *= conj(x)^count_ones(rp & (2^N - 1))
            acc += amp * IF[sp+1, rp+1]
        end
        ρ[a+1, b+1] = acc
    end
    return ρ
end

# ------------------------------ 相干态表象转换（严格 IF → Grassmann 多项式） ------------------------------
# 逐 cell 相干矩阵元（2×2 bath 矩阵、Grassmann 多项式元）乘积 + 求迹：
#   ⟨η̄⁺_k|B|η⁺_{k-1}⟩[b',b] = Σ_{s',s} (−ā⁺_k)^{s'} B[s',s][b',b] (−a⁺_{k-1})^s
#   ⟨η̄⁻_{k-1}|B†|η⁻_k⟩[b,b'] = Σ_{r,r'} (−ā⁻_{k-1})^r (B[r',r])†[b,b'] (−a⁻_k)^{r'}
# 注意：必须用完整相干矩阵元（含各 cell 独立求和的"占据不一致"项），Fock 路径版
# 只在 Berezin 积分后才与之一致，作为多项式不同（g=0 时后者不含部分占据重叠项）。
# 稀疏多项式 = Dict{mask, coeff}（mask 为腿位 bitmask，单项式取腿位升序）
function pmul(A::Dict{Int,ComplexF64}, B::Dict{Int,ComplexF64})
    out = Dict{Int,ComplexF64}()
    for (ma, ca) in A
        for (mb, cb) in B
            (ma & mb) == 0 || continue
            m = ma | mb
            out[m] = get(out, m, 0.0im) + ca * cb * merge_sign(ma, mb)
        end
    end
    return out
end
const ONEP = Dict{Int,ComplexF64}(0 => 1.0 + 0.0im)

# 单 cell 相干矩阵元：ket 侧（分支 br，out 节点 o，in 节点 i）
function cellket(Blk, o::Int, i::Int, br::Symbol)
    f = Matrix{Dict{Int,ComplexF64}}(undef, 2, 2)
    for b1 in 1:2, b2 in 1:2
        d = Dict{Int,ComplexF64}()
        for s1 in 0:1, s0 in 0:1
            c = Blk[s1+1, s0+1][b1, b2]
            iszero(c) && continue
            seq = Int[]
            s1 == 1 && push!(seq, pos(o, br, true))
            s0 == 1 && push!(seq, pos(i, br, false))
            m, s = isempty(seq) ? (0, 1) : mono_from_seq(seq)
            d[m] = get(d, m, 0.0im) + c * ((s1 + s0) % 2 == 1 ? -1 : 1) * s
        end
        f[b1, b2] = d
    end
    return f
end

# 单 cell 相干矩阵元：bra 侧（B†，out 节点 o=k−1，in 节点 i=k）
function cellbra(Blk, o::Int, i::Int)
    f = Matrix{Dict{Int,ComplexF64}}(undef, 2, 2)
    for b1 in 1:2, b2 in 1:2
        d = Dict{Int,ComplexF64}()
        for r in 0:1, r1 in 0:1
            c = conj(Blk[r1+1, r+1][b2, b1])       # (B[r1,r])†[b,b'] = conj(B[r1,r][b',b])
            iszero(c) && continue
            seq = Int[]
            r == 1 && push!(seq, pos(o, :-, true))
            r1 == 1 && push!(seq, pos(i, :-, false))
            m, s = isempty(seq) ? (0, 1) : mono_from_seq(seq)
            d[m] = get(d, m, 0.0im) + c * ((r + r1) % 2 == 1 ? -1 : 1) * s
        end
        f[b1, b2] = d
    end
    return f
end

pmmul(X, F) = [begin
        d = Dict{Int,ComplexF64}()
        for k in 1:2
            d = merge(+, d, pmul(X[i, k], F[k, j]))
        end
        d
    end for i in 1:2, j in 1:2]

function coherent_IF(N::Int, T::Real; g = GC)
    h = T / N
    B = cellB(h; g = g)
    Blk = [B[2s1+1:2s1+2, 2s0+1:2s0+2] for s1 in 0:1, s0 in 0:1]
    rb = rhob()
    I2 = Matrix{Dict{Int,ComplexF64}}(undef, 2, 2)
    for i in 1:2, j in 1:2
        I2[i, j] = i == j ? ONEP : Dict{Int,ComplexF64}()
    end
    X = I2
    for k in N:-1:1
        X = pmmul(X, cellket(Blk, k, k - 1, :+))
    end
    Y = I2
    for k in 1:N
        Y = pmmul(Y, cellbra(Blk, k - 1, k))
    end
    # I = tr(X·ρ_b·Y) = Σ_{i,k} X[i,k]·ρ_b[k,k]·Y[k,i]
    tot = Dict{Int,ComplexF64}()
    for i in 1:2, k in 1:2
        tot = merge(+, tot, pmul(pmul(X[i, k], Dict{Int,ComplexF64}(0 => rb[k, k])), Y[k, i]))
    end
    L = 4 * (N + 1)
    data = zeros(ComplexF64, 1 << L)
    for (m, c) in tot
        data[m+1] = c
    end
    # 除自由重叠：(1+ā⁺_k a⁺_{k-1})⁻¹(1+ā⁻_{k-1} a⁻_k)⁻¹ ⟸ 乘 (1−āa)(1−āa)
    one_ = 1.0 + 0.0im
    for k in 1:N
        data = gds(data, sf([(Int[], one_), ([pos(k, :+, true), pos(k - 1, :+, false)], -one_)]))
        data = gds(data, sf([(Int[], one_), ([pos(k - 1, :-, true), pos(k, :-, false)], -one_)]))
    end
    return data
end

# ------------------------------ Grassmann 收缩：IF → ρ(T) ------------------------------
# 测度 e^{−āa}（每节点每分支）× 自由动能 × 初态相干核 × 末态投影，取全积分（顶单项式）
function rho_from_IF(ifdata::Vector{ComplexF64}, N::Int, T::Real)
    h = T / N; x = exp(-im * h * EPSD)
    one_ = 1.0 + 0.0im
    P = copy(ifdata)
    for k in 1:N                                    # ⟨η̄⁺_k|A|η⁺_{k-1}⟩ = 1 + x·ā⁺_k a⁺_{k-1}
        P = gds(P, sf([(Int[], one_), ([pos(k, :+, true), pos(k - 1, :+, false)], x)]))
        P = gds(P, sf([(Int[], one_), ([pos(k - 1, :-, true), pos(k, :-, false)], conj(x))]))
    end
    # 初态 ⟨η̄⁺₀|ρ₀|η⁻₀⟩ = Σ_{s,r} (−ā⁺(0))ˢ ρ₀[s,r] (−a⁻(0))ʳ
    P = gds(P, sf([(Int[], RHO0[1, 1]),
                   ([pos(0, :+, true)], -RHO0[2, 1]),
                   ([pos(0, :-, false)], -RHO0[1, 2]),
                   ([pos(0, :+, true), pos(0, :-, false)], RHO0[2, 2])]))
    # Berezin 测度（节点 0..N，两分支）
    for j in 0:N, br in (:+, :-)
        P = gds(P, sf([(Int[], one_), ([pos(j, br, true), pos(j, br, false)], -one_)]))
    end
    # 末态投影 (−a⁺(N))ᵃ(−ā⁻(N))ᵇ 后全积分（顶单项式）。
    # 通道符号 (−1)^{a(1−b)}：费米相干态密度矩阵核的符号约定（⟨−η̄|⋯|η⟩ 型负号约定），
    # 由自检 D（g=0 时与自由杂质精确演化逐位一致，N=1..3）校准标定。
    ρ = zeros(ComplexF64, 2, 2)
    for a in 0:1, b in 0:1
        seq = Int[]
        a == 1 && push!(seq, pos(N, :+, false))
        b == 1 && push!(seq, pos(N, :-, true))
        Q = gds(P, sf([(seq, (-one_)^(a + b))]))
        ρ[a+1, b+1] = (a == 1 && b == 0 ? -1 : 1) * Q[end]
    end
    return ρ
end

# ------------------------------ QuAPI 侧 ------------------------------
function quapi_IF(N::Int, T::Real; g = GC)
    corr = fermionic_Δt(DiracDelta(EPS0; α = g^2), β = BETA, N = N, t = T, μ = 0.0)
    hmat = GrassmannHybridizationMatrix(corr)
    for (k, leg) in enumerate(hmat.basis)          # basis 顺序与 pos 约定一致性自检
        pos(leg.j - 1, leg.branch, leg.conj) == k ||
            error("QuAPI basis 顺序与脚本 pos 约定不一致")
    end
    return GrassmannIF(hmat), corr
end

# ------------------------------ 主对比 ------------------------------
function compare_case(N::Int, T::Real; g = GC)
    edata = coherent_IF(N, T; g = g)
    qif, corr = quapi_IF(N, T; g = g)
    qdata = qif.data
    L = 4 * (N + 1)
    full = (1 << L) - 1
    unused = (1 << (pos(0, :+, true) - 1)) | (1 << (pos(N, :+, false) - 1)) |
             (1 << (pos(0, :-, false) - 1)) | (1 << (pos(N, :-, true) - 1))
    # (B) 对齐的位重映射（几何对齐：链侧腿 ↔ 其来源 cell 的 QuAPI 行/列）：
    #   + 分支：链侧 ā⁺(k) 来自 cell k（节点左侧）↦ QuAPI ā⁺(k−1)（k=1..N）；a⁺ 不动
    #   − 分支（时间后向，偏移反向）：链侧 a⁻(k) 来自 cell k ↦ QuAPI a⁻(k−1)；ā⁻ 不动
    shiftmap = Dict{Int,Int}()
    for j in 0:N, br in (:+, :-)
        shiftmap[pos(j, br, false) - 1] = br === :- && j >= 1 ?
                                          pos(j - 1, br, false) - 1 : pos(j, br, false) - 1
        shiftmap[pos(j, br, true) - 1] = br === :+ && j >= 1 ?
                                         pos(j - 1, br, true) - 1 : pos(j, br, true) - 1
    end
    maxA = 0.0; maxB = 0.0; degmax = Dict{Int,Float64}()
    nmon = 0; nbadB = 0
    for m in 0:full
        (m & unused) == 0 || continue
        nmon += 1
        dA = abs(edata[m+1] - qdata[m+1])
        maxA = max(maxA, dA)
        qm = 0
        for b in 0:L-1
            ((m >> b) & 1) == 1 && (qm |= 1 << shiftmap[b])
        end
        dB = abs(edata[m+1] - qdata[qm+1])
        maxB = max(maxB, dB)
        deg = count_ones(m)
        degmax[deg] = max(get(degmax, deg, 0.0), dB)
        dB > 1e-12 && (nbadB += 1)
    end
    return (edata = edata, qdata = qdata, corr = corr, maxA = maxA, maxB = maxB,
            nmon = nmon, nbadB = nbadB, degmax = degmax)
end

# ------------------------------ 解析核（O(g²) 预告） ------------------------------
Ttri(ε, h) = abs(ε) < 1e-10 ? h^2 / 2 : (1 - im * ε * h - exp(-im * ε * h)) / ε^2
Rrect(ε, m, h) = exp(-im * ε * m * h) * 2 * (1 - cos(ε * h)) / ε^2

creim(z) = @sprintf("(%+.6e%+.6ei)", real(z), imag(z))

# ------------------------------ 主流程 ------------------------------
function main()
    println("=" ^ 86)
    println("费米 QuAPI IF 严格性验证：单模 bath 准绝热链严格 IF（Fock→相干态） vs QuAPI.Dense")
    @printf("模型: εd=%.2f ε₀=%.2f g=%.2f β=%.2f  → n=%.4f n̄=%.4f\n",
            EPSD, EPS0, GC, BETA, nocc(), nbar())
    println("=" ^ 86)

    cases = [(1, 0.5), (2, 0.9), (3, 1.2)]
    res = Dict{Tuple{Int,Float64},typeof(compare_case(1, 0.5))}()
    ok = true

    # ---------- 自检 ----------
    println("\n[自检] A: Fock 路径收缩 vs 直接 4 维严格演化（机器精度）")
    println("       B: g=0 ⇒ I_coh ≡ 1（重叠除法正确）")
    println("       D: g=0 ⇒ Grassmann 骨架全积分 vs 自由杂质解析演化（机器精度）")
    println("       E: 小 g 二次核 vs 块迹闭式 tr((B00−B11)ρB00†ᴺ)、−tr(B10B01B00ᴺ⁻¹ρB00†ᴺ)（机器精度）")
    println("       C: I_exact × 离散骨架 的 ρ(T)（信息项：组合本身是 O(δt²) 离散化，非机器精度）")
    gE = 0.01
    N3, T3 = 3, 1.2
    h3 = T3 / N3
    eE = coherent_IF(N3, T3; g = gE)
    Bk = [cellB(h3; g = gE)[2s1+1:2s1+2, 2s0+1:2s0+2] for s1 in 0:1, s0 in 0:1]
    rbd = Matrix(rhob())
    B00ᴺ = Bk[1,1]^N3; B00ᴹ = Bk[1,1]^(N3 - 1); B00dᴺ = (Bk[1,1]')^N3
    v_same_dir = tr(B00ᴺ * rbd * B00dᴺ) - tr(B00ᴹ * Bk[2,2] * rbd * B00dᴺ)
    v_cross_dir = -tr(Bk[2,1] * Bk[1,2] * B00ᴹ * rbd * B00dᴺ)
    v_same = kval(eE, pos(1, :+, true), pos(0, :+, false))
    v_cross = kval(eE, pos(3, :+, true), pos(1, :+, false))
    dE = max(abs(v_same - v_same_dir), abs(v_cross - v_cross_dir))
    for (N, T) in cases
        IFv = fock_IF(N, T)
        dA = norm(fock_rho(N, T, IFv) - exact_rho(N, T))
        e0 = coherent_IF(N, T; g = 0.0)
        nzc = count(!iszero, e0) - 1
        edata = coherent_IF(N, T)
        dC = norm(rho_from_IF(edata, N, T) - exact_rho(N, T))
        rf = [exp(-im * T * EPSD * (a - b)) * RHO0[a+1, b+1] for a in 0:1, b in 0:1]
        dD = norm(rho_from_IF(e0, N, T) - rf)
        @printf("  N=%d T=%.2f:  A=%.2e  B(非常数项数)=%d  C=%.2e  D=%.2e\n",
                N, T, dA, nzc, dC, dD)
        ok &= dA < 1e-13 && nzc == 0 && dD < 1e-12
    end
    @printf("  E (g=%.2f):  同cell %s vs 闭式 %s (|Δ|=%.1e)；跨cell %s vs 闭式 %s (|Δ|=%.1e)\n",
            gE, creim(v_same), creim(v_same_dir), abs(v_same - v_same_dir),
            creim(v_cross), creim(v_cross_dir), abs(v_cross - v_cross_dir))
    ok &= dE < 1e-12
    ok || error("自检未通过：链侧构造有误，后续对比无意义")
    println("  自检 A/B/D/E 全部通过（机器精度）。")

    # ---------- IF 逐单项式对比 ----------
    println("\n[主检验] 严格链 IF vs QuAPI GrassmannIF（逐单项式比较，链侧腿子集 2^{4N} 项）")
    println("  (A) 直接对齐（包内约定）；(B) ā 腿左移一格对齐（同 cell 系数 ↔ QuAPI 对角元）")
    for (N, T) in cases
        r = compare_case(N, T)
        res[(N, T)] = r
        @printf("  N=%d T=%.2f（%d 个链侧单项式）:  max‖Δ‖(A)=%.3e   max‖Δ‖(B)=%.3e   超 1e-12 项数(B)=%d\n",
                N, T, r.nmon, r.maxA, r.maxB, r.nbadB)
        ds = sort(collect(r.degmax))
        println("    按次数(B): " * join(["$(d)次: $(sprintf2(v))" for (d, v) in ds], "  "))
    end

    # ---------- O(g²) 核对照（小 g，N=3） ----------
    println("\n[核对照] g=0.05（O(g⁴)≈6e-9 可忽略）：二次系数与常数项（(B) 对齐，均取排序单项式约定）")
    g2 = 0.05
    N, T = 3, 1.2
    h = T / N
    r2 = compare_case(N, T; g = g2)
    println("  —— 同 cell G₊₊（链侧 ā⁺(k)a⁺(k−1) ↔ QuAPI ā⁺(k−1)a⁺(k−1) 对角元）——")
    println("   解析（O(g²) 主阶）: 链 = g²[n̄T(ε₀) − n·δt²/2]（= tr((B00−B11)ρ)），QuAPI = g²[n̄T − nT*(ε₀)]")
    println("   （n 部分仅差 g²n[T*(ε₀) − δt²/2] = O(g²δt³) 相位项；n̄ 部分一致）")
    for k in 1:N
        ex = kval(r2.edata, pos(k, :+, true), pos(k - 1, :+, false))
        qu = kval(r2.qdata, pos(k - 1, :+, true), pos(k - 1, :+, false))
        @printf("   k=%d:  链=%s  QuAPI=%s  |Δ|=%.2e\n", k, creim(ex), creim(qu), abs(ex - qu))
    end
    println("  —— 跨 cell G₊₊（链侧 ā⁺(k)a⁺(j)，分离 m 个 cell；双方 O(g²) = +n̄g²R(m)）——")
    for (k, j) in ((3, 1), (3, 0), (2, 0))
        m = k - j - 1
        ex = kval(r2.edata, pos(k, :+, true), pos(j, :+, false))
        qu = kval(r2.qdata, pos(k - 1, :+, true), pos(j, :+, false))
        @printf("   (k,j)=(%d,%d) m=%d:  链=%s  QuAPI=%s  解析=%s  |Δ|=%.2e\n",
                k, j, m, creim(ex), creim(qu), creim(nbar() * g2^2 * Rrect(EPS0, m, h)), abs(ex - qu))
    end
    println("  —— 常数项（主导差异源）——")
    Bk2 = [cellB(h; g = g2)[2s1+1:2s1+2, 2s0+1:2s0+2] for s1 in 0:1, s0 in 0:1]
    leak = 1 - real(tr((Bk2[1,1]^N) * Matrix(rhob()) * ((Bk2[1,1]^N)')))
    @printf("   链 = %s  （= tr(B00ᴺρB00†ᴺ)，比 1 低 O(N·g²δt²)=%.2e：cell 内跳回概率泄漏）\n",
            creim(r2.edata[1]), leak)
    @printf("   QuAPI = %s（严格 = 1）  |Δ|=%.2e\n", creim(r2.qdata[1]), abs(r2.edata[1] - r2.qdata[1]))
    println("  —— 同节点 / 交叉分支示例（数值对照）——")
    ex = kval(r2.edata, pos(2, :+, true), pos(2, :+, false))     # 同节点 ā⁺(2)a⁺(2)
    qu = kval(r2.qdata, pos(1, :+, true), pos(2, :+, false))
    @printf("   同节点 ā⁺(2)a⁺(2):  链=%s  QuAPI(B)=%s  |Δ|=%.2e\n", creim(ex), creim(qu), abs(ex - qu))
    ex = kval(r2.edata, pos(2, :+, true), pos(1, :-, false))     # 交叉分支
    qu = kval(r2.qdata, pos(1, :+, true), pos(1, :-, false))
    @printf("   交叉 ā⁺(2)a⁻(1):   链=%s  QuAPI(B)=%s  |Δ|=%.2e\n", creim(ex), creim(qu), abs(ex - qu))

    # ---------- ρ(T) 层面：收敛阶 ----------
    println("\n[ρ(T) 层面] QuAPI IF 收缩（包内 (A) 约定，与 quapi_rlm_rt_gt 同型）vs 严格链")
    T0 = 0.9
    prev = NaN
    for N in 1:3
        qif, _ = quapi_IF(N, T0)
        err = norm(rho_from_IF(qif.data, N, T0) - exact_rho(N, T0))
        @printf("  N=%d (δt=%.3f):  ‖ρ_QuAPI − ρ_链‖_F = %.3e%s\n", N, T0 / N, err,
                isnan(prev) ? "" : @sprintf("   (误差比 %.2f)", prev / err))
        prev = err
    end

    # ---------- 结论 ----------
    println("\n" * "=" ^ 86)
    exact_match = all(res[(N, T)].maxB < 1e-12 for (N, T) in cases)
    if exact_match
        println("结论：QuAPI.Dense 费米 IF 与严格准绝热链 IF 在相干态表象逐位一致（≈机器精度）——严格。")
    else
        mA = maximum(res[(N, T)].maxA for (N, T) in cases)
        mB = maximum(res[(N, T)].maxB for (N, T) in cases)
        @printf("结论：非逐位严格。直接对齐(A) max=%.3e；几何对齐(B) max=%.3e。差异来源：\n", mA, mB)
        println("  (i) 常数项（主导，O(N·g²δt²)）：链 IF = tr(B00ᴺρB00†ᴺ) = 1 − N·λ —— cell 内")
        println("      杂质↔bath 跳去又跳回的虚过程对节点相干插入不可见（概率泄漏 λ），QuAPI 恒为 1；")
        println("  (ii) 二次核：同分支（G₊₊ 型）O(g²) 主阶逐项一致——跨 cell 吻合至 O(g⁴)，同 cell 吻合至")
        println("      O(g²δt³)（n 部分相位 T*(ε₀) vs 无相位 δt²/2）；混合分支（G₊₋/G₋₋ 型）系数带")
        println("      相位/符号约定差；同节点系数主阶反号；")
        println("  (iii) 四次及以上单项式 O(g⁴) 的非 Gaussian 内容（链的 B 因子含 cell 内模混合，")
        println("      而 exp(ΣāΔa) 为 Wick/Gaussian 因子化）。")
        println("  ⟹ 费米杂质算符 d†/d 无法像玻色对角耦合那样在路径扇内冻结为常数，cell 内虚过程使")
        println("     离散 IF 带内在 O(g²δt²)/cell 误差（随 δt² 逐 cell 消失）；ρ(T) 层面保持一阶收敛")
        println("     （上表误差比 ≈ 2.6/1.7 → 2）。区别于玻色对角耦合情形（diagcoupling.jl：逐位一致）。")
    end
    println("=" ^ 86)
    return exact_match
end

# 小工具：格式化浮点（供按次数行使用）
sprintf2(v) = @sprintf("%.1e", v)

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
