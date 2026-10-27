# =============================================================================
# QuAPI 费米 IF 严格性验证（虚时）—— 单模费米 bath：严格准绝热环链 IF vs QuAPI.Dense
#
# 背景与 realtime 版（fermionicIF.jl）的结论：费米 QuAPI IF 非逐位严格——
#   差异为 (i) 常数项 O(N·g²δt²) 的 cell 内"跳去又跳回"概率泄漏、(ii) 二次核
#   O(g²δt³) 相位项、(iii) O(g⁴) 非 Gaussian 高次项；局标二阶（O(δt²)/cell）、
#   定总时长累积为一阶（O(δt)），ρ(T) 一阶收敛。
# 本脚本在虚时（τ 环，结构更简单：无边界腿、无窗口舍弃、链侧 2Nτ 条腿与
#   GrassmannIF(fermionic_Δτ) 的 basis 一一对应）重复同样的检验：
#
# 模型  H = ε_d d†d + ε₀ c†c + g(d†c + c†d)，单 bath 模（谱 g²δ(ε−ε₀)），μ=0
# 链    Z = Tr[(B_τ(δτ) A_τ(δτ))^{Nτ}]，δτ = β/Nτ，A_τ = e^{−δτ ε_d d†d}，
#       B_τ = expm[−δτ(ε₀ c†c + g(d†c+c†d))]   ←—— 杂质-环境耦合项的（虚时）演化算符，严格计算
# IF    节点相干插入 + 环形 bath 迹：
#       I_exact = Tr_b[Π_{k=1..Nτ}⟨η̄_k|B_τ|η_{k−1}⟩] / (Π_k(1+ā_k a_{k−1}) · Z_bath)
#       （节点指标环回：cell k 的 in 节点为 k−1 mod Nτ；除自由重叠与 bath 配分
#         Z_bath = 1+e^{−βε₀} 归一 ⟹ I(g=0)≡1，与 QuAPI IF 的归一化一致）
# 收缩骨架（与包内 test/quapi/imag.jl 的 quapi_rlm_gt 同型）：
#       Z = ∫Π_k[dā_k da_k e^{−ā_k a_k}] · Π_{k≥2}(1+y ā_k a_{k−1}) · (1−y ā_1 a_N)
#           · I ,   y = e^{−δτ ε_d}（−y 回绕 = 费米迹的反周期约定）
#
# 自检（须全部机器精度通过，否则报错）：
#   A'. Fock 环路径收缩 Z == 直接 4 维严格演化 Tr[(B_τA_τ)^{Nτ}]
#   B'. g=0 ⇒ I ≡ 1（重叠/Z_bath 归一正确）
#   D'. g=0 ⇒ 骨架全积分 == 自由杂质配分 1+e^{−βε_d}（含反周期回绕符号）
#   E'. 小 g 二次核 == 块迹闭式 tr((B00−B11)B00^{Nτ−1})/Z_bath（同 cell）、
#       ∓tr(B10·B00^{m−1}·B01·B00^{Nτ−m−1})/Z_bath（跨 cell，分离 m）
#   C'.（信息项）骨架×I_exact 的 Z：分立骨架本身是 O(δτ²) 离散化，非机器精度
#
# 主检验：I_exact vs GrassmannIF(fermionic_Δτ) 全 2^{2Nτ} 系数逐单项式比较：
#   (A) 直接对齐（节点腿 ↔ Δτ 行列）；(B) ā 腿环向左移一格（几何对齐：
#       链侧 ā(node k) 来自以节点 k 为右端的 cell ⟺ QuAPI 行 k−1，同 cell 系数
#       ↔ QuAPI 对角元；环回：ā₁ ↦ ā_N）
#
# 运行：julia --startup-file=no --project=../test performance/fermionicIFtau.jl
# =============================================================================
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf

# ------------------------------ 模型参数 ------------------------------
const EPSD = 0.6                                  # 杂质能级 ε_d
const EPS0 = 0.9                                  # 单 bath 模能量 ε₀（谱 g²δ(ε−ε₀)）
const GC   = 0.35                                 # 杂质-环境耦合 g
const BETA = 1.2                                  # 逆温度（μ=0）

nocc() = 1 / (exp(BETA * EPS0) + 1)               # bath 模占据 n
nbar() = 1 - nocc()
Zbath() = 1 + exp(-BETA * EPS0)                   # 自由 bath 配分（IF 归一用）

# ------------------------------ 双模费米代数（4 维 Fock） ------------------------------
# 基底 |n_d n_c⟩ = (d†)^{n_d}(c†)^{n_c}|0⟩，线性指标 = 2n_d + n_c + 1；d† 无 JW 串，
# c† 带 (−1)^{n_d}。
function cellBtau(dτ; g = GC)
    ddag = zeros(ComplexF64, 4, 4); ddag[3, 1] = 1; ddag[4, 2] = 1
    cdag = zeros(ComplexF64, 4, 4); cdag[2, 1] = 1; cdag[4, 3] = -1
    dd, c = ddag', cdag'
    return exp(-dτ * (EPS0 * (cdag * c) + g * (ddag * c + cdag * dd)))
end
cellAtau(dτ) = kron(exp(-dτ * EPSD * ComplexF64[0 0; 0 1]), Matrix(I, 2, 2))

# ------------------------------ 外代数（Grassmann）工具 ------------------------------
# 腿位（1 基，与 QuAPI 虚时 basis 顺序一致：节点 j 的 (a, ā) 相邻）：pos(j, ā?)
pos(j::Int, ā::Bool) = 2 * j - (ā ? 0 : 1)

function mono_from_seq(seq::Vector{Int})
    s = 0
    for i in 1:length(seq), j in (i + 1):length(seq)
        seq[i] > seq[j] && (s += 1)
    end
    mask = 0
    for p in seq; mask |= 1 << (p - 1); end
    return mask, isodd(s) ? -1 : 1
end

function merge_sign(A::Int, B::Int)
    s = 0
    while B > 0
        lb = B & -B; B -= lb
        s += count_ones(A & ~(2 * lb - 1))
    end
    return isodd(s) ? -1 : 1
end

function sf(terms)
    out = Tuple{Int,ComplexF64}[]
    for (seq, c) in terms
        c = ComplexF64(c)
        if isempty(seq)
            push!(out, (0, c))
        else
            m, s = mono_from_seq(seq)
            push!(out, (m, c * s))
        end
    end
    return out
end

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

kval(data::Vector{<:Number}, ps::Int...) =
    data[foldl((a, p) -> a | (1 << (p - 1)), ps; init = 0) + 1]

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

# 单 cell 相干矩阵元 ⟨η̄_o|B_τ|η_i⟩[b',b] = Σ_{s',s} (−ā_o)^{s'} B[s',s][b',b] (−a_i)^s
function cellket(Blk, o::Int, i::Int)
    f = Matrix{Dict{Int,ComplexF64}}(undef, 2, 2)
    for b1 in 1:2, b2 in 1:2
        d = Dict{Int,ComplexF64}()
        for s1 in 0:1, s0 in 0:1
            c = Blk[s1+1, s0+1][b1, b2]
            iszero(c) && continue
            seq = Int[]
            s1 == 1 && push!(seq, pos(o, true))
            s0 == 1 && push!(seq, pos(i, false))
            m, s = isempty(seq) ? (0, 1) : mono_from_seq(seq)
            d[m] = get(d, m, 0.0im) + c * ((s1 + s0) % 2 == 1 ? -1 : 1) * s
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

# ------------------------------ 严格链：Z、Fock 环收缩、相干 IF ------------------------------
Zdirect(N::Int; g = GC) = tr((cellBtau(BETA / N; g = g) * cellAtau(BETA / N))^N)

# 自检 A' 用的 Fock 环路径收缩（杂质占据路径枚举 + bath 块环乘求迹）
function fock_Z_tau(N::Int; g = GC)
    dτ = BETA / N
    Blk = let B = cellBtau(dτ; g = g)
        [B[2s1+1:2s1+2, 2s0+1:2s0+2] for s1 in 0:1, s0 in 0:1]
    end
    y = exp(-dτ * EPSD)
    Z = 0.0 + 0.0im
    for mask in 0:(1 << N)-1
        sk(k) = (mask >> (k - 1)) & 1            # bit k−1 = 节点 k 的占据 s_k
        M = Matrix{ComplexF64}(I, 2, 2)
        for k in N:-1:1
            M = M * Blk[sk(k)+1, sk(mod1(k - 1, N))+1]
        end
        amp = y ^ sum(sk(k) for k in 1:N)        # Π_k e^{−δτε_d s_{k−1}}（环上等价于 Π e^{−δτε_d s_k}）
        Z += amp * tr(M)
    end
    return Z
end

# 严格 IF（相干态表象，环）：Tr_b[Π⟨η̄_k|B|η_{k−1}⟩] / (Π(1+ā_k a_{k−1})·Z_bath)
function coherent_IF_tau(N::Int; g = GC)
    dτ = BETA / N
    Blk = let B = cellBtau(dτ; g = g)
        [B[2s1+1:2s1+2, 2s0+1:2s0+2] for s1 in 0:1, s0 in 0:1]
    end
    I2 = [i == j ? ONEP : Dict{Int,ComplexF64}() for i in 1:2, j in 1:2]
    X = I2
    for k in N:-1:1
        X = pmmul(X, cellket(Blk, k, mod1(k - 1, N)))
    end
    tot = merge(+, X[1, 1], X[2, 2])             # 纯 bath 迹（无 ρ_b）
    L = 2 * N
    data = zeros(ComplexF64, 1 << L)
    for (m, c) in tot
        data[m+1] = c
    end
    one_ = 1.0 + 0.0im
    for k in 1:N                                 # 除自由重叠（环回：cell 1 的 a 腿 = 节点 N）
        data = gds(data, sf([(Int[], one_), ([pos(k, true), pos(mod1(k - 1, N), false)], -one_)]))
    end
    return data ./ Zbath()                       # bath 配分归一 ⟹ I(g=0)≡1
end

# ------------------------------ 收缩骨架：IF → Z（包内 quapi_rlm_gt 同型） ------------------------------
function Z_from_IF(ifdata::Vector{<:Number}, N::Int)
    dτ = BETA / N
    y = exp(-dτ * EPSD)
    one_ = 1.0 + 0.0im
    P = ComplexF64.(ifdata)
    for k in 2:N                                  # 前向动能 1 + y·ā_k a_{k−1}
        P = gds(P, sf([(Int[], one_), ([pos(k, true), pos(k - 1, false)], y)]))
    end
    P = gds(P, sf([(Int[], one_), ([pos(1, true), pos(N, false)], -y)]))  # 反周期回绕
    for j in 1:N                                  # Berezin 测度 e^{−ā_j a_j}
        P = gds(P, sf([(Int[], one_), ([pos(j, true), pos(j, false)], -one_)]))
    end
    return P[end]                                 # 顶单项式系数
end

# ------------------------------ QuAPI 侧 ------------------------------
function quapi_IF_tau(N::Int; g = GC)
    corr = fermionic_Δτ(DiracDelta(EPS0; α = g^2), β = BETA, N = N, μ = 0.0)
    hmat = GrassmannHybridizationMatrix(corr)
    for (k, leg) in enumerate(hmat.basis)         # basis 顺序与 pos 约定一致性自检
        (leg.branch === :τ && leg.band == 1 && pos(leg.j, leg.conj) == k) ||
            error("QuAPI basis 顺序与脚本 pos 约定不一致")
    end
    return GrassmannIF(hmat)
end

# ------------------------------ 主对比 ------------------------------
function compare_tau(N::Int; g = GC)
    edata = coherent_IF_tau(N; g = g)
    qdata = quapi_IF_tau(N; g = g).data
    L = 2 * N
    # (B) 对齐：链侧 ā(node k) 来自以节点 k 为右端的 cell ↦ QuAPI 行 k−1（环回 ā₁ ↦ ā_N）
    shiftmap = Dict{Int,Int}()
    for j in 1:N
        shiftmap[pos(j, true) - 1] = pos(mod1(j - 1, N), true) - 1
        shiftmap[pos(j, false) - 1] = pos(j, false) - 1
    end
    maxA = 0.0; maxB = 0.0; degmax = Dict{Int,Float64}()
    nbadB = 0
    for m in 0:(1 << L)-1
        dA = abs(edata[m+1] - qdata[m+1])
        maxA = max(maxA, dA)
        qm = 0
        for b in 0:L-1
            ((m >> b) & 1) == 1 && (qm |= 1 << shiftmap[b])
        end
        dB = abs(edata[m+1] - qdata[qm+1])
        maxB = max(maxB, dB)
        d = count_ones(m)
        degmax[d] = max(get(degmax, d, 0.0), dB)
        dB > 1e-12 && (nbadB += 1)
    end
    return (edata = edata, qdata = qdata, maxA = maxA, maxB = maxB, nbadB = nbadB, degmax = degmax)
end

# 虚时解析核：三角核 T∓ 与整 cell 矩形核 K₁
Tm(ε, dτ) = (exp(-dτ * ε) - (1 - dτ * ε)) / ε^2
Tp(ε, dτ) = (exp(dτ * ε) - (1 + dτ * ε)) / ε^2
creim(z) = @sprintf("(%+.6e%+.6ei)", real(z), imag(z))

# ------------------------------ 主流程 ------------------------------
function main()
    println("=" ^ 86)
    println("费米 QuAPI IF 严格性验证（虚时 τ 环）：单模 bath 准绝热环链严格 IF vs QuAPI.Dense")
    @printf("模型: εd=%.2f ε₀=%.2f g=%.2f β=%.2f  → n=%.4f n̄=%.4f Z_bath=%.6f\n",
            EPSD, EPS0, GC, BETA, nocc(), nbar(), Zbath())
    println("=" ^ 86)

    cases = [2, 3, 4]
    ok = true

    # ---------- 自检 ----------
    println("\n[自检] A': Fock 环路径收缩 vs 直接 4 维严格演化（机器精度）")
    println("       B': g=0 ⇒ I ≡ 1；D': g=0 ⇒ 骨架 Z == 1+e^{−βε_d}（含反周期回绕符号）")
    println("       E': 小 g 二次核 vs 块迹闭式（机器精度）")
    println("       C':（信息项）骨架×I_exact 的 Z：分立骨架是 O(δτ²) 离散化")
    gE = 0.01
    NE = 4
    dτE = BETA / NE
    eE = coherent_IF_tau(NE; g = gE)
    Bk = let B = cellBtau(dτE; g = gE)
        [B[2s1+1:2s1+2, 2s0+1:2s0+2] for s1 in 0:1, s0 in 0:1]
    end
    B00ᴹ = Bk[1,1]^(NE - 1)
    v_same_dir = tr((Bk[1,1] - Bk[2,2]) * B00ᴹ) / Zbath()
    p, q = 4, 2                                   # 分离 m = p−q−1 = 1 个 cell（不跨回绕）
    msep = p - q - 1                              # p 与 q+1 之间的 cell 数 = msep−1
    v_cross_dir = -tr(Bk[2,1] * Bk[1,1]^(p - q - 2) * Bk[1,2] * Bk[1,1]^(NE - p + q)) / Zbath()
    v_same = kval(eE, pos(p, true), pos(p - 1, false))
    v_cross = kval(eE, pos(p, true), pos(q, false))
    dE = max(abs(v_same - v_same_dir), abs(v_cross - v_cross_dir))
    for N in cases
        dA = abs(fock_Z_tau(N) - Zdirect(N))
        e0 = coherent_IF_tau(N; g = 0.0)
        nzc = count(!iszero, e0) - 1
        dD = abs(Z_from_IF(e0, N) - (1 + exp(-BETA * EPSD)))
        dC = abs(Z_from_IF(coherent_IF_tau(N), N) - Zdirect(N) / Zbath())
        @printf("  Nτ=%d:  A'=%.2e  B'(非常数项数)=%d  D'=%.2e  C'=%.2e\n", N, dA, nzc, dD, dC)
        ok &= dA < 1e-13 && nzc == 0 && dD < 1e-12
    end
    @printf("  E' (g=%.2f, Nτ=%d):  同cell %s vs 闭式 %s (|Δ|=%.1e)；跨cell %s vs 闭式 %s (|Δ|=%.1e)\n",
            gE, NE, creim(v_same), creim(v_same_dir), abs(v_same - v_same_dir),
            creim(v_cross), creim(v_cross_dir), abs(v_cross - v_cross_dir))
    ok &= dE < 1e-12
    ok || error("自检未通过：链侧构造有误，后续对比无意义")
    println("  自检 A'/B'/D'/E' 全部通过（机器精度）。")

    # ---------- IF 逐单项式对比 ----------
    println("\n[主检验] 严格环链 IF vs QuAPI GrassmannIF（全 2^{2Nτ} 系数逐单项式）")
    println("  (A) 直接对齐（节点腿 ↔ Δτ 行/列）；(B) ā 腿环向左移一格（同 cell 系数 ↔ QuAPI 对角元）")
    for N in cases
        r = compare_tau(N)
        @printf("  Nτ=%d（%d 系数）:  max‖Δ‖(A)=%.3e   max‖Δ‖(B)=%.3e   超 1e-12 项数(B)=%d\n",
                N, 1 << (2N), r.maxA, r.maxB, r.nbadB)
        ds = sort(collect(r.degmax))
        println("    按次数(B): " * join(["$(d)次: $(@sprintf("%.1e", v))" for (d, v) in ds], "  "))
    end

    # ---------- 小 g 核对照 ----------
    println("\n[核对照] g=0.05（O(g⁴)≈6e-9 可忽略），Nτ=4（(B) 对齐，排序单项式约定）")
    g2 = 0.05
    N2 = 4
    dτ2 = BETA / N2
    r2 = compare_tau(N2; g = g2)
    Bk2 = let B = cellBtau(dτ2; g = g2)
        [B[2s1+1:2s1+2, 2s0+1:2s0+2] for s1 in 0:1, s0 in 0:1]
    end
    println("  —— 同 cell（链侧 ā_k a_{k−1} ↔ QuAPI ā_{k−1} a_{k−1} 对角元）——")
    println("   解析: 链 = tr((B00−B11)B00^{Nτ−1})/Z_bath；QuAPI(排序符号后) = g²[−n̄T₋ + nT₊]")
    println("        （T∓(ε)=∫₀^{δτ}ds∫₀^s ds' e^{∓ε(s−s')}，T₋+T₊ = 整 cell 双重积分 K₁）")
    for k in 2:N2
        ex = kval(r2.edata, pos(k, true), pos(k - 1, false))
        qu = kval(r2.qdata, pos(k - 1, true), pos(k - 1, false))
        @printf("   k=%d:  链=%s  QuAPI=%s  |Δ|=%.2e\n", k, creim(ex), creim(qu), abs(ex - qu))
    end
    println("  —— 跨 cell（分离 m 个 cell）——")
    for (p, q) in ((4, 2), (4, 1), (3, 1))
        m = p - q - 1
        ex = kval(r2.edata, pos(p, true), pos(q, false))
        qu = kval(r2.qdata, pos(p - 1, true), pos(q, false))
        @printf("   (p,q)=(%d,%d) m=%d:  链=%s  QuAPI=%s  |Δ|=%.2e\n",
                p, q, m, creim(ex), creim(qu), abs(ex - qu))
    end
    println("  —— 常数项（主导差异源）——")
    leak = 1 - real(tr(Bk2[1,1]^N2)) / Zbath()
    @printf("   链 = %s  （= tr(B00^{Nτ})/Z_bath，与 1 偏离 O(Nτ·g²δτ²)=%.2e：cell 内跳回泄漏）\n",
            creim(r2.edata[1]), leak)
    @printf("   QuAPI = %s（严格 = 1）  |Δ|=%.2e\n", creim(r2.qdata[1]), abs(r2.edata[1] - r2.qdata[1]))

    # ---------- Z 层面：收敛阶 ----------
    println("\n[Z 层面] QuAPI IF 骨架收缩（quapi_rlm_gt 同型）vs 严格链 Z/Z_bath")
    prev = NaN
    for N in 2:6
        errq = abs(Z_from_IF(quapi_IF_tau(N).data, N) - Zdirect(N) / Zbath())
        erre = abs(Z_from_IF(coherent_IF_tau(N), N) - Zdirect(N) / Zbath())
        @printf("  Nτ=%d (δτ=%.4f):  ‖Z_QuAPI−Z‖=%.3e%s   ‖Z_Iexact−Z‖=%.3e\n",
                N, BETA / N, errq, isnan(prev) ? "" : @sprintf("  (误差比 %.2f)", prev / errq), erre)
        prev = errq
    end

    # ---------- 结论 ----------
    println("\n" * "=" ^ 86)
    println("结论（虚时，与实时 fermionicIF.jl 一致）：")
    println("  1. 非逐位严格：即使虚时环结构无边界腿，严格链 IF 与 QuAPI GrassmannIF 仍有")
    println("     O(g²δτ²)/cell 的内在差异（按次数表 + 核对照）；")
    println("  2. 差异结构与实时相同：(i) 常数项泄漏 tr(B00^{Nτ})/Z_bath ≠ 1（cell 内跳去又")
    println("     跳回的虚过程对节点相干插入不可见）；(ii) 二次核主阶一致、相位项 O(g²δτ³)；")
    println("     (iii) 高次单项式 O(g⁴) 非 Gaussian 内容；")
    println("  3. 阶数定性：局标（单 cell）二阶 O(g²δτ²)，定 β 累积为一阶 O(g²β·δτ)——")
    println("     Z 层面实测误差比 1.41→1.18（渐近阶 ≈0.9，粗网格含高阶修正；且该误差由")
    println("     骨架分裂主导——Z_Iexact 与 Z_QuAPI 误差同量级，IF 差异见系数级对比）；")
    println("     这是费米杂质算符无法在路径扇内冻结的分段化误差，而非 Wick 因子化误差")
    println("     （费米的 Gaussian 迹是精确的）。")
    println("=" ^ 86)
    return all(compare_tau(N).maxB < 1e-12 for N in cases)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
