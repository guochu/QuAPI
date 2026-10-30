# =============================================================================
# exactif/brute.jl —— 有限热库多体暴力参照（定义"真值"）
#
# 用途：Ng–Reichman 严格核（连续热库）的逐层验证锚点。
# 单自旋扇区（两自旋 Anderson 的 IF 按自旋因子化，浴完全相同）。
#
# Trotter 链（Strang，杂质半步在最外）：
#   Û_tot = imp(½)·W(Δt)·imp(Δt)·W(Δt)···imp(Δt)·W(Δt)·imp(½)
# 其中 W = exp(−iΔtĤ₂)，Ĥ₂ = Σ_k E_k ĉ_k†ĉ_k + Σ_k(V_k d̂†ĉ_k + h.c.)
#（杂质 onsite ε 放进 H₁，即核的单粒子矩阵取 ε=0 —— 论文默认分割）。
# 顶点布局（每分支每步 2 个顶点，共 4N 个 Grassmann 变量，式 (S39)）：
#   正向（ket）支：α_n⁺ = s⁺_{n+}（时刻 nΔt，踢前）  n=0..N−1
#                  β_m⁺ = s⁺_{m−}（时刻 mΔt，踢后）  m=1..N
#   反向（bra）支：α_n⁻, β_m⁻ 同上。
#   W(Δt) 连接 α_n → β_{n+1}；杂质演化连接 β_m → α_m（内步全步，
#   两端 s_0 → α_0 与 β_N → s_N 为半步）。
#
# 路径值（IF 的计算基定义，无任何 Grassmann/相干态假设）：
#   I[S] = Tr_B[ 𝒲⁺(S⁺) · ρ_B · 𝒲⁻(S⁻)† ]，
#   𝒲⁺ = X[β_N⁺, α_{N−1}⁺]·X[β_{N−1}⁺, α_{N−2}⁺]···X[β_1⁺, α_0⁺]  （晚者居左）
#   𝒲⁻ = X[β_1⁻, α_0⁻]·X[β_2⁻, α_1⁻]···X[β_N⁻, α_{N−1}⁻]        （早者居左）
#   X[s,s'] ≡ ⟨s|Ŵ(Δt)|s'⟩（杂质占据态矩阵元，2^M×2^M 浴算符）。
#
# Grassmann 高斯结构（论文式 (S4)/(S7) 的数值自检，本文件不预设）：
#   I[S] =? c · Pf(Θ[S,S])，c = I[∅]，Θ[u,v] = I[{u,v}]/c。
#   由 p(η) = Π_{u<v}(1+Θ_{uv}η_uη_v)·c 展开系数即 Pfaffian 子阵。
# =============================================================================
using LinearAlgebra: eigen, exp, diagm, Hermitian, I, tr, det

# ------------------------------ Pfaffian ------------------------------ #

"""Pfaffian（Wimmer 型反对称高斯消元，含选主元；输入被破坏，用 pfaffian 留副本）。"""
function pfaffian!(A::AbstractMatrix)
    n = size(A, 1)
    n == 0 && return one(eltype(A))
    isodd(n) && return zero(eltype(A))
    pf = one(eltype(A))
    for k in 1:2:n-1
        p = k + argmax(abs.(A[k, k+1:n]))
        abs(A[k, p]) < 1e-300 && return zero(eltype(A))
        if p != k + 1
            A[k+1, :], A[p, :] = A[p, :], A[k+1, :]
            A[:, k+1], A[:, p] = A[:, p], A[:, k+1]
            pf = -pf
        end
        τ = A[k, k+1]
        pf *= τ
        # 消元更新（手算 4×4 校验：A'[3,4] = f + (cd−be)/a ⟹ Pf = af−be+cd ✓）
        for i in k+2:n, j in k+2:n
            A[i, j] -= (A[i, k+1] * A[k, j] - A[i, k] * A[k+1, j]) / τ
        end
    end
    return pf
end
pfaffian(A::AbstractMatrix) = pfaffian!(copy(A))

"""小尺寸 Pfaffian 展开参照（n ≤ 8），用于检验 pfaffian! 的正确性。"""
function pfaffian_expand(A::AbstractMatrix)
    n = size(A, 1)
    n == 0 && return one(eltype(A))
    isodd(n) && return zero(eltype(A))
    return _pf_expand(A, collect(1:n))
end
function _pf_expand(A, idx)
    length(idx) == 0 && return one(eltype(A))
    i = idx[1]
    s = zero(eltype(A))
    for (c, j) in enumerate(idx[2:end])
        A[i, j] == 0 && continue
        rest = [x for x in idx[2:end] if x != j]
        # 置换 (i, j, rest...) 相对升序的符号：把 j 移到第二位的相邻对换数 = c
        s += (-1)^(c - 1) * A[i, j] * _pf_expand(A, rest)
    end
    return s
end

# ------------------------------ 多体机构 ------------------------------ #

"""
    mb_generator(h) -> Matrix

单粒子矩阵 h（n×n）对应的多体生成元 G = Σ_ij h_ij f̂_i† f̂_j（2^n×2^n）。
基：占据数基，模 1 = 位 1（Jordan–Wigner 序，任意固定约定均可，
因只用整体矩阵元）。
"""
# Jordan–Wigner 约定：f̂_j = s_j a_j（a_j = 硬核玻色湮灭，s_j = 前 j−1 位的奇偶），
# f̂_i†f̂_j 作用在 |J⟩ 上的符号 = (−1)^(sj+si)，其中 sj 用 J 的位算、si 用
# c_j 作用后的态的位算（i ≠ j 时与 J 相同；i = j 时该项为 n̂_j 对角）。
function mb_generator(h::AbstractMatrix)
    n = size(h, 1)
    dim = 1 << n
    G = zeros(ComplexF64, dim, dim)
    for J in 0:dim-1
        for j in 1:n                       # c_j（湮灭，作用在浴/杂质模 j）
            (J >> (j - 1)) & 1 == 0 && continue
            sj = count_ones(J & ((1 << (j - 1)) - 1))   # c_j 的 JW 符号
            Jj = J ⊻ (1 << (j - 1))
            for i in 1:n                   # c_i†
                (Jj >> (i - 1)) & 1 == 1 && continue
                si = count_ones(Jj & ((1 << (i - 1)) - 1))  # c_i† 的 JW 符号
                Istate = Jj ⊻ (1 << (i - 1))
                h[i, j] == 0 && continue
                G[Istate+1, J+1] += (-1)^(sj + si) * h[i, j]
            end
        end
    end
    return G
end

"""
    single_particle_h(Es, Vs) -> (1+M)×(1+M)

杂质（位 1，ε=0）+ 浴；h[1,1+k] = V_k，h[1+k,1] = conj(V_k)，h[1+k,1+k] = E_k δ。
"""
function single_particle_h(Es::AbstractVector{<:Real}, Vs::AbstractVector{<:Number})
    M = length(Es)
    h = zeros(ComplexF64, M + 1, M + 1)
    for k in 1:M
        h[1, 1+k] = Vs[k]
        h[1+k, 1] = conj(Vs[k])
        h[1+k, 1+k] = Es[k]
    end
    return h
end

"""
    impurity_blocks(Wmb, M) -> X（4×... × 2^M×2^M 的函数闭包）

X[s,s'] = ⟨s|Ŵ_mb|s'⟩：把 Wmb 按杂质位（位 1）分块。
返回函数 (s::Int, sp::Int) -> Matrix（s,sp ∈ {0,1} 占据）。
"""
function impurity_blocks(Wmb::AbstractMatrix, M::Int)
    dimb = 1 << M
    X = Dict{Tuple{Int,Int},Matrix{ComplexF64}}()
    for s in 0:1, sp in 0:1
        blk = zeros(ComplexF64, dimb, dimb)
        for b1 in 0:dimb-1, b2 in 0:dimb-1
            # 模序：mode 1 = 杂质（bit 0），modes 2..M+1 = 浴（bits 1..M）
            blk[b1+1, b2+1] = Wmb[((b1 << 1) | s) + 1, ((b2 << 1) | sp) + 1]
        end
        X[(s, sp)] = blk
    end
    return X
end

# ------------------------------ 路径值 / Θ 提取 ------------------------------ #

"""
    if_full_tensor(Es, Vs, β, N, Δt) -> NamedTuple

**完整过程张量**：I[{s}]（每顶点一个二值 s ∈ {0,1}，共 4N 个，
(S39) 顶点序）。I[{s}] = Tr_B[𝒲⁺({s⁺})·ρ_B·𝒲̃⁻({s⁻})]，
其中每支的 W 链因子每元胞一项：X[β_m^·, α_{m−1}^·]（β 作行、α 作列）。
矩阵元视角：I[{s}] = ⟨s⁻ 轮廓|𝒯|s⁺ 轮廓⟩，𝒯 是 2N 个轮廓模上的算符。

顶点位（同 (S39)）：β_m⁺ 在 4m−1，α_{m−1}⁺ 在 4m−3；
β_m⁻ 在 4m，α_{m−1}⁻ 在 4m−2。
"""
function if_full_tensor(Es, Vs, β, N, Δt)
    M = length(Es)
    h = single_particle_h(Es, Vs)
    Wmb = exp(-im * Δt * mb_generator(h))
    X = impurity_blocks(Wmb, M)
    nb = [1 / (1 + exp(β * Es[k])) for k in 1:M]
    ρB = diagm(0 => [prod((1 - nb[k], nb[k])[(j >> (k - 1)) & 1 + 1] for k in 1:M)
                     for j in 0:(1<<M)-1])
    nvert = 4N
    val(S, p) = (S >> (p - 1)) & 1
    Ivals = Vector{ComplexF64}(undef, 1 << nvert)
    for S in 0:(1 << nvert) - 1
        # ket 支：𝒲⁺ = X[β_N⁺,α_{N−1}⁺]···X[β₂⁺,α₁⁺]·X[β₁⁺,α₀⁺]（最晚居左：
        # Û_tot 的最右因子先作用 ⟹ 最早胞元在最右）
        Wp = Matrix{ComplexF64}(I, 1 << M, 1 << M)
        for m in 1:N
            Wp = X[(val(S, 4m - 1), val(S, 4m - 3))] * Wp
        end
        # bra 支：𝒲̃⁻ = X[β₁⁻,α₀⁻]†·X[β₂⁻,α₁⁻]†···（最早居左：右侧乘法先作用最早）
        Wm = Matrix{ComplexF64}(I, 1 << M, 1 << M)
        for m in 1:N
            Wm = Wm * X[(val(S, 4m), val(S, 4m - 2))]'
        end
        Ivals[S+1] = tr(Wp * ρB * Wm)
    end
    return (; Ivals, X, ρB, Wmb, nvert)
end

"""
    if_path_values(Es, Vs, β, N, Δt) -> NamedTuple

过程张量 + 对角配对核提取：`c` = I[∅]，`Θ`[u,v] = I[{u,v}]/c。
"""
function if_path_values(Es, Vs, β, N, Δt; return_more::Bool=false)
    full = if_full_tensor(Es, Vs, β, N, Δt)
    Ivals = full.Ivals
    nvert = 4N
    c = Ivals[1]
    Θ = zeros(ComplexF64, nvert, nvert)
    for u in 1:nvert, v in u+1:nvert
        Θ[u, v] = Ivals[((1 << (u - 1)) | (1 << (v - 1))) + 1] / c
        Θ[v, u] = -Θ[u, v]
    end
    out = (; Ivals, c, Θ)
    return_more && (out = merge(out, (; full)))
    return out
end

"""
    brute_ABZN(Es, Vs, β, N, Δt)

直接由有限浴单粒子矩阵计算论文小矩阵（独立于连续机制）：
- `F1` = F(Δt) = [W]₁₁
- `Ad`，`Bd`：A[i,j] = r⁽ⁱ⁾·v⁽ʲ⁾，B[i,j] = r⁽ⁱ⁾·f·v⁽ʲ⁾，
  r⁽ⁱ⁾ = a·U_Bⁱ（行），v⁽ʲ⁾ = (U_B†)ʲ·ā（列），a = [W]₁,₂:，U_B = [W]₂:,₂:。
- `ZNoZ` = det(I − Bd)（应等于 if_path_values 的 c）。
"""
function brute_ABZN(Es, Vs, β, N, Δt)
    h = single_particle_h(Es, Vs)
    W = exp(-im * Δt * h)
    M = length(Es)
    F1 = W[1, 1]
    a = vec(W[1, 2:end])               # 列向量（M）
    UB = W[2:end, 2:end]
    r = [vec(transpose(a) * (UB^i)) for i in 0:N-1]   # 行向量：r⁽ⁱ⁾ = aᵀ U_Bⁱ（无共轭！）
    v = [vec((UB' )^j * conj.(a)) for j in 0:N-1]    # 列向量：v⁽ʲ⁾ = (U_B†)ʲ ā
    f = [1 / (1 + exp(β * Es[k])) for k in 1:M]
    Ad = [sum(r[i] .* v[j]) for i in 1:N, j in 1:N]          # 双线性 Σ_k r_k v_k
    Bd = [sum(r[i] .* (f .* v[j])) for i in 1:N, j in 1:N]
    return (; F1, a, UB, Ad, Bd, ZNoZ=det(I - Bd), r, v, W)
end

# ------------------------------ Trotter 化精确参照 ------------------------------ #

"""
    trotter_exact(Es, Vs, β, N, Δt; εs=(0.0,0.0), U=0.0, ρ0=nothing)
        -> ρ_imp（4×4，基 |0⟩,|↑⟩,|↓⟩,|↑↓⟩）

两自旋、每自旋独立浴的 Trotter 化全动力学（多体矩阵直乘，机器精度）。
ρ0：4×4 初始杂质密度矩阵（默认 |0⟩⟨0|）。
"""
function trotter_exact(Es, Vs, β, N, Δt; εs=(0.0, 0.0), U=0.0, ρ0=nothing)
    M = length(Es)
    hσ = single_particle_h(Es, Vs)
    # 模序：d↑, bath↑…, d↓, bath↓…
    h = zeros(ComplexF64, 2 * (M + 1), 2 * (M + 1))
    h[1:M+1, 1:M+1] .= hσ
    h[M+2:end, M+2:end] .= hσ
    Wmb = exp(-im * Δt * mb_generator(h))
    # 踢：对角相位 ε_↑n↑ + ε_↓n↓ + U n↑n↓（d↑ = 模 1，d↓ = 模 M+2）
    dim = 1 << (2 * (M + 1))
    Eimp = [εs[1] * ((J >> 0) & 1) + εs[2] * ((J >> (M + 1)) & 1) +
            U * ((J >> 0) & 1) * ((J >> (M + 1)) & 1) for J in 0:dim-1]
    Kfull = diagm(0 => exp.(-im * Δt .* Eimp))
    Khalf = diagm(0 => exp.(-(im * Δt / 2) .* Eimp))
    Utot = Khalf * Wmb * (Kfull * Wmb)^(N - 1) * Khalf
    # 热浴（两自旋独立）
    nb = [1 / (1 + exp(β * Es[k])) for k in 1:M]
    ρB1 = diagm(0 => [prod(((1 - nb[k], nb[k])[(j >> (k - 1)) & 1 + 1]) for k in 1:M)
                      for j in 0:(1<<M)-1])
    ρB = kron(ρB1, ρB1)
    ρ0 === nothing && (ρ0 = zeros(ComplexF64, 4, 4); ρ0[1, 1] = 1.0)
    # 初态：ρ0 ⊗ ρB，嵌入到全模序（d↑,…,d↓,…）：杂质占据 = 位0(d↑) 与位 M+1(d↓)
    # 全基 (J↑: 杂质+bath↑ 的 (M+1) 位) × (J↓): 重新排为 [d↑ | bath↑ | d↓ | bath↓]
    ρtot = _embed_impurity_bath(ρ0, ρB, M)
    ρt = Utot * ρtot * Utot'
    return _traceout_bath(ρt, M)
end

# ρ0(4×4) ⊗ ρB(4^M×4^M) → 全基（模序 d↑,bath↑(M),d↓,bath↓(M)）矩阵
function _embed_impurity_bath(ρ0, ρB, M)
    dim = 1 << (2 * (M + 1))
    out = zeros(ComplexF64, dim, dim)
    dimb = 1 << M
    for iup in 0:1, idn in 0:1
        for bup in 0:dimb-1, bdn in 0:dimb-1
            row = iup | (bup << 1) | (idn << (M + 1)) | (bdn << (M + 2))
            for jup in 0:1, jdn in 0:1
                for cup in 0:dimb-1, cdn in 0:dimb-1
                    col = jup | (cup << 1) | (jdn << (M + 1)) | (cdn << (M + 2))
                    out[row+1, col+1] = ρ0[2*iup + idn + 1, 2*jup + jdn + 1] *
                                        ρB[bup*dimb + bdn + 1, cup*dimb + cdn + 1]
                end
            end
        end
    end
    return out
end

function _traceout_bath(ρt, M)
    dimb = 1 << M
    out = zeros(ComplexF64, 4, 4)
    for iup in 0:1, idn in 0:1, jup in 0:1, jdn in 0:1
        s = 0.0im
        for bup in 0:dimb-1, bdn in 0:dimb-1
            row = iup | (bup << 1) | (idn << (M + 1)) | (bdn << (M + 2))
            col = jup | (bup << 1) | (jdn << (M + 1)) | (bdn << (M + 2))
            s += ρt[row+1, col+1]
        end
        out[2*iup + idn + 1, 2*jup + jdn + 1] = s
    end
    return out
end
