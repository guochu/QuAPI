# ============================================================================
# utility.jl —— 影响泛函的稠密表示与 contour 组装工具
# ============================================================================

# ------------------------------------------------------------------ #
# contour Δ（hybridization 矩阵）组装
#
# 节点行序：分支 [:+, :-, :τ] 逐分支、分支内时间升序；整条 contour 的 Δ
# 是分支块的块矩阵。元素类型由 correlation function 决定（虚时通常为实，
# 实/混合含虚部），这里不写死为复数。
# ------------------------------------------------------------------ #

_branches(::ImagCorrelationFunction) = [:τ]
_branches(::RealCorrelationFunction) = [:+, :-]
_branches(::AbstractMixedCorrelationFunction) = [:+, :-, :τ]

# 分支块矩阵 → 整块矩阵（类型自动提升）
function _assemble(brs::Vector{Symbol}, blocks::Matrix)
    dims = [size(blocks[i, 1], 1) for i in eachindex(brs)]
    offs = cumsum([0; dims])
    T = promote_type(map(eltype, blocks)...)
    M = zeros(T, offs[end], offs[end])
    for i in eachindex(brs), j in eachindex(brs)
        M[offs[i]+1:offs[i+1], offs[j]+1:offs[j+1]] .= blocks[i, j]
    end
    return M
end

function _contour_delta(x::AbstractCorrelationFunction)
    brs = _branches(x)
    blocks = [Matrix(branch(x; b1=b1, b2=b2)) for b1 in brs, b2 in brs]
    return _assemble(brs, blocks)
end

# Δ 的算符次序拆分 M1 + M2 = Δ：
# CorrelationMatrix 块区分两种等时算符次序——M1 取 ηⱼₖ（含等时 ηⱼₖ[1]），
# M2 取 ηₖⱼ（含等时 ηₖⱼ[1]）。普通矩阵块（实/混合 contour 的 ζ 交叉块）
# 无次序信息，按上三角 / 严格下三角拆分（等时项归入 M1）。
function _split_block(B::CorrelationMatrix)
    m = size(B, 1)
    T = eltype(B)
    A1 = zeros(T, m, m)
    A2 = zeros(T, m, m)
    for i in 1:m, j in 1:m
        if i > j
            A1[i, j] = B.ηⱼₖ[i-j+1]
        elseif i < j
            A2[i, j] = B.ηₖⱼ[j-i+1]
        else
            A1[i, i] = B.ηⱼₖ[1]
            A2[i, i] = B.ηₖⱼ[1]
        end
    end
    return A1, A2
end
_split_block(B::AbstractMatrix) = (Matrix(triu(B)), Matrix(tril(B, -1)))

function _contour_delta_split(x::AbstractCorrelationFunction)
    brs = _branches(x)
    parts = [_split_block(branch(x; b1=b1, b2=b2)) for b1 in brs, b2 in brs]
    return _assemble(brs, [p[1] for p in parts]),
           _assemble(brs, [p[2] for p in parts])
end

# Δ 行序的 (branch, node) 列表
_contour_nodes(x::ImagCorrelationFunction) =
    [(:τ, i) for i in 1:size(x, 1)]
_contour_nodes(x::RealCorrelationFunction) =
    vcat([(:+, i) for i in 1:size(x, 1)], [(:-, i) for i in 1:size(x, 1)])
_contour_nodes(x::AbstractMixedCorrelationFunction) =
    vcat([(:+, i) for i in 1:rsize(x)], [(:-, i) for i in 1:rsize(x)],
         [(:τ, i) for i in 1:isize(x)])

# ------------------------------------------------------------------ #
# GTEMPO / TEMPO 默认 ordering（order=1，单 band 单 flavor）下的链序：
#   IM  : τ_N .. τ_1                        （虚时降序）
#   RT  : (+,k) (-,k) .. (+,1) (-,1)        （时间降序，+ 先 − 后）
#   MX  : τ_Nτ .. τ_1  ++  (+,1) (-,1) .. (+,kt) (-,kt)
# ------------------------------------------------------------------ #

_if_legorder(x::ImagCorrelationFunction) =
    [(:τ, i) for i in size(x, 1):-1:1]
_if_legorder(x::RealCorrelationFunction) =
    [(b, i) for i in size(x, 1):-1:1 for b in (:+, :-)]
function _if_legorder(x::AbstractMixedCorrelationFunction)
    kt, Nτ = rsize(x), isize(x)
    return vcat([(:τ, i) for i in Nτ:-1:1],
                [(b, i) for i in 1:kt for b in (:+, :-)])
end

# ------------------------------------------------------------------ #
# 指标与稠密表示的容器类型
# ------------------------------------------------------------------ #

"""
    GrassmannContourIndex(j, branch, band, conj)

费米稠密 IF 的一条 Grassmann 腿：`j` 为实时/虚时步数，`branch` 为
`:+`、`:-` 或 `:τ`，`band` 为能带（BCS bath 含 2 个 band，其余为 1），
`conj = false` 为 `a`（湮灭），`true` 为 `ā`（共轭）。
"""
struct GrassmannContourIndex
    j::Int
    branch::Symbol
    band::Int
    conj::Bool
end

"""
    BosonicContourIndex(j, branch)

玻色稠密 IF 的一条 contour 节点：`j` 为实时/虚时步数，`branch` 为
`:+`、`:-` 或 `:τ`。
"""
struct BosonicContourIndex
    j::Int
    branch::Symbol
end

"""
    DenseFermionicIF{T<:Number}

费米影响泛函的稠密表示：`basis` 为 `L = 2n` 条腿（链序，每节点相邻的
`a`、`ā` 各一条），`data` 为 `2^L` 个 Grassmann monomial 系数（腿 1
最快）。`T` 由 correlation function 的元素类型决定。
"""
struct DenseFermionicIF{T<:Number}
    basis::Vector{GrassmannContourIndex}
    data::Vector{T}
end

"""
    DenseBosonicADT{T<:Number}

玻色影响泛函在对角（ADT）耦合下的稠密表示：`basis` 为 `n` 个 contour
节点（链序），`d` 为物理维度（耦合算符本征值矢量的长度），`data` 为
`dⁿ` 个路径权重（节点 1 最快）——指标数与 `basis` 数一致。
"""
struct DenseBosonicADT{T<:Number}
    basis::Vector{BosonicContourIndex}
    d::Int
    data::Vector{T}
end

"""
    DenseBosonicPT{T<:Number}

玻色影响泛函在非对角（PT）耦合下的稠密表示：`basis` 为 `n` 个 contour
节点（链序），`d` 为系统维度，`data` 为 `d^{2n}` 个矩阵元（指标 2p−1 =
节点 p 的 ket 路径、2p = 节点 p 的 bra 路径，指标 1 最快）——指标数是
`basis` 数的两倍。
"""
struct DenseBosonicPT{T<:Number}
    basis::Vector{BosonicContourIndex}
    d::Int
    data::Vector{T}
end

# ------------------------------------------------------------------ #
# 费米：稠密 Grassmann 张量
# ------------------------------------------------------------------ #

"""
    DenseFermionicIF(x::AbstractCorrelationFunction) -> DenseFermionicIF

费米影响泛函 IF 的**稠密 Grassmann 表示**：

```text
I[ā, a] = exp(Σ_pq ā_p Δ_pq a_q)
```

每个 contour 时间节点对应相邻两条物理腿 `(a, ā)`
（`basis` 元素 `GrassmannContourIndex(j, branch, band=1, conj)`，各维度
2），`data` 为 `2^L`（`L = 2n`，`n` 为节点数）个 Grassmann monomial
系数。腿序与 GTEMPO 三种 lattice 的默认 ordering（`order=1`）一致
（见 `_if_legorder`）：虚时降序、实时时间降序（+ 先 − 后）、混合为虚时
降序接实时升序；GTEMPO 链上的左边界辅助腿不进入 IF。

元素类型 `T` 取自 correlation function（虚时通常为实）。

腿数 `L ≤ 30`（`2^30` 个系数即内存上限）；超过时抛 `ArgumentError`。

前缀割（链前 `k` 个节点块，即腿 `1..2k`）的 plain SVD 熵等于系数态的
OSEE（协方差快速算法，`GaussianStates.osee`），可用于交叉验证；节点
`p` 在 Δ 行序（`[:+; :-; :τ]`，分支内时间升序）中的指标即
GaussianStates 配对态的模式指标（site `p` = 模式对 `[p, n+p]`）。
"""
function DenseFermionicIF(x::AbstractCorrelationFunction)
    Δ = _contour_delta(x)
    n = size(Δ, 1)
    L = 2n
    L ≤ 30 || throw(ArgumentError(
        "稠密 IF 张量阶数 L = $(L) 超过上限 30（2^L 个系数将耗尽内存）；" *
        "请减小 contour 节点数 n（当前 $(n)，上限 15），或使用协方差快速算法"))
    T = eltype(Δ)
    # 链腿序下的配对系数：C[leg(ā_p), leg(a_q)] = Δ[p, q]
    nodes = _contour_nodes(x)
    chain = _if_legorder(x)
    pos = Dict(bi => k for (k, bi) in enumerate(chain))
    C = zeros(T, L, L)
    for p in 1:n, q in 1:n
        iszero(Δ[p, q]) && continue
        C[2pos[nodes[p]], 2pos[nodes[q]]-1] = Δ[p, q]
    end
    # 逐对折叠 Π (1 + C_ij ξ_i ξ_j)：data[mask+1] = 该 monomial 的系数
    data = zeros(T, 1 << L)
    data[1] = one(T)
    for i in 1:L, j in 1:L
        c = C[i, j]
        iszero(c) && continue
        a, b = minmax(i, j)
        cc = i < j ? c : -c
        bita, bitb = 1 << (a - 1), 1 << (b - 1)
        # 目标序递减：同一因子内源 data[m2] 尚未被本因子修改
        for m2 in (1 << L - 1):-1:0
            (m2 & bita) == 0 || continue
            (m2 & bitb) == 0 || continue
            m = m2 | bita | bitb
            # 费米交换符号：ξ_a ξ_b 移过 a、b 之下的奇腿
            sgn = isodd(count_ones(m2 & (bita - 1)) + count_ones(m2 & (bitb - 1))) ?
                  -1.0 : 1.0
            data[m + 1] += cc * sgn * data[m2 + 1]
        end
    end
    basis = GrassmannContourIndex[GrassmannContourIndex(j, b, 1, c)
                                  for (b, j) in chain for c in (false, true)]
    return DenseFermionicIF(basis, data)
end

# ------------------------------------------------------------------ #
# 玻色：ADT（对角耦合）与 PT（非对角耦合）
# ------------------------------------------------------------------ #

"""
    DenseBosonicADT(x::AbstractCorrelationFunction, λ::AbstractVector{<:Real})
        -> DenseBosonicADT

玻色影响泛函在**对角（ADT）耦合**下的稠密表示。耦合形式

```text
H_int = hyb ⊗ Σ_k g_k (a†_k + a_k) ,   hyb = diag(λ)
```

在 `hyb` 的本征基下每个 contour 节点的路径值取其本征值，`data` 元素为

```text
I[v₁..v_n] = exp(Σ_pq Δ_pq v_p v_q)
```

其中 `Δ = _contour_delta(x)`（节点行序 `[:+; :-; :τ]`），`n` 为节点数；
与 brute-force QuAPI 的 `exp(dot(v, G, v))` 约定一致（含对角项与上下
三角全求和）。`λ` 为本征值矢量（自旋-玻色取 `[1, -1]` 等；顺序决定
张量指标到路径值的映射）；`basis` 为链序节点（GTEMPO/TEMPO 默认
ordering，见 `_if_legorder`），`d = length(λ)`。

复杂度 O(n² dⁿ)；`dⁿ > 2³⁰` 时抛 `ArgumentError`（d=2 即 n ≤ 30）。
"""
function DenseBosonicADT(x::AbstractCorrelationFunction,
                         λ::AbstractVector{<:Real})
    isempty(λ) && throw(ArgumentError("λ 不能为空"))
    d = length(λ)
    Δ = _contour_delta(x)
    nodes = _contour_nodes(x)
    chain = _if_legorder(x)
    n = length(chain)
    n * log2(d) ≤ 30 || throw(ArgumentError(
        "稠密 IF 张量总元素数 d^n = $(d)^$(n) 超过上限 2^30；" *
        "请减小 contour 节点数 n（当前 $(n)）或路径值个数 d（当前 $(d)）"))
    T = eltype(Δ)
    nodeidx = Dict(bi => k for (k, bi) in enumerate(nodes))
    data = ones(T, ntuple(_ -> d, n))
    # 逐对吸收 exp(Δ_pq λ_a λ_b)：I = Π_pq exp(Δ_pq v_p v_q)
    for p in 1:n, q in 1:n
        c = Δ[nodeidx[chain[p]], nodeidx[chain[q]]]
        iszero(c) && continue
        F = exp.(c .* λ .* transpose(λ))
        @inbounds for I in CartesianIndices(data)
            data[I] *= F[I[p], I[q]]
        end
    end
    basis = BosonicContourIndex[BosonicContourIndex(j, b) for (b, j) in chain]
    return DenseBosonicADT(basis, d, vec(data))
end

"""
    DenseBosonicPT(x::AbstractCorrelationFunction, hyb::AbstractMatrix)
        -> DenseBosonicPT

玻色影响泛函在**非对角（PT）耦合**下的稠密表示，`hyb` 为任意方阵
（`NonDiagonalHyb` 约定：`op*a + op'*a'`，`op = hyb`，不要求厄米）。
在 `hyb` 的本征基以外，IF 不再是路径权重的对角嵌入，而是作用在系统
路径上的算符：

```text
I = exp(H)
H = Σ_pq [ M1[p,q]·op1_p·op2_q + M2[p,q]·op2_p·op1_q ]
```

`M1 + M2 = Δ` 是 `_contour_delta(x)` 的算符次序拆分（`CorrelationMatrix`
块用 `ηⱼₖ` / `ηₖⱼ` 区分等时两种次序，见 `_contour_delta_split`），
`op1 = hyb`、`op2 = hyb'`，`op_p` 表示作用在第 `p` 个节点的系统
自由度上。对对角 `hyb` 上式退化为 ADT 路径权重的对角嵌入。这里直接
对显式构造的 `H` 做矩阵指数，**不含任何 Prony/指数展开近似**。

`data` 为 `d^{2n}` 个矩阵元：指标 `2p−1` = 节点 `p` 的 ket 路径，
`2p` = 节点 `p` 的 bra 路径，指标 1 最快；`basis` 为链序节点，
`d = size(hyb, 1)`。总元素数 `d^{2n} > 2³⁰` 时抛 `ArgumentError`
（d=2 即 n ≤ 15）。
"""
function DenseBosonicPT(x::AbstractCorrelationFunction, hyb::AbstractMatrix)
    size(hyb, 1) == size(hyb, 2) || throw(ArgumentError("hyb 必须是方阵"))
    d = size(hyb, 1)
    M1, M2 = _contour_delta_split(x)
    # 张量指标按链序（与 basis 一致）排列：把 M1/M2 从 Δ 行序重排到链序
    chain = _if_legorder(x)
    nodeidx = Dict(bi => k for (k, bi) in enumerate(_contour_nodes(x)))
    ord = [nodeidx[bi] for bi in chain]
    M1 = M1[ord, ord]
    M2 = M2[ord, ord]
    n = length(chain)
    n * log2(d) ≤ 15 || throw(ArgumentError(
        "稠密 PT 张量总元素数 d^{2n} = $(d)^$(2n) 超过上限 2^30；" *
        "请减小 contour 节点数 n（当前 $(n)）或系统维度 d（当前 $(d)）"))
    T = promote_type(eltype(M1), eltype(hyb))
    op1 = Matrix{T}(hyb)
    op2 = Matrix{T}(adjoint(hyb))
    dim = d^n
    O1 = [_op_at(op1, p, n) for p in 1:n]
    O2 = [_op_at(op2, p, n) for p in 1:n]
    # H = Σ_p O1[p]·(Σ_q M1[p,q] O2[q]) + O2[p]·(Σ_q M2[p,q] O1[q])
    H = zeros(T, dim, dim)
    for p in 1:n
        W1 = zeros(T, dim, dim)
        W2 = zeros(T, dim, dim)
        for q in 1:n
            c1 = M1[p, q]
            iszero(c1) || (W1 .+= c1 .* O2[q])
            c2 = M2[p, q]
            iszero(c2) || (W2 .+= c2 .* O1[q])
        end
        H .+= O1[p] * W1
        H .+= O2[p] * W2
    end
    F = exp(H)
    # reshape 得 (ket₁..ket_n, bra₁..bra_n)，再交错成每节点 (ket, bra)
    Ten = reshape(F, ntuple(_ -> d, 2n))
    perm = vcat([[p, n + p] for p in 1:n]...)
    basis = BosonicContourIndex[BosonicContourIndex(j, b) for (b, j) in chain]
    return DenseBosonicPT(basis, d, vec(permutedims(Ten, perm)))
end

# op 作用在 d^n 维空间的第 p 个张量因子上（节点 1 为最快指标，与 data 的
# 线性序一致）：从慢到快的因子依次为节点 n, n−1, …, 1
function _op_at(op::AbstractMatrix{T}, p::Int, n::Int) where {T}
    d = size(op, 1)
    slow = d^(n - p)      # 节点 n..p+1
    fast = d^(p - 1)      # 节点 p−1..1
    return kron(Matrix{T}(I, slow, slow), op, Matrix{T}(I, fast, fast))
end
