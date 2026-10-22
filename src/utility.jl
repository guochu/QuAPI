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

# BCS bath：2 个 flavor（band ↑=1, ↓=2），四个块的 contour 结构相同
_bands(::AbstractCorrelationFunction) = 1
_bands(::BCSCorrelationFunction) = 2
_branches(x::BCSCorrelationFunction) = _branches(x.uu)

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
# BCS：四个块的节点结构相同，取任一块
_contour_nodes(x::BCSCorrelationFunction) = _contour_nodes(x.uu)

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
_if_legorder(x::BCSCorrelationFunction) = _if_legorder(x.uu)

# ------------------------------------------------------------------ #
# 稠密表示的容器类型（腿指标见 contourindex.jl）
# ------------------------------------------------------------------ #

"""
    GrassmannIF{T<:Number}

费米影响泛函的稠密 Grassmann 表示：`basis` 为 `L` 条腿（链序，每节点
band 逐次展开、每个 band 内 `a`、`ā` 相邻），`data` 为 `2^L` 个
Grassmann monomial 系数（腿 1 最快）。`T` 由 correlation function 的
元素类型决定。
"""
struct GrassmannIF{T<:Number}
    basis::Vector{GrassmannContourIndex}
    data::Vector{T}
end

"""
    BosonicIF{T<:Number}

玻色影响泛函在对角（ADT）耦合下的稠密表示：`basis` 为 `n` 个 contour
节点（链序），`d` 为物理维度（耦合算符本征值矢量的长度），`data` 为
`dⁿ` 个路径权重（节点 1 最快）——指标数与 `basis` 数一致。
"""
struct BosonicIF{T<:Number}
    basis::Vector{BosonicContourIndex}
    d::Int
    data::Vector{T}
end

"""
    NonDiagonalBosonicIF{T<:Number}

玻色影响泛函在非对角（PT）耦合下的稠密表示：`basis` 为 `n` 个 contour
节点（链序），`d` 为系统维度，`data` 为 `d^{2n}` 个矩阵元（指标 2p−1 =
节点 p 的 ket 路径、2p = 节点 p 的 bra 路径，指标 1 最快）——指标数是
`basis` 数的两倍。
"""
struct NonDiagonalBosonicIF{T<:Number}
    basis::Vector{BosonicContourIndex}
    d::Int
    data::Vector{T}
end

# ------------------------------------------------------------------ #
# 费米：稠密 Grassmann 张量
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
# 费米稠密 Grassmann 表示的三件套（供 GrassmannIF 与外部复用）：
#   _if_legs(x)         → 链序腿列表（每节点逐 band、a/ā 相邻）
#   _if_coefficients(x) → 双线性指数矩阵 C（行 = 腿，C[i,j] 配 ξ_i ξ_j）
#   _fold_grassmann(C)  → exp(Σ C_ij ξ_i ξ_j) 的稠密系数向量
# ------------------------------------------------------------------ #

# 稠密 IF 的输入：普通 correlation function 或 BCS 包装（后者不是
# AbstractCorrelationFunction 的子类型）
const IFSource = Union{AbstractCorrelationFunction,BCSCorrelationFunction}

# 链序腿：节点（链序）× band × (a, ā)
function _if_legs(x::IFSource)
    nb = _bands(x)
    return GrassmannContourIndex[GrassmannContourIndex(j, b, band, c)
                                 for (b, j) in _if_legorder(x)
                                 for band in 1:nb for c in (false, true)]
end

# GTEMPO 腿序的稠密系数向量 ← 双线性指数矩阵（Grassmann 幂零 ⟹
# exp(Σ C_ij ξ_i ξ_j) = Π (1 + 项)，逐项折叠；含费米交换符号）
function _fold_grassmann(C::AbstractMatrix{T}) where {T<:Number}
    L = size(C, 1)
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
    return data
end

# 从稠密系数向量读单项式系数：`legs` 为该 monomial 的（奇数）腿指标
function _grassmann_coefficient(data::AbstractVector, legs::Integer...)
    L = round(Int, log2(length(data)))
    mask = (1 << L) - 1
    for k in legs
        mask &= ~(1 << (k - 1))
    end
    return data[mask + 1]
end

# IF 指数中的 Δ 项：((band₁, conj₁), (band₂, conj₂), Δ 块)
# 普通费米高斯 bath：ā Δ a
_delta_terms(x::AbstractCorrelationFunction) =
    [((1, true), (1, false), _contour_delta(x))]
# BCS：正常 ↑↑ (ā↑Δa↑)、正常 ↓↓ (a↓Δā↓)、反常 ↑↓ (ā↑Δā↓)、反常 ↓↑ (a↓Δa↑)
_delta_terms(x::BCSCorrelationFunction) =
    [((1, true), (1, false), _contour_delta(x.uu)),
     ((2, false), (2, true), _contour_delta(x.dd)),
     ((1, true), (2, true), _contour_delta(x.ud)),
     ((2, false), (1, false), _contour_delta(x.du))]

function _if_eltype(x::IFSource)
    return promote_type(map(t -> eltype(t[3]), _delta_terms(x))...)
end

# 双线性指数矩阵 C（腿序同 _if_legs），C[i, j] 为 ξ_i ξ_j 的系数
function _if_coefficients(x::IFSource)
    legs = _if_legs(x)
    T = _if_eltype(x)
    L = length(legs)
    idx = Dict(l => k for (k, l) in enumerate(legs))
    chain = _if_legorder(x)
    nodes = _contour_nodes(x)
    nodeidx = Dict(bi => k for (k, bi) in enumerate(nodes))
    C = zeros(T, L, L)
    for (t1, t2, Δ) in _delta_terms(x)
        for p in eachindex(chain), q in eachindex(chain)
            c = Δ[nodeidx[chain[p]], nodeidx[chain[q]]]
            iszero(c) && continue
            i = idx[GrassmannContourIndex(chain[p][2], chain[p][1], t1[1], t1[2])]
            j = idx[GrassmannContourIndex(chain[q][2], chain[q][1], t2[1], t2[2])]
            C[i, j] += c
        end
    end
    return C
end

"""
    GrassmannIF(x::Union{AbstractCorrelationFunction,BCSCorrelationFunction})
        -> GrassmannIF

费米影响泛函 IF 的**稠密 Grassmann 表示**：

```text
I[ā, a] = exp(Σ_pq ā_p Δ_pq a_q)
```

每个 contour 时间节点对应相邻的物理腿（`basis` 元素
`GrassmannContourIndex(j, branch, band, conj)`，各维度 2）；单 flavor
bath 每节点 2 条腿 `(a, ā)`，**BCS bath 每节点 4 条腿**（`band = 1, 2`
两个 flavor，每个 flavor 的 `(a, ā)` 相邻，加反常配对块）。`data` 为
`2^L` 个 Grassmann monomial 系数（`L = length(basis)`，腿 1 最快）。

腿序与 GTEMPO 三种 lattice 的默认 ordering（`order=1`，2 flavor 时
`bands=2` 的 `A1Ā1B1B̄1` 型）一致（见 `_if_legorder`）：虚时降序、实时
时间降序（+ 先 − 后）、混合为虚时降序接实时升序；GTEMPO 链上的左边界
辅助腿不进入 IF。

元素类型 `T` 取自 correlation function（虚时通常为实）。

腿数 `L ≤ 30`（`2^30` 个系数即内存上限）；超过时抛 `ArgumentError`。

前缀割（链前 `k` 个节点块，即腿 `1..2k·bands`）的 plain SVD 熵等于
系数态的 OSEE（协方差快速算法，`GaussianStates.osee`），可用于交叉
验证；节点 `p` 在 Δ 行序（`[:+; :-; :τ]`，分支内时间升序）中的指标即
GaussianStates 配对态的模式指标。
"""
function GrassmannIF(x::IFSource)
    basis = _if_legs(x)
    L = length(basis)
    L ≤ 30 || throw(ArgumentError(
        "稠密 IF 张量阶数 L = $(L) 超过上限 30（2^L 个系数将耗尽内存）；" *
        "请减小 contour 节点数 n 或 flavor 数（当前 L = $(L)），" *
        "或使用协方差快速算法"))
    return GrassmannIF(basis, _fold_grassmann(_if_coefficients(x)))
end

# ------------------------------------------------------------------ #
# 玻色：ADT（对角耦合）与 PT（非对角耦合）
# ------------------------------------------------------------------ #

"""
    BosonicIF(x::AbstractCorrelationFunction, λ::AbstractVector{<:Real})
        -> BosonicIF

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
function BosonicIF(x::AbstractCorrelationFunction,
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
    return BosonicIF(basis, d, vec(data))
end

"""
    NonDiagonalBosonicIF(x::AbstractCorrelationFunction, hyb::AbstractMatrix)
        -> NonDiagonalBosonicIF

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
function NonDiagonalBosonicIF(x::AbstractCorrelationFunction, hyb::AbstractMatrix)
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
    return NonDiagonalBosonicIF(basis, d, vec(permutedims(Ten, perm)))
end

# op 作用在 d^n 维空间的第 p 个张量因子上（节点 1 为最快指标，与 data 的
# 线性序一致）：从慢到快的因子依次为节点 n, n−1, …, 1
function _op_at(op::AbstractMatrix{T}, p::Int, n::Int) where {T}
    d = size(op, 1)
    slow = d^(n - p)      # 节点 n..p+1
    fast = d^(p - 1)      # 节点 p−1..1
    return kron(Matrix{T}(I, slow, slow), op, Matrix{T}(I, fast, fast))
end
