# ============================================================================
# utility.jl —— 影响泛函的稠密张量构造与 contour 组装工具
# ============================================================================

# ------------------------------------------------------------------ #
# 整条 contour 的 Δ（hybridization 配对矩阵）组装：
#   ImagCorrelationFunction          [:τ]         （N 个节点）
#   RealCorrelationFunction          [:+, :-]     （各 rsize 个节点）
#   AbstractMixedCorrelationFunction [:+, :-, :τ]（2·rsize + isize 个节点）
# 节点行序：[:+; :-; :τ] 逐分支、分支内时间升序。
# 供 GaussianStates extension（bdg_matrix 重载）与 dense_fermionicIF 复用。
# ------------------------------------------------------------------ #

_contour_delta(x::ImagCorrelationFunction) =
    Matrix{ComplexF64}(branch(x; b1=:τ, b2=:τ))

function _contour_delta(x::RealCorrelationFunction)
    N = size(x, 1)
    Δ = Matrix{ComplexF64}(undef, 2N, 2N)
    for (i, b1) in enumerate((:+, :-)), (j, b2) in enumerate((:+, :-))
        Δ[(i-1)*N+1:i*N, (j-1)*N+1:j*N] .= branch(x; b1=b1, b2=b2)
    end
    return Δ
end

function _contour_delta(x::AbstractMixedCorrelationFunction)
    Nr, Ni = rsize(x), isize(x)
    branches = (:+, :-, :τ)
    dims = (Nr, Nr, Ni)
    offs = (0, Nr, 2Nr)
    Δ = Matrix{ComplexF64}(undef, 2Nr + Ni, 2Nr + Ni)
    for (i, b1) in enumerate(branches), (j, b2) in enumerate(branches)
        Δ[offs[i]+1:offs[i]+dims[i], offs[j]+1:offs[j]+dims[j]] .=
            branch(x; b1=b1, b2=b2)
    end
    return Δ
end

# Δ 行序的 (branch, node) 列表（与 _contour_delta 的组装一致）
_contour_nodes(x::ImagCorrelationFunction) =
    [(:τ, i) for i in 1:size(x, 1)]
_contour_nodes(x::RealCorrelationFunction) =
    vcat([(:+, i) for i in 1:size(x, 1)], [(:-, i) for i in 1:size(x, 1)])
_contour_nodes(x::AbstractMixedCorrelationFunction) =
    vcat([(:+, i) for i in 1:rsize(x)], [(:-, i) for i in 1:rsize(x)],
         [(:τ, i) for i in 1:isize(x)])

# ------------------------------------------------------------------ #
# GTEMPO 默认 ordering（order=1，单 band 单 flavor）下的链序：
# 每个时间节点 → 相邻的两条 Grassmann 腿 (a, ā)，a 在前。三种 contour
# 的节点链序（相对顺序与 GTEMPO GrassmannLattice 一致，去掉其左边界
# 辅助腿与不进入 IF 的连接节点）：
#   IM  : τ_N .. τ_1                        （虚时降序）
#   RT  : (+,k) (-,k) .. (+,1) (-,1)        （时间降序，+ 先 − 后）
#   MX  : τ_Nτ .. τ_1  ++  (+,1) (-,1) .. (+,kt) (-,kt)
#         （虚时降序接实时升序；与 GTEMPO/osee 联合回归的 bond 对应一致）
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

"""
    dense_fermionicIF(x::AbstractCorrelationFunction) -> Array{ComplexF64}

费米影响泛函 IF 的**稠密 Grassmann 张量**：

```text
I[ā, a] = exp(Σ_pq ā_p Δ_pq a_q)
```

每个 contour 时间节点对应相邻两条物理腿 `(a, ā)`（各维度 2），返回
`2×2×…×2`（`L = 2n` 阶，`n` 为节点数）张量，元素为 Grassmann monomial
系数。腿序与 GTEMPO 三种 lattice 的默认 ordering（`order=1`）一致
（见 [`_if_legorder`](@ref)）：虚时降序、实时时间降序（+ 先 − 后）、
混合为虚时降序接实时升序；GTEMPO 链上的左边界辅助腿不进入 IF，不在
本张量中。

张量阶数 `L ≤ 30`（`2^30` 个复系数即内存上限）；超过时抛 `ArgumentError`。

前缀割（链前 `k` 个节点块，即腿 `1..2k`）的 plain SVD 熵等于系数态的
OSEE（协方差快速算法，`GaussianStates.osee`），可用于交叉验证；节点
`p` 在 `_contour_delta` 行序（`[:+; :-; :τ]`，分支内时间升序）中的指标
即 GaussianStates 配对态的模式指标（site `p` = 模式对 `[p, n+p]`）。
"""
function dense_fermionicIF(x::AbstractCorrelationFunction)
    Δ = _contour_delta(x)
    n = size(Δ, 1)
    L = 2n
    L ≤ 30 || throw(ArgumentError(
        "稠密 IF 张量阶数 L = $(L) 超过上限 30（2^L 个复系数将耗尽内存）；" *
        "请减小 contour 节点数 n（当前 $(n)，上限 15），或使用协方差快速算法"))
    # 链腿序下的配对系数：C[leg(ā_p), leg(a_q)] = Δ[p, q]
    nodes = _contour_nodes(x)
    chain = _if_legorder(x)
    pos = Dict(bi => k for (k, bi) in enumerate(chain))
    C = zeros(ComplexF64, L, L)
    for p in 1:n, q in 1:n
        iszero(Δ[p, q]) && continue
        C[2pos[nodes[p]], 2pos[nodes[q]]-1] = Δ[p, q]
    end
    # 逐对折叠 Π (1 + C_ij ξ_i ξ_j)：T[mask+1] = 该 monomial 的系数
    T = zeros(ComplexF64, 1 << L)
    T[1] = 1.0
    for i in 1:L, j in 1:L
        c = C[i, j]
        iszero(c) && continue
        a, b = minmax(i, j)
        cc = i < j ? c : -c
        bita, bitb = 1 << (a - 1), 1 << (b - 1)
        # 目标序递减：同一因子内源 T[m2] 尚未被本因子修改
        for m2 in (1 << L - 1):-1:0
            (m2 & bita) == 0 || continue
            (m2 & bitb) == 0 || continue
            m = m2 | bita | bitb
            # 费米交换符号：ξ_a ξ_b 移过 a、b 之下的奇腿
            sgn = isodd(count_ones(m2 & (bita - 1)) + count_ones(m2 & (bitb - 1))) ?
                  -1.0 : 1.0
            T[m + 1] += cc * sgn * T[m2 + 1]
        end
    end
    return reshape(T, ntuple(_ -> 2, L))
end

"""
    dense_bosonicIF(x::AbstractCorrelationFunction, λ::AbstractVector{<:Real})
        -> Array{ComplexF64}

玻色影响泛函 IF 的**稠密路径张量**。耦合形式

```text
H_int = hyb ⊗ Σ_k g_k (a†_k + a_k)
```

（`hyb` 作用于杂质自由度，自旋-玻色模型取 `hyb = σ_z`）。对 bath 精确
积分后 IF 只依赖杂质路径——在 `hyb` 的本征基下每个 contour 节点的路径值
取其本征值。`λ` 即该**本征值矢量**（厄米矩阵的本征值为实，自旋-玻色取
`[1, -1]` 等；矢量顺序决定张量指标到路径值的映射），张量元素为

```text
I[v₁..v_n] = exp(Σ_pq Δ_pq v_p v_q)
```

其中 `Δ = _contour_delta(x)` 为整条 contour 的 hybridization 矩阵（节点
行序 `[:+; :-; :τ]`），`n` 为节点数；这与 brute-force QuAPI 的
`exp(dot(v, G, v))` 约定一致（含对角项与上下三角全求和）。返回
`d×d×…×d`（`n` 阶，`d = length(λ)`）张量，维度顺序与
`dense_fermionicIF` 相同（GTEMPO 默认 ordering，见 `_if_legorder`）。

复杂度 O(n² dⁿ)；`dⁿ > 2³⁰` 时抛 `ArgumentError`（自旋-玻色 d=2 即
n ≤ 30）。
"""
function dense_bosonicIF(x::AbstractCorrelationFunction, λ::AbstractVector{<:Real})
    isempty(λ) && throw(ArgumentError("λ 不能为空"))
    d = length(λ)
    Δ = _contour_delta(x)
    nodes = _contour_nodes(x)
    chain = _if_legorder(x)
    n = length(chain)
    # 内存上限：总元素数 d^n ≤ 2^30（用 n·log2(d) 防整数溢出）
    n * log2(d) ≤ 30 || throw(ArgumentError(
        "稠密 IF 张量总元素数 d^n = $(d)^$(n) 超过上限 2^30；" *
        "请减小 contour 节点数 n（当前 $(n)）或路径值个数 d（当前 $(d)）"))
    nodeidx = Dict(bi => k for (k, bi) in enumerate(nodes))
    T = ones(ComplexF64, ntuple(_ -> d, n))
    # 逐对吸收 exp(Δ_pq λ_a λ_b)：I = Π_pq exp(Δ_pq v_p v_q)
    for p in 1:n, q in 1:n
        c = Δ[nodeidx[chain[p]], nodeidx[chain[q]]]
        iszero(c) && continue
        F = exp.(c .* λ .* transpose(λ))
        @inbounds for I in CartesianIndices(T)
            T[I] *= F[I[p], I[q]]
        end
    end
    return T
end
