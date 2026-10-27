# ------------------------------------------------------------------ #
# contour Δ（hybridization 矩阵）组装
#
# 节点行序：分支 [:+, :-, :τ] 逐分支、分支内时间升序；整条 contour 的 Δ
# 是分支块的块矩阵。元素类型由 correlation function 决定（虚时通常为实，
# 实/混合含虚部），这里不写死为复数。
# ------------------------------------------------------------------ #

# 按统一的 contour 节点序（contour_ordering）把分支块重排成扁平矩阵。
# blocks[(b1,b2)] 是 branch(x; b1=b1, b2=b2) 的矩阵（分支内指标 1..N_b）。
function _assemble_ordered(nodes, blocks::AbstractDict)
    T = promote_type(map(eltype, values(blocks))...)
    n = length(nodes)
    M = zeros(T, n, n)
    for (p, (b1, j1)) in enumerate(nodes), (q, (b2, j2)) in enumerate(nodes)
        M[p, q] = blocks[(b1, b2)][j1, j2]
    end
    return M
end

function _contour_delta(x::AbstractCorrelationFunction)
    brs = branches(x)
    blocks = Dict((b1, b2) => Matrix(branch(x; b1=b1, b2=b2))
                  for b1 in brs, b2 in brs)
    return _assemble_ordered(contour_ordering(x), blocks)
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
    brs = branches(x)
    parts = Dict((b1, b2) => _split_block(branch(x; b1=b1, b2=b2))
                 for b1 in brs, b2 in brs)
    nodes = contour_ordering(x)
    M1 = _assemble_ordered(nodes, Dict(k => v[1] for (k, v) in parts))
    M2 = _assemble_ordered(nodes, Dict(k => v[2] for (k, v) in parts))
    return M1, M2
end

# ------------------------------------------------------------------ #
# 稠密 IF 的腿基
# ------------------------------------------------------------------ #

# 稠密 IF 的输入：普通 correlation function 或 BCS 包装（后者不是
# AbstractCorrelationFunction 的子类型）
const IFSource = Union{AbstractCorrelationFunction,BCSCorrelationFunction}

# flavor（band）个数：BCS bath 有 2 个
_bands(::AbstractCorrelationFunction) = 1
_bands(::BCSCorrelationFunction) = 2

"""
    contour_ordering(x::Union{AbstractCorrelationFunction,BCSCorrelationFunction})
        -> Vector{Tuple{Symbol,Int}}

contour 节点的**统一序**，返回 `(branch, j)` 列表。它同时定义
Grassmann/玻色 IF 张量的轴序与 `_contour_delta` 的行列序：

```text
IM : (:τ,1) (:τ,2) .. (:τ,Nτ)                    （虚时升序）
RT : (+,1) (-,1) (+,2) (-,2) ..                  （时间升序，+ 先 − 后）
MX : (:τ,1) .. (:τ,Nτ)  ++  (+,1) (-,1) ..       （虚时升序接实时升序）
```

费米多 band 时，`basis` 在每个节点内按 band 1,2,… 展开，每个 band 内
`a`、`ā` 相邻（`ā` 紧跟 `a`）。
"""
contour_ordering(x::ImagCorrelationFunction) =
    [(:τ, i) for i in 1:size(x, 1)]
contour_ordering(x::RealCorrelationFunction) =
    [(b, i) for i in 1:size(x, 1) for b in (:+, :-)]
function contour_ordering(x::AbstractMixedCorrelationFunction)
    kt, Nτ = rsize(x), isize(x)
    return vcat([(:τ, i) for i in 1:Nτ],
                [(b, i) for i in 1:kt for b in (:+, :-)])
end
contour_ordering(x::BCSCorrelationFunction) = contour_ordering(x.uu)

"""
    basis(x::Union{AbstractCorrelationFunction,BCSCorrelationFunction})
        -> Vector{GrassmannContourIndex}

稠密 Grassmann 影响泛函的腿基：节点按 [`contour_ordering`](@ref) 排列，
每个节点内逐 band（1, 2, …）展开，每个 band 内 `a`、`ā` 相邻（`ā` 紧跟 `a`）。
BCS bath 每节点 4 条腿（两个 band）。

单 flavor 的腿数为 `2n`（n 为 contour 节点数），BCS 为 `4n`。
"""
function basis(x::IFSource)
    nb = _bands(x)
    return GrassmannContourIndex[GrassmannContourIndex(j, b, band, c)
                                 for (b, j) in contour_ordering(x)
                                 for band in 1:nb for c in (false, true)]
end

# ------------------------------------------------------------------ #
# 稠密 IF 之间的乘法（Grassmann / 玻色）
#
# 两个 IF 的 basis 可以不同：共享的腿按 Grassmann 乘法融合
# （1·ξ = ξ，ξ·ξ = 0，重排符号计入），只出现在一边的腿成为外腿；
# 玻色情形为共享节点上的逐点（Hadamard）乘积。
# ------------------------------------------------------------------ #

# 并基：bA 在前，bB 独有的元素按 bB 的顺序追加；返回并基以及 bB 各元素在
# 并基中的位置。
function _union_basis(bA, bB)
    ubasis = collect(bA)
    where = Dict(l => k for (k, l) in enumerate(ubasis))
    posB = Vector{Int}(undef, length(bB))
    for (k, l) in enumerate(bB)
        p = get(where, l, 0)
        if p == 0
            push!(ubasis, l)
            p = length(ubasis)
            where[l] = p
        end
        posB[k] = p
    end
    return ubasis, posB
end
