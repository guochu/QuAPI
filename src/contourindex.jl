# ============================================================================
# contourindex.jl —— 稠密 IF 的腿指标，以及其在 correlation function 上的访问
# ============================================================================

"""
    AbstractContourIndex

稠密影响泛函的腿指标抽象类型：定位一条 Grassmann / 路径腿在 contour 上的
位置（时间步 `j` 与分支 `branch`）。
"""
abstract type AbstractContourIndex end

"""
    GrassmannContourIndex(j, branch, band, conj)

费米稠密 IF 的一条 Grassmann 腿：`j` 为实时/虚时步数，`branch` 为
`:+`、`:-` 或 `:τ`，`band` 为能带（BCS bath 含 2 个 band，其余为 1；
BCS 下 band 奇数 = 上自旋，偶数 = 下自旋），`conj = false` 为 `a`
（湮灭），`true` 为 `ā`（共轭）。
"""
struct GrassmannContourIndex <: AbstractContourIndex
    j::Int
    branch::Symbol
    band::Int
    conj::Bool
end

"""
    BosonicContourIndex(j, branch)

玻色稠密 IF 的一条路径腿：`j` 为实时/虚时步数，`branch` 为
`:+`、`:-` 或 `:τ`。
"""
struct BosonicContourIndex <: AbstractContourIndex
    j::Int
    branch::Symbol
end

"腿指标所在的分支（`:+`、`:-` 或 `:τ`）。"
branch(x::AbstractContourIndex) = x.branch

# ------------------------------------------------------------------ #
# 用腿指标访问 hybridization function：
#     x[i::AbstractContourIndex, j::AbstractContourIndex] ==
#         index(x, i.j, j.j; b1 = branch(i), b2 = branch(j))
# ------------------------------------------------------------------ #

"""
    x[i::GrassmannContourIndex, j::GrassmannContourIndex]

读取费米 hybridization function 的矩阵元，等价于
`index(x, i.j, j.j; b1=branch(i), b2=branch(j))`。

普通（单 flavor）费米 bath 的 hybridization function 只有**带内耦合**，
因此要求两个指标的 `band` 相等（BCS bath 见下一个方法）。
"""
function Base.getindex(x::AbstractCorrelationFunction,
                       i::GrassmannContourIndex, j::GrassmannContourIndex)
    i.band == j.band || throw(ArgumentError(
        "普通费米 bath 的 hybridization function 只有带内耦合：" *
        "两个 GrassmannContourIndex 的 band 必须相等（$(i.band) vs $(j.band)）"))
    return index(x, i.j, j.j; b1=branch(i), b2=branch(j))
end

"""
    x[i::GrassmannContourIndex, j::GrassmannContourIndex]   (BCS bath)

BCS bath 的 `band` 按 flavor 解释：奇数 = 上自旋（`band = 1`），偶数 =
下自旋（`band = 2`），允许带间（反常）块。行为等价于
`index(x[σᵢ, σⱼ], i.j, j.j; b1=branch(i), b2=branch(j))`。
"""
function Base.getindex(x::BCSCorrelationFunction,
                       i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.band ≥ 1) || throw(ArgumentError("band 必须 ≥ 1（当前 $(i.band)）"))
    (j.band ≥ 1) || throw(ArgumentError("band 必须 ≥ 1（当前 $(j.band)）"))
    σᵢ = isodd(i.band) ? :↑ : :↓
    σⱼ = isodd(j.band) ? :↑ : :↓
    return index(x[σᵢ, σⱼ], i.j, j.j; b1=branch(i), b2=branch(j))
end

"""
    x[i::BosonicContourIndex, j::BosonicContourIndex]

读取玻色 hybridization function 的矩阵元，等价于
`index(x, i.j, j.j; b1=branch(i), b2=branch(j))`（玻色 bath 无 band）。
"""
function Base.getindex(x::AbstractCorrelationFunction,
                       i::BosonicContourIndex, j::BosonicContourIndex)
    return index(x, i.j, j.j; b1=branch(i), b2=branch(j))
end
