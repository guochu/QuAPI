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
# - 普通费米 bath：按 conj 组合给出 Wick 收缩（a-ā 取 Δ，ā-a 取 -Δ，
#   同 conj 为零）；
# - BCS bath：按 band 选择自旋块；
# - 玻色 bath：直接取 Δ 矩阵元。
# ------------------------------------------------------------------ #

"""
    x[i::GrassmannContourIndex, j::GrassmannContourIndex]

Bath contraction of a normal (non-BCS) fermionic bath between the two
Grassmann legs `i` and `j`, selected by the `conj` flags:

- `a-ā` (`i.conj == false`, `j.conj == true`):

      x[i, j] = Δ[i.j, j.j] ;

- `ā-a` (`i.conj == true`, `j.conj == false`):

      x[i, j] = -Δ[i.j, j.j] ;

- `a-a` and `ā-ā`: zero.

Δ is the `branch` block carried by the two indices; the legs may belong to
different bands (no band matching is required).
"""
function Base.getindex(x::AbstractCorrelationFunction,
                       i::GrassmannContourIndex, j::GrassmannContourIndex)
    D = branch(x; b1=branch(i), b2=branch(j))
    if !i.conj && j.conj          # a-ā
        return D[i.j, j.j]
    elseif i.conj && !j.conj      # ā-a
        return -D[i.j, j.j]
    else                          # a-a, ā-ā
        return zero(eltype(D))
    end
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
