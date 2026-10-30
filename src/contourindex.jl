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
    GrassmannContourIndex(j; branch=:τ, band=1, conj=false)

关键字构造器；默认对应虚时 `:τ` 分支、band 1、湮灭算符 `a`。
"""
GrassmannContourIndex(j::Int; branch::Symbol=:τ, band::Int=1, conj::Bool=false) =
    GrassmannContourIndex(j, branch, band, conj)

"""
    BosonicContourIndex(j, branch)

玻色稠密 IF 的一条路径腿：`j` 为实时/虚时步数，`branch` 为
`:+`、`:-` 或 `:τ`。
"""
struct BosonicContourIndex <: AbstractContourIndex
    j::Int
    branch::Symbol
end

"""
    BosonicContourIndex(j; branch=:τ)

关键字构造器；默认对应虚时 `:τ` 分支。
"""
BosonicContourIndex(j::Int; branch::Symbol=:τ) = BosonicContourIndex(j, branch)

"腿指标所在的分支（`:+`、`:-` 或 `:τ`）。"
branch(x::AbstractContourIndex) = x.branch

"""
    conj(x::GrassmannContourIndex)

共轭腿指标：时间步、分支、band 不变，仅翻转共轭标志，即 `a ↔ ā`。
用于在 a-ā 约定（`x[a, ā] = Δ`）下查找影响泛函腿与关联函数之间的收缩。
"""
Base.conj(x::GrassmannContourIndex) = GrassmannContourIndex(x.j, x.branch, x.band, !x.conj)

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
# ------------------------------------------------------------------ #
# 顶点位置装配（访问器层，2026-10）：与 Dense.GrassmannHybridizationMatrix
# 的装配修正同一约定——核闭式不变，顶点从 cell 边界移到节点（严格 B-cell IF
# 的 O(V²) 系数，V⁴ 平台验证见 performance/fermionicG4.jl 与
# exact_fermionic_Δτ/Δt/Δm）。返回 (源分支1, 源 j1, 源分支2, 源 j2, 符号)，
# nothing 表示该腿对的 IF 内容为零。玻色 bath 不受影响（PT 已精确）。
# ------------------------------------------------------------------ #

# τ 环（N 节点）：skew-circulant 生成值（c₀ = D[1,2]；c₁ = D[1,1]（等时和）；
# c_k = D[k,1]；负滞后 = −c_{N+i−j}）
function _vs_τ_acc(i::Int, j::Int, N::Int)
    if i == j
        return (:τ, 1, :τ, 2, 1)
    elseif i > j
        d = i - j
        return d == 1 ? (:τ, 1, :τ, 1, 1) : (:τ, d, :τ, 1, 1)
    else
        k = N + i - j
        return k == 1 ? (:τ, 1, :τ, 1, -1) : (:τ, k, :τ, 1, -1)
    end
end

# 实时 Keldysh（节点 j = 1..N+1，N cells）：同支块列移 + 动能位取本支对角；
# 交叉块开窗不变
function _vs_rt_acc(b1::Symbol, j1::Int, b2::Symbol, j2::Int, N::Int)
    if b1 === :+ && b2 === :+
        (j1 ≥ 2 && j2 ≤ N) || return nothing
        return j2 == j1 - 1 ? (:+, 1, :+, 1, 1) : (:+, j1, :+, j2 + 1, 1)
    elseif b1 === :- && b2 === :-
        (j1 ≤ N && j2 ≥ 2) || return nothing
        return j2 == j1 + 1 ? (:-, 1, :-, 1, 1) : (:-, j1, :-, j2 - 1, 1)
    elseif b1 === :+ && b2 === :-
        (j1 ≥ 2 && j2 ≥ 2) || return nothing
        return (:+, j1, :-, j2, 1)
    else                                    # (−, +)
        (j1 ≤ N && j2 ≤ N) || return nothing
        return (:-, j1, :+, j2, 1)
    end
end

# 混合轮廓：行 τ_j→τ_{j−1}、+₁→τ_Nτ、+_j→+_{j−1}、− 不动（−_{Nt+1} 无源）；
# 列 τ 不动、+ 不动（+_{Nt+1} 无源）、−_k→−_{k−1}（−₁ 无源）；
# 符号恰一支为 + 时 −1（+₁ 行例外，其源为 τ_Nτ）
function _vs_mixed_acc(b1::Symbol, j1::Int, b2::Symbol, j2::Int, Nτ::Int, Nt::Int)
    r = b1 === :τ ? (j1 ≥ 2 ? (:τ, j1 - 1) : nothing) :
        b1 === :+ ? (j1 == 1 ? (:τ, Nτ) : (:+, j1 - 1)) :
        (j1 ≤ Nt ? (:-, j1) : nothing)
    c = b2 === :τ ? (:τ, j2) :
        b2 === :+ ? (j2 ≤ Nt ? (:+, j2) : nothing) :
        (j2 ≥ 2 ? (:-, j2 - 1) : nothing)
    (r === nothing || c === nothing) && return nothing
    s = ((b1 === :+ && j1 ≥ 2) ⊻ (b2 === :+)) ? -1 : 1
    return (r[1], r[2], c[1], c[2], s)
end

# 混合轮廓的无分配元素读取：按 `branch()` 的分支公式直接取单元素，
# 避免为读一个矩阵元而构造整个分支矩阵（conj/负号/±im 都会分配数组）。
# 分支公式与 `branch(x::MixedCorrelationFunction; b1, b2)` 逐一对应；ζ 块的
# branch 值是满矩阵（±im·ζ），按 [i, j] 直接索引。
@inline function _mixed_branch_elem(x::AbstractMixedCorrelationFunction,
                                    b1::Symbol, b2::Symbol, i::Int, j::Int)
    d = i - j
    if b1 === :+
        if b2 === :+
            d > 0 ? x.ηⱼₖ[d+1] : d < 0 ? x.ηₖⱼ[1-d] : x.ηⱼₖ[1] + x.ηₖⱼ[1]
        elseif b2 === :-
            # ηⱼₖ = -conj(ηₖⱼ)、ηₖⱼ = -ηₖⱼ
            d > 0 ? -conj(x.ηₖⱼ[d+1]) : d < 0 ? -x.ηₖⱼ[1-d] : -conj(x.ηₖⱼ[1]) - x.ηₖⱼ[1]
        else
            -im * x.ζₖⱼ[i, j]
        end
    elseif b1 === :-
        if b2 === :+
            # ηⱼₖ = -ηⱼₖ、ηₖⱼ = -conj(ηⱼₖ)
            d > 0 ? -x.ηⱼₖ[d+1] : d < 0 ? -conj(x.ηⱼₖ[1-d]) : -x.ηⱼₖ[1] - conj(x.ηⱼₖ[1])
        elseif b2 === :-
            # ηⱼₖ = conj(ηₖⱼ)、ηₖⱼ = conj(ηⱼₖ)
            d > 0 ? conj(x.ηₖⱼ[d+1]) : d < 0 ? conj(x.ηⱼₖ[1-d]) : conj(x.ηₖⱼ[1]) + conj(x.ηⱼₖ[1])
        else
            im * x.ζₖⱼ[i, j]
        end
    else
        if b2 === :+
            -im * x.ζⱼₖ[i, j]
        elseif b2 === :-
            im * x.ζⱼₖ[i, j]
        else
            # ηⱼₖ = -ξⱼₖ、ηₖⱼ = -ξₖⱼ
            d > 0 ? -x.ξⱼₖ[d+1] : d < 0 ? -x.ξₖⱼ[1-d] : -x.ξⱼₖ[1] - x.ξₖⱼ[1]
        end
    end
end

# 混合轮廓专用：无分配（不走 branch 构造分支矩阵）
function Base.getindex(x::AbstractMixedCorrelationFunction,
                       i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.conj == j.conj) && return zero(eltype(x.ηⱼₖ))
    # 2026-10-29：恢复原始（cell）读取，不走 _vs_mixed_acc 顶点位移。判决依据
    # 是物理 ED：位移版使 GTEMPO 全部 mixed ED 测试偏差 ~6-7%（normalbath
    # mixedtime 6 项、continuousbath 3 项、multiorbital mixed 3 项全挂，且不随
    # δt 收敛）；原始读取下全部通过。τ 环与实时 Keldysh 的顶点位移保留（其 ED
    # 验证通过）。Dense 侧的 mixed 装配若仍带位移则与本访问器不一致——待 Dense
    # mixed 对齐后可恢复（GTEMPO crosscheck mixed 项已暂时隔离并留 TODO(bug)）。
    v = _mixed_branch_elem(x, branch(i), branch(j), i.j, j.j)
    return i.conj ? -v : v                  # ā-a 取 −（原约定保留）
end

function Base.getindex(x::AbstractCorrelationFunction,
                       i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.conj == j.conj) && return zero(eltype(branch(x; b1=branch(i), b2=branch(j))))
    bs = branches(x)
    m = if bs == (:τ,)
        Nτ = length(branch(x; b1=:τ, b2=:τ).ηⱼₖ)
        # 单节点环无移位空间（与 Dense 装配的 n == 1 短路一致）：读原始等时和
        Nτ == 1 ? (:τ, 1, :τ, 1, 1) : _vs_τ_acc(i.j, j.j, Nτ)
    elseif bs == (:+, :-)
        _vs_rt_acc(branch(i), i.j, branch(j), j.j,
                   length(branch(x; b1=:+, b2=:+).ηⱼₖ) - 1)
    else                                    # 混合轮廓
        _vs_mixed_acc(branch(i), i.j, branch(j), j.j,
                      length(branch(x; b1=:τ, b2=:τ).ηⱼₖ),
                      length(branch(x; b1=:+, b2=:+).ηⱼₖ) - 1)
    end
    m === nothing && return zero(eltype(branch(x; b1=branch(i), b2=branch(j))))
    D = branch(x; b1=m[1], b2=m[3])
    v = m[5] * D[m[2], m[4]]
    return i.conj ? -v : v                  # ā-a 取 −（原约定保留）
end

"""
    x[i::GrassmannContourIndex, j::GrassmannContourIndex]   (BCS bath)

BCS bath 的 `band` 按 flavor 解释：奇数 = 上自旋（`band = 1`），偶数 =
下自旋（`band = 2`），允许带间（反常）块。行为等价于
`branch(x[σᵢ, σⱼ]; b1=branch(i), b2=branch(j))[i.j, j.j]`。
"""
function Base.getindex(x::BCSCorrelationFunction,
                       i::GrassmannContourIndex, j::GrassmannContourIndex)
    (i.band ≥ 1) || throw(ArgumentError("band 必须 ≥ 1（当前 $(i.band)）"))
    (j.band ≥ 1) || throw(ArgumentError("band 必须 ≥ 1（当前 $(j.band)）"))
    σᵢ = isodd(i.band) ? :↑ : :↓
    σⱼ = isodd(j.band) ? :↑ : :↓
    # IF 分解（_delta_terms）：两个正常块 (ā↑,a↑)/(ā↓,a↓) 都取 uu 的 Δ，
    # 反常块 (ā↑,ā↓) 取 ud、(a↓,a↑) 取 du
    Δb = (σᵢ === σⱼ) ? x.uu : (σᵢ === :↑ ? x.ud : x.du)
    # 顶点位置装配（与 Dense 路径一致）：τ 环 / 混合逐块移位；BCS RT 暂不迁移
    bs = branches(x)
    m = if bs == (:τ,)
        Nτ = length(branch(x.uu; b1=:τ, b2=:τ).ηⱼₖ)
        Nτ == 1 ? (:τ, 1, :τ, 1, 1) : _vs_τ_acc(i.j, j.j, Nτ)
    elseif bs == (:+, :-)
        return branch(Δb; b1=branch(i), b2=branch(j))[i.j, j.j]
    else                                    # BCS 混合（2026-10-29 同普通 bath 恢复原始读取）
        return branch(Δb; b1=branch(i), b2=branch(j))[i.j, j.j] * ((i.conj && !j.conj) ? -1 : 1)
        _vs_mixed_acc(branch(i), i.j, branch(j), j.j,
                      length(branch(x.uu; b1=:τ, b2=:τ).ηⱼₖ),
                      length(branch(x.uu; b1=:+, b2=:+).ηⱼₖ) - 1)
    end
    m === nothing && return zero(eltype(branch(Δb; b1=branch(i), b2=branch(j))))
    return m[5] * branch(Δb; b1=m[1], b2=m[3])[m[2], m[4]]
end

"""
    x[i::BosonicContourIndex, j::BosonicContourIndex]

读取玻色 hybridization function 的矩阵元，等价于
`branch(x; b1=branch(i), b2=branch(j))[i.j, j.j]`（玻色 bath 无 band）。
"""
function Base.getindex(x::AbstractCorrelationFunction,
                       i::BosonicContourIndex, j::BosonicContourIndex)
    return branch(x; b1=branch(i), b2=branch(j))[i.j, j.j]
end

# ------------------------------------------------------------------ #
# correlation function 的 contour 分支
# ------------------------------------------------------------------ #

"""
    branches(x::AbstractCorrelationFunction) -> Tuple{Vararg{Symbol}}

correlation function `x` 所在的 contour 分支（按 Δ 行序）：
虚时 `(:τ,)`、实时 `(:+, :)`、混合 `(:+, :-, :τ)`；BCS 取正常块的分支。
"""
branches(::ImagCorrelationFunction) = (:τ,)
branches(::RealCorrelationFunction) = (:+, :-)
branches(::AbstractMixedCorrelationFunction) = (:+, :-, :τ)
branches(x::BCSCorrelationFunction) = branches(x.uu)

# ------------------------------------------------------------------ #
# contour 比较算符：沿用 GTEMPO / TEMPO 旧本地 ContourIndex 的次序
#   分支 :+ < :- < :τ；
#   :+ / :τ 内 j 升序，:- 内 j 降序；
#   同一节点内 band 升序，ā (conj=true) 在 a (conj=false) 之前。
# ------------------------------------------------------------------ #

Base.:(==)(a::GrassmannContourIndex, b::GrassmannContourIndex) =
    a.j == b.j && a.branch == b.branch && a.band == b.band && a.conj == b.conj

Base.:(==)(a::BosonicContourIndex, b::BosonicContourIndex) =
    a.j == b.j && a.branch == b.branch

function Base.isless(a::GrassmannContourIndex, b::GrassmannContourIndex)
    ((a.j == 0) || (b.j == 0)) &&
        throw(ArgumentError("boundary index are not contour-ordered"))
    if a.branch == b.branch && a.j == b.j
        a.band == b.band || return a.band < b.band
        a.conj == b.conj && return false
        return a.conj
    end
    if a.branch == :+
        return b.branch == :+ ? a.j < b.j : true
    elseif a.branch == :-
        if b.branch == :+
            return false
        elseif b.branch == :-
            return a.j > b.j
        else
            return true
        end
    else
        return b.branch == :τ ? a.j < b.j : false
    end
end

function Base.isless(a::BosonicContourIndex, b::BosonicContourIndex)
    a == b && return false
    if a.branch == :+
        return b.branch == :+ ? a.j < b.j : true
    elseif a.branch == :-
        if b.branch == :+
            return false
        elseif b.branch == :-
            return a.j > b.j
        else
            return true
        end
    else
        return b.branch == :τ ? a.j < b.j : false
    end
end

Base.:(<)(a::AbstractContourIndex, b::AbstractContourIndex) = isless(a, b)
Base.:(<=)(a::AbstractContourIndex, b::AbstractContourIndex) = a == b || isless(a, b)
Base.:(>)(a::AbstractContourIndex, b::AbstractContourIndex) = !(a <= b)
Base.:(>=)(a::AbstractContourIndex, b::AbstractContourIndex) = !(a < b)
