# ------------------------------------------------------------------ #
# 稠密表示的容器类型（腿指标见 contourindex.jl）
# ------------------------------------------------------------------ #

"""
    GrassmannADT{T<:Number}

费米影响泛函的稠密 Grassmann 表示：`basis` 为 `L` 条腿（链序，每节点
band 逐次展开、每个 band 内 `a`、`ā` 相邻），`data` 为 `2^L` 个
Grassmann monomial 系数（腿 1 最快）。`T` 由 correlation function 的
元素类型决定。用 [`GrassmannIF`](@ref) 从 correlation function 构造。
"""
struct GrassmannADT{T<:Number}
    basis::Vector{GrassmannContourIndex}
    data::Vector{T}
end

# ------------------------------------------------------------------ #
# 费米：稠密 Grassmann 张量
# ------------------------------------------------------------------ #

# ------------------------------------------------------------------ #
# 稠密 IF 的两步构造：先算 hybridization matrix，再构造 IF
# ------------------------------------------------------------------ #

"""
    GrassmannHybridizationMatrix{T<:Number}

费米稠密 IF 的杂交矩阵：`basis` 为链序腿列表（见 [`basis`](@ref)），
`Δ` 为这些腿上的双线性系数矩阵（`Δ[i,j]` 是 ξ_i ξ_j 的系数；普通
bath 给 `ā Δ a`，BCS 还含反常 `āā`、`aa` 块）。

用 `GrassmannHybridizationMatrix(corr)` 从 correlation function 计算，
之后交给 [`GrassmannIF`](@ref) 折叠成稠密 IF。
"""
struct GrassmannHybridizationMatrix{T<:Number}
    basis::Vector{GrassmannContourIndex}
    Δ::Matrix{T}
end

# --------------------- 顶点位置装配（费米 PT 修正） --------------------- #
# fermionic_Δτ/Δt 的核闭式（胞元积分）逐位等于严格 B-cell IF 的 O(V²) 系数，
# 但旧装配把顶点放在 cell 边界（滞后偏一格），对严格方案的残差因此停在 O(V²)
# （V 减半残差比值 4，performance/fermionicG4.jl 实测）。此处把顶点移到节点
# ——与玻色 non-diagonal PT 在 BosonicIF 里做的分支规则修复属同一层级的
# 装配修正，核数值完全复用：
#
#   虚时 τ 环（生成值重组为 skew-circulant，N 为节点数）：
#     c₀ = ηₖⱼ[2]，c₁ = ηⱼₖ[1]+ηₖⱼ[1]（等时和），c_k = ηⱼₖ[k]（k ≥ 2）
#     C[āᵢ, aⱼ] = i > j ? c_{i−j} : i < j ? −c_{N+i−j} : c₀
#   实时 Keldysh（节点 0..N，腿位 pos(j,br,ā) = 4j + (br==:+ ? 0:2) + (ā?2:1)）：
#     (ā⁺,a⁻)/(ā⁻,a⁺) 交叉块：值不变，仅开窗（ā⁺₀ 行、a⁻₀ 列 / a⁺_N 列、
#       ā⁻_N 行置零——顶点窗口 = cell 1..N）
#     (ā⁺,a⁺) 块：列移一位（aⱼ→aⱼ₊₁）；动能位 (j, j−1) 取 (ā⁻,a⁻) 块对角
#     (ā⁻,a⁻) 块：镜像（列移 aⱼ→aⱼ₋₁；动能位 (j−1, j) 取 (ā⁺,a⁺) 块对角）
#
# 精度：核内容对严格方案（exact_fermionic_Δτ/Δt 参照）的残差按 V⁴ 标度
# （比值≈16）——与玻色 non-diagonal PT 同级（O(V²) 系数精确 + V⁴ 截断）。
# 整体 Z/ρ 另带一个 O(V²) 归一化常数（"跳去又跳回"泄漏 = exact_fermionic_*
# 的 prefactor；迹归一化的观测量中自动消去）。BCS τ 轮廓按块同构移位（严格
# BCS 链验证待补）；BCS 实时与混合轮廓暂沿用旧装配。

function _shift_τblock!(D::AbstractMatrix, C::AbstractMatrix, n::Int, rowpos, colpos)
    # 单个 (ā,a) 型块（或反常 (ā,ā)/(a,a) 块）的顶点移位：rowpos(i)/colpos(j)
    # 给出块内 (āᵢ, aⱼ) 的腿位；生成值从块的 Toeplitz 读取（滞后偏一格重组）
    n == 1 && return
    c = Vector{eltype(C)}(undef, n - 1)
    c[1] = C[rowpos(1), colpos(1)]              # c₁ = 块对角（等时和）
    for k in 2:(n - 1)
        c[k] = C[rowpos(k), colpos(1)]          # c_k = 滞后 k−1（节点 k → 1）
    end
    c₀ = C[rowpos(1), colpos(2)]                # c₀ = 滞后 −1（节点 1 → 2）
    for i in 1:n, j in 1:n
        D[rowpos(i), colpos(j)] = i == j ? c₀ : i > j ? c[i-j] : -c[n+i-j]
    end
    return
end

function _vertex_shifted_τ(C::AbstractMatrix, n::Int)
    D = copy(C)
    _shift_τblock!(D, C, n, i -> 2i, j -> 2j - 1)
    return Matrix(D)
end

function _vertex_shifted_bcs_τ(C::AbstractMatrix, n::Int)
    # BCS τ（每节点 a↑,ā↑,a↓,ā↓ 四腿）：正常块 (ā↑,a↑)/(ā↓,a↓) 与反常块
    # (ā↑,ā↓)/(a↓,a↑) 逐块做同构移位（顶点统一在节点；BCS 严格链验证待补）
    D = copy(C)
    _shift_τblock!(D, C, n, i -> 4i - 2, j -> 4j - 3)
    _shift_τblock!(D, C, n, i -> 4i, j -> 4j - 1)
    _shift_τblock!(D, C, n, i -> 4i - 2, j -> 4j)
    _shift_τblock!(D, C, n, i -> 4i - 1, j -> 4j - 3)
    return Matrix(D)
end

function _vertex_shifted_rt(C::AbstractMatrix{T}, N::Int) where {T}
    pos(j0, br, ā) = 4 * j0 + (br === :+ ? 0 : 2) + (ā ? 2 : 1)
    diag₊₊ = C[pos(0, :+, true), pos(0, :+, false)]   # (+,+) 块对角（Toeplitz 常值）
    diag₋₋ = C[pos(0, :-, true), pos(0, :-, false)]   # (−,−) 块对角
    D = zeros(T, 4 * (N + 1), 4 * (N + 1))
    for j1 in 0:N, j2 in 0:N
        j1 ≥ 1 && j2 ≥ 1 &&                       # (ā⁺,a⁻)：值不变，开窗（行 ā⁺₀ / 列 a⁻₀）
            (D[pos(j1, :+, true), pos(j2, :-, false)] = C[pos(j1, :+, true), pos(j2, :-, false)])
        j2 ≤ N - 1 && j1 ≤ N - 1 &&               # (ā⁻,a⁺)：值不变，开窗（列 a⁺_N / 行 ā⁻_N）
            (D[pos(j1, :-, true), pos(j2, :+, false)] = C[pos(j1, :-, true), pos(j2, :+, false)])
        if j1 ≥ 1 && j2 ≤ N - 1                 # (ā⁺,a⁺)：列移一位，动能位取本支对角
            D[pos(j1, :+, true), pos(j2, :+, false)] =
                j2 == j1 - 1 ? diag₊₊ : C[pos(j1, :+, true), pos(j2 + 1, :+, false)]
        end
        if j1 ≤ N - 1 && j2 ≥ 1                 # (ā⁻,a⁻)：镜像（动能位取本支对角）
            D[pos(j1, :-, true), pos(j2, :-, false)] =
                j2 == j1 + 1 ? diag₋₋ : C[pos(j1, :-, true), pos(j2 - 1, :-, false)]
        end
    end
    return D
end

function _vertex_shifted_mixed(C::AbstractMatrix, Nτ::Int, Nt::Int)
    # 混合轮廓顶点移位（数值标定 vs exact_fermionic_Δm，Richardson O(V²)）：
    #   行源 r：τ_j→τ_{j−1}（j≥2；τ₁ 无源）；+₁→τ_Nτ；+_j→+_{j−1}（j≥2）；
    #          −_j→−_j（j ≤ Nt；−_{Nt+1} 无源）
    #   列源 c：τ_k→τ_k；+_k→+_k（k ≤ Nt；+_{Nt+1} 无源）；−_k→−_{k−1}（k≥2；−₁ 无源）
    #   符号 s：行源、列源恰一个在 + 支时取 −1（+₁ 行的源是 τ_Nτ，行侧不算 + 支）
    ipos(b, j, ā) = 2 * _mxnode(b, j, Nτ) - (ā ? 0 : 1)
    D = zero(C)
    for b1 in (:τ, :+, :-), j1 in _mjran(b1, Nτ, Nt),
        b2 in (:τ, :+, :-), j2 in _mjran(b2, Nτ, Nt)
        r = b1 === :τ ? (j1 ≥ 2 ? (:τ, j1 - 1) : nothing) :
            b1 === :+ ? (j1 == 1 ? (:τ, Nτ) : (:+, j1 - 1)) :
            (j1 ≤ Nt ? (:-, j1) : nothing)
        c = b2 === :τ ? (:τ, j2) :
            b2 === :+ ? (j2 ≤ Nt ? (:+, j2) : nothing) :
            (j2 ≥ 2 ? (:-, j2 - 1) : nothing)
        (r === nothing || c === nothing) && continue
        # 符号：恰一支为 + 时取 −1；+₁ 行例外（其源为 τ_Nτ 行，跨支读不翻符号）
        s = ((b1 === :+ && j1 ≥ 2) ⊻ (b2 === :+)) ? -1 : 1
        D[ipos(b1, j1, true), ipos(b2, j2, false)] =
            s * C[ipos(r[1], r[2], true), ipos(c[1], c[2], false)]
    end
    return Matrix(D)
end

# 混合轮廓链位：τ_j（j=1..Nτ）在前，随后按时间 jj=1..Nt+1 的 (+,−) 对
function _mxnode(b::Symbol, j::Int, Nτ::Int)
    b === :τ && return j
    return Nτ + 2 * (j - 1) + (b === :+ ? 1 : 2)
end
_mjran(b::Symbol, Nτ::Int, Nt::Int) = b === :τ ? (1:Nτ) : (1:Nt+1)

# 行/列源与符号（单 flavor 混合轮廓；BCS 混合逐 band 同构复用）
function _mx_maps(b1::Symbol, j1::Int, b2::Symbol, j2::Int, Nτ::Int, Nt::Int)
    r = b1 === :τ ? (j1 ≥ 2 ? (:τ, j1 - 1) : nothing) :
        b1 === :+ ? (j1 == 1 ? (:τ, Nτ) : (:+, j1 - 1)) :
        (j1 ≤ Nt ? (:-, j1) : nothing)
    c = b2 === :τ ? (:τ, j2) :
        b2 === :+ ? (j2 ≤ Nt ? (:+, j2) : nothing) :
        (j2 ≥ 2 ? (:-, j2 - 1) : nothing)
    s = ((b1 === :+ && j1 ≥ 2) ⊻ (b2 === :+)) ? -1 : 1
    return r, c, s
end

function _vertex_shifted_bcs_mixed(C::AbstractMatrix, Nτ::Int, Nt::Int)
    # BCS 混合（每节点 a↑,ā↑,a↓,ā↓ 四腿）：正常块 (ā↑,a↑)/(ā↓,a↓) 与反常块
    # (ā↑,ā↓)/(a↓,a↑) 按单 flavor 规则逐块同构移位（严格 BCS 链验证待补）
    ipos(b, j, band, ā) = 4 * (_mxnode(b, j, Nτ) - 1) + (band - 1) * 2 + (ā ? 2 : 1)
    D = zero(C)
    blocks = (((1, true, 1, false), (1, true, 1, false)),   # (ā↑,a↑)=uu
              ((2, true, 2, false), (2, true, 2, false)),   # (ā↓,a↓)=dd
              ((1, true, 2, true), (1, true, 2, true)),     # (ā↑,ā↓)=ud
              ((2, false, 1, false), (2, false, 1, false))) # (a↓,a↑)=du
    for b1 in (:τ, :+, :-), j1 in _mjran(b1, Nτ, Nt),
        b2 in (:τ, :+, :-), j2 in _mjran(b2, Nτ, Nt)
        r, c, s = _mx_maps(b1, j1, b2, j2, Nτ, Nt)
        (r === nothing || c === nothing) && continue
        for (dst, src) in blocks
            D[ipos(b1, j1, dst[1], dst[2]), ipos(b2, j2, dst[3], dst[4])] =
                s * C[ipos(r[1], r[2], src[1], src[2]), ipos(c[1], c[2], src[3], src[4])]
        end
    end
    return Matrix(D)
end

function _vertex_shifted_IF(C::AbstractMatrix, bs::Vector{GrassmannContourIndex})
    brs = branch.(bs)
    single = all(b -> b.band == 1, bs)
    mixed = any(==(:τ), brs) && any(b -> b === :+ || b === :-, brs)
    if all(==(:τ), brs) && single
        return _vertex_shifted_τ(C, length(bs) ÷ 2)
    elseif all(==(:τ), brs)                     # BCS τ：逐块移位（严格链验证待补）
        return _vertex_shifted_bcs_τ(C, length(bs) ÷ 4)
    elseif !mixed && single
        return _vertex_shifted_rt(C, length(bs) ÷ 4 - 1)
    elseif mixed && single
        Nτ = count(==(:τ), brs) ÷ 2
        Nt = (length(brs) - count(==(:τ), brs)) ÷ 4 - 1
        return _vertex_shifted_mixed(C, Nτ, Nt)
    elseif mixed                                 # BCS 混合：逐块同构移位（严格链验证待补）
        Nτ = count(b -> branch(b) === :τ && b.band == 1, bs) ÷ 2
        Nt = (length(bs) - count(==(:τ), brs)) ÷ 8 - 1
        return _vertex_shifted_bcs_mixed(C, Nτ, Nt)
    end
    return C                                    # BCS RT 暂沿用旧装配
end

GrassmannHybridizationMatrix(x::IFSource) =
    GrassmannHybridizationMatrix(basis(x), _vertex_shifted_IF(_if_coefficients(x), basis(x)))

# 稠密系数向量 ← 双线性指数矩阵（Grassmann 幂零 ⟹
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

"""
    x[legs::GrassmannContourIndex...]

稠密 Grassmann IF 的单项式系数。`x.data` 是
`ψ = Σ_m x.data[m+1] ξ^m`（`ξ^m = ∏_{k∈m} ξ_k`，腿序同 `basis`）的系数
向量；本方法按 contour 腿在 `basis` 中定位后，返回单项式
`ξ_leg₁⋯ξ_legᵣ` 的系数：

    x[ξ₁]        → ξ₁ 的系数
    x[ξ₁, ξ₂]    → ξ₁ξ₂ 的系数
    x[]          → 常数项（1 的系数）

参数顺序与 `basis` 顺序不一致时按 Grassmann 反对易计入重排符号，例如
`x[ξ₂, ξ₁] == -x[ξ₁, ξ₂]`；腿重复（`ξ² = 0`）时返回零。腿不在 `basis`
中时抛 `ArgumentError`。
"""
function Base.getindex(x::GrassmannADT, legs::GrassmannContourIndex...)
    pos = Vector{Int}(undef, length(legs))
    for (q, leg) in enumerate(legs)
        k = findfirst(==(leg), x.basis)
        k === nothing && throw(ArgumentError("leg $leg 不在 IF 的 basis 中"))
        pos[q] = k
    end
    allunique(pos) || return zero(eltype(x.data))
    mask = 0
    for k in pos
        mask |= 1 << (k - 1)
    end
    inv = 0
    for a in 1:length(pos), b in (a + 1):length(pos)
        inv += pos[a] > pos[b]
    end
    c = x.data[mask + 1]
    return isodd(inv) ? -c : c
end

# 把 contour 腿映射为 basis 中的位置；缺失/重复时报错
function _leg_positions(basis, legs)
    pos = Vector{Int}(undef, length(legs))
    for (q, leg) in enumerate(legs)
        k = findfirst(==(leg), basis)
        k === nothing && throw(ArgumentError("leg $leg 不在 basis 中"))
        pos[q] = k
    end
    allunique(pos) || throw(ArgumentError("legs 不能重复"))
    return pos
end

"""
    complementindex(basis::AbstractVector{<:AbstractContourIndex},
                    legs::AbstractContourIndex...)

返回 `basis` 中不在 `legs` 里的腿（保持 `basis` 顺序），费米和玻色腿都适用。
与 `getindex` 配合可读取带插入的互补单项式系数：

    x[complementindex(x.basis, ξ₁, ξ₃)...]

`legs` 重复或不在 `basis` 中时报错。
"""
function complementindex(basis::AbstractVector{<:AbstractContourIndex},
                         legs::AbstractContourIndex...)
    _leg_positions(basis, legs)
    return filter(l -> !(l in legs), basis)
end

"""
    margin(x::GrassmannADT, legs::GrassmannContourIndex...) -> GrassmannADT

对未列出的腿的占据状态求和（Grassmann 版边缘化），返回只含 `legs` 的
[`GrassmannADT`](@ref)；`legs` 的顺序即结果 `basis` 的顺序，`margin(x)`
返回只含常数项的 IF。重复/不在 `basis` 中的腿报错。
"""
function margin(x::GrassmannADT, legs::GrassmannContourIndex...)
    pos = _leg_positions(x.basis, legs)
    L, r = length(x.basis), length(pos)
    out = zeros(eltype(x.data), 1 << r)
    for m in 0:(1 << L) - 1
        mr = 0
        for (q, k) in enumerate(pos)
            mr |= ((m >> (k - 1)) & 1) << (q - 1)
        end
        out[mr + 1] += x.data[m + 1]
    end
    return GrassmannADT(GrassmannContourIndex[legs...], out)
end

# IF 指数中的 Δ 项：((band₁, conj₁), (band₂, conj₂), Δ 块)
# 普通费米高斯 bath：ā Δ a
_delta_terms(x::AbstractCorrelationFunction) =
    [((1, true), (1, false), _contour_delta(x))]
# BCS：正常 ↑↑ (ā↑Δa↑)、正常 ↓↓ (ā↓Δa↓，自旋对称下与 ↑↑ 同块)、
# 反常 ↑↓ (ā↑Δā↓)、反常 ↓↑ (a↓Δa↑)
_delta_terms(x::BCSCorrelationFunction) =
    [((1, true), (1, false), _contour_delta(x.uu)),
     ((2, true), (2, false), _contour_delta(x.uu)),
     ((1, true), (2, true), _contour_delta(x.ud)),
     ((2, false), (1, false), _contour_delta(x.du))]

function _if_eltype(x::IFSource)
    return promote_type(map(t -> eltype(t[3]), _delta_terms(x))...)
end

# 双线性指数矩阵 C（腿序即 basis(x)），C[i, j] 为 ξ_i ξ_j 的系数。
# Δ 的行列序与 contour_ordering 一致，故节点 p ↔ Δ 的第 p 行/列。
function _if_coefficients(x::IFSource)
    legs = basis(x)
    T = _if_eltype(x)
    L = length(legs)
    idx = Dict(l => k for (k, l) in enumerate(legs))
    chain = contour_ordering(x)
    C = zeros(T, L, L)
    for (t1, t2, Δ) in _delta_terms(x)
        for (p, (b1, j1)) in enumerate(chain), (q, (b2, j2)) in enumerate(chain)
            c = Δ[p, q]
            iszero(c) && continue
            i = idx[GrassmannContourIndex(j1, b1, t1[1], t1[2])]
            j = idx[GrassmannContourIndex(j2, b2, t2[1], t2[2])]
            C[i, j] += c
        end
    end
    return C
end

"""
    GrassmannIF(x::Union{AbstractCorrelationFunction,BCSCorrelationFunction})
        -> GrassmannADT

费米影响泛函 IF 的**稠密 Grassmann 表示**：

```text
I[ā, a] = exp(Σ_pq ā_p Δ_pq a_q)
```

每个 contour 时间节点对应相邻的物理腿（`basis` 元素
`GrassmannContourIndex(j, branch, band, conj)`，各维度 2）；单 flavor
bath 每节点 2 条腿 `(a, ā)`，**BCS bath 每节点 4 条腿**（`band = 1, 2`
两个 flavor，每个 flavor 的 `(a, ā)` 相邻，加反常配对块）。`data` 为
`2^L` 个 Grassmann monomial 系数（`L = length(basis)`，腿 1 最快）。

节点按 `contour_ordering` 排列：虚时/实时时间升序（实时每时刻 `+` 先 `−` 后），
混合为虚时升序接实时升序；每个节点内按 band 1,2,…，band 内 `(a, ā)`。
下游 lattice 的 ordering 若不同，比较时对稠密张量做一次 `permutedims` 即可；
GTEMPO 链上的左边界辅助腿不进入 IF。

元素类型 `T` 取自 correlation function（虚时通常为实）。

腿数满足 `2^L ≤ MAXTENSORDIM`（`2^30` 个系数即内存上限）；超过时抛
`ArgumentError`。

前缀割（链前 `k` 个节点块，即腿 `1..2k·bands`）的 plain SVD 熵等于
系数态的 OSEE（协方差快速算法，`QuAPIGaussianStatesExt`：
`FermionicGaussianState(corr)` + `entropy`），可用于交叉验证；节点 `p`
在 Δ 行序（`[:+; :-; :τ]`，分支内时间升序）中的指标即 GaussianStates
系数态的模式指标。
"""
function GrassmannIF(h::GrassmannHybridizationMatrix)
    L = length(h.basis)
    L ≤ log2(MAXTENSORDIM) || throw(ArgumentError(
        "稠密 IF 张量阶数 L = $(L) 超过上限 2^$(Int(log2(MAXTENSORDIM)))（2^L 个系数将耗尽内存）；" *
        "请减小 contour 节点数 n 或 flavor 数（当前 L = $(L)），" *
        "或使用协方差快速算法"))
    return GrassmannADT(h.basis, _fold_grassmann(h.Δ))
end

GrassmannIF(x::IFSource) = GrassmannIF(GrassmannHybridizationMatrix(x))

"""
    *(A::GrassmannADT, B::GrassmannADT) -> GrassmannADT

两个稠密 Grassmann IF 的乘法，允许 `basis` 不同：共享腿按 Grassmann
乘法融合（`1·ξ = ξ`、`ξ·ξ = 0`，重排符号计入），只在一边出现的腿成为
外腿。结果 `basis` 为 `A.basis` 后接 `B` 中不在 `A.basis` 里的腿；
`data` 按该顺序展开（腿 1 最快）。
"""
function Base.:*(A::GrassmannADT, B::GrassmannADT)
    ubasis, posB = _union_basis(A.basis, B.basis)
    LA, LB, L = length(A.basis), length(B.basis), length(ubasis)
    T = promote_type(eltype(A.data), eltype(B.data))
    data = zeros(T, 1 << L)
    for mA in 0:(1 << LA) - 1
        a = A.data[mA + 1]
        iszero(a) && continue
        for mB in 0:(1 << LB) - 1
            b = B.data[mB + 1]
            iszero(b) && continue
            bitsB = 0
            ok = true
            for k in 1:LB
                ((mB >> (k - 1)) & 1) == 0 && continue
                p = posB[k]
                if ((mA >> (p - 1)) & 1) == 1 || ((bitsB >> (p - 1)) & 1) == 1
                    ok = false
                    break
                end
                bitsB |= 1 << (p - 1)
            end
            ok || continue
            # A 的变量在前（并基序升序）；B 的变量按 B 的 basis 序，统计
            # 把它们重排成并基升序所需的逆序数
            inv = 0
            for k in 1:LB
                ((mB >> (k - 1)) & 1) == 0 && continue
                p = posB[k]
                inv += count_ones(mA & ~((1 << p) - 1))
                for l in 1:(k - 1)
                    ((mB >> (l - 1)) & 1) != 0 && posB[l] > p && (inv += 1)
                end
            end
            data[(mA | bitsB) + 1] += isodd(inv) ? -a * b : a * b
        end
    end
    return GrassmannADT(ubasis, data)
end

"""
    Zvalue(x::GrassmannADT)

稠密 Grassmann IF 的 Berezin 积分：全腿单项式的系数（`x.data` 的末项）。
"""
Zvalue(x::GrassmannADT) = x.data[end]
