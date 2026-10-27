"""
    BosonicADT{T<:Number}

玻色影响泛函在对角（ADT）耦合下的稠密表示：`basis` 为 `n` 个 contour
节点（链序），`d` 为物理维度（耦合算符本征值矢量的长度），`data` 为
`dⁿ` 个路径权重（节点 1 最快）——指标数与 `basis` 数一致。用
[`BosonicIF`](@ref) 从 correlation function 构造。
"""
struct BosonicADT{T<:Number}
    basis::Vector{BosonicContourIndex}
    d::Int
    data::Vector{T}
end

"""
    BosonicPT{T<:Number}

玻色影响泛函在非对角（PT）耦合下的稠密表示：`basis` 为 `n` 个 contour
节点（链序），`d` 为系统维度，`data` 为 `d^{2n}` 个矩阵元（指标 2p−1 =
节点 p 的 ket 路径、2p = 节点 p 的 bra 路径，指标 1 最快）——指标数是
`basis` 数的两倍。用 [`BosonicIF`](@ref) 从 correlation function 构造。
"""
struct BosonicPT{T<:Number}
    basis::Vector{BosonicContourIndex}
    d::Int
    data::Vector{T}
end

# ------------------------------------------------------------------ #
# 玻色：ADT（对角耦合）与 PT（非对角耦合）
# ------------------------------------------------------------------ #

"""
    BosonicSplitDelta

contour 杂交矩阵 `Δ = M1 + M2` 的**等时算符次序拆分**（`CorrelationMatrix`
块用 `ηⱼₖ`/`ηₖⱼ` 区分等时两种次序），仅供 [`BosonicIF`](@ref) 的 PT 构造。

它与 BdG 的 `(h, g)` 块（粒子-空穴 / 配对）含义无关，故用独立类型包装：
不会被误当作两个矩阵 splat 进 `bosonic_bdg_matrix(h, g)`。
"""
struct BosonicSplitDelta{T<:Number}
    M1::Matrix{T}
    M2::Matrix{T}
end

function BosonicSplitDelta(M1::AbstractMatrix, M2::AbstractMatrix)
    T = promote_type(eltype(M1), eltype(M2))
    return BosonicSplitDelta{T}(Matrix{T}(M1), Matrix{T}(M2))
end

"""
    BosonicHybridizationMatrix

玻色稠密 IF 的杂交矩阵：`basis` 为链序节点，`Δ` 为 contour 杂交矩阵
（节点行序 `[:+; :-; :τ]`）。

- `split=false`（默认）：`Δ = M1 + M2`（普通矩阵），供 [`BosonicIF`](@ref) 的
  ADT 构造；
- `split=true`：`Δ = ` [`BosonicSplitDelta`](@ref)`(M1, M2)` 为算符次序拆分，
  供 PT 构造。
"""
struct BosonicHybridizationMatrix{T}
    basis::Vector{BosonicContourIndex}
    Δ::T
end

function BosonicHybridizationMatrix(x::AbstractCorrelationFunction; split::Bool=false)
    chain = contour_ordering(x)
    bs = BosonicContourIndex[BosonicContourIndex(j, br) for (br, j) in chain]
    Δ = split ? BosonicSplitDelta(_contour_delta_split(x)...) : _contour_delta(x)
    return BosonicHybridizationMatrix(bs, Δ)
end

"""
    BosonicIF(h::BosonicHybridizationMatrix, λ::AbstractVector{<:Real})
        -> BosonicADT

玻色影响泛函在**对角（ADT）耦合**下的稠密表示。耦合形式

```text
H_int = hyb ⊗ Σ_k g_k (a†_k + a_k) ,   hyb = diag(λ)
```

在 `hyb` 的本征基下每个 contour 节点的路径值取其本征值，`data` 元素为

```text
I[v₁..v_n] = exp(Σ_pq Δ_pq v_p v_q)
```

其中 `Δ = h.Δ`（`split=false` 构造的杂交矩阵，节点行序 `[:+; :-; :τ]`），
`n` 为节点数；与 brute-force QuAPI 的 `exp(dot(v, G, v))` 约定一致（含
对角项与上下三角全求和）。`λ` 为本征值矢量（自旋-玻色取 `[1, -1]` 等；
顺序决定张量指标到路径值的映射）；`basis` 为链序节点（见
`contour_ordering`），`d = length(λ)`。

复杂度 O(n² dⁿ)；`dⁿ > MAXTENSORDIM`（`2³⁰`）时抛 `ArgumentError`
（d=2 即 n ≤ 30）。
"""
function BosonicIF(h::BosonicHybridizationMatrix,
                         λ::AbstractVector{<:Real})
    isempty(λ) && throw(ArgumentError("λ 不能为空"))
    h.Δ isa AbstractMatrix || throw(ArgumentError(
        "该 BosonicHybridizationMatrix 存的是 PT 的 BosonicSplitDelta(M1, M2)；" *
        "ADT 构造需要 split=false"))
    d = length(λ)
    Δ = h.Δ
    chain = h.basis
    n = length(chain)
    n * log2(d) ≤ log2(MAXTENSORDIM) || throw(ArgumentError(
        "稠密 IF 张量总元素数 d^n = $(d)^$(n) 超过上限 MAXTENSORDIM = 2³⁰；" *
        "请减小 contour 节点数 n（当前 $(n)）或路径值个数 d（当前 $(d)）"))
    T = eltype(Δ)
    data = ones(T, ntuple(_ -> d, n))
    # 逐对吸收 exp(Δ_pq λ_a λ_b)：I = Π_pq exp(Δ_pq v_p v_q)
    # Δ 的行列序与 h.basis（contour_ordering）一致 ⟹ 节点 p ↔ Δ[p, :]
    for p in 1:n, q in 1:n
        c = Δ[p, q]
        iszero(c) && continue
        F = exp.(c .* λ .* transpose(λ))
        @inbounds for I in CartesianIndices(data)
            data[I] *= F[I[p], I[q]]
        end
    end
    return BosonicADT(chain, d, vec(data))
end

# 便捷入口：一步从 correlation function 构造（内部先算杂交矩阵）
BosonicIF(x::AbstractCorrelationFunction, λ::AbstractVector{<:Real}) =
    BosonicIF(BosonicHybridizationMatrix(x), λ)

"""
    BosonicIF(h::BosonicHybridizationMatrix, hyb::AbstractMatrix)
        -> BosonicPT

玻色影响泛函在**非对角（PT）耦合**下的稠密表示，`hyb` 为任意方阵
（`NonDiagonalHyb` 约定：`op*a + op'*a'`，`op = hyb`，不要求厄米）。
在 `hyb` 的本征基以外，IF 不再是路径权重的对角嵌入，而是作用在系统
路径上的算符：

```text
I = exp(H)
```

`H` 为 contour 节点对上的双线性型。**精确性**（对单模与双模 bath 的
严格 IF——精确 cell 演化算符 + 偏迹环境——数值验证）：对角 `hyb`、
`op1 = op2`（厄米/实对称）及 `[op1,op2] = 0` 的情形 IF 严格为
Gaussian，上式**机器精度成立**；一般 `op1 ≠ op2` 时严格 IF 含 O(V⁴)
非高斯累积量（Magnus 对易子项），上式与严格 IF 一致到 **O(V²)**
（生成元逐元素），全链残差按 **V⁴** 标度。

设 `op1 = hyb`、`op2 = hyb'`。节点 ν 的**行腿矩阵**（张量指标为
(out, in)，元素语言：+ 支为 `op[a,m]`、− 支为 `conj(op[b,mm])`、
τ 支为 `op[a,m]`）：

```text
X1[ν] = op1ᵀ（+ 分支）  op1'（− 分支）  op1（τ 分支）
X2[ν] = op2ᵀ（+ 分支）  op2'（− 分支）  op2（τ 分支）
```

`M1 + M2 = Δ` 是杂交矩阵的等时算符次序拆分（`CorrelationMatrix` 块用
`ηⱼₖ` / `ηₖⱼ` 区分等时两种次序，见 `_contour_delta_split`；因此 PT
需要 `split=true` 构造的 `BosonicHybridizationMatrix`）。`H` 按块组装
（`j_ν` 为节点 ν 的 contour 步数；纯实时/纯虚时轮廓各对均适用；混合
轮廓中 ± 之间的对用纯 RT 规则）：

- **同支不同节点**（`(+,+)`、`(−,−)`、`(τ,τ)`，ν 靠后 / μ 靠前）：
  `M1[ν,μ]·X1[ν]·X2[μ] + M2[μ,ν]·X2[ν]·X1[μ]`（`−` 分支 M1↔M2 互换）；
- **同节点等时**（矩阵积作用于单节点，`(op·op')` 型双插入）：`+` 分支
  `M1[νν]·X2·X1 + M2[νν]·X1·X2`；`−` 与 `τ` 分支 M1↔M2 互换；
- **跨支对**（`(+ν, −μ)`，含同 cell 等时，从 + 侧处理）：**求和核配
  同算符**——`(M1+M2)[ν,μ]·X1[ν]·X1[μ] + (M1+M2)[μ,ν]·X2[ν]·X2[μ]`；
- **τ×± 交叉对**（混合轮廓，从 τ 侧处理）：物理收缩 `b@τ ↔ b†@±`
  给出**混合配对**——核留原槽位，± 侧算符标签与跨支相反：
  `(M1+M2)[ν,μ]·X1[ν]·X̃[μ] + (M1+M2)[μ,ν]·X2[ν]·X̃'[μ]`，其中
  `X̃[μ] = X2[μ]`（μ 为 + 节点）或 `X1[μ]`（μ 为 − 节点）。

对 `op1 = op2`（实对称/厄米 `hyb`）或对角 `hyb`，以上全部退化为
`H = Σ M1·op1·op2 + M2·op2·op1`，后者进一步退化为 ADT 路径权重的
对角嵌入。这里直接对显式构造的 `H` 做矩阵指数，**不含任何 Prony/
指数展开近似**。

`data` 为 `d^{2n}` 个矩阵元：指标 `2p−1` = 节点 `p` 的 ket 路径，
`2p` = 节点 `p` 的 bra 路径，指标 1 最快；`basis` 为链序节点，
`d = size(hyb, 1)`。总元素数 `d^{2n} > MAXTENSORDIM`（`2³⁰`）时抛
`ArgumentError`（d=2 即 n ≤ 15）。
"""
function BosonicIF(h::BosonicHybridizationMatrix, hyb::AbstractMatrix)
    size(hyb, 1) == size(hyb, 2) || throw(ArgumentError("hyb 必须是方阵"))
    h.Δ isa BosonicSplitDelta || throw(ArgumentError(
        "PT 构造需要 split=true 的 BosonicHybridizationMatrix（存 BosonicSplitDelta(M1, M2)）"))
    d = size(hyb, 1)
    M1, M2 = h.Δ.M1, h.Δ.M2
    chain = h.basis
    # M1/M2 已按 contour_ordering（= h.basis）组装，无需再重排
    n = length(chain)
    2n * log2(d) ≤ log2(MAXTENSORDIM) || throw(ArgumentError(
        "稠密 PT 张量总元素数 d^{2n} = $(d)^$(2n) 超过上限 MAXTENSORDIM = 2³⁰；" *
        "请减小 contour 节点数 n（当前 $(n)）或系统维度 d（当前 $(d)）"))
    T = promote_type(eltype(M1), eltype(hyb))
    op1 = Matrix{T}(hyb)
    op2 = Matrix{T}(adjoint(hyb))
    # 行腿矩阵：+ 支（ket）的元素为 op[a,m] ⇒ 矩阵 opᵀ；− 支（bra）的
    # 元素为 conj(op[b,mm]) ⇒ 矩阵 op'；τ 支无 bra 结构，用普通矩阵
    x1 = [l.branch === :+ ? Matrix{T}(transpose(op1)) :
          l.branch === :- ? Matrix{T}(adjoint(op1)) : op1 for l in chain]
    x2 = [l.branch === :+ ? Matrix{T}(transpose(op2)) :
          l.branch === :- ? Matrix{T}(adjoint(op2)) : op2 for l in chain]
    dim = d^n
    X1 = [_op_at(x1[p], p, n) for p in 1:n]
    X2 = [_op_at(x2[p], p, n) for p in 1:n]
    H = zeros(T, dim, dim)
    mixed = any(l.branch === :τ for l in chain)
    for p in 1:n, q in 1:n
        bp, bq = chain[p].branch, chain[q].branch
        if (bp === :+ && bq === :-) || (bp === :- && bq === :+)
            # 跨支对：统一从 + 侧处理一次。求和核（完整跨支块）配同算符（纯 RT 规则）
            bp === :+ || continue
            c = M1[p, q] + M2[p, q]           # (+,−) 块完整值
            iszero(c) || (H .+= c .* (X1[p] * X1[q]))
            c = M1[q, p] + M2[q, p]           # (−,+) 块完整值
            iszero(c) || (H .+= c .* (X2[p] * X2[q]))
        elseif bp === :τ && bq !== :τ
            # τ×± 交叉对（MX，从 τ 侧处理一次）：物理收缩 b@τ ↔ b†@± ⟹ 混合配对，
            # τ 侧放 X1/X2（同 M1 槽位），± 侧取与 RT 跨支相反的算符标签
            # （(+,−) 块配 X2@+/X1@−，(−,+) 块配 X1@+/X2@−）
            c = M1[p, q] + M2[p, q]
            iszero(c) || (H .+= c .* (X1[p] * (bq === :+ ? X2[q] : X1[q])))
            c = M1[q, p] + M2[q, p]
            iszero(c) || (H .+= c .* (X2[p] * (bq === :+ ? X1[q] : X2[q])))
        elseif bq === :τ && bp !== :τ
            continue                              # 已从 τ 侧处理
        elseif bp === bq
            minus = bp === :-                 # − 支非对角：M1↔M2 互换（τ 同 + 支）
            flipd = bp !== :+                 # 同节点等时：+ 支 M1↔X2·X1；− / τ 支 M1↔X1·X2
            jp, jq = chain[p].j, chain[q].j
            if jp > jq                        # 同支不同节点（ν 靠后）
                c = M1[p, q]
                iszero(c) || (H .+= c .* (minus ? X2[p] * X1[q] : X1[p] * X2[q]))
            elseif jp < jq
                c = M2[p, q]
                iszero(c) || (H .+= c .* (minus ? X1[q] * X2[p] : X2[q] * X1[p]))
            else                              # 同节点等时：矩阵积于单节点
                c = M1[p, p]
                iszero(c) || (H .+= c .* (flipd ? X1[p] * X2[p] : X2[p] * X1[p]))
                c = M2[p, p]
                iszero(c) || (H .+= c .* (flipd ? X2[p] * X1[p] : X1[p] * X2[p]))
            end
        end
    end
    F = exp(H)
    # reshape 得 (ket₁..ket_n, bra₁..bra_n)，再交错成每节点 (ket, bra)
    Ten = reshape(F, ntuple(_ -> d, 2n))
    perm = vcat([[p, n + p] for p in 1:n]...)
    return BosonicPT(chain, d, vec(permutedims(Ten, perm)))
end

# 便捷入口：一步从 correlation function 构造（内部先算拆分后的杂交矩阵）
BosonicIF(x::AbstractCorrelationFunction, hyb::AbstractMatrix) =
    BosonicIF(BosonicHybridizationMatrix(x; split=true), hyb)

# op 作用在 d^n 维空间的第 p 个张量因子上（节点 1 为最快指标，与 data 的
# 线性序一致）：从慢到快的因子依次为节点 n, n−1, …, 1
function _op_at(op::AbstractMatrix{T}, p::Int, n::Int) where {T}
    d = size(op, 1)
    slow = d^(n - p)      # 节点 n..p+1
    fast = d^(p - 1)      # 节点 p−1..1
    return kron(Matrix{T}(I, slow, slow), op, Matrix{T}(I, fast, fast))
end

"""
    *(A::BosonicADT, B::BosonicADT) -> BosonicADT

ADT 路径权重的逐点（Hadamard）乘积，允许节点集合不同：共享节点上的
指标对齐相乘，独有节点取外积。
"""
function Base.:*(A::BosonicADT, B::BosonicADT)
    A.d == B.d || throw(DimensionMismatch(
        "path dimension mismatch: $(A.d) vs $(B.d)"))
    ubasis, posB = _union_basis(A.basis, B.basis)
    nA, nB, n = length(A.basis), length(B.basis), length(ubasis)
    T = promote_type(eltype(A.data), eltype(B.data))
    data = zeros(T, ntuple(_ -> A.d, n))
    for I in CartesianIndices(data)
        ia, stride = 0, 1
        for p in 1:nA
            ia += (I[p] - 1) * stride
            stride *= A.d
        end
        ib, stride = 0, 1
        for q in 1:nB
            ib += (I[posB[q]] - 1) * stride
            stride *= B.d
        end
        data[I] = A.data[ia + 1] * B.data[ib + 1]
    end
    return BosonicADT(ubasis, A.d, vec(data))
end

"""
    Zvalue(x::BosonicADT)

ADT 路径权重之和（对角耦合下的配分函数 / 迹）。
"""
Zvalue(x::BosonicADT) = sum(x.data)

"""
    margin(A::BosonicADT, legs::BosonicContourIndex...) -> BosonicADT

对未列出的节点求和（边缘化），返回只含 `legs` 的 ADT；`legs` 的顺序即
结果 `basis` 的顺序。`margin(A)` 返回不含节点的 IF，其唯一元素为
`Zvalue(A)`。
"""
function margin(A::BosonicADT, legs::BosonicContourIndex...)
    pos = _leg_positions(A.basis, legs)
    d, n, r = A.d, length(A.basis), length(pos)
    T = reshape(A.data, ntuple(_ -> d, n))
    out = zeros(eltype(A.data), ntuple(_ -> d, r))
    for I in CartesianIndices(T)
        J = CartesianIndex(ntuple(q -> I[pos[q]], r))
        out[J] += T[I]
    end
    return BosonicADT(BosonicContourIndex[legs...], d, vec(out))
end


# PT 的矩阵视图：行 = 各节点 ket（节点 1 最快），列 = 各节点 bra
function _pt_matrix(P::BosonicPT)
    d, n = P.d, length(P.basis)
    n == 0 && return reshape([P.data[1]], 1, 1)
    T = reshape(P.data, ntuple(_ -> d, 2n))
    perm = vcat(collect(1:2:2n), collect(2:2:2n))   # (ket₁..ketₙ, bra₁..braₙ)
    return reshape(permutedims(T, perm), d^n, d^n)
end

function _pt_from_matrix(M::AbstractMatrix, basis, d::Int)
    n = length(basis)
    n == 0 && return BosonicPT(collect(basis), d, [M[1]])
    T = reshape(M, ntuple(_ -> d, 2n))              # (ket₁..ketₙ, bra₁..braₙ)
    perm = vcat([[p, n + p] for p in 1:n]...)       # → (ket₁, bra₁, ket₂, bra₂, …)
    return BosonicPT(collect(basis), d, vec(permutedims(T, perm)))
end

# 把 basis 为 ubasis 子集的 PT 提升到 ubasis：缺失节点取恒等超算符
function _pt_matrix_on(P::BosonicPT, ubasis, posP)
    M = _pt_matrix(P)
    n, nP, d = length(ubasis), length(P.basis), P.d
    n == nP && return M
    nP == 0 && return M[1] .* Matrix{eltype(M)}(I, d^n, d^n)
    dims = ntuple(_ -> d, n)
    lin = LinearIndices(dims)
    linP = LinearIndices(ntuple(_ -> d, nP))
    out = zeros(eltype(M), d^n, d^n)
    for I in CartesianIndices(dims), J in CartesianIndices(dims)
        ok = true
        for p in 1:n
            (p in posP || I[p] == J[p]) || (ok = false; break)
        end
        ok || continue
        IP = CartesianIndex(ntuple(q -> I[posP[q]], nP))
        JP = CartesianIndex(ntuple(q -> J[posP[q]], nP))
        out[lin[I], lin[J]] = M[linP[IP], linP[J]]
    end
    return out
end

"""
    *(A::BosonicPT, B::BosonicPT) -> BosonicPT

PT 视为矩阵（行 = 各节点 ket、列 = 各节点 bra）的矩阵乘积
`reshape(A) * reshape(B)`；basis 不同时缺失节点按恒等超算符补齐，
结果定义在并基上。
"""
function Base.:*(A::BosonicPT, B::BosonicPT)
    A.d == B.d || throw(DimensionMismatch(
        "path dimension mismatch: $(A.d) vs $(B.d)"))
    ubasis, posB = _union_basis(A.basis, B.basis)
    posA = collect(1:length(A.basis))
    return _pt_from_matrix(_pt_matrix_on(A, ubasis, posA) *
                           _pt_matrix_on(B, ubasis, posB), ubasis, A.d)
end

"""
    Zvalue(x::BosonicPT)

PT 的闭合收缩：每个时间点的 ket（输出）与下一个时间点的 bra（输入）
收缩，边界按 contour 类别闭合（与 TEMPO 的 `integrate` 一致）：

- 虚时：`bra_p = ket_{p+1}`，首尾 cyclic（`bra_n = ket_1`）；
- 实时：最早时刻 `bra₊ = bra₋`，相邻时刻 `bra_{±,i} = ket_{±,i-1}`，
  最晚时刻 `ket₊ = ket₋`；
- 混合：τ 支 `bra_p = ket_{p+1}` 链式收缩，再与实时支按 Kadanoff
  轮廓的端点连接，实时支同实时情形。

`basis` 的分支/时间标签决定收缩方式（与 `basis` 的内部顺序无关）。
"""
function Zvalue(x::BosonicPT)
    d, n = x.d, length(x.basis)
    data = reshape(x.data, ntuple(_ -> d, 2n))
    cons = _pt_constraints(x)
    s = zero(eltype(x.data))
    for I in CartesianIndices(data)
        ok = true
        for (a, b) in cons
            I[a] == I[b] || (ok = false; break)
        end
        ok && (s += data[I])
    end
    return s
end


# PT 闭合收缩的指标等式约束（指标 2p−1 = ket/输出，2p = bra/输入）
function _pt_constraints(x::BosonicPT)
    pos = Dict(b => k for (k, b) in enumerate(x.basis))
    ket(p) = 2p - 1
    bra(p) = 2p
    cons = Tuple{Int,Int}[]
    brs = unique(b.branch for b in x.basis)
    if brs == [:τ]                                  # 虚时
        n = length(x.basis)
        for p in 1:n-1
            push!(cons, (bra(p), ket(p + 1)))
        end
        push!(cons, (bra(n), ket(1)))
    else
        js = sort(unique(b.j for b in x.basis if b.branch != :τ))
        jmin, jmax = first(js), last(js)
        pplus, pminus = pos[BosonicContourIndex(jmin, :+)],
                        pos[BosonicContourIndex(jmin, :-)]
        if :τ ∉ brs                                 # 实时
            push!(cons, (bra(pplus), bra(pminus)))
        else                                        # 混合
            τs = [pos[b] for b in x.basis if b.branch == :τ]
            for k in 1:length(τs)-1
                push!(cons, (bra(τs[k]), ket(τs[k + 1])))
            end
            push!(cons, (bra(τs[end]), bra(pminus)))
            push!(cons, (bra(pplus), ket(τs[1])))
        end
        for (k, j) in enumerate(js)
            p1, p2 = pos[BosonicContourIndex(j, :+)],
                      pos[BosonicContourIndex(j, :-)]
            if k > 1
                jp = js[k - 1]
                push!(cons, (bra(p1), ket(pos[BosonicContourIndex(jp, :+)])))
                push!(cons, (bra(p2), ket(pos[BosonicContourIndex(jp, :-)])))
            end
            k == length(js) && push!(cons, (ket(p1), ket(p2)))
        end
    end
    return cons
end
