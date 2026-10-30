"""
    abstract type AbstractRealCorrelationFunction <: AbstractCorrelationFunction

Abstract supertype of real-time hybridization representations.  `branch`,
`branches`, `contour_ordering` are defined uniformly on this level; the
concrete `RealCorrelationFunction` stores the two diagonal Keldysh blocks
(see its docstring for the unified branch convention).
"""
abstract type AbstractRealCorrelationFunction <: AbstractCorrelationFunction end

branches(::AbstractRealCorrelationFunction) = (:+, :-)

"""
    struct RealCorrelationFunction{A} <: AbstractRealCorrelationFunction

Real-time hybridization function on the closed Keldysh contour, storing only
the two CROSS branch blocks — `G₊₋` (the (+,−) block, the lesser-type
correlator Δ^<) and `G₋₊` (the (−,+) block, the greater-type Δ^>) — the two
physically fundamental Keldysh functions — both of the SAME parametric
matrix type `A < AbstractMatrix{ComplexF64}`.  With `A =
CorrelationMatrix` (Toeplitz-vector storage) the pair carries the four
independent Toeplitz vectors (a, b, c, d) of a general complex-valued
spectrum; the DIAGONAL blocks are recovered by `_G₊₊(G₊₋, G₋₊)` /
`_G₋₋(G₊₋, G₋₊)`.

Unified branch convention (statistics independent; identical to
`MixedCorrelationFunction` on the Kadanoff contour — the CROSS blocks carry
the Keldysh relative minus signs inside the blocks):

    (+,+): (a, b)                (+,−): (−c, −b)
    (−,+): (−a, −d)              (−,−): (c, d)

For `CorrelationMatrix` storage the diagonal blocks are derived by `branch`
exactly — pure sign combinations of the stored halves, no conjugation, hence
valid for arbitrarily complex-valued spectra.  The sharing was audited at
μ ≠ 0 (real and complex spectra) against the historical four-block formulas:
it is exact (structural identity of the construction, not a thermal one; see
the `cross-block sharing at μ ≠ 0` testset in test/realtime.jl).

For a general (non-Toeplitz) matrix type `A` the cross branches return the
stored blocks directly; the diagonal branches are derived by the sector rule
(triangles from the partner blocks, diagonal by convention — see `_G₊₊` /
`_G₋₋`), mirroring the general-matrix handling of the previous
diagonal-pair storage.

`i, j` run over the `N+1` nodes tᵢ = i·δt of each branch.

The exponent of the discretized Feynman-Vernon influence functional built
from these blocks depends on the statistics of the bath (the minus sign is
absorbed into Δ, see `AbstractCorrelationFunction`).  With branch labels
`b₁, b₂ ∈ (:+, :-)` and `i, j` running over the `N+1` nodes of the
respective branch,

- fermionic bath, with Grassmann fields `a` (annihilation) and `ā`
  (conjugate):

      Σ_{b₁b₂} Σᵢⱼ āᵢ^{b₁} Δᵢⱼ^{b₁b₂} aⱼ^{b₂} ;

- bosonic bath with the coupling `s ⊗ B† + s† ⊗ B` (the PT case of TEMPO,
  also covering non-Hermitian `s`, e.g. a ladder operator):

      Σ_{b₁b₂} Σᵢⱼ s†ᵢ^{b₁} Δᵢⱼ^{b₁b₂} sⱼ^{b₂} ,

  where `s^{b}` is the system field on branch `b`.  For Hermitian `s`
  (`s† = s`) the exponent degenerates to the single-operator (ADT) form.
"""
struct RealCorrelationFunction{A<:AbstractMatrix{ComplexF64}} <: AbstractRealCorrelationFunction
    G₊₋::A
    G₋₊::A
end

"""
    RealCorrelationFunction(G₊₊, G₋₋)

Two-matrix constructor (auto-generated): stores the diagonal pair directly;
both arguments must be of the same matrix type (single type parameter).

    RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋) -> RealCorrelationFunction

Legacy four-block constructor (any matrix type): validates the unified
convention — for `CorrelationMatrix` inputs at the vector level (exact), for
general matrices at the sector level (strict triangles pairwise + the
Keldysh diagonal identity) — and keeps only the independent CROSS pair.
"""
function RealCorrelationFunction(G₊₊::CorrelationMatrix{ComplexF64},
                                 G₊₋::CorrelationMatrix{ComplexF64},
                                 G₋₊::CorrelationMatrix{ComplexF64},
                                 G₋₋::CorrelationMatrix{ComplexF64})
    ok = isapprox(G₊₋.ηⱼₖ, -G₋₋.ηⱼₖ; atol=1.0e-12) &&
         isapprox(G₊₋.ηₖⱼ, -G₊₊.ηₖⱼ; atol=1.0e-12) &&
         isapprox(G₋₊.ηⱼₖ, -G₊₊.ηⱼₖ; atol=1.0e-12) &&
         isapprox(G₋₊.ηₖⱼ, -G₋₋.ηₖⱼ; atol=1.0e-12)
    ok || throw(ArgumentError("cross blocks do not follow the unified convention " *
                              "(G₊₋ = CM(−G₋₋.ηⱼₖ, −G₊₊.ηₖⱼ), G₋₊ = CM(−G₊₊.ηⱼₖ, −G₋₋.ηₖⱼ))"))
    return RealCorrelationFunction(G₊₋, G₋₊)
end

function RealCorrelationFunction(G₊₊::AbstractMatrix{ComplexF64},
                                 G₊₋::AbstractMatrix{ComplexF64},
                                 G₋₊::AbstractMatrix{ComplexF64},
                                 G₋₋::AbstractMatrix{ComplexF64})
    # sector self-consistency: the four blocks must share their strict
    # triangles pairwise (and the diagonals must satisfy the Keldysh
    # identity).  The cross-diagonal contact split itself is fixed by the
    # convention of `_G₊₊` / `_G₋₋`.
    n = size(G₊₊, 1)
    (size(G₊₋, 1) == n && size(G₋₊, 1) == n && size(G₋₋, 1) == n) ||
        throw(DimensionMismatch("all four blocks must be $n×$n"))
    ok = all(G₊₋[i, j] ≈ -G₋₋[i, j] for i in 1:n, j in 1:n if i > j) &&
         all(G₊₋[i, j] ≈ -G₊₊[i, j] for i in 1:n, j in 1:n if i < j) &&
         all(G₋₊[i, j] ≈ -G₊₊[i, j] for i in 1:n, j in 1:n if i > j) &&
         all(G₋₊[i, j] ≈ -G₋₋[i, j] for i in 1:n, j in 1:n if i < j) &&
         all((G₊₊[i, i] + G₋₋[i, i] + G₊₋[i, i] + G₋₊[i, i]) ≈ 0 for i in 1:n)
    ok || throw(ArgumentError("blocks are not sector-consistent (unified " *
                              "convention): strict triangles must pairwise coincide " *
                              "up to sign and the diagonals must satisfy " *
                              "G₊₊ + G₋₋ + G₊₋ + G₋₊ = 0"))
    return RealCorrelationFunction(G₊₋, G₋₊)
end

# —— 从交叉对推导对角块的两个函数（branch 与 Mixed 共用） ——————————
# 统一约定下 G₊₋ = CM(−c, −b)、G₋₊ = CM(−a, −d)，故
#   G₊₊ = CM(a, b) = CM(−G₋₊.ηⱼₖ, −G₊₋.ηₖⱼ)
#   G₋₋ = CM(c, d) = CM(−G₊₋.ηⱼₖ, −G₋₊.ηₖⱼ)
# —— 纯符号组合、无共轭，对任意复值谱精确。

"""
    _G₊₊(G₊₋, G₋₊) -> G₊₊ block

Time-ordered (+,+) block from the stored cross pair.  Exact (pure sign
combinations) for `CorrelationMatrix` inputs; for general matrices the
sector rule applies: lower triangle from −G₋₊, upper triangle and diagonal
from −G₊₋ (the diagonal convention mirrors the general-matrix cross
derivation of the previous diagonal-pair storage).
"""
_G₊₊(G₊₋::CorrelationMatrix{ComplexF64}, G₋₊::CorrelationMatrix{ComplexF64}) =
    CorrelationMatrix(.-G₋₊.ηⱼₖ, .-G₊₋.ηₖⱼ)
function _G₊₊(G₊₋::AbstractMatrix{ComplexF64}, G₋₊::AbstractMatrix{ComplexF64})
    n = size(G₊₋, 1)
    M = zeros(ComplexF64, n, n)
    for i in 1:n, j in 1:n
        M[i, j] = i > j ? -G₋₊[i, j] : -G₊₋[i, j]
    end
    return M
end

"""
    _G₋₋(G₊₋, G₋₊) -> G₋₋ block

Anti-time-ordered (−,−) block from the stored cross pair; same exactness
and sector-rule conventions as `_G₊₊`.
"""
_G₋₋(G₊₋::CorrelationMatrix{ComplexF64}, G₋₊::CorrelationMatrix{ComplexF64}) =
    CorrelationMatrix(.-G₊₋.ηⱼₖ, .-G₋₊.ηₖⱼ)
function _G₋₋(G₊₋::AbstractMatrix{ComplexF64}, G₋₊::AbstractMatrix{ComplexF64})
    n = size(G₊₋, 1)
    M = zeros(ComplexF64, n, n)
    for i in 1:n, j in 1:n
        M[i, j] = i > j ? -G₊₋[i, j] : -G₋₊[i, j]
    end
    return M
end

# transpose 交换两交叉块（[xᵀ]^{+−} = (x^{−+})ᵀ）
Base.transpose(x::RealCorrelationFunction) = RealCorrelationFunction(transpose(x.G₋₊), transpose(x.G₊₋))
Base.iszero(x::RealCorrelationFunction) = iszero(x.G₊₋) && iszero(x.G₋₊)

Base.size(x::RealCorrelationFunction, i::Int...) = size(x.G₊₋, i...)

function Base.show(io::IO, ::MIME"text/plain", A::RealCorrelationFunction)
    print(io, "Real Correlation Function [$(size(A.G₊₋, 1))] (2-block cross)")
end
Base.:-(x::RealCorrelationFunction) = RealCorrelationFunction(-x.G₊₋, -x.G₋₊)
Base.:+(A::RealCorrelationFunction, B::RealCorrelationFunction) = RealCorrelationFunction(A.G₊₋ + B.G₊₋, A.G₋₊ + B.G₋₊)

# branch for Toeplitz-vector (CorrelationMatrix) storage: cross blocks derived
# exactly by pure sign combinations of the stored halves (no conjugation —
# valid for arbitrarily complex-valued spectra)
function branch(x::RealCorrelationFunction{CorrelationMatrix{ComplexF64}}, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    if b1 == :+
        if b2 == :+
            return _G₊₊(x.G₊₋, x.G₋₊)
        else
            return x.G₊₋
        end
    else
        if b2 == :+
            return x.G₋₊
        else
            return _G₋₋(x.G₊₋, x.G₋₊)
        end
    end
end

# fallback for general (non-Toeplitz) matrix storage: the cross branches
# return the stored blocks; the diagonal branches are derived by the sector
# rule through _G₊₊ / _G₋₋ (triangles from the partner cross blocks,
# diagonal by convention).  The four Keldysh blocks of any correlator share
# their strict triangles pairwise, so the cross pair determines all four
# triangles; only the per-node CONTACT SPLIT of the two diagonal diagonals
# (n extra degrees of freedom) is not carried, fixed by the convention
# exactly as in the previous diagonal-pair storage.  CorrelationMatrix
# storage carries the contacts exactly — for Toeplitz data prefer it.
function branch(x::RealCorrelationFunction, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    if b1 == :+
        return b2 == :+ ? _G₊₊(x.G₊₋, x.G₋₊) : x.G₊₋
    else
        return b2 == :+ ? x.G₋₊ : _G₋₋(x.G₊₋, x.G₋₊)
    end
end

branch(x::RealCorrelationFunction; b1::Symbol, b2::Symbol) = branch(x, b1, b2)

# equality/comparison on the stored diagonal pair (equivalent to block-wise
# comparison: the cross blocks are functions of the pair)
Base.:(==)(x::RealCorrelationFunction, y::RealCorrelationFunction) = (x.G₊₋ == y.G₊₋) && (x.G₋₊ == y.G₋₊)
Base.isapprox(x::RealCorrelationFunction, y::RealCorrelationFunction; kwargs...) =
    isapprox(x.G₊₋, y.G₊₋; kwargs...) && isapprox(x.G₋₊, y.G₋₊; kwargs...)
