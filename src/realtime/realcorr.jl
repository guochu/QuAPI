"""
    abstract type AbstractRealCorrelationFunction <: AbstractCorrelationFunction

Abstract supertype of the real-time hybridization representations: the compact
two-vector `RealCorrelationFunction` and the full four-block
`GenericRealCorrelationFunction`.  `branch`, `branches`, `contour_ordering`,
`transpose`, arithmetic and comparisons are defined uniformly on this level.

Both representations follow ONE unified branch convention (the same as
`MixedCorrelationFunction` on the Kadanoff contour; see `changes.md`): with
`a = ηⱼₖ` (forward lags) and `b = ηₖⱼ` (backward lags), the four Keldysh
blocks are

    (+,+): (a, b)                (+,−): (−c, −b)
    (−,+): (−a, −d)              (−,−): (c, d)

with c = ηⱼₖ₂, d = ηₖⱼ₂ (for real-valued spectra c = conj(b), d = conj(a)),
i.e. the CROSS branches carry the Keldysh relative minus signs inside the
blocks and the derivation uses only negation.  This is statistics
independent — fermionic and bosonic hybridizations store the same data and
derive the blocks identically.

`branch(x; b1, b2)` with `b1, b2 ∈ (:+, :-)` selects the block; `i, j` run
over the `N+1` nodes tᵢ = i·δt of each branch.
"""
abstract type AbstractRealCorrelationFunction <: AbstractCorrelationFunction end

branches(::AbstractRealCorrelationFunction) = (:+, :-)

"""
    struct GenericRealCorrelationFunction{...} <: AbstractRealCorrelationFunction

Full four-block real-time hybridization function on the closed Keldysh
contour, storing the four branch sub-blocks `G₊₊`, `G₊₋`, `G₋₊`, `G₋₋`
(usually `CorrelationMatrix` instances) explicitly.  The blocks follow the
unified branch convention documented on `AbstractRealCorrelationFunction`
(cross branches carry the Keldysh minus signs).  This representation carries
the four independent Toeplitz vectors of a general (possibly complex-valued)
spectrum and is the form to use whenever the cross-branch vectors cannot be
derived by conjugation of the primary pair (e.g. the BCS anomalous components
with a complex gap phase: there `conj` would conjugate the spectral weight
itself, which breaks the derivation by a factor `e^{2iφ}`).

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
struct GenericRealCorrelationFunction{A<:AbstractMatrix{ComplexF64}, B<:AbstractMatrix{ComplexF64}, C<:AbstractMatrix{ComplexF64}, D<:AbstractMatrix{ComplexF64}} <: AbstractRealCorrelationFunction
    G₊₊::A
    G₊₋::B
    G₋₊::C
    G₋₋::D
end

Base.transpose(x::GenericRealCorrelationFunction) = GenericRealCorrelationFunction(transpose(x.G₊₊), transpose(x.G₋₊), transpose(x.G₊₋), transpose(x.G₋₋))
Base.iszero(x::GenericRealCorrelationFunction) = iszero(x.G₊₊) && iszero(x.G₊₋) && iszero(x.G₋₊) && iszero(x.G₋₋)

Base.size(x::GenericRealCorrelationFunction, i::Int...) = size(x.G₊₊, i...)

function Base.show(io::IO, ::MIME"text/plain", A::GenericRealCorrelationFunction)
    print(io, "Real Correlation Function [$(size(A.G₊₊, 1))] (generic)")
end
Base.:-(x::GenericRealCorrelationFunction) = GenericRealCorrelationFunction(-x.G₊₊, -x.G₊₋, -x.G₋₊, -x.G₋₋)
Base.:+(A::GenericRealCorrelationFunction, B::GenericRealCorrelationFunction) = GenericRealCorrelationFunction(A.G₊₊ + B.G₊₊, A.G₊₋ + B.G₊₋, A.G₋₊ + B.G₋₊, A.G₋₋ + B.G₋₋)
branch(x::GenericRealCorrelationFunction, f1::Bool, f2::Bool) = ifelse(f1, ifelse(f2, x.G₊₊, x.G₊₋), ifelse(f2, x.G₋₊, x.G₋₋))
function branch(x::GenericRealCorrelationFunction, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    if b1 == :+
        if b2 == :+
            return x.G₊₊
        else
            return x.G₊₋
        end
    else
        if b2 == :+
            return x.G₋₊
        else
            return x.G₋₋
        end
    end
end
branch(x::GenericRealCorrelationFunction; b1::Symbol, b2::Symbol) = branch(x, b1, b2)

# =============================================================================
# Compact two-vector representation
# =============================================================================
"""
    struct RealCorrelationFunction <: AbstractRealCorrelationFunction

Compact real-time hybridization function storing the four independent
Toeplitz vectors as plain `Vector{ComplexF64}` (mirroring the storage of
`MixedCorrelationFunction` on the Kadanoff contour):

- `ηⱼₖ`  — (+,+) forward lags,  thermal weight g₁,
- `ηₖⱼ`  — (+,+) backward lags, thermal weight g₂,
- `ηⱼₖ₂` — (+,−) forward lags,  thermal weight g₂,
- `ηₖⱼ₂` — (−,+) backward lags, thermal weight g₁.

The four branch blocks are derived in `branch` following the unified
(statistics-independent) convention — using only negation, **no conjugation**:

    (+,+): (a, b)                (+,−): (−c, −b)
    (−,+): (−a, −d)              (−,−): (c, d)

with `a = ηⱼₖ`, `b = ηₖⱼ`, `c = ηⱼₖ₂`, `d = ηₖⱼ₂`.  For real-valued spectra
`c = conj(b)` and `d = conj(a)` hold exactly (cell-kernel conjugation
symmetry with real thermal weights), so the two-vector constructor derives
them by conjugation; for complex-valued spectra — the BCS anomalous coherence
weights — conjugation would conjugate the spectral weight itself (breaking
the blocks by a factor `e^{2iφ}`, the very same trap as in `bcs_Δm`), and
the four-vector constructor must be used with `ηⱼₖ₂`/`ηₖⱼ₂` computed by
their own quadrature (`fermionic_Δt` does this unconditionally).

`fermionic_Δt` / `bosonic_Δt` return this compact form.
"""
struct RealCorrelationFunction <: AbstractRealCorrelationFunction
    ηⱼₖ::Vector{ComplexF64}
    ηₖⱼ::Vector{ComplexF64}
    ηⱼₖ₂::Vector{ComplexF64}
    ηₖⱼ₂::Vector{ComplexF64}
end

"""
    RealCorrelationFunction(ηⱼₖ, ηₖⱼ) -> RealCorrelationFunction

Two-vector constructor for **real-valued spectra**: derives the secondary
pair by conjugation (`ηⱼₖ₂ = conj(ηₖⱼ)`, `ηₖⱼ₂ = conj(ηⱼₖ)`), which is exact
when the spectral weights are real.
"""
function RealCorrelationFunction(ηⱼₖ::AbstractVector{<:Number}, ηₖⱼ::AbstractVector{<:Number})
    a = collect(ComplexF64, ηⱼₖ)
    b = collect(ComplexF64, ηₖⱼ)
    return RealCorrelationFunction(a, b, conj.(b), conj.(a))
end

"""
    RealCorrelationFunction(ηⱼₖ, ηₖⱼ, ηⱼₖ₂, ηₖⱼ₂) -> RealCorrelationFunction

Four-vector constructor (complex-spectrum safe): takes the secondary pair
explicitly, computed with the spectral weight un-conjugated.
"""
RealCorrelationFunction(ηⱼₖ::AbstractVector{<:Number}, ηₖⱼ::AbstractVector{<:Number},
                        ηⱼₖ₂::AbstractVector{<:Number}, ηₖⱼ₂::AbstractVector{<:Number}) =
    RealCorrelationFunction(collect(ComplexF64, ηⱼₖ), collect(ComplexF64, ηₖⱼ),
                            collect(ComplexF64, ηⱼₖ₂), collect(ComplexF64, ηₖⱼ₂))

Base.transpose(x::RealCorrelationFunction) = RealCorrelationFunction(x.ηₖⱼ, x.ηⱼₖ, x.ηₖⱼ₂, x.ηⱼₖ₂)
Base.iszero(x::RealCorrelationFunction) = iszero(x.ηⱼₖ) && iszero(x.ηₖⱼ) && iszero(x.ηⱼₖ₂) && iszero(x.ηₖⱼ₂)

Base.size(x::RealCorrelationFunction) = (length(x.ηⱼₖ), length(x.ηⱼₖ))
Base.size(x::RealCorrelationFunction, i::Int) = size(x)[i]

function Base.show(io::IO, ::MIME"text/plain", A::RealCorrelationFunction)
    print(io, "Real Correlation Function [$(length(A.ηⱼₖ))] (compact)")
end
Base.:-(x::RealCorrelationFunction) = RealCorrelationFunction(-x.ηⱼₖ, -x.ηₖⱼ, -x.ηⱼₖ₂, -x.ηₖⱼ₂)
Base.:+(A::RealCorrelationFunction, B::RealCorrelationFunction) =
    RealCorrelationFunction(A.ηⱼₖ + B.ηⱼₖ, A.ηₖⱼ + B.ηₖⱼ, A.ηⱼₖ₂ + B.ηⱼₖ₂, A.ηₖⱼ₂ + B.ηₖⱼ₂)

# unified branch derivation (statistics independent; cross branches carry the
# Keldysh minus signs — same convention as MixedCorrelationFunction).  Only
# negation is used, so the derivation is exact for complex-valued spectra as
# long as ηⱼₖ₂/ηₖⱼ₂ were supplied un-conjugated.
function branch(x::RealCorrelationFunction, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    a, b, c, d = x.ηⱼₖ, x.ηₖⱼ, x.ηⱼₖ₂, x.ηₖⱼ₂
    if b1 == :+
        if b2 == :+
            return CorrelationMatrix(a, b)
        else
            return CorrelationMatrix(.-c, .-b)
        end
    else
        if b2 == :+
            return CorrelationMatrix(.-a, .-d)
        else
            return CorrelationMatrix(c, d)
        end
    end
end

branch(x::RealCorrelationFunction; b1::Symbol, b2::Symbol) = branch(x, b1, b2)
function branch(x::RealCorrelationFunction, f1::Bool, f2::Bool)
    return branch(x, f1 ? :+ : :-, f2 ? :+ : :-)
end

# uniform block-wise equality / comparison across real representations
function Base.:(==)(x::AbstractRealCorrelationFunction, y::AbstractRealCorrelationFunction)
    return all(branch(x, b1, b2) == branch(y, b1, b2) for b1 in (:+, :-) for b2 in (:+, :-))
end
function Base.isapprox(x::AbstractRealCorrelationFunction, y::AbstractRealCorrelationFunction; kwargs...)
    return all(isapprox(branch(x, b1, b2), branch(y, b1, b2); kwargs...)
               for b1 in (:+, :-) for b2 in (:+, :-))
end

"""
    RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋) -> GenericRealCorrelationFunction

Legacy four-block constructor: builds a `GenericRealCorrelationFunction`
from the four branch matrices (unified convention — cross branches carry the
Keldysh minus signs).
"""
RealCorrelationFunction(G₊₊::AbstractMatrix{ComplexF64}, G₊₋::AbstractMatrix{ComplexF64}, G₋₊::AbstractMatrix{ComplexF64}, G₋₋::AbstractMatrix{ComplexF64}) =
    GenericRealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
