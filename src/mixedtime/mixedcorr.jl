"""
    abstract type AbstractMixedCorrelationFunction <: AbstractCorrelationFunction end

Abstract supertype of correlation functions defined on the L-shaped
(Kadanoff-Baym) contour, which combines an imaginary-time branch `:τ` of
`Nτ` nodes and a real-time Keldysh part with branches `:+`/`:-` of `Nt+1`
nodes each.  See `MixedCorrelationFunction` for the concrete implementation.
"""
abstract type AbstractMixedCorrelationFunction <: AbstractCorrelationFunction end

"""
	struct MixedCorrelationFunction

Correlations on the L-shaped (Kadanoff-Baym) contour.

The contour consists of an imaginary-time branch τ ∈ [0, β) discretized into
`Nτ` steps of δτ = β/Nτ, followed by the closed Keldysh part t ∈ [0, t] with
`Nt+1` nodes per branch (δt = t/Nt).  Owing to the double contour integral the
hybridization is stored as nine blocks, grouped in three pairs:

- real-time blocks `ηⱼₖ`/`ηₖⱼ` (length `Nt+1`), entering the (++), (+-), (-+),
  (--) Keldysh combinations;
- imaginary-time blocks `ξⱼₖ`/`ξₖⱼ` (length `Nτ`), entering the (τ,τ) block;
- mixed-time blocks `ζⱼₖ` (size `Nτ × (Nt+1)`) and `ζₖⱼ` (size `(Nt+1) × Nτ`),
  entering the cross (+/-,τ) combinations.

`branch(x; b1, b2)` with `b1, b2 ∈ (:+, :-, :τ)` returns the corresponding
sub-block as a matrix (the nine combinations G¹¹ ... G³³); `index` gives the
individual matrix element.

The exponent of the discretized Feynman-Vernon influence functional built
from these blocks depends on the statistics of the bath (the minus sign is
absorbed into Δ, see `AbstractCorrelationFunction`).  With branch labels
`b₁, b₂ ∈ (:+, :-, :τ)` and `i, j` running over the nodes of the respective
branch (`Nt+1` for `:+`/`:-`, `Nτ` for `:τ`),

- fermionic bath, with Grassmann fields `a` (annihilation) and `ā`
  (conjugate):

      Σ_{b₁b₂} Σᵢⱼ āᵢ^{b₁} Δᵢⱼ^{b₁b₂} aⱼ^{b₂} ;

- bosonic bath with the coupling `s ⊗ B† + s† ⊗ B` (the PT case of TEMPO,
  also covering non-Hermitian `s`, e.g. a ladder operator):

      Σ_{b₁b₂} Σᵢⱼ s†ᵢ^{b₁} Δᵢⱼ^{b₁b₂} sⱼ^{b₂} ,

  where `s^{b}` is the system field on branch `b`.  For Hermitian `s`
  (`s† = s`) the exponent degenerates to the single-operator (ADT) form

      Σ_{b₁b₂} Σᵢⱼ sᵢ^{b₁} Δᵢⱼ^{b₁b₂} sⱼ^{b₂} .
"""
struct MixedCorrelationFunction <: AbstractMixedCorrelationFunction
    # real time
    ηⱼₖ::Vector{ComplexF64}
    ηₖⱼ::Vector{ComplexF64}
    # imag time
    ξⱼₖ::Vector{ComplexF64}
    ξₖⱼ::Vector{ComplexF64}
    # mix time
    ζⱼₖ::Matrix{ComplexF64}
    ζₖⱼ::Matrix{ComplexF64}
end

function Base.show(io::IO, ::MIME"text/plain", A::AbstractMixedCorrelationFunction)
    print(io, "Mixed Correlation Function [$(isize(A))+$(rsize(A))]")
end

Base.transpose(x::MixedCorrelationFunction) = MixedCorrelationFunction(x.ηₖⱼ, x.ηⱼₖ, x.ξₖⱼ, x.ξⱼₖ, transpose(x.ζₖⱼ), transpose(x.ζⱼₖ))
Base.iszero(x::MixedCorrelationFunction) = iszero(x.ηⱼₖ) && iszero(x.ηₖⱼ) && iszero(x.ξⱼₖ) && iszero(x.ξₖⱼ) && iszero(x.ζⱼₖ) && iszero(x.ζₖⱼ)
Base.:-(x::MixedCorrelationFunction) = MixedCorrelationFunction(-x.ηⱼₖ, -x.ηₖⱼ, -x.ξⱼₖ, -x.ξₖⱼ, -x.ζⱼₖ, -x.ζₖⱼ)
Base.:+(A::MixedCorrelationFunction, B::MixedCorrelationFunction) = MixedCorrelationFunction(A.ηⱼₖ + B.ηⱼₖ, A.ηₖⱼ + B.ηₖⱼ, A.ξⱼₖ + B.ξⱼₖ, A.ξₖⱼ + B.ξₖⱼ, A.ζⱼₖ + B.ζⱼₖ, A.ζₖⱼ + B.ζₖⱼ)
# branch(x::MixedCorrelationFunction, f1::Symbol, f2::Symbol) = ifelse(f1, ifelse(f2, x.G₊₊, x.G₊₋), ifelse(f2, x.G₋₊, x.G₋₋))

isize(x::MixedCorrelationFunction) = length(x.ξⱼₖ)
rsize(x::MixedCorrelationFunction) = length(x.ηⱼₖ)


function branch(x::MixedCorrelationFunction; b1::Symbol, b2::Symbol)
    @boundscheck begin
        (b1 in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
        (b2 in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
    end
    if b1 == :+
        if b2 == :+
            return CorrelationMatrix(x.ηⱼₖ, x.ηₖⱼ)
        elseif b2 == :-
            ηₖⱼ = -(x.ηₖⱼ)
            ηⱼₖ = conj(ηₖⱼ)          # 等时也为共轭对：{b·S₋, b·S₊}，和 = b·f²
            return CorrelationMatrix(ηⱼₖ, ηₖⱼ)
        else
            return -im*x.ζₖⱼ
        end
    elseif b1 == :-
        if b2 == :+
            ηⱼₖ = -x.ηⱼₖ
            ηₖⱼ = -conj(x.ηⱼₖ)       # 等时也为共轭对：{bb·S₋, bb·S₊}，和 = bb·f²
            return CorrelationMatrix(ηⱼₖ, ηₖⱼ)
        elseif b2 == :-
            ηⱼₖ = conj(x.ηₖⱼ)            # 等时为共轭对：{−b·S₋, −bb·S₊}（与 (++) 的
            ηₖⱼ = conj(x.ηⱼₖ)            # {−bb·S₋, −b·S₊} 相配，ket=bra 路径权重为实）
            return CorrelationMatrix(ηⱼₖ, ηₖⱼ)
        else
            return im*x.ζₖⱼ
        end
    else
        if b2 == :+
            return -im*x.ζⱼₖ
        elseif b2 == :-
            return im*x.ζⱼₖ
        else
            ηⱼₖ = -x.ξⱼₖ
            ηₖⱼ = -x.ξₖⱼ
            return CorrelationMatrix(ηⱼₖ, ηₖⱼ)
        end
    end
end

function Base.:(==)(x::MixedCorrelationFunction, y::MixedCorrelationFunction)
    return (x.ηⱼₖ == y.ηⱼₖ) && (x.ηₖⱼ == y.ηₖⱼ) && (x.ξⱼₖ == y.ξⱼₖ) && (x.ξₖⱼ == y.ξₖⱼ) && (x.ζⱼₖ == y.ζⱼₖ) && (x.ζₖⱼ == y.ζₖⱼ)
end
function Base.isapprox(x::MixedCorrelationFunction, y::MixedCorrelationFunction; kwargs...)
    return isapprox(x.ηⱼₖ, y.ηⱼₖ; kwargs...) && isapprox(x.ηₖⱼ, y.ηₖⱼ; kwargs...) && isapprox(x.ξⱼₖ, y.ξⱼₖ; kwargs...) && isapprox(
                    x.ξₖⱼ, y.ξₖⱼ; kwargs...) && isapprox(x.ζⱼₖ, y.ζⱼₖ; kwargs...) && isapprox(x.ζₖⱼ, y.ζₖⱼ; kwargs...)
end 
