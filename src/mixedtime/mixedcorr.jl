"""
    abstract type AbstractMixedCorrelationFunction <: AbstractCorrelationFunction end

Abstract supertype of correlation functions defined on the L-shaped
(Kadanoff-Baym) contour, which combines an imaginary-time branch `:τ` of
`Nτ` nodes and a real-time Keldysh part with branches `:+`/`:-` of `Nt+1`
nodes each.  See `MixedCorrelationFunction` for the concrete implementation.
"""
abstract type AbstractMixedCorrelationFunction <: AbstractCorrelationFunction end

"""
    struct MixedCorrelationFunction{A} <: AbstractMixedCorrelationFunction

Correlations on the L-shaped (Kadanoff-Baym) contour, with the same minimal
block storage as `RealCorrelationFunction`:

- `G₊₋::A`, `G₋₊::A` (parametric matrix type `A`, both fields identical
  type) — the two CROSS real-time Keldysh blocks (lesser/greater-type
  correlators, the physically fundamental pair).  The diagonal blocks are
  recovered by `_G₊₊(G₊₋, G₋₊)` / `_G₋₋(G₊₋, G₋₊)` (pure sign
  combinations, exact for arbitrarily complex-valued spectra):

      (+,−) = G₊₋                      (−,+) = G₋₊
      (+,+) = CM(−G₋₊.ηⱼₖ, −G₊₋.ηₖⱼ)   (−,−) = CM(−G₊₋.ηⱼₖ, −G₋₊.ηₖⱼ)

- `G₊τ::Matrix{ComplexF64}` — the (+,τ) cross block, `(Nt+1)×Nτ`;
- `Gτ₊::Matrix{ComplexF64}` — the (τ,+) cross block, `Nτ×(Nt+1)`;
  the two cross blocks are stored as independent inputs (the ζ blocks of the
  construction; they are independent in general — the fermionic
  antisymmetry ζⱼₖ = −ζₖⱼᵀ holds only at μ = 0).  The (−,τ)/(τ,−)
  variants are their negations.
- `Gττ::CorrelationMatrix` — the (τ,τ) block, `Nτ×Nτ` Toeplitz.

`branch(x; b1, b2)` with `b1, b2 ∈ (:+, :-, :τ)` returns the corresponding
sub-block; the nine combinations share the IF-exponent semantics of
`AbstractCorrelationFunction`.  Construct as
`MixedCorrelationFunction(G₊₋, G₋₊, G₊τ, Gτ₊, Gττ)` (auto-generated
constructor; the two real-time cross blocks must share one matrix type).
"""
struct MixedCorrelationFunction{A<:AbstractMatrix{ComplexF64}} <: AbstractMixedCorrelationFunction
    G₊₋::A
    G₋₊::A
    G₊τ::Matrix{ComplexF64}
    Gτ₊::Matrix{ComplexF64}
    Gττ::CorrelationMatrix{ComplexF64}
end


function Base.show(io::IO, ::MIME"text/plain", A::AbstractMixedCorrelationFunction)
    print(io, "Mixed Correlation Function [$(isize(A))+$(rsize(A))]")
end

Base.transpose(x::MixedCorrelationFunction) =
    MixedCorrelationFunction(transpose(x.G₋₊), transpose(x.G₊₋),
                             Matrix(transpose(x.Gτ₊)), Matrix(transpose(x.G₊τ)),
                             transpose(x.Gττ))
Base.iszero(x::MixedCorrelationFunction) = iszero(x.G₊₋) && iszero(x.G₋₊) &&
    iszero(x.G₊τ) && iszero(x.Gτ₊) && iszero(x.Gττ)
Base.:-(x::MixedCorrelationFunction) =
    MixedCorrelationFunction(-x.G₊₋, -x.G₋₊, -x.G₊τ, -x.Gτ₊, -x.Gττ)
Base.:+(A::MixedCorrelationFunction, B::MixedCorrelationFunction) =
    MixedCorrelationFunction(A.G₊₋ + B.G₊₋, A.G₋₊ + B.G₋₊,
                             A.G₊τ + B.G₊τ, A.Gτ₊ + B.Gτ₊, A.Gττ + B.Gττ)

isize(x::MixedCorrelationFunction) = size(x.Gττ, 1)
rsize(x::MixedCorrelationFunction) = size(x.G₊₋, 1)

function branch(x::MixedCorrelationFunction{A}; b1::Symbol, b2::Symbol) where {A<:CorrelationMatrix{ComplexF64}}
    @boundscheck begin
        (b1 in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
        (b2 in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
    end
    if b1 == :+
        if b2 == :+
            return _G₊₊(x.G₊₋, x.G₋₊)
        elseif b2 == :-
            return x.G₊₋
        else
            return x.G₊τ
        end
    elseif b1 == :-
        if b2 == :+
            return x.G₋₊
        elseif b2 == :-
            return _G₋₋(x.G₊₋, x.G₋₊)
        else
            return .-(x.G₊τ)
        end
    else
        if b2 == :+
            return x.Gτ₊
        elseif b2 == :-
            return .-(x.Gτ₊)
        else
            return x.Gττ
        end
    end
end

# general (non-Toeplitz) real-time storage: real-time cross blocks derived by
# the sector rule (triangles from the partner blocks, diagonal from the
# same-branch partner — same convention as RealCorrelationFunction)
function branch(x::MixedCorrelationFunction; b1::Symbol, b2::Symbol)
    @boundscheck begin
        (b1 in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
        (b2 in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
    end
    if b1 == :τ
        if b2 == :+
            return x.Gτ₊
        elseif b2 == :-
            return .-(x.Gτ₊)
        else
            return x.Gττ
        end
    elseif b2 == :τ
        return b1 == :+ ? x.G₊τ : .-(x.G₊τ)
    end
    if b1 == :+
        return b2 == :+ ? _G₊₊(x.G₊₋, x.G₋₊) : x.G₊₋
    else
        return b2 == :+ ? x.G₋₊ : _G₋₋(x.G₊₋, x.G₋₊)
    end
end

function Base.:(==)(x::MixedCorrelationFunction, y::MixedCorrelationFunction)
    return (x.G₊₋ == y.G₊₋) && (x.G₋₊ == y.G₋₊) && (x.G₊τ == y.G₊τ) &&
           (x.Gτ₊ == y.Gτ₊) && (x.Gττ == y.Gττ)
end
function Base.isapprox(x::MixedCorrelationFunction, y::MixedCorrelationFunction; kwargs...)
    return isapprox(x.G₊₋, y.G₊₋; kwargs...) && isapprox(x.G₋₊, y.G₋₊; kwargs...) &&
           isapprox(x.G₊τ, y.G₊τ; kwargs...) && isapprox(x.Gτ₊, y.Gτ₊; kwargs...) &&
           isapprox(x.Gττ, y.Gττ; kwargs...)
end