"""
    struct RealCorrelationFunction{...} <: AbstractCorrelationFunction

Real-time hybridization function on the closed Keldysh contour.

The Keldysh contour runs forward (`+`) from t = 0 to t = t and backward (`-`)
back to t = 0; each branch carries N+1 nodes tᵢ = i·δt (i = 0..N, δt = t/N).
The four branch sub-blocks are stored as (N+1)×(N+1) matrices (usually
`CorrelationMatrix` instances):

- `G₊₊`:  time-ordered (+/+),
- `G₊₋`:  lesser-like (+/−),
- `G₋₊`:  greater-like (−/+),
- `G₋₋`:  anti-time-ordered (−/−),

and are selected by `branch(x; b1, b2)` with `b1, b2 ∈ (:+, :-)`.

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
  (`s† = s`) the exponent degenerates to the single-operator (ADT) form

      Σ_{b₁b₂} Σᵢⱼ sᵢ^{b₁} Δᵢⱼ^{b₁b₂} sⱼ^{b₂} .
"""
struct RealCorrelationFunction{A<:AbstractMatrix{ComplexF64}, B<:AbstractMatrix{ComplexF64}, C<:AbstractMatrix{ComplexF64}, D<:AbstractMatrix{ComplexF64}} <: AbstractCorrelationFunction
    G₊₊::A
    G₊₋::B
    G₋₊::C
    G₋₋::D
end

Base.transpose(x::RealCorrelationFunction) = RealCorrelationFunction(transpose(x.G₊₊), transpose(x.G₋₊), transpose(x.G₊₋), transpose(x.G₋₋))
Base.iszero(x::RealCorrelationFunction) = iszero(x.G₊₊) && iszero(x.G₊₋) && iszero(x.G₋₊) && iszero(x.G₋₋)

Base.size(x::RealCorrelationFunction, i::Int...) = size(x.G₊₊, i...)

function Base.show(io::IO, ::MIME"text/plain", A::RealCorrelationFunction)
    print(io, "Real Correlation Function [$(size(A.G₊₊, 1))]")
end
Base.:-(x::RealCorrelationFunction) = RealCorrelationFunction(-x.G₊₊, -x.G₊₋, -x.G₋₊, -x.G₋₋)
Base.:+(A::RealCorrelationFunction, B::RealCorrelationFunction) = RealCorrelationFunction(A.G₊₊ + B.G₊₊, A.G₊₋ + B.G₊₋, A.G₋₊ + B.G₋₊, A.G₋₋ + B.G₋₋)
branch(x::RealCorrelationFunction, f1::Bool, f2::Bool) = ifelse(f1, ifelse(f2, x.G₊₊, x.G₊₋), ifelse(f2, x.G₋₊, x.G₋₋))
function branch(x::RealCorrelationFunction, b1::Symbol, b2::Symbol)
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
branch(x::RealCorrelationFunction; b1::Symbol, b2::Symbol) = branch(x, b1, b2)

Base.:(==)(x::RealCorrelationFunction, y::RealCorrelationFunction) = (x.G₊₊ == y.G₊₊) && (x.G₊₋ == y.G₊₋) && (x.G₋₊ == y.G₋₊) && (x.G₋₋ == y.G₋₋)
function Base.isapprox(x::RealCorrelationFunction, y::RealCorrelationFunction; kwargs...)
    return isapprox(x.G₊₊, y.G₊₊; kwargs...) && isapprox(x.G₊₋, y.G₊₋; kwargs...) && isapprox(x.G₋₊, y.G₋₊; kwargs...) && isapprox(x.G₋₋, y.G₋₋; kwargs...)
end 
