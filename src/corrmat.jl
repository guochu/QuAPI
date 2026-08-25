"""
    struct CorrelationMatrix{T <: Number} <: AbstractMatrix{T}

Storage of a time-translationally invariant correlation matrix.

The correlation function of a stationary environment depends only on the time
difference |i - j|, so only the two Toeplitz vectors are stored:
- `ηⱼₖ[d]` = Δ[i, j] for i ≥ j, with `d = i - j + 1` (equal time at `d = 1`);
- `ηₖⱼ[d]` = Δ[i, j] for i < j, with `d = j - i + 1`.

The two equal-time entries `ηⱼₖ[1]` and `ηₖⱼ[1]` are stored separately on
purpose, because they correspond to the two orderings of the environment
operators in the influence-functional exponent at equal time.  For a
discretized (Keldysh) contour the exponent contains a term
Σᵢⱼ sᵢᵘ Δᵢⱼᵘᵛ sⱼᵛ (spin indices suppressed), with `sᵢᵘ` the system
coupling operator on branch `u` at node `i` and `Δᵢⱼᵘᵛ` the corresponding
bath correlation.  For the equal-time node `i = j` the product
`sᵢᵘ Δᵢᵢᵘᵛ sᵢᵛ` contains two environment operators at the same instant, and
whether the two orderings can be combined depends on their commutation:
- *Diagonal coupling* (`sᵢᵘ = sᵢᵛ = ŝ`, a single Hermitian operator such as
  σ̂z): the operator commutes with itself, the two orderings give identical
  contributions that add directly, `Δ[i, i] = ηⱼₖ[1] + ηₖⱼ[1]`.
- *Off-diagonal coupling* (e.g. `ŝ = σ̂₊, σ̂₋` or a vector coupling `ŝ·B̂`):
  the two operators at equal time do not commute, the two orderings give
  different contributions to the exponent, and `ηⱼₖ[1]` and `ηₖⱼ[1]` must be
  tracked separately.

See Appendix C of `docs/GuoChen2026a.pdf` for the explicit treatment of the
equal-time correlation matrix.  Accordingly the element `A[i, j]` reads
`ηⱼₖ[i-j+1]` (i > j), `ηₖⱼ[j-i+1]` (i < j) and `ηⱼₖ[1] + ηₖⱼ[1]` (i = j).

The type is a thin wrapper around the two vectors: arithmetic (`+`, `-`,
`transpose`) and `iszero` act on the stored data directly, while `==` and
`isapprox` require elementwise agreement of both Toeplitz vectors, treating
the two equal-time parts as physically distinct.
"""
struct CorrelationMatrix{T <: Number} <: AbstractMatrix{T}
    ηⱼₖ::Vector{T}  # j >= k
    ηₖⱼ::Vector{T}    
end

Base.size(x::CorrelationMatrix, i::Int) = ifelse(i <= 2, length(x.ηⱼₖ), 1)
Base.size(x::CorrelationMatrix) = (length(x.ηⱼₖ), length(x.ηⱼₖ))
Base.iszero(x::CorrelationMatrix) = iszero(x.ηⱼₖ) && iszero(x.ηₖⱼ)
Base.:(==)(x::CorrelationMatrix, y::CorrelationMatrix) = x.ηⱼₖ == y.ηⱼₖ && x.ηₖⱼ == y.ηₖⱼ

function Base.isapprox(x::CorrelationMatrix, y::CorrelationMatrix; kwargs...)
    return isapprox(x.ηⱼₖ, y.ηⱼₖ; kwargs...) && isapprox(x.ηₖⱼ, y.ηₖⱼ; kwargs...)
end

Base.similar(a::CorrelationMatrix, dims::Tuple{Int, Int}) = Matrix{eltype(a)}(undef, dims)
Base.similar(a::CorrelationMatrix, T::Type, dims::Tuple{Int, Int}) = Matrix{T}(undef, dims)

function Base.getindex(A::CorrelationMatrix, i::Int, j::Int)
    L = size(A, 1)
    ((1 <= i <= L) && (1 <= j <= L)) || throw(BoundsError())
    if (i > j)
        A.ηⱼₖ[i-j+1]
    elseif (i < j)
        A.ηₖⱼ[j-i+1]
    else
        A.ηⱼₖ[1] + A.ηₖⱼ[1]
    end
end
# function Base.getindex(A::CorrelationMatrix, i::Int)
#     if i == 0
#         return A.ηⱼₖ[1] + A.ηₖⱼ[1]
#     elseif i > 0
#         return A.ηⱼₖ[i+1]
#     else
#         return A.ηₖⱼ[i-1]
#     end
# end

function Base.:+(x::CorrelationMatrix, y::CorrelationMatrix)
    (size(x) == size(y)) || throw(DimensionMismatch("correlation matrix size mismatch"))
    return CorrelationMatrix(x.ηⱼₖ + y.ηⱼₖ, x.ηₖⱼ + y.ηₖⱼ)
end
Base.:(-)(x::CorrelationMatrix) = CorrelationMatrix(-x.ηⱼₖ, -x.ηₖⱼ)
Base.:-(x::CorrelationMatrix, y::CorrelationMatrix) = x + (-y)
Base.transpose(x::CorrelationMatrix) = CorrelationMatrix(x.ηₖⱼ, x.ηⱼₖ)