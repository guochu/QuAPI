"""
    struct ImagCorrelationFunction{M <: AbstractMatrix} <: AbstractCorrelationFunction

Imaginary-time hybridization function Δ(τ) on the Matsubara contour
τ ∈ [0, β), discretized into N steps of δτ = β/N.

Only the single branch `:τ` exists.  The object wraps a matrix (typically a
`CorrelationMatrix` storing the Toeplitz data) of size (N, N), whose
`(i, j)` element is the hybridization Δ[i, j] between contour nodes
`τᵢ = (i-1)δτ` and `τⱼ = (j-1)δτ`, entering the discretized influence
functional as `exp[Σᵢⱼ āᵢ Δᵢⱼ aⱼ]` (minus sign absorbed, see
`AbstractCorrelationFunction`).
"""
struct ImagCorrelationFunction{M<:AbstractMatrix} <: AbstractCorrelationFunction
    data::M
end

Base.size(x::ImagCorrelationFunction, i::Int...) = size(x.data, i...)

function Base.show(io::IO, ::MIME"text/plain", A::ImagCorrelationFunction)
    print(io, "Imaginary Correlation Function [$(size(A.data, 1))]")
end
Base.:(-)(x::ImagCorrelationFunction) = ImagCorrelationFunction(-x.data)
Base.:+(A::ImagCorrelationFunction, B::ImagCorrelationFunction) = ImagCorrelationFunction(A.data + B.data)
Base.:(==)(x::ImagCorrelationFunction, y::ImagCorrelationFunction) = x.data == y.data
Base.transpose(x::ImagCorrelationFunction) = ImagCorrelationFunction(transpose(x.data))
Base.iszero(x::ImagCorrelationFunction) = iszero(x.data)

Base.isapprox(x::ImagCorrelationFunction, y::ImagCorrelationFunction; kwargs...) = isapprox(x.data, y.data; kwargs...)
function index(x::ImagCorrelationFunction, i::Int, j::Int; b1::Symbol=:τ, b2::Symbol=:τ)
    ((b1 == :τ) && (b2 == :τ)) || throw(ArgumentError("branch must be :τ"))
    x.data[i, j]
end 
function branch(x::ImagCorrelationFunction, b1::Symbol, b2::Symbol)
    ((b1 == :τ) && (b2 == :τ)) || throw(ArgumentError("branch must be :τ"))
    return x.data
end 
branch(x::ImagCorrelationFunction; b1::Symbol=:τ, b2::Symbol=:τ) = branch(x, b1, b2)