"""
    struct ImagCorrelationFunction{M <: AbstractMatrix} <: AbstractCorrelationFunction

Imaginary-time hybridization function Δ(τ) on the Matsubara contour
τ ∈ [0, β), discretized into N steps of δτ = β/N.

Only the single branch `:τ` exists.  The object wraps a matrix (typically a
`CorrelationMatrix` storing the Toeplitz data) of size (N, N), whose
`(i, j)` element is the hybridization Δ[i, j] between contour nodes
`τᵢ = (i-1)δτ` and `τⱼ = (j-1)δτ`.

The exponent of the discretized Feynman-Vernon influence functional built
from Δ depends on the statistics of the bath (the minus sign is absorbed
into Δ, see `AbstractCorrelationFunction`).

For a fermionic bath, with Grassmann fields `a` (annihilation) and `ā`
(conjugate), the exponent is

    Σᵢⱼ āᵢ Δᵢⱼ aⱼ .

For a bosonic bath, the system couples to the bath operators `B†` and `B`
through an operator `s` and its Hermitian conjugate,

    H_int = s ⊗ B† + s† ⊗ B ,

which is the PT (process-tensor) case of TEMPO and also covers
non-Hermitian `s` (e.g. a ladder operator).  The exponent is then

    Σᵢⱼ s†ᵢ Δᵢⱼ sⱼ ,

where `sᵢ = s(τᵢ)`.  For Hermitian `s` (`s† = s`) the two couplings
coincide and the exponent degenerates to the single-operator (ADT) form

    Σᵢⱼ sᵢ Δᵢⱼ sⱼ ;

for diagonal `s` this is exactly the exponent of the standard
imaginary-time QuAPI/TEMPO ADT framework.
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