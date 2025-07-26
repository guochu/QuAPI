struct ImagCorrelationFunction{M<:AbstractMatrix{Float64}} <: AbstractCorrelationFunction
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
index(x::ImagCorrelationFunction, i::Int, j::Int) = x.data[i, j]
function branch(x::ImagCorrelationFunction; b1::Symbol=:τ, b2::Symbol=:τ)
    ((b1 == :τ) && (b2 == :τ)) || throw(ArgumentError("branch must be :τ"))
    return x.data
end 
