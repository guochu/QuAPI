struct RealCorrelationFunction{A<:AbstractMatrix{ComplexF64}, B<:AbstractMatrix{ComplexF64}, C<:AbstractMatrix{ComplexF64}, D<:AbstractMatrix{ComplexF64}} <: AbstractCorrelationFunction
    G₊₊::A
    G₊₋::B
    G₋₊::C
    G₋₋::D
end

Base.transpose(x::RealCorrelationFunction) = RealCorrelationFunction(transpose(x.G₊₊), transpose(x.G₋₊), transpose(x.G₊₋), transpose(x.G₋₋))


Base.size(x::RealCorrelationFunction, i::Int...) = size(x.G₊₊, i...)

function Base.show(io::IO, ::MIME"text/plain", A::RealCorrelationFunction)
    print(io, "Real Correlation Function [$(size(A.G₊₊, 1))]")
end
index(x::RealCorrelationFunction, i::Int, j::Int; b1::Symbol, b2::Symbol) = branch(x, b1, b2)[i, j]

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
