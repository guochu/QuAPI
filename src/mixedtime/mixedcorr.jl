abstract type AbstractMixedCorrelationFunction <: AbstractCorrelationFunction end
"""
	struct MixedCorrelationFunction

Correlations on the L-shaped contour,
there are 9 of them due to the double integral
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


function index(A::MixedCorrelationFunction, i::Int, j::Int; b1::Symbol, b2::Symbol)
    b₁, b₂ = b1, b2
    @boundscheck begin
        (b₁ in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
        (b₂ in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))

        if b₁ == :τ
            (1 <= i <= isize(A)) || throw(BoundsError(1:isize(A), i))
        else
            (1 <= i <= rsize(A)) || throw(BoundsError(1:rsize(A), i))
        end
        if b₂ == :τ
            (1 <= j <= isize(A)) || throw(BoundsError(1:isize(A), j))
        else
            (1 <= j <= rsize(A)) || throw(BoundsError(1:rsize(A), j))
        end
    end
    # here b₁ and b₂ are branches
    if (b₁ == :+ && b₂ == :+) # G¹¹=G⁺⁺
        if (i == j)
            A.ηⱼₖ[1] + A.ηₖⱼ[1]
        elseif (i > j)
            A.ηⱼₖ[i-j+1]
        else
            A.ηₖⱼ[j-i+1]
        end
    elseif (b₁ == :+ && b₂ == :-) # G¹²=G⁺⁻
        if (i == j)
            -2*A.ηₖⱼ[1]
        elseif (i > j)
            -A.ηₖⱼ[i-j+1]'
        else
            -A.ηₖⱼ[j-i+1]
        end
    elseif (b₁ == :+ && b₂ == :τ) # G¹³
        -im*A.ζₖⱼ[i,j]
    elseif (b₁ == :- && b₂ == :+) # G²¹=G⁻⁺
        if (i == j)
            -2*A.ηⱼₖ[1]
        elseif (i > j)
            -A.ηⱼₖ[i-j+1]
        else
            -A.ηⱼₖ[j-i+1]'
        end
    elseif (b₁ == :- && b₂ == :-) # G²²=G⁻⁻
        if (i == j)
            A.ηⱼₖ[1]+A.ηₖⱼ[1]
        elseif (i > j)
            A.ηₖⱼ[i-j+1]'
        else
            A.ηⱼₖ[j-i+1]'
        end
    elseif (b₁ == :- && b₂ == :τ) # G²³
        im*A.ζₖⱼ[i,j]
    elseif (b₁ == :τ && b₂ == :+) # G³¹
        -im*A.ζⱼₖ[i,j]
    elseif (b₁ == :τ && b₂ == :-) # G³²
        im*A.ζⱼₖ[i,j]
    elseif (b₁ == :τ && b₂ == :τ) # G³³=G(τ)
        if (i == j)
            -(A.ξⱼₖ[1]+A.ξₖⱼ[1])
        elseif (i > j)
            -A.ξⱼₖ[i-j+1]
        else
            -A.ξₖⱼ[j-i+1]
        end
    end
end

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
            ηⱼₖ = conj(ηₖⱼ)
            ηⱼₖ[1] = ηₖⱼ[1]
            return CorrelationMatrix(ηⱼₖ, ηₖⱼ)
        else
            return -im*x.ζₖⱼ
        end
    elseif b1 == :-
        if b2 == :+
            ηⱼₖ = -x.ηⱼₖ
            ηₖⱼ = -conj(x.ηₖⱼ)
            ηₖⱼ[1] = ηⱼₖ[1]
            return CorrelationMatrix(ηⱼₖ, ηₖⱼ)
        elseif b2 == :-
            ηⱼₖ = conj(x.ηⱼₖ)
            ηⱼₖ[1] = x.ηⱼₖ[1]
            ηₖⱼ = conj(x.ηₖⱼ)
            ηₖⱼ[1] = x.ηₖⱼ[1]
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
