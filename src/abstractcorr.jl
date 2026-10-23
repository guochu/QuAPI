"""
    abstract type AbstractCorrelationFunction end

Abstract supertype of all correlation-function containers used to build the
Feynman-Vernon influence functional of QuAPI.

Concrete subtypes store the discretized hybridization function Δ (the bath
correlation function) on the relevant contour.  Following the package-wide
convention, the minus sign of the influence-functional exponent is absorbed
into Δ, i.e. the discretized influence functional reads

    I[ψ̄, ψ] = exp[Σᵢⱼ ψ̄ᵢ Δᵢⱼ ψⱼ].

All subtypes must implement `size` and `branch`:

- `size(x)`:  the number of contour nodes along each branch;
- `branch(x; b1, b2)`:  the branch sub-block as an `AbstractMatrix`.
"""
abstract type AbstractCorrelationFunction end


Base.size(x::AbstractCorrelationFunction, i::Int...) = error("size not defined for correlation type $(typeof(x))")

"""
    branch(x::AbstractCorrelationFunction; b1::Symbol, b2::Symbol)

Return the branch sub-block Δ[b1,b2] of the hybridization function as an
`AbstractMatrix`; its `(i, j)` element is the matrix element Δᵢⱼ.

Valid branch labels are `:τ` for imaginary time, `:+`/`:-` for real time and
all three for mixed time; unknown labels raise an `ArgumentError`.
"""
branch(x::AbstractCorrelationFunction; b1::Symbol=:τ, b2::Symbol=:τ) = error("branch not defined for correlation type $(typeof(x))")


Base.:-(A::AbstractCorrelationFunction, B::AbstractCorrelationFunction) = A + (-B)