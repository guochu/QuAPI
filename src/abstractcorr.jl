"""
    abstract type AbstractCorrelationFunction end

Abstract supertype of all correlation-function containers used to build the
Feynman-Vernon influence functional of QuAPI.

Concrete subtypes store the discretized hybridization function Δ (the bath
correlation function) on the relevant contour.  Following the package-wide
convention, the minus sign of the influence-functional exponent is absorbed
into Δ, i.e. the discretized influence functional reads

    I[ψ̄, ψ] = exp[Σᵢⱼ ψ̄ᵢ Δᵢⱼ ψⱼ].

All subtypes must implement `size`, `index` and `branch`:

- `size(x)`:  the number of contour nodes along each branch;
- `index(x, i, j; b1, b2)`:  the matrix element Δ[b1,b2][i, j];
- `branch(x; b1, b2)`:  the branch sub-block as an `AbstractMatrix`.
"""
abstract type AbstractCorrelationFunction end


Base.size(x::AbstractCorrelationFunction, i::Int...) = error("size not defined for correlation type $(typeof(x))")

"""
    index(x::AbstractCorrelationFunction, i::Int, j::Int; b1::Symbol, b2::Symbol)

Matrix element Δᵢⱼ of the discretized hybridization function on the branch
pair `(b1, b2)`, i.e. `branch(x; b1, b2)[i, j]`.

Branch labels depend on the concrete contour:
- imaginary time:  only `:τ`;
- real time (closed Keldysh):  `:+` (forward) and `:-` (backward);
- mixed time (L-shaped Kadanoff-Baym):  `:+`, `:-` and `:τ`.

`i` and `j` are 1-based discrete contour node indices along their respective
branches.  Unknown branch labels raise an `ArgumentError`.
"""
index(x::AbstractCorrelationFunction, i::Int, j::Int) = error("index not defined for correlation type $(typeof(x))")

"""
    branch(x::AbstractCorrelationFunction; b1::Symbol, b2::Symbol)

Return the branch sub-block Δ[b1,b2] of the hybridization function as an
`AbstractMatrix` whose `(i, j)` element is `index(x, i, j; b1, b2)`.

Valid branch labels are `:τ` for imaginary time, `:+`/`:-` for real time and
all three for mixed time; unknown labels raise an `ArgumentError`.
"""
branch(x::AbstractCorrelationFunction; b1::Symbol=:τ, b2::Symbol=:τ) = error("branch not defined for correlation type $(typeof(x))")


Base.:-(A::AbstractCorrelationFunction, B::AbstractCorrelationFunction) = A + (-B)