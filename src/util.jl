_g₁(β::Real, μ::Real, ε::Float64) = 1 - fermidirac(β, μ, ε)
_g₂(β::Real, μ::Real, ε::Float64) = fermidirac(β, μ, ε)

const tol = 1.0e-6