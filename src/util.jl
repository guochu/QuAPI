_f₁(β::Real, μ::Real, ε::Float64) = 1 - fermidirac(β, μ, ε)
_f₂(β::Real, μ::Real, ε::Float64) = fermidirac(β, μ, ε)

_g₁(β::Real, μ::Real, ε::Float64) = 1 + boseeinstein(β, μ, ε)
_g₂(β::Real, μ::Real, ε::Float64) = boseeinstein(β, μ, ε)


const QuAPI_tol = 1.0e-6