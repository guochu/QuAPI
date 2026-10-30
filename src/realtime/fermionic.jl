"""
    Δt(bath::AbstractFermionicNormalBath; N::Int, t::Real)
    Δt(bath::AbstractBosonicNormalBath; N::Int, t::Real)
    Δt(bath::AbstractBCSBath; N::Int, t::Real)

Discretized real-time hybridization function on the closed Keldysh contour,
covering the total evolution time t with `N` steps δt = t/N.  Returns a
`RealCorrelationFunction` with the four branch sub-blocks G₊₊, G₊₋, G₋₊,
G₋₋, or a `BCSCorrelationFunction` for a superconducting bath.

Convenience wrapper that dispatches to `fermionic_Δt` (fermionic baths),
`bosonic_Δt` (bosonic baths) and `bcs_Δt` (BCS baths), reading `β` and `μ`
from the bath object.  The minus sign of the influence-functional exponent is
absorbed into Δ.
"""
Δt(bath::AbstractFermionicNormalBath; N::Int, t::Real) = fermionic_Δt(bath.spectrum, β=bath.β, N=N, t=t, μ=bath.μ)
fermionic_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, μ::Real=0) = fermionic_Δt(f, β, t, N, μ)
fermionic_Δt(f::AbstractBoundedFunction, disperse::Function; β::Real, N::Int, t::Real, μ::Real=0) = fermionic_Δt(f, disperse, β, t, N, μ)

"""
    fermionic_Δt(bath::AbstractFermionicNormalBath; N::Int, t::Real)
    fermionic_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, μ::Real=0)
    fermionic_Δt(f::AbstractBoundedFunction, disperse::Function; β::Real, N::Int, t::Real, μ::Real=0)

Fermionic hybridization function Δᵢⱼᵘᵛ on the real-time (closed Keldysh)
contour.

The Feynman-Vernon influence functional has the form

    I[ā, a] = e^{Σᵤᵥ Σᵢⱼ āᵢᵘ Δᵢⱼᵘᵛ aⱼᵛ},

where the minus sign is absorbed into the definition of Δ compared with the
usual convention; `i`, `j` are discrete real-time step indices (nodes
tᵢ = i·δt, δt = t/N) and `u`, `v ∈ (:+, :-)` label the forward/backward
branch.  Returns a `RealCorrelationFunction` whose four blocks are selected by
`branch(x; b1, b2)`.

Arguments:
- `f`:  the bath spectrum function A(ε) as an `AbstractBoundedFunction`;
- `disperse`:  optional dispersion applied to the integration variable,
  `ε -> disperse(ε)` (default `identity`);
- `β`:  the inverse temperature;
- `μ`:  the chemical potential;
- `t`:  the total evolution time;
- `N`:  number of discrete real-time steps, such that δt = t/N.

The `bath` convenience form reads `β` and `μ` from the bath object.
"""
fermionic_Δt(f::AbstractBoundedFunction, β::Real, t::Real, N::Int, μ::Real) = fermionic_Δt(f, identity, β, t, N, μ)

# function fermionic_Δt(f::AbstractBoundedFunction, β::Real, t::Real, N::Int, μ::Real)
#     β = convert(Float64, β)
#     μ = convert(Float64, μ)
#     δt = t / N
#     g₁(ε) = _f₁(β, μ, ε); g₂(ε) = _f₂(β, μ, ε)
#     fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
#     fⱼⱼ = _fⱼⱼ_r(f, δt)

#     ### G₊₊
#     ηⱼₖ = zeros(ComplexF64, N+1)
#     ηₖⱼ = zeros(ComplexF64, N+1)

#     ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
#     for i = 1:N
#         ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
#     end
#     ηₖⱼ[1] = quadgkwrapper(fⱼⱼ' * g₂)
#     for i = 1:N
#         ηₖⱼ[i+1] = quadgkwrapper(fⱼₖ(i)' * g₂)
#     end
#     G₊₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)   

#     ### G₊₋
#     ηⱼₖ = zeros(ComplexF64, N+1)
#     ηₖⱼ = zeros(ComplexF64, N+1)

#     ηⱼₖ[1] = quadgkwrapper(fⱼⱼ * g₂)
#     for i = 1:N
#         ηⱼₖ[i+1] = quadgkwrapper(fⱼₖ(i) * g₂)
#     end
#     ηₖⱼ[1] = quadgkwrapper(fⱼⱼ' * g₂)
#     for i = 1:N
#         ηₖⱼ[i+1] = ηⱼₖ[i+1]'
#     end
#     G₊₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

#     ### G₋₊
#     ηⱼₖ = zeros(ComplexF64, N+1)
#     ηₖⱼ = zeros(ComplexF64, N+1)

#     ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
#     for i = 1:N
#         ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i)*g₁)
#     end
#     ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₁)
#     for i = 1:N
#         ηₖⱼ[i+1] = ηⱼₖ[i+1]'
#     end
#     G₋₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

#     ### G₋₋
#     ηⱼₖ = zeros(ComplexF64, N+1)
#     ηₖⱼ = zeros(ComplexF64, N+1)

#     ηⱼₖ[1] = quadgkwrapper(fⱼⱼ * g₂)
#     for i = 1:N
#         ηⱼₖ[i+1] = quadgkwrapper(fⱼₖ(i) * g₂)
#     end
#     ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₁)
#     for i = 1:N
#         ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(i)' * g₁)
#     end
#     G₋₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)
#     return RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
# end

function fermionic_Δt(f::AbstractBoundedFunction, disperse::Function, β::Real, t::Real, N::Int, μ::Real)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    δt = t / N
    g₁(ε) = fermidirac(β, μ - disperse(ε)); g₂(ε) = fermidirac(β, disperse(ε) - μ)
    fⱼₖ(Δk) = _fⱼₖ_r_disperse(f, disperse, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r_disperse(f, disperse, δt)
    fₖₖ = _fₖₖ_r_disperse(f, disperse, δt)

    # The four branch blocks share four distinct Toeplitz vectors:
    #   G₊₊ = CM(a, b)   G₊₋ = CM(c, b)   G₋₊ = CM(a, d)   G₋₋ = CM(c, d)
    # computed once here (previously each was quadruplicated, doubling the
    # quadgk cost).  Backward-sector entries use the negative-lag cell kernel
    # fⱼₖ(-i) with the spectrum un-conjugated — the physical correlator carries
    # the same weight at forward and backward lags; conjugation is only used
    # for the branch DERIVATION, which is exact for real-valued spectra.
    a = zeros(ComplexF64, N+1)
    b = zeros(ComplexF64, N+1)
    a[1] = quadgkwrapper(-fⱼⱼ * g₁)
    b[1] = quadgkwrapper(fₖₖ * g₂)
    for i = 1:N
        a[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
        b[i+1] = quadgkwrapper(fⱼₖ(-i) * g₂)
    end
    # secondary pair computed by its own quadrature with the spectrum
    # un-conjugated — complex-spectrum safe (BCS anomalous coherence weights),
    # and for real-valued spectra bit-identical to conj(b)/conj(a)
    ηⱼₖ₂ = zeros(ComplexF64, N+1)
    ηₖⱼ₂ = zeros(ComplexF64, N+1)
    ηⱼₖ₂[1] = quadgkwrapper(fⱼⱼ * g₂)
    ηₖⱼ₂[1] = quadgkwrapper(-fₖₖ * g₁)
    for i = 1:N
        ηⱼₖ₂[i+1] = quadgkwrapper(fⱼₖ(i) * g₂)
        ηₖⱼ₂[i+1] = quadgkwrapper(-fⱼₖ(-i) * g₁)
    end
    return RealCorrelationFunction(a, b, ηⱼₖ₂, ηₖⱼ₂)
end
