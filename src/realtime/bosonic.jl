Δt(bath::AbstractBosonicNormalBath; N::Int, t::Real) = bosonic_Δt(bath.spectrum, β=bath.β, N=N, t=t, μ=bath.μ)
bosonic_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, μ::Real=0) = bosonic_Δt(f, β, N, t/N, μ)
bosonic_Δt(f::AbstractBoundedFunction, disperse::Function; β::Real, N::Int, t::Real, μ::Real=0) = bosonic_Δt(f, disperse, β, N, t/N, μ)

"""
    bosonic_Δt(bath::AbstractBosonicNormalBath; N::Int, t::Real)
    bosonic_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, μ::Real=0)
    bosonic_Δt(f::AbstractBoundedFunction, disperse::Function; β::Real, N::Int, t::Real, μ::Real=0)

Bosonic hybridization function Δᵢⱼᵘᵛ on the real-time (closed Keldysh)
contour.

Same convention as `fermionic_Δt`: the influence functional reads
I = e^{Σᵤᵥ Σᵢⱼ sᵢᵘ Δᵢⱼᵘᵛ sⱼᵛ} (spin indices suppressed) with the minus sign
absorbed into Δ; `i`, `j` are discrete real-time step indices (nodes
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

Note:  for μ = 0 the Bose factor carries a 1/ε singularity at ε = 0, and the
small-ε branch of the integrand diverges as `1/(β ε)`, which is `Inf` exactly
at ε = 0.  Spectra with support touching ε = 0 (a flat band, an ohmic-like
edge, or a delta peak at zero frequency) can therefore make the quadrature
raise a `DomainError`; shifting the band edge slightly above 0 avoids the
instability.
"""
bosonic_Δt(f::AbstractBoundedFunction, β::Real, N::Int, δt::Real, μ::Real) = bosonic_Δt(f, identity, β, N, δt, μ)

# function bosonic_Δt(f::AbstractBoundedFunction, β::Real, N::Int, δt::Real, μ::Real)
#     β = convert(Float64, β)
#     μ = convert(Float64, μ)
#     g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
#     fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
#     fⱼⱼ = _fⱼⱼ_r(f, δt)

#     ### G₊₊
#     ηⱼₖ = zeros(ComplexF64, N+1)
#     ηₖⱼ = zeros(ComplexF64, N+1)

#     ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
#     for i = 1:N
#         ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
#     end
#     ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₂)
#     for i = 1:N
#         ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(i)' * g₂)
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

#     ηⱼₖ[1] = quadgkwrapper(fⱼⱼ * g₁)
#     for i = 1:N
#         ηⱼₖ[i+1] = quadgkwrapper(fⱼₖ(i) * g₁)
#     end
#     ηₖⱼ[1] = quadgkwrapper(fⱼⱼ' * g₁)
#     for i = 1:N
#         ηₖⱼ[i+1] = ηⱼₖ[i+1]'
#     end
#     G₋₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

#     ### G₋₋
#     ηⱼₖ = zeros(ComplexF64, N+1)
#     ηₖⱼ = zeros(ComplexF64, N+1)

#     ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₂)
#     for i = 1:N
#         ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₂)
#     end
#     ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₁)
#     for i = 1:N
#         ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(i)' * g₁)
#     end
#     G₋₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)
#     return RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
# end

function bosonic_Δt(f::AbstractBoundedFunction, disperse::Function, β::Real, N::Int, δt::Real, μ::Real)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    # g₁(ε) = _g₁(β, μ, disperse(ε)); g₂(ε) = _g₂(β, μ, disperse(ε))
    fⱼₖ(Δk) = _bosonic_fⱼₖ_r_disperse(f, disperse, β, μ, Δk, δt)
    fₖⱼ(Δk) = _bosonic_fₖⱼ_r_disperse(f, disperse, β, μ, Δk, δt)
    fⱼⱼ = _bosonic_fⱼⱼ_r_disperse(f, disperse, β, μ, δt)
    fₖₖ = _bosonic_fₖₖ_r_disperse(f, disperse, β, μ, δt)
    # the equal-time entries of G₊₋/G₋₊ are stored as the two ordering pieces
    # (conjugate pair of the gⱼⱼ/gₖₖ kernels), matching the legacy `bosonic_Δt_deprecated`
    # storage convention; only for diagonal coupling do the two parts add
    fⱼⱼg₂ = f * (ε -> boseeinstein(β, μ, disperse(ε)) * _gⱼⱼ_r(δt)(disperse(ε)))
    fₖₖg₁ = f * (ε -> (1 + boseeinstein(β, μ, disperse(ε))) * _gₖₖ_r(δt)(disperse(ε)))

    ### G₊₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i))
    end
    ηₖⱼ[1] = quadgkwrapper(-fₖₖ)
    for i = 1:N
        ηₖⱼ[i+1] = quadgkwrapper(-fₖⱼ(-i))
    end
    G₊₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)   

    ### G₊₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(fⱼⱼg₂)
    ηₖⱼ[1] = quadgkwrapper(fₖₖ)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(fₖⱼ(i))
    end
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₊₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(fⱼⱼ)
    ηₖⱼ[1] = quadgkwrapper(fₖₖg₁)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(fⱼₖ(i))
    end
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₋₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(-fₖₖ')
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(-fₖⱼ(i))
    end
    ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ')
    for i = 1:N
        ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(-i))
    end
    G₋₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)
    return RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
end

function _bosonic_fⱼₖ_r_disperse(f::AbstractBoundedFunction, disperse, β, μ, Δk::Int, δt)
    g = _bosonic_gⱼₖ_r(β, μ, Δk, δt)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′  
end
_bosonic_fⱼₖ_r(f::AbstractBoundedFunction, β, μ, Δk::Int, δt) = f * _bosonic_gⱼₖ_r(β, μ, Δk, δt)
function _bosonic_gⱼₖ_r(β, μ, Δk::Int, δt)
    return ε -> begin
        if abs(ε) > get_quapi_tol()
            return 2*(1 + boseeinstein(β, μ, ε))*exp(-im*ε*Δk*δt)*(1-cos(ε*δt))/ε^2
        else
            if μ == zero(μ)
                return exp(-im*ε*Δk*δt)*δt^2/(β*ε)
            else
                return (1 + boseeinstein(β, μ, ε))*exp(-im*ε*Δk*δt)*δt^2
            end
        end
    end
end

function _bosonic_fₖⱼ_r_disperse(f::AbstractBoundedFunction, disperse, β, μ, Δk::Int, δt)
    g = _bosonic_gₖⱼ_r(β, μ, Δk, δt)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′  
end
_bosonic_fₖⱼ_r(f::AbstractBoundedFunction, β, μ, Δk::Int, δt) = f * _bosonic_gₖⱼ_r(β, μ, Δk, δt)
function _bosonic_gₖⱼ_r(β, μ, Δk::Int, δt)
    return ε -> begin
        if abs(ε) > get_quapi_tol()
            return 2*boseeinstein(β, μ, ε)*exp(-im*ε*Δk*δt)*(1-cos(ε*δt))/ε^2
        else
            if μ == zero(μ)
                return exp(-im*ε*Δk*δt)*δt^2/(β*ε)
            else
                return boseeinstein(β, μ, ε) * exp(-im*ε*Δk*δt)*δt^2
            end
        end
    end
end

function _bosonic_fⱼⱼ_r_disperse(f::AbstractBoundedFunction, disperse, β, μ, δt)
    g = _bosonic_gⱼⱼ_r(β, μ, δt)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′  
end
_bosonic_fⱼⱼ_r(f::AbstractBoundedFunction, β, μ, δt) = f * _bosonic_gⱼⱼ_r(β, μ, δt)
function _bosonic_gⱼⱼ_r(β, μ, δt)
    return ε -> begin
        if abs(ε) > get_quapi_tol()
            return (1 + boseeinstein(β, μ, ε))*((1-im*ε*δt)-exp(-im*ε*δt))/ε^2
        else
            if μ == zero(μ)
                return 0.5*δt^2/(β*ε)
            else
                return (1 + boseeinstein(β, μ, ε))*0.5*δt^2
            end
        end
    end
end

function _bosonic_fₖₖ_r_disperse(f::AbstractBoundedFunction, disperse, β, μ, δt)
    g = _bosonic_gₖₖ_r(β, μ, δt)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′      
end
_bosonic_fₖₖ_r(f::AbstractBoundedFunction, β, μ, δt) = f * _bosonic_gₖₖ_r(β, μ, δt)
function _bosonic_gₖₖ_r(β, μ, δt)
    return ε -> begin
        if abs(ε) > get_quapi_tol()
            return boseeinstein(β, μ, ε)*((1+im*ε*δt)-exp(im*ε*δt))/ε^2
        else
            if μ == zero(μ)
                return 0.5*δt^2/(β*ε)
            else
                return boseeinstein(β, μ, ε)*0.5*δt^2
            end
        end
    end
end
