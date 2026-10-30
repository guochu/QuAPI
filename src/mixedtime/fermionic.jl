"""
    Δm(bath::AbstractFermionicNormalBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ)
    Δm(bath::AbstractBosonicNormalBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ)
    Δm(bath::AbstractBCSBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ)

Discretized hybridization function on the L-shaped (Kadanoff-Baym) contour:
an imaginary-time branch of `Nτ` steps δτ = β/Nτ followed by the closed Keldysh
part of total time t with `Nt` steps δt = t/Nt.  Returns a
`MixedCorrelationFunction` (or `BCSCorrelationFunction` for a superconducting
bath) with the nine branch combinations G¹¹ ... G³³.

Convenience wrapper that dispatches to `fermionic_Δm` (fermionic baths),
`bosonic_Δm` (bosonic baths) and `bcs_Δm` (BCS baths), reading `β` and `μ`
from the bath object.  The minus sign of the influence-functional exponent is
absorbed into Δ.
"""
Δm(bath::AbstractFermionicNormalBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = fermionic_Δm(bath.spectrum, β=bath.β, μ=bath.μ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)


"""
    fermionic_Δm(bath::AbstractFermionicNormalBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ)
    fermionic_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
    fermionic_Δm(f::AbstractBoundedFunction, disperse::Function; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)

Fermionic hybridization function Δᵢⱼᵘᵛ on the Kadanoff-Baym (L-shaped)
contour.

The Feynman-Vernon influence functional has the form

    I[ā, a] = e^{Σᵤᵥ Σᵢⱼ āᵢᵘ Δᵢⱼᵘᵛ aⱼᵛ},

where the minus sign is absorbed into the definition of Δ compared with the
usual convention; `i`, `j` are discrete step indices and `u`, `v ∈ (:+, :-, :τ)`
label the branch: `:+`/`:-` the forward/backward real-time branch (Nt+1 nodes
tᵢ = i·δt, δt = t/Nt), `:τ` the imaginary-time branch (Nτ nodes
τᵢ = (i-1)δτ, δτ = β/Nτ).  Returns a `MixedCorrelationFunction`.

Arguments:
- `f`:  the bath spectrum function A(ε) as an `AbstractBoundedFunction`;
- `disperse`:  optional dispersion applied to the integration variable,
  `ε -> disperse(ε)` (default `identity`);
- `β`:  the inverse temperature;
- `μ`:  the chemical potential;
- `Nτ`:  number of imaginary-time steps, δτ = β/Nτ;
- `t`:  the total real evolution time;
- `Nt`:  number of real-time steps, δt = t/Nt;
- `δτ`:  the imaginary-time step size (defaults to β/Nτ).

The `bath` convenience form reads `β` and `μ` from the bath object.
"""
fermionic_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ) = fermionic_Δm(
                f, identity, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)

# function fermionic_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
#     δt = t / Nt
#     g₁(ε) = _f₁(β, μ, ε)
#     g₂(ε) = _f₂(β, μ, ε)
#     # real time
#     fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
#     fⱼⱼ = _fⱼⱼ_r(f, δt)
#     fₖₖ = _fₖₖ_r(f, δt)
#     # imaginary time
#     # hⱼₖ(Δk::Int) = _fⱼₖ_i(f, Δk, δτ)
#     # hⱼⱼ = _fⱼⱼ_i(f, δτ)
#     # hₖₖ = _fₖₖ_i(f, δτ)
#     hⱼₖ(Δk::Int) = _fermionic_fⱼₖ_i(f, β, μ, Δk, δτ)
#     hₖⱼ(Δk::Int) = _fermionic_fₖⱼ_i(f, β, μ, Δk, δτ)
#     hⱼⱼ = _fermionic_fⱼⱼ_i(f, β, μ, δτ)
#     hₖₖ = _fermionic_fₖₖ_i(f, β, μ, δτ)

#     # mixed time
#     # lⱼₖ(j::Int, k::Int) = _lⱼₖ(f, j, k, δt, δτ)
#     # lₖⱼ(k::Int, j::Int) = _lₖⱼ(f, k, j, δt, δτ)
#     lⱼₖ(j::Int, k::Int) = _fermionic_lⱼₖ(f, β, μ, j, k, δt, δτ)
#     lₖⱼ(k::Int, j::Int) = _fermionic_lₖⱼ(f, β, μ, k, j, δt, δτ)

#     N, M = Nt, Nτ

#     # real time
#     ηⱼₖ = zeros(ComplexF64, Nt+1)
#     ηₖⱼ = zeros(ComplexF64, Nt+1)

#     ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
#     for i = 1:Nt
#         ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
#     end
#     ηₖⱼ[1] = quadgkwrapper(fₖₖ * g₂)
#     for i = 1:Nt
#         ηₖⱼ[i+1] = quadgkwrapper(fⱼₖ(-i) * g₂)
#     end

#     # imag time
#     ξⱼₖ = zeros(ComplexF64, M) # j >= k
#     ξₖⱼ = zeros(ComplexF64, M)

#     ξⱼₖ[1] = quadgkwrapper(-hⱼⱼ)
#     for k in 2:M
#         ξⱼₖ[k] = quadgkwrapper(-hⱼₖ(k-1))
#     end
    
#     ξₖⱼ[1] = quadgkwrapper(-hₖₖ)
#     for k in 2:M
#         ξₖⱼ[k] = quadgkwrapper(-hₖⱼ(k-1))
#     end

#     # mix time
#     ζⱼₖ = zeros(ComplexF64, M, N+1)
#     ζₖⱼ = zeros(ComplexF64, N+1, M)

#     for j = 1:M, k = 1:N+1
#         ζⱼₖ[j,k] = quadgkwrapper(lⱼₖ(j-1,k-1))
#         ζₖⱼ[k,j] = quadgkwrapper(lₖⱼ(k-1,j-1))
#     end
#     MixedCorrelationFunction(ηⱼₖ,ηₖⱼ,ξⱼₖ,ξₖⱼ,ζⱼₖ,ζₖⱼ)
# end

function fermionic_Δm(f::AbstractBoundedFunction, disperse::Function; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
    δt = t / Nt
    g₁(ε) = fermidirac(β, μ - disperse(ε))
    g₂(ε) = fermidirac(β, disperse(ε) - μ)
    # real time
    fⱼₖ(Δk) = _fⱼₖ_r_disperse(f, disperse, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r_disperse(f, disperse, δt)
    fₖₖ = _fₖₖ_r_disperse(f, disperse, δt)
    # imaginary time
    # hⱼₖ(Δk::Int) = _fⱼₖ_i(f, Δk, δτ)
    # hⱼⱼ = _fⱼⱼ_i(f, δτ)
    # hₖₖ = _fₖₖ_i(f, δτ)
    hⱼₖ(Δk::Int) = _fermionic_fⱼₖ_i_disperse(f, disperse, β, μ, Δk, δτ)
    hₖⱼ(Δk::Int) = _fermionic_fₖⱼ_i_disperse(f, disperse, β, μ, Δk, δτ)
    hⱼⱼ = _fermionic_fⱼⱼ_i_disperse(f, disperse, β, μ, δτ)
    hₖₖ = _fermionic_fₖₖ_i_disperse(f, disperse, β, μ, δτ)

    # mixed time
    # lⱼₖ(j::Int, k::Int) = _lⱼₖ(f, j, k, δt, δτ)
    # lₖⱼ(k::Int, j::Int) = _lₖⱼ(f, k, j, δt, δτ)
    lⱼₖ(j::Int, k::Int) = _fermionic_lⱼₖ_disperse(f, disperse, β, μ, j, k, δt, δτ)
    lₖⱼ(k::Int, j::Int) = _fermionic_lₖⱼ_disperse(f, disperse, β, μ, k, j, δt, δτ)

    N, M = Nt, Nτ

    # real time
    ηⱼₖ = zeros(ComplexF64, Nt+1)
    ηₖⱼ = zeros(ComplexF64, Nt+1)
    ηⱼₖ₂ = zeros(ComplexF64, Nt+1)
    ηₖⱼ₂ = zeros(ComplexF64, Nt+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
    for i = 1:Nt
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(fₖₖ * g₂)
    for i = 1:Nt
        ηₖⱼ[i+1] = quadgkwrapper(fⱼₖ(-i) * g₂)
    end
    # secondary pair (un-conjugated spectrum; equals conj of the primary pair
    # for real-valued spectra)
    ηⱼₖ₂[1] = quadgkwrapper(fⱼⱼ * g₂)
    for i = 1:Nt
        ηⱼₖ₂[i+1] = quadgkwrapper(fⱼₖ(i) * g₂)
    end
    ηₖⱼ₂[1] = quadgkwrapper(-fₖₖ * g₁)
    for i = 1:Nt
        ηₖⱼ₂[i+1] = quadgkwrapper(-fⱼₖ(-i) * g₁)
    end

    # imag time
    ξⱼₖ = zeros(ComplexF64, M) # j >= k
    ξₖⱼ = zeros(ComplexF64, M)

    ξⱼₖ[1] = quadgkwrapper(-hⱼⱼ)
    for k in 2:M
        ξⱼₖ[k] = quadgkwrapper(-hⱼₖ(k-1))
    end
    
    ξₖⱼ[1] = quadgkwrapper(-hₖₖ)
    for k in 2:M
        ξₖⱼ[k] = quadgkwrapper(-hₖⱼ(k-1))
    end

    # mix time
    ζⱼₖ = zeros(ComplexF64, M, N+1)
    ζₖⱼ = zeros(ComplexF64, N+1, M)

    for j = 1:M, k = 1:N+1
        ζⱼₖ[j,k] = quadgkwrapper(lⱼₖ(j-1,k-1))
        ζₖⱼ[k,j] = quadgkwrapper(lₖⱼ(k-1,j-1))
    end
    MixedCorrelationFunction(ηⱼₖ,ηₖⱼ,ηⱼₖ₂,ηₖⱼ₂,ξⱼₖ,ξₖⱼ,ζⱼₖ,ζₖⱼ)
end

# function _fermionic_lⱼₖ(f::AbstractBoundedFunction, β, μ, j::Int, k::Int, δt, δτ)
#     # g(ε) = ifelse(abs(ε) > QuAPI_tol, im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/ε^2, exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt)
#     function g(ε)
#         if β >= μ
#             x = exp(-safe_mult(β, ε-μ))+1
#             if abs(ε) > QuAPI_tol
#                 -im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/(ε^2 * x)
#             else
#                 -exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt / x
#             end
#         else 
#             x = exp(safe_mult(β, ε-μ)) + 1
#             y = exp(safe_mult(β-j*δτ, ε-μ)) * exp(-j*δτ*μ)
#             if abs(ε) > QuAPI_tol
#                 -im*y*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/(ε^2 * x)
#             else
#                 -y*exp(im*ε*k*δt)*δτ*δt / x
#             end
#         end
#     end
#     return f * g
# end

# function _fermionic_lₖⱼ(f::AbstractBoundedFunction, β, μ, k::Int, j::Int, δt, δτ)
#     # g(ε) = ifelse(abs(ε) > QuAPI_tol, im*exp(ε*j*δτ)*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/ε^2, exp(ε*j*δτ)*exp(-im*ε*k*δt)*δτ*δt)
#     function g(ε)
#         if β >= μ
#             x = exp(-safe_mult(β, ε-μ))+1
#             y = exp(-safe_mult(β - j*δτ, ε-μ)) * exp(j*δτ*μ)
#             if abs(ε) > QuAPI_tol
#                 im*y*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/(ε^2 * x)
#             else
#                 y*exp(-im*ε*k*δt)*δτ*δt / x
#             end
#         else
#             x = exp(safe_mult(β, ε-μ))+1
#             y = exp(j*δτ * ε)

#             if abs(ε) > QuAPI_tol
#                 im*y*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/(ε^2 * x)
#             else
#                 y*exp(-im*ε*k*δt)*δτ*δt / x
#             end
#         end
#     end
#     return f * g
# end

function _fermionic_lⱼₖ_disperse(f::AbstractBoundedFunction, disperse, β, μ, j::Int, k::Int, δt, δτ)
    g = _fermionic_gⱼₖ(β, μ, j, k, δt, δτ)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′  
end
_fermionic_lⱼₖ(f::AbstractBoundedFunction, β, μ, j::Int, k::Int, δt, δτ) = f * _fermionic_gⱼₖ(β, μ, j, k, δt, δτ)
function _fermionic_gⱼₖ(β, μ, j::Int, k::Int, δt, δτ)
    return ε -> begin
        if β >= μ
            x = exp(-safe_mult(β, ε-μ))+1
            if abs(ε) > get_quapi_tol()
                -im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/(ε^2 * x)
            else
                -exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt / x
            end
        else 
            x = exp(safe_mult(β, ε-μ)) + 1
            y = exp(safe_mult(β-j*δτ, ε-μ)) * exp(-j*δτ*μ)
            if abs(ε) > get_quapi_tol()
                -im*y*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/(ε^2 * x)
            else
                -y*exp(im*ε*k*δt)*δτ*δt / x
            end
        end
    end
end

function _fermionic_lₖⱼ_disperse(f::AbstractBoundedFunction, disperse, β, μ, k::Int, j::Int, δt, δτ)
    g = _fermionic_gₖⱼ(β, μ, k, j, δt, δτ)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′      
end
_fermionic_lₖⱼ(f::AbstractBoundedFunction, β, μ, k::Int, j::Int, δt, δτ) = f * _fermionic_gₖⱼ(β, μ, k, j, δt, δτ)
function _fermionic_gₖⱼ(β, μ, k::Int, j::Int, δt, δτ)
    return ε -> begin
        if β >= μ
            x = exp(-safe_mult(β, ε-μ))+1
            y = exp(-safe_mult(β - j*δτ, ε-μ)) * exp(j*δτ*μ)
            if abs(ε) > get_quapi_tol()
                im*y*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/(ε^2 * x)
            else
                y*exp(-im*ε*k*δt)*δτ*δt / x
            end
        else
            x = exp(safe_mult(β, ε-μ))+1
            y = exp(j*δτ * ε)

            if abs(ε) > get_quapi_tol()
                im*y*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/(ε^2 * x)
            else
                y*exp(-im*ε*k*δt)*δτ*δt / x
            end
        end
    end
end
