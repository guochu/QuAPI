Δm(bath::AbstractFermionicBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = fermionic_Δm(bath.spectrum, β=bath.β, μ=bath.μ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)


"""
    fermionic_Δm(f, β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0)

Fermionic hybiriization function Δᵢⱼᵘᵛ on the Kadanoff-Baym contour

The Feynman-Vernon influence functional has the form 
        I[ā, a] = e^{ΣᵤᵥΣᵢⱼāᵢᵘΔᵢⱼᵘᵛaⱼᵛ},
where we have absorbed the minus sign into the definition of Δ compared
to the usually used convention
Here i, j are discrete time step indices, u, v = ±, τ denotes the branch
labels (+ means forward branch, - means backward branch, τ means imaginary
time branch)

f: the spectrum function
β: the inverse temperature
Nτ: number of discrete real time steps
μ: the chemical potential
t: the total evolution time
Nt: number of discrete real time steps,
such that we have δt = t/Nt, δτ = β/Nτ
"""
function fermionic_Δm(f::AbstractSpectrumFunction; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
    δt = t / Nt
    g₁(ε) = _f₁(β, μ, ε); g₂(ε) = _f₂(β, μ, ε)
    # real time
    fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r(f, δt)
    fₖₖ = _fₖₖ_r(f, δt)
    # imaginary time
    hⱼₖ(Δk::Int) = _fⱼₖ_i(f, Δk, δτ)
    hⱼⱼ = _fⱼⱼ_i(f, δτ)
    hₖₖ = _fₖₖ_i(f, δτ)
    # mixed time
    lⱼₖ(j::Int, k::Int) = _lⱼₖ(f, j, k, δt, δτ)
    lₖⱼ(k::Int, j::Int) = _lₖⱼ(f, k, j, δt, δτ)
    N, M = Nt, Nτ

    # real time
    ηⱼₖ = zeros(ComplexF64, Nt+1)
    ηₖⱼ = zeros(ComplexF64, Nt+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
    for i = 1:Nt
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(fₖₖ * g₂)
    for i = 1:Nt
        ηₖⱼ[i+1] = quadgkwrapper(fⱼₖ(-i) * g₂)
    end

    # imag time
    ξⱼₖ = zeros(ComplexF64, M) # j >= k
    ξₖⱼ = zeros(Float64, M)

    ξⱼₖ[1] = quadgkwrapper(-hⱼⱼ * g₁)
    for k = 2:M
        ξⱼₖ[k] = quadgkwrapper(-hⱼₖ(k-1) * g₁)
    end
    
    ξₖⱼ[1] = quadgkwrapper(hₖₖ * g₂)
    for k = 2:M
        ξₖⱼ[k] = quadgkwrapper(hⱼₖ(1-k) * g₂)
    end

    # mix time
    ζⱼₖ = zeros(ComplexF64, M, N+1)
    ζₖⱼ = zeros(ComplexF64, N+1, M)

    for j = 1:M, k = 1:N+1
        ζⱼₖ[j,k] = quadgkwrapper(-lⱼₖ(j-1,k-1) * g₁)
        ζₖⱼ[k,j] = quadgkwrapper(lₖⱼ(k-1,j-1) * g₂)
    end
    MixedCorrelationFunction(ηⱼₖ,ηₖⱼ,ξⱼₖ,ξₖⱼ,ζⱼₖ,ζₖⱼ)
end

# # mixed part lⱼₖ for G31 and lₖⱼ for G13
# function _lⱼₖ(f, j::Int, k::Int, ε::Float64, δt, δτ)
#     if (abs(ε) > tol)
#         im*f(ε)/ε^2*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)
#     else
#         f(ε)*exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt
#     end
# end

# function _lₖⱼ(f, k::Int, j::Int, ε::Float64, δt, δτ)
#     if (abs(ε) > tol)
#         im*f(ε)/ε^2*exp(ε*j*δτ)*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)
#     else
#         f(ε)*exp(ε*j*δτ)*exp(-im*ε*k*δt)*δτ*δt
#     end
# end

function _lⱼₖ(f::AbstractBoundedFunction, j::Int, k::Int, δt, δτ)
    g(ε) = ifelse(abs(ε) > tol, im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/ε^2, exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt)
    return f * g
end

function _lₖⱼ(f::AbstractBoundedFunction, k::Int, j::Int, δt, δτ)
    g(ε) = ifelse(abs(ε) > tol, im*exp(ε*j*δτ)*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/ε^2, exp(ε*j*δτ)*exp(-im*ε*k*δt)*δτ*δt)
    return f * g
end
