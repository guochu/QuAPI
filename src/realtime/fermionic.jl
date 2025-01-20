Δt(bath::AbstractFermionicBath; N::Int, t::Real) = fermionic_Δt(bath.spectrum, β=bath.β, N=N, t=t, μ=bath.μ)
fermionic_Δt(f::AbstractSpectrumFunction; β::Real, N::Int, t::Real, μ::Real=0) = fermionic_Δt(f, β, t, N, μ)

"""
    fermionic_Δt(f, β::Real, N::Int, δt::Real, μ::Real)

Fermionic hybiriization function Δᵢⱼᵘᵛ on the real-time axis

The Feynman-Vernon influence functional has the form 
        I[ā, a] = e^{ΣᵤᵥΣᵢⱼāᵢᵘΔᵢⱼᵘᵛaⱼᵛ},
where we have absorbed the minus sign into the definition of Δ compared
to the usually used convention
Here i, j are discrete time step indices, u, v = ± denotes the branch
labels (+ means forward branch, - means backward branch)

f: the spectrum function
β: the inverse temperature
μ: the chemical potential
t: the total evolution time
N: number of discrete real time steps,
such that we have δt = t/Nt
"""
function fermionic_Δt(f::AbstractSpectrumFunction, β::Real, t::Real, N::Int, μ::Real)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    δt = t / N
    g₁(ε) = _f₁(β, μ, ε); g₂(ε) = _f₂(β, μ, ε)
    fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r(f, δt)

    ### G₊₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(fⱼⱼ' * g₂)
    for i = 1:N
        ηₖⱼ[i+1] = quadgkwrapper(fⱼₖ(i)' * g₂)
    end
    G₊₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)   

    ### G₊₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(fⱼⱼ * g₂)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(fⱼₖ(i) * g₂)
    end
    ηₖⱼ[1] = quadgkwrapper(fⱼⱼ' * g₂)
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₊₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i)*g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₁)
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₋₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(fⱼⱼ * g₂)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(fⱼₖ(i) * g₂)
    end
    ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₁)
    for i = 1:N
        ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(i)' * g₁)
    end
    G₋₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)
    return RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
end


# # assume j > k
# function _fⱼₖ_r(f, Δk::Int, ε::Float64, δt)
#     if (abs(ε) > tol)
#         2*f(ε)/ε^2*exp(-im*ε*Δk*δt)*(1-cos(ε*δt))
#     else
#         f(ε)*exp(-im*ε*Δk*δt)*δt^2
#     end
# end

# function _fⱼⱼ_r(f, ε::Float64, δt)
#     if (abs(ε) > tol)
#         f(ε)/ε^2*((1-im*ε*δt)-exp(-im*ε*δt))
#     else
#         0.5*f(ε)*δt^2
#     end
# end

# function _fₖₖ_r(f, ε::Float64, δt)
#     if (abs(ε) > tol)
#         f(ε)/ε^2*((1+im*ε*δt)-exp(im*ε*δt))
#     else
#         0.5*f(ε)*δt^2
#     end
# end

function _fⱼₖ_r(f::AbstractBoundedFunction, Δk::Int, δt)
    g(ε) = ifelse(abs(ε) > tol, 2exp(-im*ε*Δk*δt)*(1-cos(ε*δt))/ε^2, exp(-im*ε*Δk*δt)*δt^2)
    return f * g
end

function _fⱼⱼ_r(f::AbstractBoundedFunction, δt)
    g(ε) = ifelse(abs(ε) > tol, ((1-im*ε*δt)-exp(-im*ε*δt))/ε^2, 0.5*δt^2)
    return f * g
end

function _fₖₖ_r(f::AbstractBoundedFunction, δt)
    g(ε) = ifelse(abs(ε) > tol, ((1+im*ε*δt)-exp(im*ε*δt))/ε^2, 0.5*δt^2)
    return f * g
end
