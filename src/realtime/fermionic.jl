Δt(bath::AbstractFermionicNormalBath; N::Int, t::Real) = fermionic_Δt(bath.spectrum, β=bath.β, N=N, t=t, μ=bath.μ)
fermionic_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, μ::Real=0) = fermionic_Δt(f, β, t, N, μ)
fermionic_Δt(f::AbstractBoundedFunction, disperse::Function; β::Real, N::Int, t::Real, μ::Real=0) = fermionic_Δt(f, disperse, β, t, N, μ)

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
    g₁(ε) = _f₁(β, μ, disperse(ε)); g₂(ε) = _f₂(β, μ, disperse(ε))
    fⱼₖ(Δk) = _fⱼₖ_r_disperse(f, disperse, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r_disperse(f, disperse, δt)

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
