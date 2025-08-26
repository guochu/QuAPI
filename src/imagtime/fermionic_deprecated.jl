Δτ2(bath::AbstractFermionicNormalBath; N::Int, δτ::Real=bath.β/N) = fermionic_Δτ2(bath.spectrum, β=bath.β, N=N, μ=bath.μ, δτ=δτ)
fermionic_Δτ2(f::AbstractBoundedFunction; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) = fermionic_Δτ2(f, β, N, μ, δτ)



"""
    fermionic_Δτ(f, β::Real, N::Int, μ::Real)

Fermionic hybiriization function Δᵢⱼ on the imaginary-time axis

The Feynman-Vernon influence functional has the form 
        I[ā, a] = e^{ΣᵢⱼāᵢΔᵢⱼaⱼ},
where we have absorbed the minus sign into the definition of Δ compared
to the usually used convention, i, j are discrete time step indices

f: the spectrum function
β: the inverse temperature
μ: the chemical potential
N: number of discrete imaginary time steps, 
such that δτ = β/N
"""
function fermionic_Δτ2(f0::AbstractBoundedFunction, β::Real, N::Int, μ::Real, δτ::Real=β / N)
    # f′, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    f = spectrumshift(f0, μ)
    # δτ = β / N
    
    g₁(ϵ) = _f₁(β, 0., ϵ)
    g₂(ϵ) = _f₂(β, 0., ϵ)
    fⱼₖ(Δk::Int) = _fⱼₖ_i(f, Δk, δτ)
    fₖⱼ(Δk::Int) = _fₖⱼ_i(f, Δk, δτ)
    fⱼⱼ = _fⱼⱼ_i(f, δτ)
    fₖₖ = _fₖₖ_i(f, δτ)

    # j >= k
    L = N
    ηⱼₖ = zeros(Float64, L)
    ηⱼₖ[1] = quadgkwrapper(fⱼⱼ * g₁)
    for k = 1:L-1
        ηⱼₖ[k+1] = quadgkwrapper(fⱼₖ(k) * g₁)
    end

    ηₖⱼ = zeros(Float64, L)
    ηₖⱼ[1] = quadgkwrapper(fₖₖ * (ε->-g₂(ε)))
    for k = 1:L-1
        ηₖⱼ[k+1] = quadgkwrapper(fₖⱼ(k) * (ε -> -g₂(ε)))
    end
    ImagCorrelationFunction(CorrelationMatrix{Float64}(ηⱼₖ, ηₖⱼ))
end