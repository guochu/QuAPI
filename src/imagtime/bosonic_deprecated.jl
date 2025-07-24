Δτ2(bath::AbstractBosonicBath; N::Int, δτ::Real=bath.β/N) = bosonic_Δτ2(bath.spectrum, β=bath.β, N=N, μ=bath.μ, δτ=δτ)
bosonic_Δτ2(f::AbstractBoundedFunction; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) = bosonic_Δτ2(f, β, N, μ, δτ)

"""
    bosonic_Δτ(f, β::Real, N::Int)

f: the spectrum function
"""
function bosonic_Δτ2(f0::AbstractBoundedFunction, β::Real, N::Int, μ::Real, δτ::Real=β / N)
    # f′, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    # lb -= μ
    # ub -= μ
    # f = bounded(ϵ->f′(ϵ + μ), lb, ub) 
    f = spectrumshift(f0, μ)
    # δτ = β / N
    
    g₁(ϵ) = _g₁(β, 0., ϵ)
    g₂(ϵ) = _g₂(β, 0., ϵ)
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
    ηₖⱼ[1] = quadgkwrapper(fₖₖ * g₂)
    for k = 1:L-1
        ηₖⱼ[k+1] = quadgkwrapper(fₖⱼ(k) * g₂)
    end
    ImagCorrelationFunction(CorrelationMatrix{Float64}(ηⱼₖ, ηₖⱼ))
end