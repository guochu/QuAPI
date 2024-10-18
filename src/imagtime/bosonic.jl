Cτ(bath::AbstractBosonicBath; N::Int, δτ::Real=bath.β/N) = bosonic_Cτ(bath.spectrum, β=bath.β, N=N, μ=bath.μ, δτ=δτ)
bosonic_Cτ(f::SpectrumFunction; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) = bosonic_Cτ(f, β, N, μ, δτ)

"""
    bosonic_Cτ(f, β::Real, N::Int)

f is the spectrum function
"""
function bosonic_Cτ(f0::SpectrumFunction, β::Real, N::Int, μ::Real, δτ::Real=β / N)
    f′, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    f(ϵ) = f′(ϵ + μ)
    lb -= μ
    ub -= μ
    # δτ = β / N
    
    g₁(ϵ) = _g₁(β, 0., ϵ)
    g₂(ϵ) = _g₂(β, 0., ϵ)
    fⱼₖ(Δk::Int, ε::Float64) = _fⱼₖ_i(f, Δk, ε, δτ)
    fₖⱼ(Δk::Int, ε::Float64) = _fₖⱼ_i(f, Δk, ε, δτ)
    fⱼⱼ(ε::Float64) = _fⱼⱼ_i(f, ε, δτ)
    fₖₖ(ε::Float64) = _fₖₖ_i(f, ε, δτ)

    # j >= k
    L = N
    ηⱼₖ = zeros(Float64, L)
    ηⱼₖ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε), lb, ub)[1]
    for k = 1:L-1
        ηⱼₖ[k+1] = quadgk(ε -> g₁(ε)*fⱼₖ(k,ε), lb, ub)[1]
    end

    ηₖⱼ = zeros(Float64, L)
    ηₖⱼ[1] = quadgk(ε -> g₂(ε)*fₖₖ(ε), lb, ub)[1]
    for k = 1:L-1
        ηₖⱼ[k+1] = quadgk(ε -> g₂(ε)*fₖⱼ(k,ε), lb, ub)[1]
    end
    ImagCorrelationFunction(CorrelationMatrix{Float64}(ηⱼₖ, ηₖⱼ))
end