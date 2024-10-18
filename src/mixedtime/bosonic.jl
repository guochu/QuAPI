Cm(bath::AbstractBosonicBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = bosonic_Cm(bath.spectrum, β=bath.β, μ=bath.μ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)
function bosonic_Cm(f0::SpectrumFunction; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
	f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    δt = t / Nt
    g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
    # real time
    fⱼₖ(Δk, ε) = _fⱼₖ_r(f, Δk, ε, δt)
    fⱼⱼ(ε) = _fⱼⱼ_r(f, ε, δt)
    fₖₖ(ε) = _fₖₖ_r(f, ε, δt)
    # imaginary time
    hⱼₖ(Δk::Int, ε::Float64) = _fⱼₖ_i(f, Δk, ε, δτ)
    hⱼⱼ(ε::Float64) = _fⱼⱼ_i(f, ε, δτ)
    hₖₖ(ε::Float64) = _fₖₖ_i(f, ε, δτ)
    # mixed time
    lⱼₖ(j::Int, k::Int, ε::Float64) = _lⱼₖ(f, j, k, ε, δt, δτ)
    lₖⱼ(k::Int, j::Int, ε::Float64) = _lₖⱼ(f, k, j, ε, δt, δτ)
    N, M = Nt, Nτ

    # real time
    ηⱼₖ = zeros(ComplexF64, Nt+1)
    ηₖⱼ = zeros(ComplexF64, Nt+1)

    ηⱼₖ[1] = quadgk(ε -> g₁(ε)*fⱼⱼ(ε), lb, ub)[1]
    for i = 1:Nt
        ηⱼₖ[i+1] = quadgk(ε -> g₁(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> g₂(ε)*fₖₖ(ε), lb, ub)[1]
    for i = 1:Nt
        ηₖⱼ[i+1] = quadgk(ε -> g₂(ε)*fⱼₖ(-i,ε), lb, ub)[1]
    end

    # imag time
    ξⱼₖ = zeros(ComplexF64, M) # j >= k
    ξₖⱼ = zeros(Float64, M)

    ξⱼₖ[1] = quadgk(ε -> g₁(ε)*hⱼⱼ(ε), lb, ub)[1]
    for k = 2:M
        ξⱼₖ[k] = quadgk(ε -> g₁(ε)*hⱼₖ(k-1,ε), lb, ub)[1]
    end
    
    ξₖⱼ[1] = quadgk(ε -> g₂(ε)*hₖₖ(ε), lb, ub)[1]
    for k = 2:M
        ξₖⱼ[k] = quadgk(ε -> g₂(ε)*hⱼₖ(1-k,ε), lb, ub)[1]
    end

    # mix time
    ζⱼₖ = zeros(ComplexF64, M, N+1)
    ζₖⱼ = zeros(ComplexF64, N+1, M)

    for j = 1:M, k = 1:N+1
        ζⱼₖ[j,k] = quadgk(ε -> g₁(ε)*lⱼₖ(j-1,k-1,ε), lb, ub)[1]
        ζₖⱼ[k,j] = quadgk(ε -> g₂(ε)*lₖⱼ(k-1,j-1,ε), lb, ub)[1]
    end
    MixedCorrelationFunction(ηⱼₖ,ηₖⱼ,ξⱼₖ,ξₖⱼ,ζⱼₖ,ζₖⱼ)
end