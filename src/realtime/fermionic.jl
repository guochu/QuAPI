Ct(bath::AbstractFermionicBath; N::Int, t::Real) = fermionic_Ct(bath.spectrum, β=bath.β, N=N, t=t, μ=bath.μ)
fermionic_Ct(f::SpectrumFunction; β::Real, N::Int, t::Real, μ::Real=0) = fermionic_Ct(f, β, N, t/N, μ)

function fermionic_Ct(f0::SpectrumFunction, β::Real, N::Int, δt::Real, μ::Real)
    f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
    fⱼₖ(Δk, ε) = _fⱼₖ_r(f, Δk, ε, δt)
    fⱼⱼ(ε) = _fⱼⱼ_r(f, ε, δt)

    ### G₊₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε), lb, ub)[1] 
    for i = 1:N
        ηⱼₖ[i+1] = quadgk(ε -> -g₁(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> g₂(ε)*fⱼⱼ(ε)', lb, ub)[1]
    for i = 1:N
        ηₖⱼ[i+1] = quadgk(ε -> g₂(ε)*fⱼₖ(i,ε)', lb, ub)[1]
    end
    G₊₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)   

    ### G₊₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgk(ε -> g₂(ε)*fⱼⱼ(ε), lb, ub)[1] 
    for i = 1:N
        ηⱼₖ[i+1] = quadgk(ε -> g₂(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> g₂(ε)*fⱼⱼ(ε)', lb, ub)[1]
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₊₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε), lb, ub)[1] 
    for i = 1:N
        ηⱼₖ[i+1] = quadgk(ε -> -g₁(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε)', lb, ub)[1]
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₋₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgk(ε -> g₂(ε)*fⱼⱼ(ε), lb, ub)[1] 
    for i = 1:N
        ηⱼₖ[i+1] = quadgk(ε -> g₂(ε)*fⱼₖ(i,ε), lb, ub)[1]
    end
    ηₖⱼ[1] = quadgk(ε -> -g₁(ε)*fⱼⱼ(ε)', lb, ub)[1]
    for i = 1:N
        ηₖⱼ[i+1] = quadgk(ε -> -g₁(ε)*fⱼₖ(i,ε)', lb, ub)[1]
    end
    G₋₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)
    return RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
end


# assume j > k
function _fⱼₖ_r(f, Δk::Int, ε::Float64, δt)
    if (abs(ε) > tol)
        2*f(ε)/ε^2*exp(-im*ε*Δk*δt)*(1-cos(ε*δt))
    else
        f(ε)*exp(-im*ε*Δk*δt)*δt^2
    end
end

function _fⱼⱼ_r(f, ε::Float64, δt)
    if (abs(ε) > tol)
        f(ε)/ε^2*((1-im*ε*δt)-exp(-im*ε*δt))
    else
        0.5*f(ε)*δt^2
    end
end

function _fₖₖ_r(f, ε::Float64, δt)
    if (abs(ε) > tol)
        f(ε)/ε^2*((1+im*ε*δt)-exp(im*ε*δt))
    else
        0.5*f(ε)*δt^2
    end
end
