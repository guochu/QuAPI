Ct(bath::AbstractFermionicBath; N::Int, t::Real) = fermionic_Ct(bath.spectrum, β=bath.β, N=N, t=t, μ=bath.μ)
fermionic_Ct(f::SpectrumFunction; β::Real, N::Int, t::Real, μ::Real=0) = fermionic_Ct(f, β, N, t/N, μ)

function fermionic_Ct(f::SpectrumFunction, β::Real, N::Int, δt::Real, μ::Real)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
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
