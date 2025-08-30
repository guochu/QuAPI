Δt(bath::AbstractBosonicNormalBath; N::Int, t::Real) = bosonic_Δt(bath.spectrum, β=bath.β, N=N, t=t, μ=bath.μ)
bosonic_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, μ::Real=0) = bosonic_Δt(f, β, N, t/N, μ)
bosonic_Δt(f::AbstractBoundedFunction, disperse::Function; β::Real, N::Int, t::Real, μ::Real=0) = bosonic_Δt(f, disperse, β, N, t/N, μ)


function bosonic_Δt(f::AbstractBoundedFunction, β::Real, N::Int, δt::Real, μ::Real)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
    fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r(f, δt)

    ### G₊₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₂)
    for i = 1:N
        ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(i)' * g₂)
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

    ηⱼₖ[1] = quadgkwrapper(fⱼⱼ * g₁)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(fⱼⱼ' * g₁)
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₋₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₂)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₂)
    end
    ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₁)
    for i = 1:N
        ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(i)' * g₁)
    end
    G₋₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)
    return RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
end

function bosonic_Δt(f::AbstractBoundedFunction, disperse::Function, β::Real, N::Int, δt::Real, μ::Real)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    g₁(ε) = _g₁(β, μ, disperse(ε)); g₂(ε) = _g₂(β, μ, disperse(ε))
    fⱼₖ(Δk) = _fⱼₖ_r_disperse(f, disperse, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r_disperse(f, disperse, δt)

    ### G₊₊
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₂)
    for i = 1:N
        ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(i)' * g₂)
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

    ηⱼₖ[1] = quadgkwrapper(fⱼⱼ * g₁)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(fⱼⱼ' * g₁)
    for i = 1:N
        ηₖⱼ[i+1] = ηⱼₖ[i+1]'
    end
    G₋₊ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)  

    ### G₋₋
    ηⱼₖ = zeros(ComplexF64, N+1)
    ηₖⱼ = zeros(ComplexF64, N+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₂)
    for i = 1:N
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₂)
    end
    ηₖⱼ[1] = quadgkwrapper(-fⱼⱼ' * g₁)
    for i = 1:N
        ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(i)' * g₁)
    end
    G₋₋ = CorrelationMatrix(ηⱼₖ, ηₖⱼ)
    return RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋)
end
