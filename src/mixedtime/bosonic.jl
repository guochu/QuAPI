Δm(bath::AbstractBosonicNormalBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = bosonic_Δm(bath.spectrum, β=bath.β, μ=bath.μ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)
function bosonic_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
    δt = t / Nt
    g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
    # real time
    fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r(f, δt)
    fₖₖ = _fₖₖ_r(f, δt)
    # imaginary time
    # hⱼₖ(Δk::Int) = _fⱼₖ_i(f, Δk, δτ)
    # hⱼⱼ = _fⱼⱼ_i(f, δτ)
    # hₖₖ = _fₖₖ_i(f, δτ)
    hⱼₖ(Δk::Int) = _bosonic_fⱼₖ_i(f, β, μ, Δk, δτ)
    hₖⱼ(Δk::Int) = _bosonic_fₖⱼ_i(f, β, μ, Δk, δτ)
    hⱼⱼ = _bosonic_fⱼⱼ_i(f, β, μ, δτ)
    hₖₖ = _bosonic_fₖₖ_i(f, β, μ, δτ)

    # mixed time
    # lⱼₖ(j::Int, k::Int) = _lⱼₖ(f, j, k, δt, δτ)
    # lₖⱼ(k::Int, j::Int) = _lₖⱼ(f, k, j, δt, δτ)
    lⱼₖ(j::Int, k::Int) = _bosonic_lⱼₖ(f, β, μ, j, k, δt, δτ)
    lₖⱼ(k::Int, j::Int) = _bosonic_lₖⱼ(f, β, μ, k, j, δt, δτ)

    N, M = Nt, Nτ

    # real time
    ηⱼₖ = zeros(ComplexF64, Nt+1)
    ηₖⱼ = zeros(ComplexF64, Nt+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
    for i = 1:Nt
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(-fₖₖ * g₂)
    for i = 1:Nt
        ηₖⱼ[i+1] = quadgkwrapper(-fⱼₖ(-i) * g₂)
    end

    # imag time
    ξⱼₖ = zeros(ComplexF64, M) # j >= k
    ξₖⱼ = zeros(ComplexF64, M)

    ξⱼₖ[1] = quadgkwrapper(-hⱼⱼ)
    for k in 2:M
        ξⱼₖ[k] = quadgkwrapper(-hⱼₖ(k-1))
    end
    
    ξₖⱼ[1] = quadgkwrapper(-hₖₖ)
    for k in 2:M
        ξₖⱼ[k] = quadgkwrapper(-hₖⱼ(k-1))
    end

    # mix time
    ζⱼₖ = zeros(ComplexF64, M, N+1)
    ζₖⱼ = zeros(ComplexF64, N+1, M)

    for j = 1:M, k = 1:N+1
        ζⱼₖ[j,k] = quadgkwrapper(lⱼₖ(j-1,k-1))
        ζₖⱼ[k,j] = quadgkwrapper(lₖⱼ(k-1,j-1))
    end
    MixedCorrelationFunction(ηⱼₖ,ηₖⱼ,ξⱼₖ,ξₖⱼ,ζⱼₖ,ζₖⱼ)
end


function _bosonic_lⱼₖ(f::AbstractBoundedFunction, β, μ, j::Int, k::Int, δt, δτ)
    function g(ε)
        x = 1-exp(-safe_mult(β, ε-μ))
        if abs(ε) > QuAPI_tol
            -im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/(ε^2 * x)
        else
            -exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt / x
        end
    end
    return f * g
end

function _bosonic_lₖⱼ(f::AbstractBoundedFunction, β, μ, k::Int, j::Int, δt, δτ)
    function g(ε)
        x = 1 - exp(-safe_mult(β, ε-μ))
        y = exp(-safe_mult(β - j*δτ, ε-μ)) * exp(j*δτ*μ)
        if abs(ε) > QuAPI_tol
            -im*y*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/(ε^2 * x)
        else
            -y*exp(-im*ε*k*δt)*δτ*δt / x
        end
    end
    return f * g
end