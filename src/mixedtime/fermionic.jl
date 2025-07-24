Δm(bath::AbstractFermionicBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = fermionic_Δm(bath.spectrum, β=bath.β, μ=bath.μ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)


"""
    fermionic_Δm(f, β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0)

Fermionic hybiriization function Δᵢⱼᵘᵛ on the Kadanoff-Baym contour

The Feynman-Vernon influence functional has the form 
        I[ā, a] = e^{ΣᵤᵥΣᵢⱼāᵢᵘΔᵢⱼᵘᵛaⱼᵛ},
where we have absorbed the minus sign into the definition of Δ compared
to the usually used convention
Here i, j are discrete time step indices, u, v = ±, τ denotes the branch
labels (+ means forward branch, - means backward branch, τ means imaginary
time branch)

f: the spectrum function
β: the inverse temperature
Nτ: number of discrete real time steps
μ: the chemical potential
t: the total evolution time
Nt: number of discrete real time steps,
such that we have δt = t/Nt, δτ = β/Nτ
"""
function fermionic_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
    δt = t / Nt
    g₁(ε) = _f₁(β, μ, ε)
    g₂(ε) = _f₂(β, μ, ε)
    # real time
    fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r(f, δt)
    fₖₖ = _fₖₖ_r(f, δt)
    # imaginary time
    # hⱼₖ(Δk::Int) = _fⱼₖ_i(f, Δk, δτ)
    # hⱼⱼ = _fⱼⱼ_i(f, δτ)
    # hₖₖ = _fₖₖ_i(f, δτ)
    hⱼₖ(Δk::Int) = _fermionic_fⱼₖ_i(f, β, μ, Δk, δτ)
    hₖⱼ(Δk::Int) = _fermionic_fₖⱼ_i(f, β, μ, Δk, δτ)
    hⱼⱼ = _fermionic_fⱼⱼ_i(f, β, μ, δτ)
    hₖₖ = _fermionic_fₖₖ_i(f, β, μ, δτ)

    # mixed time
    # lⱼₖ(j::Int, k::Int) = _lⱼₖ(f, j, k, δt, δτ)
    # lₖⱼ(k::Int, j::Int) = _lₖⱼ(f, k, j, δt, δτ)
    lⱼₖ(j::Int, k::Int) = _fermionic_lⱼₖ(f, β, μ, j, k, δt, δτ)
    lₖⱼ(k::Int, j::Int) = _fermionic_lₖⱼ(f, β, μ, k, j, δt, δτ)

    N, M = Nt, Nτ

    # real time
    ηⱼₖ = zeros(ComplexF64, Nt+1)
    ηₖⱼ = zeros(ComplexF64, Nt+1)

    ηⱼₖ[1] = quadgkwrapper(-fⱼⱼ * g₁)
    for i = 1:Nt
        ηⱼₖ[i+1] = quadgkwrapper(-fⱼₖ(i) * g₁)
    end
    ηₖⱼ[1] = quadgkwrapper(fₖₖ * g₂)
    for i = 1:Nt
        ηₖⱼ[i+1] = quadgkwrapper(fⱼₖ(-i) * g₂)
    end

    # imag time
    ξⱼₖ = zeros(ComplexF64, M) # j >= k
    ξₖⱼ = zeros(Float64, M)

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

function _fermionic_lⱼₖ(f::AbstractBoundedFunction, β, μ, j::Int, k::Int, δt, δτ)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/ε^2, exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt)
    function g(ε)
        if β >= μ
            x = exp(-safe_mult(β, ε-μ))+1
            if abs(ε) > QuAPI_tol
                -im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/(ε^2 * x)
            else
                -exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt / x
            end
        else 
            x = exp(safe_mult(β, ε-μ)) + 1
            y = exp(safe_mult(β-j*δτ, ε-μ)) * exp(-j*δτ*μ)
            if abs(ε) > QuAPI_tol
                -im*y*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/(ε^2 * x)
            else
                -y*exp(im*ε*k*δt)*δτ*δt / x
            end
        end
    end
    return f * g
end

function _fermionic_lₖⱼ(f::AbstractBoundedFunction, β, μ, k::Int, j::Int, δt, δτ)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, im*exp(ε*j*δτ)*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/ε^2, exp(ε*j*δτ)*exp(-im*ε*k*δt)*δτ*δt)
    function g(ε)
        if β >= μ
            x = exp(-safe_mult(β, ε-μ))+1
            y = exp(-safe_mult(β - j*δτ, ε-μ)) * exp(j*δτ*μ)
            if abs(ε) > QuAPI_tol
                im*y*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/(ε^2 * x)
            else
                y*exp(-im*ε*k*δt)*δτ*δt / x
            end
        else
            x = exp(safe_mult(β, ε-μ))+1
            y = exp(j*δτ * ε)

            if abs(ε) > QuAPI_tol
                im*y*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/(ε^2 * x)
            else
                y*exp(-im*ε*k*δt)*δτ*δt / x
            end
        end
    end
    return f * g
end