Δτ(bath::AbstractFermionicBath; N::Int, δτ::Real=bath.β/N) = fermionic_Δτ(bath.spectrum, β=bath.β, N=N, μ=bath.μ, δτ=δτ)
fermionic_Δτ(f::AbstractSpectrumFunction; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) = fermionic_Δτ(f, β, N, μ, δτ)



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
function fermionic_Δτ(f0::AbstractSpectrumFunction, β::Real, N::Int, μ::Real, δτ::Real=β / N)
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

# # from now on j > k
# function _fⱼₖ_i(f, Δk::Int, ε::Float64, δτ)
#     if (abs(ε) > tol)
#         -2f(ε)/ε^2*exp(-Δk*δτ*ε)*(1-cosh(δτ*ε))
#     else
#         f(ε)*exp(-Δk*δτ*ε)*δτ^2
#     end
# end

# function _fₖⱼ_i(f, Δk::Int, ε::Float64, δτ)
#     if (abs(ε) > tol)
#         -2f(ε)/ε^2*exp(Δk*δτ*ε)*(1-cosh(δτ*ε))
#     else
#         f(ε)*exp(Δk*δτ*ε)*δτ^2
#     end
# end

# function _fⱼⱼ_i(f, ε::Float64, δτ)
#     if (abs(ε) > tol)
#         f(ε)/ε^2*(exp(-δτ*ε)-(1-δτ*ε))
#     else
#         0.5f(ε)*δτ^2
#     end
# end

# function _fₖₖ_i(f, ε::Float64, δτ)
#     if (abs(ε) > tol)
#         f(ε)/ε^2*(exp(δτ*ε)-(1+δτ*ε))
#     else
#         0.5*f(ε)*δτ^2/2
#     end
# end

function _fⱼₖ_i(f::AbstractBoundedFunction, Δk::Int, δτ)
    g(ε) = ifelse(abs(ε) > tol, -2exp(-Δk*δτ*ε)*(1-cosh(δτ*ε))/ε^2, exp(-Δk*δτ*ε)*δτ^2)
    return f * g
end

function _fₖⱼ_i(f::AbstractBoundedFunction, Δk::Int, δτ)
    g(ε) = ifelse(abs(ε) > tol, -2exp(Δk*δτ*ε)*(1-cosh(δτ*ε))/ε^2, exp(Δk*δτ*ε)*δτ^2)
    return f * g
end

function _fⱼⱼ_i(f::AbstractBoundedFunction, δτ)
    g(ε) = ifelse(abs(ε) > tol, (exp(-δτ*ε)-(1-δτ*ε))/ε^2, 0.5*δτ^2)
    return f * g
end

function _fₖₖ_i(f::AbstractBoundedFunction, δτ)
    g(ε) = ifelse(abs(ε) > tol, (exp(δτ*ε)-(1+δτ*ε))/ε^2, 0.5*δτ^2/2)
    return f * g
end

function Δiw_to_Δτ(Δiw::AbstractVector{<:Number}; β::Real, N::Int)
    iseven(length(Δiw)) || throw(ArgumentError("even number of frequencies expected"))
    δτ = β/N
    nmax = div(length(Δiw), 2) - 1

    ηⱼₖ = zeros(ComplexF64, N)
    ηₖⱼ = zeros(ComplexF64, N)

    for (i, ωₙ) in enumerate(ifrequencies(β, nmax)) 
        a = -(2/(β*ωₙ^2)) * Δiw[i] * (1 - cos(ωₙ * δτ))
        ηⱼₖ[1] += a / 2
        ηₖⱼ[1] += a / 2
        for k in 1:N-1
            ηⱼₖ[k+1] += a * exp(-im*ωₙ*k*δτ)
            ηₖⱼ[k+1] += a * exp(im*ωₙ*k*δτ)
        end
    end
    for i in 1:N
        (abs(imag(ηⱼₖ[i])) < 1.0e-8) || error("imag part of ηⱼₖ[$i] is nonzero")
        (abs(imag(ηₖⱼ[i])) < 1.0e-8) || error("imag part of ηₖⱼ[$i] is nonzero")
    end

    ImagCorrelationFunction(CorrelationMatrix{Float64}(real(ηⱼₖ), real(ηₖⱼ)))
end