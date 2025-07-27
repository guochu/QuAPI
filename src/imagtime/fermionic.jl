Δτ(bath::AbstractFermionicBath; N::Int, δτ::Real=bath.β/N) = fermionic_Δτ(bath.spectrum, β=bath.β, N=N, μ=bath.μ, δτ=δτ)
fermionic_Δτ(f::AbstractBoundedFunction; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) = fermionic_Δτ(f, β, N, μ, δτ)



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
function fermionic_Δτ(f0::AbstractBoundedFunction, β::Real, N::Int, μ::Real, δτ::Real=β / N)
    # f′, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
    β = convert(Float64, β)
    μ = convert(Float64, μ)
    f = spectrumshift(f0, μ)
    # δτ = β / N
    
    fⱼₖ(Δk::Int) = _fermionic_fⱼₖ_i(f, β, 0, Δk, δτ)
    fₖⱼ(Δk::Int) = _fermionic_fₖⱼ_i(f, β, 0, Δk, δτ)
    fⱼⱼ = _fermionic_fⱼⱼ_i(f, β, 0, δτ)
    fₖₖ = _fermionic_fₖₖ_i(f, β, 0, δτ)

    # j >= k
    L = N
    tmp = quadgkwrapper(fⱼⱼ)
    T = typeof(tmp)
    ηⱼₖ = zeros(T, L)
    ηⱼₖ[1] = tmp
    for k = 1:L-1
        ηⱼₖ[k+1] = quadgkwrapper(fⱼₖ(k))
    end

    ηₖⱼ = zeros(T, L)
    ηₖⱼ[1] = quadgkwrapper(fₖₖ)
    for k = 1:L-1
        ηₖⱼ[k+1] = quadgkwrapper(fₖⱼ(k))
    end
    ImagCorrelationFunction(CorrelationMatrix{T}(ηⱼₖ, ηₖⱼ))
end



function _fermionic_fⱼₖ_i(f::AbstractBoundedFunction, β, μ, Δk::Int, δτ)
    function g(ε)
        if ε >= μ
            x = exp(-safe_mult(β, ε-μ))+1
            if abs(ε) > QuAPI_tol
                -2exp(-Δk*δτ*ε)*(1-cosh(δτ*ε))/(ε^2 * x)
            else
                exp(-Δk*δτ*ε)*δτ^2 / x
            end
        else 
            x = exp(safe_mult(β, ε-μ)) + 1
            y = exp(safe_mult(β-Δk*δτ, ε-μ)) * exp(-Δk*δτ*μ)
            if abs(ε) > QuAPI_tol
                -2y*(1-cosh(δτ*ε))/(ε^2 * x)
            else
                y*δτ^2 / x
            end            
        end
    end
    return f * g
end

function _fermionic_fⱼⱼ_i(f::AbstractBoundedFunction, β, μ, δτ)
    function g(ε)
        if ε >= μ
            x = exp(-safe_mult(β, ε-μ))+1
            if abs(ε) > QuAPI_tol
                (exp(-δτ*ε)-(1-δτ*ε))/(ε^2 * x)
            else
                0.5*δτ^2 / x
            end 
        else
            x = exp(safe_mult(β, ε-μ))
            if abs(ε) > QuAPI_tol
                x*(exp(-δτ*ε)-(1-δτ*ε))/(ε^2 * (x+1))
            else
                0.5*δτ^2*x / (x+1)
            end 
        end
    end
    return f * g
end

function _fermionic_fₖⱼ_i(f::AbstractBoundedFunction, β, μ, Δk::Int, δτ)
    function g(ε)
        if ε >= μ
            x = exp(-safe_mult(β, ε-μ))+1
            y = exp(-safe_mult(β - Δk*δτ, ε-μ)) * exp(Δk*δτ*μ)
            if abs(ε) > QuAPI_tol
                2y*(1-cosh(δτ*ε))/(ε^2 * x)
            else
                -y*δτ^2 / x
            end 
        else
            x = exp(safe_mult(β, ε-μ))+1
            y = exp(Δk*δτ * ε)
            if abs(ε) > QuAPI_tol
                2y*(1-cosh(δτ*ε))/(ε^2 * x)
            else
                -y*δτ^2 / x
            end 

        end
    end
    return f * g
end

function _fermionic_fₖₖ_i(f::AbstractBoundedFunction, β, μ, δτ)
    function g(ε)
        if ε >= μ
            x = exp(-safe_mult(β, ε-μ))
            if abs(ε) > QuAPI_tol
                -(exp(δτ*ε)-(1+δτ*ε)) * x / (ε^2 * (x+1))
            else
                -0.5*δτ^2 * x / 2 * (x+1)
            end 
        else
            x = exp(safe_mult(β, ε-μ))
            if abs(ε) > QuAPI_tol
                -(exp(δτ*ε)-(1+δτ*ε)) / (ε^2 * (x+1))
            else
                -0.5*δτ^2 / 2 * (x+1)
            end 
        end
    end
    return f * g
end
