function _fⱼₖ_i(f::AbstractBoundedFunction, Δk::Int, δτ)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, -2exp(-Δk*δτ*ε)*(1-cosh(δτ*ε))/ε^2, exp(-Δk*δτ*ε)*δτ^2)
    function g(ε)
        if abs(ε) > QuAPI_tol
            -2exp(-Δk*δτ*ε)*(1-cosh(δτ*ε))/ε^2
        else
            exp(-Δk*δτ*ε)*δτ^2
        end
    end
    return f * g
end


# this function is unstable for large β!!!
function _fₖⱼ_i(f::AbstractBoundedFunction, Δk::Int, δτ)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, -2exp(Δk*δτ*ε)*(1-cosh(δτ*ε))/ε^2, exp(Δk*δτ*ε)*δτ^2)
    function g(ε)
        if abs(ε) > QuAPI_tol
            -2exp(Δk*δτ*ε)*(1-cosh(δτ*ε))/ε^2
        else
            exp(Δk*δτ*ε)*δτ^2
        end
    end
    return f * g
end

function _fⱼⱼ_i(f::AbstractBoundedFunction, δτ)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, (exp(-δτ*ε)-(1-δτ*ε))/ε^2, 0.5*δτ^2)
    function g(ε)
        if abs(ε) > QuAPI_tol
            (exp(-δτ*ε)-(1-δτ*ε))/ε^2
        else
            0.5*δτ^2
        end
    end
    return f * g
end

function _fₖₖ_i(f::AbstractBoundedFunction, δτ)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, (exp(δτ*ε)-(1+δτ*ε))/ε^2, 0.5*δτ^2/2)
    function g(ε)
        if abs(ε) > QuAPI_tol
            (exp(δτ*ε)-(1+δτ*ε))/ε^2
        else
            0.5*δτ^2/2
        end
    end
    return f * g
end

function safe_mult(β::Real, ϵ::Real)
    if (β == Inf) && (ϵ == 0)
        return zero(ϵ)
    end
    return β * ϵ
end