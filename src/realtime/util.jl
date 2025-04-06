function _fⱼₖ_r(f::AbstractBoundedFunction, Δk::Int, δt)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, 2exp(-im*ε*Δk*δt)*(1-cos(ε*δt))/ε^2, exp(-im*ε*Δk*δt)*δt^2)
    function g(ε)
        if abs(ε) > QuAPI_tol
            2exp(-im*ε*Δk*δt)*(1-cos(ε*δt))/ε^2
        else
            exp(-im*ε*Δk*δt)*δt^2
        end
    end
    return f * g
end

function _fⱼⱼ_r(f::AbstractBoundedFunction, δt)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, ((1-im*ε*δt)-exp(-im*ε*δt))/ε^2, 0.5*δt^2)
    function g(ε)
        if abs(ε) > QuAPI_tol
            ((1-im*ε*δt)-exp(-im*ε*δt))/ε^2
        else
            0.5*δt^2
        end
    end
    return f * g
end

function _fₖₖ_r(f::AbstractBoundedFunction, δt)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, ((1+im*ε*δt)-exp(im*ε*δt))/ε^2, 0.5*δt^2)
    function g(ε)
        if abs(ε) > QuAPI_tol
            ((1+im*ε*δt)-exp(im*ε*δt))/ε^2
        else
            0.5*δt^2
        end
    end
    return f * g
end
