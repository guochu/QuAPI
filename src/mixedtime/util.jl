function _lⱼₖ(f::AbstractBoundedFunction, j::Int, k::Int, δt, δτ)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/ε^2, exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt)
    function g(ε)
        if abs(ε) > QuAPI_tol
            im*exp(-ε*j*δτ)*exp(im*ε*k*δt)*(exp(-ε*δτ)-1)*(exp(im*ε*δt)-1)/ε^2
        else
            exp(-ε*j*δτ)*exp(im*ε*k*δt)*δτ*δt
        end
    end
    return f * g
end

function _lₖⱼ(f::AbstractBoundedFunction, k::Int, j::Int, δt, δτ)
    # g(ε) = ifelse(abs(ε) > QuAPI_tol, im*exp(ε*j*δτ)*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/ε^2, exp(ε*j*δτ)*exp(-im*ε*k*δt)*δτ*δt)
    function g(ε)
        if abs(ε) > QuAPI_tol
            im*exp(ε*j*δτ)*exp(-im*ε*k*δt)*(exp(ε*δτ)-1)*(exp(-im*ε*δt)-1)/ε^2
        else
            exp(ε*j*δτ)*exp(-im*ε*k*δt)*δτ*δt
        end
    end
    return f * g
end