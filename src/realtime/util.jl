# function _fⱼₖ_r(f::AbstractBoundedFunction, Δk::Int, δt)
#     # g(ε) = ifelse(abs(ε) > QuAPI_tol, 2exp(-im*ε*Δk*δt)*(1-cos(ε*δt))/ε^2, exp(-im*ε*Δk*δt)*δt^2)
#     function g(ε)
#         if abs(ε) > QuAPI_tol
#             2exp(-im*ε*Δk*δt)*(1-cos(ε*δt))/ε^2
#         else
#             exp(-im*ε*Δk*δt)*δt^2
#         end
#     end
#     return f * g
# end

# function _fⱼⱼ_r(f::AbstractBoundedFunction, δt)
#     # g(ε) = ifelse(abs(ε) > QuAPI_tol, ((1-im*ε*δt)-exp(-im*ε*δt))/ε^2, 0.5*δt^2)
#     function g(ε)
#         if abs(ε) > QuAPI_tol
#             ((1-im*ε*δt)-exp(-im*ε*δt))/ε^2
#         else
#             0.5*δt^2
#         end
#     end
#     return f * g
# end

# function _fₖₖ_r(f::AbstractBoundedFunction, δt)
#     # g(ε) = ifelse(abs(ε) > QuAPI_tol, ((1+im*ε*δt)-exp(im*ε*δt))/ε^2, 0.5*δt^2)
#     function g(ε)
#         if abs(ε) > QuAPI_tol
#             ((1+im*ε*δt)-exp(im*ε*δt))/ε^2
#         else
#             0.5*δt^2
#         end
#     end
#     return f * g
# end


function _fⱼₖ_r_disperse(f::AbstractBoundedFunction, disperse, Δk::Int, δt)
    g = _gⱼₖ_r(Δk, δt)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′  
end
_fⱼₖ_r(f::AbstractBoundedFunction, Δk::Int, δt) = f * _gⱼₖ_r(Δk, δt)
function _gⱼₖ_r(Δk::Int, δt)
    return ε -> begin
        if abs(ε) > QuAPI_tol
            2exp(-im*ε*Δk*δt)*(1-cos(ε*δt))/ε^2
        else
            exp(-im*ε*Δk*δt)*δt^2
        end
    end
end

function _fⱼⱼ_r_disperse(f::AbstractBoundedFunction, disperse, δt)
    g = _gⱼⱼ_r(δt)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′  
end
_fⱼⱼ_r(f::AbstractBoundedFunction, δt) = f * _gⱼⱼ_r(δt)
function _gⱼⱼ_r(δt)
    return ε -> begin
        if abs(ε) > QuAPI_tol
            ((1-im*ε*δt)-exp(-im*ε*δt))/ε^2
        else
            0.5*δt^2
        end
    end
end

function _fₖₖ_r_disperse(f::AbstractBoundedFunction, disperse, δt)
    g = _gₖₖ_r(δt)
    g′ = ϵ -> g(disperse(ϵ))
    return f * g′      
end
_fₖₖ_r(f::AbstractBoundedFunction, δt) = f * _gₖₖ_r(δt)
function _gₖₖ_r(δt)
    return ε -> begin
        if abs(ε) > QuAPI_tol
            ((1+im*ε*δt)-exp(im*ε*δt))/ε^2
        else
            0.5*δt^2
        end
    end
end