bcs_Δτ(bath::AbstractFermionicBath, Δ::Real; N::Int, δτ::Real=bath.β/N) = bcs_Δτ(bath.spectrum, Δ, β=bath.β, N=N, μ=bath.μ, δτ=δτ)


function bcs_Δτ(f::AbstractBoundedFunction, Δ::Real; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) 

end


