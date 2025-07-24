bcs_Δt(bath::AbstractFermionicBath, Δ::Real; N::Int, t::Real) = bcs_Δt(bath.spectrum, Δ, β=bath.β, μ=bath.μ, N=N, t=t)


function bcs_Δt(f::AbstractBoundedFunction, Δ::Real; β::Real, N::Int, t::Real, μ::Real=0) 

end
