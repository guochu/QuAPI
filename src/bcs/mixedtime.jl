bcs_Δm(bath::AbstractFermionicBath, Δ::Real; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = bcs_Δm(bath.spectrum, Δ, β=bath.β, μ=bath.μ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)

function bcs_Δm(f::AbstractBoundedFunction, Δ::Real; β::Real, Nτ::Int, t::Real, Nt::Int, μ::Real=0, δτ::Real=β/Nτ)
	
end