bcs_Δτ(bath::AbstractFermionicBath, Δ::Real; N::Int, δτ::Real=bath.β/N) = bcs_Δτ(bath.spectrum, Δ, β=bath.β, N=N, μ=bath.μ, δτ=δτ)


function bcs_Δτ(f::AbstractBoundedFunction, Δ::Real; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) 
	fu2, fuv, fv2 = get_all_fs(f, Δ)
	Δτ_u2 = Δτ(fu2, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_uv = Δτ(fuv, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_v2 = Δτ(fv2, β=β, N=N, μ=μ, δτ=δτ)
	c11 = Δτ_uv + transpose(Δτ_uv)
	c22 = c11
	c12 = Δτ_u2 + Δτ_v2
	c21 = Δτ_u2 + Δτ_v2
	return BCSCorrelationFunction(c11, c12, c21, c22)
end


