bcs_Δτ(bath::AbstractFermionicBath, Δ::Real; N::Int, δτ::Real=bath.β/N) = bcs_Δτ(bath.spectrum, Δ, β=bath.β, N=N, μ=bath.μ, δτ=δτ)


function bcs_Δτ(f::AbstractBoundedFunction, Δ::Real; β::Real, N::Int, μ::Real=0, δτ::Real=β/N) 
	fu2, fuv, fv2 = get_all_fs(f, Δ)
	Δτ_uu = Δτ(fu2, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_uv = Δτ(fuv, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_vv = Δτ(fv2, β=β, N=N, μ=μ, δτ=δτ)
	return _bcs_corr(Δτ_uu, Δτ_uv, Δτ_vv)
end


function _bcs_corr(Δ_uu, Δ_uv, Δ_vv)
	cc = Δ_uv + transpose(Δ_uv)
	nn = cc
	cn = Δ_uu - transpose(Δ_vv)
	nc = Δ_vv - transpose(Δ_uu)
	return 	BCSCorrelationFunction(cc, cn, nc, nn)
end