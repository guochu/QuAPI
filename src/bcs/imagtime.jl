Δτ(bath::AbstractBCSBath; N::Int, δτ::Real=bath.β/N) = bcs_Δτ(bath.spectrum, β=bath.β, N=N, μ=bath.μ, Δ=bath.Δ, δτ=δτ)


function bcs_Δτ(f::AbstractBoundedFunction; β::Real, N::Int, μ::Real=0, Δ::Number=0, δτ::Real=β/N) 
	fu2, fuv, fvu, fv2 = get_all_fs(f, Δ)
	Δτ_uu = fermionic_Δτ(fu2, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_uv = fermionic_Δτ(fuv, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_vu = fermionic_Δτ(fvu, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_vv = fermionic_Δτ(fv2, β=β, N=N, μ=μ, δτ=δτ)
	return _bcs_corr(Δτ_uu, Δτ_uv, Δτ_vu, Δτ_vv)
end


function _bcs_corr(Δ_uu, Δ_uv, Δ_vu, Δ_vv)
	cc = Δ_uv + transpose(Δ_uv)
	nn = Δ_vu + transpose(Δ_vu)
	cn = Δ_uu - transpose(Δ_vv)
	nc = Δ_vv - transpose(Δ_uu)
	return 	BCSCorrelationFunction(cc, cn, nc, nn)
end