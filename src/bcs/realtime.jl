Δt(bath::AbstractBCSBath; N::Int, t::Real) = bcs_Δt(bath.spectrum, β=bath.β, μ=bath.μ, Δ=bath.Δ, N=N, t=t)


function bcs_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, Δ::Number=0, μ::Real=0) 
	fu2, fuv, fvu, fv2 = get_all_fs(f, Δ)
	disperse = ϵ -> _dispersion(ϵ, Δ)
	Δt_uu = fermionic_Δt(fu2, disperse, β=β, N=N, t=t, μ=μ)
	Δt_uv = fermionic_Δt(fuv, disperse, β=β, N=N, t=t, μ=μ)
	Δt_vu = fermionic_Δt(fvu, disperse, β=β, N=N, t=t, μ=μ)
	Δt_vv = fermionic_Δt(fv2, disperse, β=β, N=N, t=t, μ=μ)
	return _bcs_corr(Δt_uu, Δt_uv, Δt_vu, Δt_vv)
end
