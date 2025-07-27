Δm(bath::AbstractBCSBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = bcs_Δm(bath.spectrum, β=bath.β, μ=bath.μ, Δ=bath.Δ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)

function bcs_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int, Δ::Number=0, μ::Real=0, δτ::Real=β/Nτ)
	fu2, fuv, fvu, fv2 = get_all_fs(f, Δ)
	Δτ_uu = fermionic_Δm(fu2, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)
	Δτ_uv = fermionic_Δm(fuv, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)
	Δτ_vu = fermionic_Δm(fvu, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)
	Δτ_vv = fermionic_Δm(fv2, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)
	return _bcs_corr(Δτ_uu, Δτ_uv, Δτ_vu, Δτ_vv)
end