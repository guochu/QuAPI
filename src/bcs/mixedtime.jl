Δm(bath::AbstractBCSBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ) = bcs_Δm(bath.spectrum, β=bath.β, μ=bath.μ, Δ=bath.Δ, Nτ=Nτ, t=t, Nt=Nt, δτ=δτ)

"""
    bcs_Δm(bath::AbstractBCSBath; Nτ::Int, t::Real, Nt::Int, δτ::Real=bath.β/Nτ)
    bcs_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int, Δ::Number=0, μ::Real=0, δτ::Real=β/Nτ)

Hybridization function of a BCS (superconducting) bath on the L-shaped
(Kadanoff-Baym) contour: an imaginary-time branch of `Nτ` steps δτ = β/Nτ
followed by the closed Keldysh part of total time t with `Nt` steps δt = t/Nt.

Same construction as `bcs_Δτ`/`bcs_Δt`: the coherence factors (u², uv, vu, v²)
are applied to the normal spectrum and `fermionic_Δm` is evaluated with the
dispersion `ε -> sqrt(ε² + |Δ|²)`.  Returns a `BCSCorrelationFunction` with
the normal blocks `uu`/`dd` and the anomalous blocks `ud`/`du`.

Arguments:
- `f`:  the normal-state bath spectrum A(ε) as an `AbstractBoundedFunction`;
- `β`:  the inverse temperature;
- `Nτ`:  number of imaginary-time steps, δτ = β/Nτ;
- `t`:  the total real evolution time;
- `Nt`:  number of real-time steps, δt = t/Nt;
- `Δ`:  the pairing gap (real or complex), default 0;
- `μ`:  the chemical potential;
- `δτ`:  the imaginary-time step size (defaults to β/Nτ).

The `bath` convenience form reads `β`, `μ` and `Δ` from the bath object.
"""
function bcs_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int, Δ::Number=0, μ::Real=0, δτ::Real=β/Nτ)
	fu2, fuv, fvu, fv2 = get_all_fs(f, Δ)
	disperse = ϵ -> _dispersion(ϵ, Δ)
	Δτ_uu = fermionic_Δm(fu2, disperse, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)
	Δτ_uv = fermionic_Δm(fuv, disperse, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)
	Δτ_vu = fermionic_Δm(fvu, disperse, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)
	Δτ_vv = fermionic_Δm(fv2, disperse, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, δτ=δτ)
	return _bcs_corr(Δτ_uu, Δτ_uv, Δτ_vu, Δτ_vv)
end