Δt(bath::AbstractBCSBath; N::Int, t::Real) = bcs_Δt(bath.spectrum, β=bath.β, μ=bath.μ, Δ=bath.Δ, N=N, t=t)


"""
    bcs_Δt(bath::AbstractBCSBath; N::Int, t::Real)
    bcs_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, Δ::Number=0, μ::Real=0)

Real-time hybridization function of a BCS (superconducting) bath on the closed
Keldysh contour covering total time t with `N` steps δt = t/N.

Same construction as `bcs_Δτ`: the coherence factors (u², uv, vu, v²) are
applied to the normal spectrum and `fermionic_Δt` is evaluated with the
dispersion `ε -> sqrt(ε² + |Δ|²)`.  Returns a `BCSCorrelationFunction` with
the normal blocks `uu`/`dd` and the anomalous blocks `ud`/`du`.

Arguments:
- `f`:  the normal-state bath spectrum A(ε) as an `AbstractBoundedFunction`;
- `β`:  the inverse temperature;
- `N`:  number of real-time steps, δt = t/N;
- `t`:  the total evolution time;
- `Δ`:  the pairing gap (real or complex), default 0;
- `μ`:  the chemical potential.

The `bath` convenience form reads `β`, `μ` and `Δ` from the bath object.
"""
function bcs_Δt(f::AbstractBoundedFunction; β::Real, N::Int, t::Real, Δ::Number=0, μ::Real=0) 
	fu2, fuv, fvu, fv2 = get_all_fs(f, Δ)
	disperse = ϵ -> _dispersion(ϵ, Δ)
	Δt_uu = fermionic_Δt(fu2, disperse, β=β, N=N, t=t, μ=μ)
	Δt_uv = fermionic_Δt(fuv, disperse, β=β, N=N, t=t, μ=μ)
	Δt_vu = fermionic_Δt(fvu, disperse, β=β, N=N, t=t, μ=μ)
	Δt_vv = fermionic_Δt(fv2, disperse, β=β, N=N, t=t, μ=μ)
	return _bcs_corr(Δt_uu, Δt_uv, Δt_vu, Δt_vv)
end
