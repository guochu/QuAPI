Δτ(bath::AbstractBCSBath; N::Int, δτ::Real=bath.β/N) = bcs_Δτ(bath.spectrum, β=bath.β, N=N, μ=bath.μ, Δ=bath.Δ, δτ=δτ)


"""
    bcs_Δτ(bath::AbstractBCSBath; N::Int, δτ::Real=bath.β/N)
    bcs_Δτ(f::AbstractBoundedFunction; β::Real, N::Int, μ::Real=0, Δ::Number=0, δτ::Real=β/N)

Imaginary-time hybridization function of a BCS (superconducting) bath on a
grid of `N` steps δτ = β/N.

Following `src/bcs/bcscorr.jl` the four blocks are obtained by applying the
Bogoliubov coherence factors (u², uv, vu, v²) to the normal spectrum and
evaluating `fermionic_Δτ` with the dispersion `ε -> sqrt(ε² + |Δ|²)`.  Returns
a `BCSCorrelationFunction` with the normal blocks `uu`/`dd` and the anomalous
blocks `ud`/`du`, entering the influence-functional exponent

    ā↑ Δ↑↑ a↑ + ā↑ Δ↑↓ ā↓ + a↓ Δ↓↑ a↑ + a↓ Δ↓↓ ā↓.

Arguments:
- `f`:  the normal-state bath spectrum A(ε) as an `AbstractBoundedFunction`;
- `β`:  the inverse temperature;
- `N`:  number of imaginary-time steps, δτ = β/N;
- `μ`:  the chemical potential;
- `Δ`:  the pairing gap (real or complex), default 0;
- `δτ`:  the step size (defaults to β/N).

The `bath` convenience form reads `β`, `μ` and `Δ` from the bath object.
"""
function bcs_Δτ(f::AbstractBoundedFunction; β::Real, N::Int, μ::Real=0, Δ::Number=0, δτ::Real=β/N) 
	fu2, fuv, fvu, fv2 = get_all_fs(f, Δ)
	disperse = ϵ -> _dispersion(ϵ, Δ)
	Δτ_uu = fermionic_Δτ(fu2, disperse, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_uv = fermionic_Δτ(fuv, disperse, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_vu = fermionic_Δτ(fvu, disperse, β=β, N=N, μ=μ, δτ=δτ)
	Δτ_vv = fermionic_Δτ(fv2, disperse, β=β, N=N, μ=μ, δτ=δτ)
	return _bcs_corr(Δτ_uu, Δτ_uv, Δτ_vu, Δτ_vv)
end


function _bcs_corr(Δ_uu, Δ_uv, Δ_vu, Δ_vv)
	cc = Δ_uv + transpose(Δ_uv)
	nn = Δ_vu + transpose(Δ_vu)
	cn = Δ_uu - transpose(Δ_vv)
	nc = Δ_vv - transpose(Δ_uu)
	return 	BCSCorrelationFunction(cc, cn, nc, nn)
end