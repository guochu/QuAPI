struct InfiniteImagCorrelationCache{B <: Vacuum}
	ηⱼₖ::Vector{Float64}
	ηₖⱼ::Vector{Float64}
	bath::B
	δτ::Float64
end


current_size(x::InfiniteImagCorrelationCache) = length(x.ηⱼₖ)
current_temperature(x::InfiniteImagCorrelationCache) = current_size(x) * x.δτ

InfiniteImagCorrelationCache(bath::Vacuum, δτ::Real) = InfiniteImagCorrelationCache(Float64[], Float64[], bath, convert(Float64, δτ))

"""
    infinite_Δτ(bath::Vacuum; δt::Real, atol::Real, β::Real)

Translationally invariant zero-temperature imaginary time hybridization function
"""
function infinite_Δτ(bath::Vacuum; δτ::Real, kwargs...)
    corr = InfiniteImagCorrelationCache(bath, δτ)
    compute!(corr; kwargs...)
    return corr
end

function compute!(x::InfiniteImagCorrelationCache; atol::Real=1.0e-6, β::Real=200)
	maxiter = round(Int, β / x.δτ)
	atol = convert(Float64, atol)
	while current_size(x) <= maxiter
		if (current_size(x) > 1) && (max(abs(x.ηⱼₖ[end]), abs(x.ηₖⱼ[end])) <= atol)
			break
		end
		compute_next!(x)
	end
	(current_size(x) > maxiter) && @warn "correlation can not decay to less than $atol for β=$(β), error=$(max(abs(x.ηⱼₖ[end]), abs(x.ηₖⱼ[end])))"
	return current_size(x)
end


function compute_next!(x::InfiniteImagCorrelationCache{<:FermionicVacuum})
	f0 = x.bath.f
	β, μ = x.bath.β, x.bath.μ
    f = spectrumshift(f0, μ)
    δτ = x.δτ

    # g₁(ϵ) = _f₁(β, 0., ϵ)
    # g₂(ϵ) = _f₂(β, 0., ϵ)
    # fⱼₖ(Δk::Int) = _fⱼₖ_i(f, Δk, δτ)
    # fₖⱼ(Δk::Int) = _fₖⱼ_i(f, Δk, δτ)
    # fⱼⱼ = _fⱼⱼ_i(f, δτ)
    # fₖₖ = _fₖₖ_i(f, δτ)

    fⱼₖ(Δk::Int) = _fermionic_fⱼₖ_i(f, β, Δk, δτ)
    fₖⱼ(Δk::Int) = _fermionic_fₖⱼ_i(f, β, Δk, δτ)
    fⱼⱼ = _fermionic_fⱼⱼ_i(f, β, δτ)
    fₖₖ = _fermionic_fₖₖ_i(f, β, δτ)


    # j >= k
    N = current_size(x)
    if N == 0
    	ηⱼₖ_new = quadgkwrapper(fⱼⱼ)
    	ηₖⱼ_new = quadgkwrapper(fₖₖ)
    else
    	ηⱼₖ_new = quadgkwrapper(fⱼₖ(N))
    	ηₖⱼ_new = quadgkwrapper(fₖⱼ(N))
    end
    push!(x.ηⱼₖ, ηⱼₖ_new)
    push!(x.ηₖⱼ, ηₖⱼ_new)
    return x
end

function compute_next!(x::InfiniteImagCorrelationCache{<:BosonicVacuum})
	f0 = x.bath.f
	β, μ = x.bath.β, x.bath.μ
    f = spectrumshift(f0, μ)
    δτ = x.δτ

    # g₁(ϵ) = _f₁(β, 0., ϵ)
    # g₂(ϵ) = _f₂(β, 0., ϵ)
    # fⱼₖ(Δk::Int) = _fⱼₖ_i(f, Δk, δτ)
    # fₖⱼ(Δk::Int) = _fₖⱼ_i(f, Δk, δτ)
    # fⱼⱼ = _fⱼⱼ_i(f, δτ)
    # fₖₖ = _fₖₖ_i(f, δτ)

    fⱼₖ(Δk::Int) = _bosonic_fⱼₖ_i(f, β, Δk, δτ)
    fₖⱼ(Δk::Int) = _bosonic_fₖⱼ_i(f, β, Δk, δτ)
    fⱼⱼ = _bosonic_fⱼⱼ_i(f, β, δτ)
    fₖₖ = _bosonic_fₖₖ_i(f, β, δτ)


    # j >= k
    N = current_size(x)
    if N == 0
    	ηⱼₖ_new = quadgkwrapper(fⱼⱼ)
    	ηₖⱼ_new = quadgkwrapper(fₖₖ)
    else
    	ηⱼₖ_new = quadgkwrapper(fⱼₖ(N))
    	ηₖⱼ_new = quadgkwrapper(fₖⱼ(N))
    end
    push!(x.ηⱼₖ, ηⱼₖ_new)
    push!(x.ηₖⱼ, ηₖⱼ_new)
    return x
end
