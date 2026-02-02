struct InfiniteRealCorrelationCache{B <: AbstractBath}
	ηⱼₖ::Vector{ComplexF64}
	ηₖⱼ::Vector{ComplexF64}
	bath::B
	δt::Float64
end

current_size(x::InfiniteRealCorrelationCache) = length(x.ηⱼₖ)
current_time(x::InfiniteRealCorrelationCache) = current_size(x) * x.δt

InfiniteRealCorrelationCache(bath::AbstractBath, δt::Real) = InfiniteRealCorrelationCache(ComplexF64[], ComplexF64[], bath, convert(Float64, δt))

function branch(x::InfiniteRealCorrelationCache, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    local ηⱼₖ, ηₖⱼ
	if b1 == :+
		if b2 == :+
			ηⱼₖ, ηₖⱼ = x.ηⱼₖ, x.ηₖⱼ
		else
			ηⱼₖ = conj(x.ηₖⱼ)
			ηₖⱼ = x.ηₖⱼ
			ηⱼₖ[1] = x.ηₖⱼ[1]
			ηₖⱼ[1] = x.ηₖⱼ[1]
		end
	else
		if b2 == :+
			ηⱼₖ = x.ηⱼₖ
			ηₖⱼ = conj(x.ηⱼₖ)
			ηⱼₖ[1] = ηₖⱼ[1] = x.ηⱼₖ[1]
		else
			ηⱼₖ = conj(x.ηₖⱼ)
			ηₖⱼ = conj(x.ηⱼₖ)
			ηⱼₖ[1] = x.ηₖⱼ[1]
			ηₖⱼ[1] = x.ηⱼₖ[1]
		end
	end	
	return CorrelationMatrix(ηⱼₖ, ηₖⱼ)
end

"""
	infinite_Δt(bath::AbstractBath; δt::Real, atol::Real, t::Real)

Translationally invariant real time hybridization function
"""
function infinite_Δt(bath::AbstractBath; δt::Real, kwargs...)
	corr = InfiniteRealCorrelationCache(bath, δt)
	compute!(corr; kwargs...)
	return corr
end

function compute!(x::InfiniteRealCorrelationCache; atol::Real=1.0e-6, t::Real=200)
	maxiter = round(Int, t / x.δt)
	atol = convert(Float64, atol)
	while current_size(x) <= maxiter
		if (current_size(x) > 1) && (max(abs(x.ηⱼₖ[end]), abs(x.ηₖⱼ[end])) <= atol)
			break
		end
		compute_next!(x)
	end
	(current_size(x) > maxiter) && @warn "correlation can not decay to less than $atol for t=$(t), error=$(max(abs(x.ηⱼₖ[end]), abs(x.ηₖⱼ[end])))"
	return current_size(x)
end

function compute_next!(x::InfiniteRealCorrelationCache{<:AbstractFermionicNormalBath})
	f = x.bath.f
	# f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
	β, μ, δt = x.bath.β, x.bath.μ, x.δt

    g₁(ε) = _f₁(β, μ, ε); g₂(ε) = _f₂(β, μ, ε)
    # real time
    fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r(f, δt)
    fₖₖ = _fₖₖ_r(f, δt)

    # j >= k
    N = current_size(x)
    if N == 0
    	ηⱼₖ_new = quadgkwrapper(-fⱼⱼ * g₁)
    	ηₖⱼ_new = quadgkwrapper(fₖₖ * g₂)
    else
    	ηⱼₖ_new = quadgkwrapper(-fⱼₖ(N) * g₁)
    	ηₖⱼ_new = quadgkwrapper(fⱼₖ(-N) * g₂)
    end
    push!(x.ηⱼₖ, ηⱼₖ_new)
    push!(x.ηₖⱼ, ηₖⱼ_new)
    return x
end

function compute_next!(x::InfiniteRealCorrelationCache{<:AbstractBosonicNormalBath})
	f = x.bath.f
	# f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
	β, μ, δt = x.bath.β, x.bath.μ, x.δt

    g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
    # real time
    fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
    fⱼⱼ = _fⱼⱼ_r(f, δt)
    # fₖₖ = _fₖₖ_r(f, δt)

    # j >= k
    N = current_size(x)
    if N == 0
    	ηⱼₖ_new = quadgkwrapper(fⱼⱼ * g₁)
    	ηₖⱼ_new = quadgkwrapper(fⱼⱼ' * g₂)
    else
    	ηⱼₖ_new = quadgkwrapper(fⱼₖ(N) * g₁)
    	ηₖⱼ_new = quadgkwrapper(fⱼₖ(N)' * g₂)
    end
    push!(x.ηⱼₖ, ηⱼₖ_new)
    push!(x.ηₖⱼ, ηₖⱼ_new)
    return x
end