"""
    struct InfiniteRealCorrelationCache{B <: AbstractBath}

Incrementally grown, translationally invariant real-time hybridization
function on the closed Keldysh contour.

The Toeplitz vectors `ηⱼₖ`/`ηₖⱼ` are extended one node at a time by
`compute_next!` until the newest element falls below `atol` (or a maximum
iteration count is reached); `compute!` drives the growth and `current_size`
/`current_time` report the current number of nodes and the covered real time.
`branch(x; b1, b2)` with `b1, b2 ∈ (:+, :-)` returns the corresponding Keldysh
sub-block as a `CorrelationMatrix`.
"""
mutable struct InfiniteRealCorrelationCache{B <: AbstractBath}
	ηⱼₖ::Vector{ComplexF64}
	ηₖⱼ::Vector{ComplexF64}
	# 跨支等时的两个专属 Bose 加权积分（玻色；平移不变只消节点指标、不消分支
	# 指标——它们与 ++ 支的 J[1]=∫−fⱼⱼ/K[1]=∫−fₖₖ 是不同权重×不同三角原语的
	# 新函数，必须单独积分存储）：
	#   eqⱼₖ_pm = ∫ fⱼⱼ·g₂   （G₊₋ 等时 ηⱼₖ[1]，g₂ = boseeinstein）
	#   eqₖₖ_mp = ∫ fₖₖ·g₁   （G₋₊ 等时 ηₖⱼ[1]，g₁ = 1 + boseeinstein）
	eqⱼₖ_pm::ComplexF64
	eqₖₖ_mp::ComplexF64
	bath::B
	δt::Float64
end

"""
    current_size(x::InfiniteRealCorrelationCache)

Number of nodes of the incrementally grown Toeplitz vectors.
"""
current_size(x::InfiniteRealCorrelationCache) = length(x.ηⱼₖ)

"""
    current_time(x::InfiniteRealCorrelationCache)

Covered real-time window, `current_size(x) * x.δt`.
"""
current_time(x::InfiniteRealCorrelationCache) = current_size(x) * x.δt

InfiniteRealCorrelationCache(bath::AbstractBath, δt::Real) = InfiniteRealCorrelationCache(ComplexF64[], ComplexF64[], ComplexF64(NaN, NaN), ComplexF64(NaN, NaN), bath, convert(Float64, δt))

# 四支 Keldysh 块：与有限格 bosonic_Δt / fermionic_Δt 逐项同构（返回副本，
# 不与 cache 共享存储）。语义（J₁ = ηⱼₖ⁺⁺[1]、K₁ = ηₖⱼ⁺⁺[1]）：
#   费米（单一三角原语对 fⱼₖ/fⱼₖ′，权重 g₁/g₂ 按支挂）：
#     (+,+): (a, b)                    (+,−): (conj(b), b)
#     (−,+): (a, conj(a))              (−,−): (conj(b), conj(a))
#   玻色（跨支等时为专属积分，见 eqⱼₖ_pm/eqₖₖ_mp）：
#     (+,+): (a, b)
#     (+,−): (eqⱼₖ_pm, −b)             （间隔 = conj(b[2:end]) / −b[2:end]）
#     (−,+): (−a, eqₖₖ_mp)             （间隔 = −a[2:end] / −conj(a[2:end])）
#     (−,−): (conj(b), conj(a))
function branch(x::InfiniteRealCorrelationCache, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    return _branch(x, Val(b1), Val(b2))
end

function _branch(x::InfiniteRealCorrelationCache{<:AbstractFermionicNormalBath},
                 ::Val{b1}, ::Val{b2}) where {b1, b2}
    a, b = x.ηⱼₖ, x.ηₖⱼ
    if (b1, b2) == (:+, :+)
        return CorrelationMatrix(copy(a), copy(b))
    elseif (b1, b2) == (:+, :-)
        return CorrelationMatrix(conj.(b), copy(b))
    elseif (b1, b2) == (:-, :+)
        return CorrelationMatrix(copy(a), conj.(a))
    else
        return CorrelationMatrix(conj.(b), conj.(a))
    end
end

function _branch(x::InfiniteRealCorrelationCache{<:AbstractBosonicNormalBath},
                 ::Val{b1}, ::Val{b2}) where {b1, b2}
    (isnan(real(x.eqⱼₖ_pm)) || isnan(real(x.eqₖₖ_mp))) &&
        error("跨支等时积分未计算：先调用 compute_next!/compute!")
    a, b = x.ηⱼₖ, x.ηₖⱼ
    ηⱼₖ = Vector{ComplexF64}(undef, length(a))
    ηₖⱼ = Vector{ComplexF64}(undef, length(a))
    if (b1, b2) == (:+, :+)
        ηⱼₖ .= a; ηₖⱼ .= b
    elseif (b1, b2) == (:+, :-)
        ηⱼₖ[1] = x.eqⱼₖ_pm;          ηⱼₖ[2:end] .= .-conj.(b[2:end])
        ηₖⱼ .= .-b
    elseif (b1, b2) == (:-, :+)
        ηⱼₖ .= .-a
        ηₖⱼ[1] = x.eqₖₖ_mp;          ηₖⱼ[2:end] .= .-conj.(a[2:end])
    else
        ηⱼₖ[1] = conj(b[1]);         ηⱼₖ[2:end] .= conj.(b[2:end])
        ηₖⱼ[1] = conj(a[1]);         ηₖⱼ[2:end] .= conj.(a[2:end])
    end
    return CorrelationMatrix(ηⱼₖ, ηₖⱼ)
end

"""
    infinite_Δt(bath::AbstractBath; δt::Real, atol::Real=1.0e-6, t::Real=200)

Translationally invariant real-time hybridization function on the closed
Keldysh contour, discretized on a grid of step `δt`.

The Toeplitz vectors are grown node by node until the newest element drops
below `atol`; `t` only sets the maximum real-time window
(`maxiter = round(Int, t/δt)`) and a warning is emitted if the correlation
does not decay to `atol` within it.

Returns an `InfiniteRealCorrelationCache`, whose `ηⱼₖ`/`ηₖⱼ` vectors can be
combined into the Keldysh branch blocks via `branch(x; b1, b2)` once the decay
criterion is met.

Arguments:
- `bath`:  a fermionic or bosonic bath (e.g. a `NormalBath` or `Vacuum`);
- `δt`:  the real-time step size;
- `atol`:  absolute tolerance for the newest `η` element (default 1.0e-6);
- `t`:  maximum real-time window for the growth (default 200).
"""
function infinite_Δt(bath::AbstractBath; δt::Real, kwargs...)
	corr = InfiniteRealCorrelationCache(bath, δt)
	compute!(corr; kwargs...)
	return corr
end

"""
    compute!(x::InfiniteRealCorrelationCache; atol::Real=1.0e-6, t::Real=200)

Grow the Toeplitz vectors by repeated `compute_next!` until the newest
element falls below `atol` or the window `t` is exhausted.  Returns the
number of nodes.
"""
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

"""
    compute_next!(x::InfiniteRealCorrelationCache)

Append one node to the Toeplitz vectors.
"""
function compute_next!(x::InfiniteRealCorrelationCache{<:AbstractFermionicNormalBath})
	f = x.bath.f
	# f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
	β, μ, δt = x.bath.β, x.bath.μ, x.δt

    g₁(ε) = fermidirac(β, μ - ε); g₂(ε) = fermidirac(β, ε - μ)
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

# function compute_next!(x::InfiniteRealCorrelationCache{<:AbstractBosonicNormalBath})
# 	f = x.bath.f
# 	# f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
# 	β, μ, δt = x.bath.β, x.bath.μ, x.δt

#     g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
#     # real time
#     fⱼₖ(Δk) = _fⱼₖ_r(f, Δk, δt)
#     fⱼⱼ = _fⱼⱼ_r(f, δt)
#     # fₖₖ = _fₖₖ_r(f, δt)

#     # j >= k
#     N = current_size(x)
#     if N == 0
#     	ηⱼₖ_new = quadgkwrapper(-fⱼⱼ * g₁)
#     	ηₖⱼ_new = quadgkwrapper(-fⱼⱼ' * g₂)
#     else
#     	ηⱼₖ_new = quadgkwrapper(-fⱼₖ(N) * g₁)
#     	ηₖⱼ_new = quadgkwrapper(-fⱼₖ(N)' * g₂)
#     end
#     push!(x.ηⱼₖ, ηⱼₖ_new)
#     push!(x.ηₖⱼ, ηₖⱼ_new)
#     return x
# end

function compute_next!(x::InfiniteRealCorrelationCache{<:AbstractBosonicNormalBath})
	f = x.bath.f
	# f, lb, ub = f0.f, lowerbound(f0), upperbound(f0)
	β, μ, δt = x.bath.β, x.bath.μ, x.δt

    # g₁(ε) = _g₁(β, μ, ε); g₂(ε) = _g₂(β, μ, ε)
    # real time
    fⱼₖ(Δk) = _bosonic_fⱼₖ_r(f, β, μ, Δk, δt)
    fₖⱼ(Δk) = _bosonic_fₖⱼ_r(f, β, μ, Δk, δt)
    fⱼⱼ = _bosonic_fⱼⱼ_r(f, β, μ, δt)
    fₖₖ = _bosonic_fₖₖ_r(f, β, μ, δt)
    # fₖₖ = _fₖₖ_r(f, δt)

    # j >= k
    N = current_size(x)
    if N == 0
    	ηⱼₖ_new = quadgkwrapper(-fⱼⱼ)
    	ηₖⱼ_new = quadgkwrapper(-fₖₖ)
    	# 跨支等时的两个专属积分（只依赖 δt，N==0 时算一次；与有限格
    	# bosonic_Δt 第 116–117 行的 fⱼⱼg₂/fₖₖg₁ 同定义）
    	x.eqⱼₖ_pm = quadgkwrapper(f * (ε -> boseeinstein(β, μ, ε) * _gⱼⱼ_r(δt)(ε)))
    	x.eqₖₖ_mp = quadgkwrapper(f * (ε -> (1 + boseeinstein(β, μ, ε)) * _gₖₖ_r(δt)(ε)))
    else
    	ηⱼₖ_new = quadgkwrapper(-fⱼₖ(N))
    	ηₖⱼ_new = quadgkwrapper(-fₖⱼ(-N))
    end
    push!(x.ηⱼₖ, ηⱼₖ_new)
    push!(x.ηₖⱼ, ηₖⱼ_new)
    return x
end
