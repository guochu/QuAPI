
_dispersion(ϵ::Real, Δ::Number) = sqrt(ϵ^2 + abs2(Δ))

function _u(ϵ::Real, Δ::Number)
	(Δ == zero(Δ)) && return (ϵ >= 0) ? one(ϵ) : zero(ϵ)
	return sqrt((1 + ϵ/_dispersion(ϵ, Δ))/2)
end 
function _v(ϵ::Real, Δ::Real)
	(Δ == zero(Δ)) && return (ϵ >= 0) ? zero(ϵ) : one(ϵ)
	return sqrt((1 - ϵ/_dispersion(ϵ, Δ))/2)
end 

function _v(ϵ::Real, Δ::Complex)
	(Δ == zero(Δ)) && return (ϵ >= 0) ? zero(ϵ) : one(ϵ)
	phi = angle(Δ)
	return sqrt((1 - ϵ/_dispersion(ϵ, Δ))/2) * exp(im*phi)
end 


function _u2(ϵ, Δ) 
	x = _u(ϵ, Δ)
	return conj(x) * x
end
function _v2(ϵ, Δ)
	x = _v(ϵ, Δ)
	return conj(x) * x
end 
function _uv(ϵ, Δ)
	return conj(_u(ϵ, Δ)) * _v(ϵ, Δ)
end 
function _vu(ϵ, Δ)
	return conj(_v(ϵ, Δ)) * _u(ϵ, Δ)
end

struct BCSCorrelationFunction{C <: AbstractCorrelationFunction}
	cc::C
	cn::C
	nc::C
	nn::C
end

function Base.getindex(x::BCSCorrelationFunction, i::Int, j::Int)
	(1 <= i <= 2) || throw(BoundsError(1:2, i))
	(1 <= j <= 2) || throw(BoundsError(1:2, j))
	return ifelse(i==1, ifelse(j==1, x.cc, x.cn), ifelse(j==1, x.nc, x.nn))
end
Base.getindex(x::BCSCorrelationFunction, c1::Bool, c2::Bool) = ifelse(c1, ifelse(c2, x.cc, x.cn), ifelse(c2, x.nc, x.nn))



_mult_f(x::AbstractBoundedFunction, f) = similar(x, w -> f(w) * x.f(w)) 
_mult_f(x::DiracDelta, f) = x * f

function get_all_fs(f::AbstractBoundedFunction, Δ::Number)
	__u2(ϵ) = _u2(ϵ, Δ)
	__v2(ϵ) = _v2(ϵ, Δ)
	__uv(ϵ) = _uv(ϵ, Δ)
	__vu(ϵ) = _vu(ϵ, Δ)
	return _mult_f(f, __u2), _mult_f(f, __uv), _mult_f(f, __vu), _mult_f(f, __v2)
end
