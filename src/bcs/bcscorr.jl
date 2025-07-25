
_u(ϵ, Δ) = (1 + ϵ/sqrt(ϵ^2 + Δ^2))/2
_v(ϵ, Δ) = (1 - ϵ/sqrt(ϵ^2 + Δ^2))/2

_u2(ϵ, Δ) = _u(ϵ, Δ)^2
_v2(ϵ, Δ) = _v(ϵ, Δ)^2
_uv(ϵ, Δ) = _u(ϵ, Δ) * _v(ϵ, Δ)

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

function get_all_fs(f::AbstractBoundedFunction, Δ::Real)
	__u2(ϵ) = _u2(ϵ, Δ)
	__v2(ϵ) = _v2(ϵ, Δ)
	__uv(ϵ) = _uv(ϵ, Δ)
	return _mult_f(f, __u2), _mult_f(f, __uv), _mult_f(f, __v2)
end
