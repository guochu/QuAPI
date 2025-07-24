
_u(ϵ, Δ) = (1 + ϵ/sqrt(ϵ^2 + Δ^2))/2
_v(ϵ, Δ) = (1 - ϵ/sqrt(ϵ^2 + Δ^2))/2

struct BCSCorrelationFunction{C <: AbstractCorrelationFunction}
	c11::C
	c12::C
	c21::C
	c22::C
end

function Base.getindex(x::BCSCorrelationFunction, i::Int, j::Int)
	(1 <= i <= 2) || throw(BoundsError(1:2, i))
	(1 <= j <= 2) || throw(BoundsError(1:2, j))
	return ifelse(i==1, ifelse(j==1, x.c11, x.c12), ifelse(j==1, x.c21, x.c22))
end
Base.getindex(x::BCSCorrelationFunction, c1::Bool, c2::Bool) = ifelse(c1, ifelse(c2, x.c11, x.c12), ifelse(c2, x.c21, x.c22))



_mult_f(x::AbstractBoundedFunction, f) = similar(x, w -> f(w) * x.f(w)) 
_mult_f(x::DiracDelta, f) = x * f
