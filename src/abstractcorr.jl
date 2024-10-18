abstract type AbstractCorrelationFunction end


Base.size(x::AbstractCorrelationFunction, i::Int...) = error("size not defined for correlation type $(typeof(x))")
index(x::AbstractCorrelationFunction, i::Int, j::Int) = error("index not defined for correlation type $(typeof(x))")
branch(x::AbstractCorrelationFunction; b1::Symbol=τ, b2::Symbol=:τ) = error("branch not defined for correlation type $(typeof(x))")