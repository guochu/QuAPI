module QuAPI

# Our convention is to absorb the negative sign in the exponent of the IF into the hybridization matrix!!!
export CorrelationMatrix
export AbstractCorrelationFunction, index, branch, branches
# imaginary time
export ImagCorrelationFunction, Δτ, fermionic_Δτ, bosonic_Δτ, Δiw_to_Δτ
# real time
export RealCorrelationFunction, Δt, fermionic_Δt, bosonic_Δt
# mixed time (L-shaped Kadanoff-Byam contour)
export AbstractMixedCorrelationFunction, MixedCorrelationFunction, Δm, fermionic_Δm, bosonic_Δm

# infinite time hybridization function
export infinite_Δτ, infinite_Δt

# bcs
export BCSCorrelationFunction, bcs_Δτ, bcs_Δt, bcs_Δm

# contour indices
export AbstractContourIndex, GrassmannContourIndex, BosonicContourIndex

using Base: @boundscheck
using ImpurityModelBase
using LinearAlgebra



include("util.jl")
include("corrmat.jl")

include("abstractcorr.jl")
include("imagtime/imagtime.jl")
include("realtime/realtime.jl")
include("mixedtime/mixedtime.jl")

# infinite time hybridization function
include("infinitetime/imagtime.jl")
include("infinitetime/realtime.jl")

# BCS
include("bcs/bcs.jl")

# contour branches present in a correlation function
"""
    branches(x::AbstractCorrelationFunction) -> Tuple{Vararg{Symbol}}

correlation function `x` 所在的 contour 分支（按 Δ 行序）：
虚时 `(:τ,)`、实时 `(:+, :)`、混合 `(:+, :-, :τ)`；BCS 取正常块的分支。
"""
branches(::ImagCorrelationFunction) = (:τ,)
branches(::RealCorrelationFunction) = (:+, :-)
branches(::AbstractMixedCorrelationFunction) = (:+, :-, :τ)
branches(x::BCSCorrelationFunction) = branches(x.uu)

# 腿指标（AbstractContourIndex）与稠密 IF（QuAPI.Dense）
include("contourindex.jl")
include("models.jl")

end