module QuAPI

# Our convention is to absorb the negative sign in the exponent of the IF into the hybridization matrix!!!
export CorrelationMatrix
export AbstractCorrelationFunction, branch, branches
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

# 费米 normal bath 二次型链的严格 Gaussian IF（谱 → 精确 cell 振幅 + Schur）
export ExactGrassmannIF, exact_fermionic_Δτ, exact_fermionic_Δt

# numerical tolerance
export get_quapi_tol, set_quapi_tol!

using Base: @boundscheck
using ImpurityModelBase
using LinearAlgebra



include("utility.jl")
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

# 腿指标（AbstractContourIndex）与稠密 IF（QuAPI.Dense）
include("contourindex.jl")
include("dense/dense.jl")

# 费米 normal bath 的严格 Gaussian IF（虚时/实时，离散谱）
include("exactfermionic.jl")

end