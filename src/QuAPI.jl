module QuAPI

# Our convention is to absorb the negative sign in the exponent of the IF into the hybridization matrix!!!
export CorrelationMatrix
export AbstractCorrelationFunction, index, branch
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

using Base: @boundscheck
using ImpurityModelBase



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

end