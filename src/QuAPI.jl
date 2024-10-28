module QuAPI

# Our convention is to absorb the negative sign in the exponent of the IF into the hybridization matrix!!!
export CorrelationMatrix
export AbstractCorrelationFunction, index, branch
# imaginary time
export ImagCorrelationFunction, Cτ, fermionic_Cτ, bosonic_Cτ, Δiw_to_Cτ
# real time
export RealCorrelationFunction, Ct, fermionic_Ct, bosonic_Ct
# mixed time (L-shaped Kadanoff-Byam contour)
export AbstractMixedCorrelationFunction, MixedCorrelationFunction, Cm, fermionic_Cm, bosonic_Cm

using Base: @boundscheck
using ImpurityModelBase



include("util.jl")
include("corrmat.jl")

include("abstractcorr.jl")
include("imagtime/imagtime.jl")
include("realtime/realtime.jl")
include("mixedtime/mixedtime.jl")


end