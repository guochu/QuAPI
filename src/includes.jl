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
