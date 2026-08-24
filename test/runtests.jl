using Test, LinearAlgebra

include("../src/includes.jl")


### fermionic Δτ

include("imagtime.jl")
include("realtime.jl")
include("mixedtime.jl")

### imaginary-time QuAPI

include("quapi/imagquapi.jl")
include("quapi/bcsimagquapi.jl")

### real-time QuAPI

include("quapi/realquapi.jl")
include("quapi/bcsrealquapi.jl")

### mixed-time (Kadanoff contour) QuAPI

include("quapi/mixedquapi.jl")

include("bcs.jl")
