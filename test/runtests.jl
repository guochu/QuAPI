using Test, LinearAlgebra

using QuAPI
using ImpurityModelBase
# deprecated implementations, kept for cross-validation tests
using QuAPI: fermionic_Δτ_deprecated, bosonic_Δτ_deprecated, bosonic_Δt_deprecated, fermionic_Δm_deprecated, bosonic_Δm_deprecated, isize, rsize


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
