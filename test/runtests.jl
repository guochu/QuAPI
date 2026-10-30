using Test, LinearAlgebra

# 加载 GaussianStates 触发 QuAPI 的 extension（Fermionic/BosonicGaussianState
# 的 correlation function 方法）；显式导入测试用到的名字
using GaussianStates
using GaussianStates: FermionicGaussianState, BosonicGaussianState,
                      entropy, covariance, nmodes, ispure
using QuAPI
using QuAPI.Dense
using ImpurityModelBase
# deprecated implementations, kept for cross-validation tests
using QuAPI: fermionic_Δτ_deprecated, bosonic_Δτ_deprecated, bosonic_Δt_deprecated, fermionic_Δm_deprecated, bosonic_Δm_deprecated, isize, rsize


### fermionic Δτ

include("imagtime.jl")
include("realtime.jl")
include("mixedtime.jl")

### GaussianStates extension: IF coefficient-state OSEE

include("osee.jl")

### imaginary-time QuAPI

include("quapi/imag.jl")
include("quapi/bcsimag.jl")

### real-time QuAPI

include("quapi/real.jl")
include("quapi/bcsreal.jl")
include("quapi/bcsmixed.jl")

### exact fermionic IF (normal bath, quadratic chain)

include("quapi/exactfermionic.jl")

### mixed-time (Kadanoff contour) QuAPI

include("quapi/mixed.jl")

include("bcs.jl")
