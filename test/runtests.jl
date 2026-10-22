using Test, LinearAlgebra

# 加载 GaussianStates 触发 QuAPI 的 extension（fermionic_/bosonic_bdg_matrix
# 的 correlation 方法）；显式导入非导出名 covariance 供测试使用
using GaussianStates
using GaussianStates: FermionicGaussianState, BosonicGaussianState,
                      osee, entropy, covariance, nmodes,
                      fermionic_bdg_matrix, bosonic_bdg_matrix
using QuAPI
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

include("quapi/imagquapi.jl")
include("quapi/bcsimagquapi.jl")

### real-time QuAPI

include("quapi/realquapi.jl")
include("quapi/bcsrealquapi.jl")

### mixed-time (Kadanoff contour) QuAPI

include("quapi/mixedquapi.jl")

include("bcs.jl")
