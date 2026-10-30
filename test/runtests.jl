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

### correlation-function API tests (Δτ/Δt/Δm, branch conventions, alignment)

include("api/runtests.jl")

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

### BCS block alignment tests (Δm ↔ Δτ/Δt)

include("api/bcs.jl")

### dense (QuAPI.Dense) vs QuAPI path-sum consistency

include("dense/imagtime_dense_vs_quapi.jl")
include("dense/realtime_dense_vs_quapi.jl")

### archive regression (performance benchmarks vs stored CSVs)

include("quapi/archive_regression.jl")


