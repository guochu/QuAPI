# ============================================================================
# dense.jl —— 影响泛函的稠密表示与 contour 组装工具（QuAPI.Dense）
# ============================================================================

module Dense

using ..QuAPI: AbstractCorrelationFunction, AbstractRealCorrelationFunction,
               BCSCorrelationFunction,
               ImagCorrelationFunction, RealCorrelationFunction,
               AbstractMixedCorrelationFunction, CorrelationMatrix,
               AbstractContourIndex, GrassmannContourIndex, BosonicContourIndex,
               branch, branches, isize, rsize
using LinearAlgebra: I, kron, triu, tril

export GrassmannADT, BosonicADT, BosonicPT,
       GrassmannHybridizationMatrix, BosonicHybridizationMatrix,
       GrassmannIF, BosonicIF, Zvalue, margin, basis, complementindex,
       contour_ordering

"稠密 IF 张量的最大元素数（2³⁰）。"
const MAXTENSORDIM = 2^30

include("common.jl")
include("fermion.jl")
include("boson.jl")

end # module Dense
