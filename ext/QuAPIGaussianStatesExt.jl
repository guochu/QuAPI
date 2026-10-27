# ============================================================================
# QuAPIGaussianStatesExt —— QuAPI × GaussianStates 扩展
#
# 由 hybridization function 直接构造影响泛函（IF）系数态的高斯态表示：
#
#   I[ā, a] = exp[Σ_pq ā_p Δ_pq a_q]     （Δ 取 QuAPI 约定：指数负号已吸收）
#
# → 系数态 = doubled 空间的纯配对高斯向量，经 `fermionic_bdg_matrix` /
#   `bosonic_bdg_matrix` 重载给出 GaussianStates 的 bdg 输入；OSEE 一律在
#   构造出的 GaussianState 上调用 `osee` 计算（协方差快速算法，与
#   GrassmannMPS bond 熵同一形式体系，见 GTEMPO/osee/osee_report.md 与
#   TEMPO/osee/osee_report.md）。
#
# 两个重载按统计性质区分（CorrelationFunction 类型本身不分玻色/费米，
# 由调用的函数名承担统计信息）：
#   fermionic_bdg_matrix(x)：W = [[0, −Δᵀ], [Δ, 0]]（复反对称，任意 contour）。
#     系数态 exp(ΣΔ c̄†c†)|0⟩ 是 Grassmann 有限和，恒良定义（含复 Δ）。
#   bosonic_bdg_matrix(x)：只支持**实 Δ**（虚时 Matsubara contour；
#     实时/混合 contour 的复 Δ 对应复协方差，快速算法失效，报错）。
#     直接配对振幅 exp(+Δ) 奇异值 >1 ⟹ I−W² 不定、公式数值失败；取
#     M = exp(−Δ)（奇异值 = e^{−Re λ(Δ)} < 1），W = [[0, M],[Mᵀ, 0]]——
#     与 exp(+Δ) 只差每模的粒子-空穴反射，时间割的 Schmidt 谱相同、熵不变
#     （TEMPO/osee 报告 §5.2；几何分布 p_n ∝ t^{2n} 在 t ↔ 1/t 下不变）。
#
# Julia 1.10 下 extension 新名字不并入父包命名空间，故本扩展只给
# GaussianStates 的既有函数加方法，不引入新名字。
# ============================================================================
module QuAPIGaussianStatesExt

using QuAPI
using QuAPI: AbstractCorrelationFunction
using QuAPI.Dense: _contour_delta, BosonicHybridizationMatrix, GrassmannHybridizationMatrix
using GaussianStates
using GaussianStates: FermionicGaussianState
using LinearAlgebra

# ------------------------------------------------------------------ #
# fermionic_bdg_matrix 重载：IF 系数态的反对称配对矩阵
#（contour Δ 组装在 src/models.jl（QuAPI.Dense）的 _contour_delta，两种实现共用）
# ------------------------------------------------------------------ #

"""
    fermionic_bdg_matrix(x::AbstractCorrelationFunction) -> Matrix{ComplexF64}

费米 hybridization function `x` 的 IF 系数态配对矩阵（GaussianStates 模式序
`[a₁..aₙ; ā₁..āₙ]`，`W[ā_p, a_q] = Δ[p,q]` 再反称化）：

```text
W = [[0, −Δᵀ], [Δ, 0]]     （左下块 = Δ）
```

系数态 `|Ψ⟩ = exp(Σ_pq Δ_pq c̄†_p c†_q)|0⟩` 是 Grassmann 有限和，对任意
contour（虚时/实时/混合，含复 Δ）恒良定义。OSEE：

```julia
s = FermionicGaussianState(fermionic_bdg_matrix(x))
n = nmodes(s) ÷ 2                              # contour 节点数
osee(s, vcat(sites, n .+ sites))               # site p ↔ 模式对 [p, n+p]
```
"""
function GaussianStates.fermionic_bdg_matrix(x::AbstractCorrelationFunction)
    # 单 flavor：与稠密 IF（GrassmannHybridizationMatrix，顶点移位后的约定）
    # 保持同一 Δ —— 从腿矩阵按节点提取 [C[āᵢ, aⱼ]]，确保与 GrassmannIF 张量
    # 的 OSEE 逐位一致；多 band（BCS）/混合 contour 沿用 _contour_delta 原序。
    C = Matrix(GrassmannHybridizationMatrix(x).Δ)
    n2 = size(C, 1)
    if n2 % 2 == 0 && all(b -> b.band == 1, GrassmannHybridizationMatrix(x).basis)
        n = n2 ÷ 2
        Δ = [C[2i, 2j-1] for i in 1:n, j in 1:n]
    else
        Δ = _contour_delta(x)
        n = size(Δ, 1)
    end
    W = zeros(ComplexF64, 2n, 2n)
    W[n+1:2n, 1:n] .= Δ
    return W - transpose(W)
end

# ------------------------------------------------------------------ #
# bosonic_bdg_matrix 重载：IF 系数态的（非正定）对称配对矩阵
# ------------------------------------------------------------------ #

"""
    bosonic_bdg_matrix(x::AbstractCorrelationFunction) -> Matrix{Float64}
    bosonic_bdg_matrix(h::BosonicHybridizationMatrix) -> Matrix{Float64}

玻色 hybridization function `x` 的 IF 系数态配对矩阵（d→∞ 高斯极限，
TEMPO/osee 报告 §5.2）。**只支持实 Δ**（虚时 Matsubara contour）；实时/
混合 contour 的复 Δ 会使系数态协方差复化、快速算法失效，调用报错。

```text
M = exp(−Δ) ,      W = [[0, M], [Mᵀ, 0]]    （实对称，谱 ⊂ (−1,1)）
```

玻色系数向量 `exp(Σ Δ b†b)` 向量化后为 `exp(Σ e^{+Δ} c†d†)|0⟩`——直接配对
振幅 `exp(+Δ)` 奇异值 > 1，`I − W²` 不定，协方差公式数值失败；取
`M = exp(−Δ)`（与 `exp(+Δ)` 只差每模粒子–空穴反射，时间割 Schmidt 谱相同，
熵不变：几何分布 `p_n ∝ t^{2n}` 在 `t ↔ 1/t` 下不变）。`M` 的奇异值
`e^{−Re λ_k(Δ)} < 1`（热 bath 谱实部为正）。OSEE：

```julia
s = BosonicGaussianState(bosonic_bdg_matrix(x))     # isbdg = false
osee(s, vcat(sites, n .+ sites))                    # site p ↔ 模式对 [p, n+p]
```

`BosonicHybridizationMatrix` 重载仅接受 `split=false`（`Δ` 为普通矩阵）。
`split=true` 存的是等时算符次序拆分 `BosonicSplitDelta(M1, M2)`——它与 BdG 的
`(h, g)` 块含义无关，不能拆成两个矩阵传给 `bosonic_bdg_matrix(h, g)`，故直接报错。
"""
function GaussianStates.bosonic_bdg_matrix(h::BosonicHybridizationMatrix)
    Δ = h.Δ
    Δ isa AbstractMatrix || throw(ArgumentError(
        "该 BosonicHybridizationMatrix 由 split=true 构造，Δ = BosonicSplitDelta(M1, M2) " *
        "是等时算符次序拆分，不是 BdG 的 (h, g) 块；玻色 OSEE 请用 " *
        "BosonicHybridizationMatrix(x)（split=false）或直接 bosonic_bdg_matrix(x)"))
    if maximum(abs, imag.(Δ)) > 1e-12 * max(1.0, maximum(abs, real.(Δ)))
        throw(ArgumentError(
            "玻色 OSEE 只支持实 Δ（虚时 Matsubara contour）；当前 Δ 含非平凡" *
            "虚部（实时/混合 contour），系数态协方差会复化，快速算法失效"))
    end
    M = exp(-real.(Δ))                                # 矩阵指数
    n = size(M, 1)
    W = zeros(Float64, 2n, 2n)
    W[1:n, n+1:2n] .= M
    W[n+1:2n, 1:n] .= transpose(M)
    return W
end

GaussianStates.bosonic_bdg_matrix(x::AbstractCorrelationFunction) =
    GaussianStates.bosonic_bdg_matrix(BosonicHybridizationMatrix(x))

end # module
