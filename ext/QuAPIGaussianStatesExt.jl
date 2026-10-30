# ============================================================================
# QuAPIGaussianStatesExt —— QuAPI × GaussianStates 扩展
#
# 由 hybridization function 直接构造影响泛函（IF）系数态的高斯态表示：
#
#   I[ā, a] = exp[Σ_pq ā_p Δ_pq a_q]     （Δ 取 QuAPI 约定：指数负号已吸收）
#
# → 系数态 = doubled 空间的纯配对高斯向量 |Ψ⟩ = exp(Σ_pq Δ_pq ā†_p a†_q)|0⟩，
#   按 GaussianStates v7「态只存实协方差」的约定表示：
#
#   fermionic：FermionicGaussianState(corr)  →  Γ = i⟨γγ⟩   （Majorana，4n×4n）
#   bosonic：  BosonicGaussianState(corr)    →  Σ = ½⟨{R,R}⟩ （(x,p)，4n×4n）
#
#   OSEE = entropy(s, modes)（GaussianStates v7 接口）：doubled 模式序
#   [a₁..aₙ; ā₁..āₙ]（模 1..n = ket contour 节点、n+1..2n = bra），
#   **site p 的割 = 模式对 [p, n+p]**，即 entropy(s, vcat(sites, n .+ sites))。
#
# 配对 → 协方差（一般路径，无 pairing 特殊分支）：系数态是一般配对高斯态
# |Ψ⟩ = exp(Σ C_ab β†_a β†_b)|0⟩（β†_a 为各腿的产生算符；BCS 的 ξξ/ξ̄ξ̄ 反常
# 块天然含内）。构造 = 湮灭条件 ξ = b − W b†（W = C − Cᵀ 一般反对称核）的
# 闭式解：G − conj(W)·G·conj(W) = I、F = −W·Gᵀ（复线性 Kronecker，O(L³)）。
# 对奇异 W（零耦合）良定义；单 flavor 情形与 GaussianStates 的 `purify`
#（Choi 路线）熵逐位互证；BCS 情形以稠密 IF 张量的前缀割 SVD 熵为锚
#（test/osee.jl）。
# 两个重载按统计性质区分（CorrelationFunction 类型本身不分玻色/费米，
# 由构造器名承担统计信息）：
#   fermionic：Δ[ā_p, a_q] 直接配对（复 Δ、任意 contour 均可）。
#   bosonic：只支持**实 Δ**（虚时 Matsubara contour；实时/混合 contour 的
#     复 Δ 对应复协方差，快速算法失效，报错）。直接配对振幅 exp(+Δ) 奇异值
#     > 1 ⟹ I−W² 不定、公式数值失败；取 M = exp(−Δ)（奇异值 = e^{−Re λ(Δ)}
#     < 1）——与 exp(+Δ) 只差每模的粒子-空穴反射，时间割的 Schmidt 谱相同、
#     熵不变（TEMPO/osee 报告 §5.2；几何分布 p_n ∝ t^{2n} 在 t ↔ 1/t 下不变）。
#
# Julia 1.10 下 extension 新名字不并入父包命名空间，故本扩展只给
# GaussianStates 的既有名字（两个 GaussianState 构造器）加方法。
# ============================================================================
module QuAPIGaussianStatesExt

using QuAPI
using QuAPI: AbstractCorrelationFunction, BCSCorrelationFunction
using QuAPI.Dense: GrassmannHybridizationMatrix, BosonicHybridizationMatrix
using GaussianStates
using GaussianStates: FermionicGaussianState, BosonicGaussianState
using LinearAlgebra

# ------------------------------------------------------------------ #
# fermionic：IF 系数态（实 Majorana 协方差；单 flavor 与 BCS 多 band 通用）
# ------------------------------------------------------------------ #

"""
    FermionicGaussianState(x::AbstractCorrelationFunction) -> FermionicGaussianState
    FermionicGaussianState(h::GrassmannHybridizationMatrix) -> FermionicGaussianState

费米 hybridization function `x` 的 IF 系数态
`|Ψ⟩ = exp(Σ_ab C_ab β†_a β†_b)|0⟩`（Grassmann 有限和；任意 contour、
任意 band 数——含 BCS 反常块——恒良定义），以 v7「只存实协方差」的
`FermionicGaussianState` 表示：`Γ = i⟨γγ⟩` 实反对称，`nmodes(s) = L`
（L = 腿数：单 flavor 2n、BCS 4n）。

一般构造 = 湮灭条件闭式（无 pairing 特殊分支）：反对称核 W = C − Cᵀ，
ξ = b − W b† 湮灭态 ⟹ G − conj(W)·G·conj(W) = I、F = −W·Gᵀ（复线性
Kronecker 解，O(L³)）；奇异 W（零耦合）良定义。单 flavor 情形与
`purify`（Choi 路线）熵逐位互证；BCS 情形以稠密 IF 张量的前缀割 SVD
熵为锚（test/osee.jl）。

OSEE 割约定：模式 = 腿重排为「全部 a 腿（contour 序），全部 ā 腿」；
单 flavor n 节点即 site p ↔ [p, n+p]（同旧约定）：

```julia
s = FermionicGaussianState(corr)
n = nmodes(s) ÷ 2                              # a 腿（ket 半区）模式数
entropy(s, vcat(sites, n .+ sites))            # site p ↔ 模式对 [p, n+p]
```
"""
function GaussianStates.FermionicGaussianState(h::GrassmannHybridizationMatrix)
    C = Matrix(h.Δ)
    L = size(C, 1)
    # 一般路径（无 pairing 特殊分支；任意 band——BCS 反常块含内；任意奇异度）：
    # 系数态 |Ψ⟩ = exp(Σ C_ab β†_a β†_b)|0⟩ 的湮灭子 ξ_a = β_a − Σ_b W_ab β†_b
    # （W = C − Cᵀ）恒湮灭该态 ⟹ 湮灭子 Nambu 矩阵 V = [I, −W]（L×2L，恒满行秩，
    # 零耦合 W = 0 亦然）。Nambu 协方差 = V 行空间投影子（purify 内部的同款
    # 湮灭子→投影子数学，F_utilities/FermionicGaussianStates 讲义 §5），再经
    # GaussianStates 的 to_majorana 变实 Majorana 协方差。
    W = C .- transpose(C)
    V = hcat(Matrix{ComplexF64}(I, L, L), -W)
    B = Matrix(qr(adjoint(V)).Q)                   # 2L×L，B†B = I，列空间 = V 行空间
    𝒞 = B * adjoint(B)                             # Nambu 投影子（2L×2L）
    # 模式约定：前 L÷2 = 全部 a 腿（contour 序），后 L÷2 = 全部 ā 腿
    # （单 flavor：site p ↔ [p, n+p]；BCS：L = 4n，2n 模/半区）
    alpleg = findall(l -> l.conj == false, h.basis)
    order = vcat(alpleg, setdiff(1:L, alpleg))     # 腿序 → 模式序
    majorder = vcat([[2p-1, 2p] for p in order]...)# 模式序 → Majorana 序
    Γ = GaussianStates.to_majorana(𝒞)
    Γ = Γ[majorder, majorder]
    return FermionicGaussianState(Γ)
end

GaussianStates.FermionicGaussianState(x::AbstractCorrelationFunction) =
    FermionicGaussianState(GrassmannHybridizationMatrix(x))
GaussianStates.FermionicGaussianState(x::BCSCorrelationFunction) =
    FermionicGaussianState(GrassmannHybridizationMatrix(x))

# 纯 ā–a 配对系数态的二点函数（模式序 [a₁..aₙ; ā₁..āₙ]）。Δ = U·diag(s)·V†
#（SVD）⟹ 配对模式占据 ε_k = s_k²/(1+s_k²)，
#   ⟨a†a⟩ = V ε V†，⟨ā†ā⟩ = Ū ε Uᵀ，⟨āa⟩ = −U s̃ V†（s̃ = s/(1+s²)）
#（GaussianStates v1 `_correlations_pairing` 同式；对稠密系数向量机器精度仲裁）
function _pairing_correlations(Δ::AbstractMatrix{<:Number})
    n = size(Δ, 1)
    Fs = svd(Matrix{ComplexF64}(Δ))
    ε = Fs.S .^ 2 ./ (1 .+ Fs.S .^ 2)                  # 配对模式占据
    s̃ = Fs.S ./ (1 .+ Fs.S .^ 2)
    A = Fs.V * Diagonal(ε) * Fs.V'                     # ⟨a†a⟩
    B = conj(Fs.U) * Diagonal(ε) * transpose(Fs.U)     # ⟨ā†ā⟩
    Cc = -(Fs.U * Diagonal(s̃) * Fs.V')                # ⟨āa⟩
    G = zeros(ComplexF64, 2n, 2n)                      # ⟨b†_i b_j⟩
    F = zeros(ComplexF64, 2n, 2n)                      # ⟨b_i b_j⟩
    G[1:n, 1:n] .= A
    G[n+1:2n, n+1:2n] .= B
    F[1:n, n+1:2n] .= .-transpose(Cc)                  # ⟨a_q ā_p⟩ = −⟨ā_p a_q⟩
    F[n+1:2n, 1:n] .= Cc
    return G, F
end

# (G, F) → Majorana 协方差 Γ_ab = i⟨γ_aγ_b⟩（γ_{2m−1} = b_m + b_m†、
# γ_{2m} = i(b_m† − b_m)），实反对称——b_m b_n† 型反厄米组合的期望纯虚，
# 任意复 (G, F) 下恒实。
function _majorana_covariance(G::AbstractMatrix{<:Number}, F::AbstractMatrix{<:Number})
    n = size(G, 1)
    M = zeros(ComplexF64, 2n, 2n)                      # ⟨γ_aγ_b⟩
    for i in 1:n, j in 1:n
        δ = ifelse(i == j, 1.0, 0.0)
        M[2i-1, 2j-1] = F[i,j] + δ - G[j,i] + G[i,j] + conj(F[j,i])
        M[2i-1, 2j]   = 1im * (δ - G[j,i] - G[i,j] - F[i,j] - conj(F[i,j]))
        M[2i, 2j-1]   = 1im * (G[i,j] + G[j,i] + F[j,i] + conj(F[j,i]) - δ)
        M[2i, 2j]     = conj(F[i,j]) + G[i,j] + δ - G[j,i] - F[i,j]
    end
    Γ = real.(1im .* M)
    return (Γ .- transpose(Γ)) ./ 2                    # 数值反对称化
end

# ------------------------------------------------------------------ #
# bosonic：IF 系数态（实 (x,p) 协方差；d→∞ 高斯极限）
# ------------------------------------------------------------------ #

"""
    BosonicGaussianState(x::AbstractCorrelationFunction) -> BosonicGaussianState
    BosonicGaussianState(h::BosonicHybridizationMatrix) -> BosonicGaussianState

玻色 hybridization function `x` 的 IF 系数态（d→∞ 高斯极限，TEMPO/osee
报告 §5.2），配对振幅 `M = exp(−Δ)`、模式序 `[a₁..aₙ; ā₁..āₙ]` 的纯配对
向量，以 v7「只存实协方差」的 `BosonicGaussianState` 表示：
`Σ = ½⟨{R,R}⟩` 为 `4n×4n` 实对称正定（辛本征值全 ½，纯态）。

**只支持实 Δ**（虚时 Matsubara contour）；实时/混合 contour 的复 Δ 会使
系数态协方差复化、快速算法失效，调用报 `ArgumentError`。

`M = exp(−Δ)`（奇异值 `e^{−λ_k(Δ)} < 1`）：玻色系数向量 `exp(Σ Δ b†b)`
向量化的直接配对振幅 `exp(+Δ)` 奇异值 > 1，`I − W²` 不定、协方差公式数值
失败；`exp(−Δ)` 与之只差每模粒子–空穴反射，时间割 Schmidt 谱相同、熵不变
（几何分布 `p_n ∝ t^{2n}` 在 `t ↔ 1/t` 下不变）。OSEE：

```julia
s = BosonicGaussianState(corr)
n = nmodes(s) ÷ 2                              # contour 节点数
entropy(s, vcat(sites, n .+ sites))            # site p ↔ 模式对 [p, n+p]
```

`BosonicHybridizationMatrix` 重载仅接受 `split=false`（`Δ` 为普通矩阵）。
`split=true` 存的是等时算符次序拆分 `BosonicSplitDelta(M1, M2)`——它不是
配对振幅，不能拆成两个矩阵传入，故直接报错。
"""
function GaussianStates.BosonicGaussianState(hm::BosonicHybridizationMatrix)
    Δ = hm.Δ
    Δ isa AbstractMatrix || throw(ArgumentError(
        "该 BosonicHybridizationMatrix 由 split=true 构造，Δ = BosonicSplitDelta(M1, M2) " *
        "是等时算符次序拆分，不是配对振幅；玻色 OSEE 请用 " *
        "BosonicHybridizationMatrix(x)（split=false）或直接 BosonicGaussianState(x)"))
    if maximum(abs, imag.(Δ)) > 1e-12 * max(1.0, maximum(abs, real.(Δ)))
        throw(ArgumentError(
            "玻色 OSEE 只支持实 Δ（虚时 Matsubara contour）；当前 Δ 含非平凡" *
            "虚部（实时/混合 contour），系数态协方差会复化，快速算法失效"))
    end
    M = exp(-real.(Δ))                                # 配对振幅（矩阵指数）
    n = size(M, 1)
    W = zeros(2n, 2n)                                 # [[0, M], [Mᵀ, 0]]（实对称）
    W[1:n, n+1:2n] .= M
    W[n+1:2n, 1:n] .= transpose(M)
    return BosonicGaussianState(_bosonic_pairing_covariance(W))
end

GaussianStates.BosonicGaussianState(x::AbstractCorrelationFunction) =
    BosonicGaussianState(BosonicHybridizationMatrix(x))

# 纯配对系数向量 |Ψ⟩ = exp(½b†Wb†)|0⟩ 的 (x,p) 协方差（W 实对称、奇异值
# < 1）。湮灭条件 (b_i − Σ_j W_ij b†_j)|Ψ⟩ = 0 给出
#   G = (I − WWᵀ)⁻¹ = ⟨bb†⟩，K = W Gᵀ = ⟨bb⟩，⟨b†b⟩ = Gᵀ − I
#   σ_xx = G + K − ½I，σ_pp = G − K − ½I，σ_xp = 0（实 W）
#（GaussianStates v1 `_bosonic_gamma_pairing` 同式；n=1 与稠密 Fock Choi
# 仲裁至机器精度，辛本征值全 ½）
function _bosonic_pairing_covariance(W::AbstractMatrix{<:Real})
    N = size(W, 1)
    Id = Matrix{Float64}(I, N, N)
    G = (Id - W * transpose(W)) \ Id
    K = W * transpose(G)
    Σ = zeros(Float64, 2N, 2N)                        # (x₁,p₁,…,x_N,p_N)
    for i in 1:N, j in 1:N
        δ = ifelse(i == j, 0.5, 0.0)
        Σ[2i-1, 2j-1] = G[i,j] + K[i,j] - δ
        Σ[2i, 2j]     = G[i,j] - K[i,j] - δ
    end
    return (Σ .+ transpose(Σ)) ./ 2
end

end # module
