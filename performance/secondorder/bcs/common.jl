# =============================================================================
# secondorder/bcs/common.jl —— 费米 BCS 浴 Toulouse 模型与严格参照
#
# 模型：H = Σ_σ ϵ_d d_σ†d_σ + Σ_kσ ε_k c_kσ†c_kσ + Σ_kσ V_k (d_σ†c_kσ + h.c.)
#         + Δ Σ_k (c_k↑†c_k↓† + c_k↓c_k↑)          （μ=0，与包内 QuAPI BCS 测试一致）
#
# 杂质两个自旋味通过正常部分与浴耦合，配对只在浴内（IF 指数
# ā↑Δ↑↑a↑ + ā↑Δ↑↓ā↓ + a↓Δ↓↑a↑ + a↓Δ↓↓ā↓ 的 Nambu 重排 ψ = (a↑, ā↓)）。
#
# 严格参照（手写 BdG，二次型）：
#   Nambu 单粒子基 Φ = (d↑, c₁↑..c_L↑ ; d̄↓†, c̄₁↓†..c̄_L↓†)，EOM 给出
#   𝓗 = [[h₊, B],[B', −h₊ᵀ]]（h₊: ϵ_d/ε_k/V_k 正常块；B = Δδ_kl 粒子←空穴；
#   B' = −B*）。传播子 S(t) = e^{−i𝓗t}，热关联 M(0) = V/(1+e^{βΛ})V⁻¹，
#   原始关联 M(t) = S(t)M(0)。G^R = −iθS（与高斯态无关）、Matsubara/平衡
#   G> 全部由谱分解闭式给出；开发期与 ImpurityModelBase 的 cdm 演化对拍。
# =============================================================================
using QuAPI, ImpurityModelBase, LinearAlgebra, Printf

"""
    bcs_toulouse(; β, ϵd, Γ, wlo, whi, Lmodes, Δ, μ) -> NamedTuple

BCS 浴 Toulouse 模型（U=0、自旋二味）：离散 BCS 浴 + 单粒子矩阵与严格参照数据。
"""
function bcs_toulouse(; β=4.0, ϵd=1.0, Γ=0.15, wlo=-3.0, whi=3.0, Lmodes=48,
                      Δ=0.5, μ=0.0)
    J = bounded(w -> Γ, wlo, whi)
    freqs = collect(range(wlo, whi; length=Lmodes + 1))
    bath = discretebcsbath(freqs, J; β=β, μ=μ, Δ=Δ)
    ws, αs = frequencies(bath), spectrumvalues(bath)
    L = length(ws)
    return (; β, ϵd, Γ, μ, Δ, bath, ws, αs, L,
            Vs=sqrt.(αs))
end

"""
    bcs_bdg(mdl) -> 𝓗 (2(L+1)×2(L+1))

Nambu 基 Φ = (d↑, c₁↑..c_L↑ ; d̄↓†, c̄₁↓†..c̄_L↓†) 的 BdG 单粒子矩阵：
𝓗 = [[h₊, B],[B*, −h₊ᵀ]]（Hermitian），B[k+1,k+1] = Δ（浴内配对，杂质无配对）。
EOM：i∂_t c_k↓† = −ε_k c_k↓† − V_k d̄↓† + Δ* c_k↑。
"""
function bcs_bdg(mdl)
    L = mdl.L
    n = L + 1
    Δ = mdl.Δ
    h = zeros(ComplexF64, n, n)
    h[1, 1] = mdl.ϵd
    for k in 1:L
        h[1, 1+k] = mdl.Vs[k]
        h[1+k, 1] = mdl.Vs[k]
        h[1+k, 1+k] = mdl.ws[k]
    end
    B = zeros(ComplexF64, n, n)
    for k in 1:L
        B[1+k, 1+k] = Δ
    end
    H = zeros(ComplexF64, 2n, 2n)
    H[1:n, 1:n] = h
    H[1:n, n+1:2n] = B
    H[n+1:2n, 1:n] = conj(B)
    H[n+1:2n, n+1:2n] = -transpose(h)
    return H
end

"BdG 谱分解：𝓗 = V·Λ·V⁻¹（含左逆），返回 (V, λ, Vinv)。"
function bcs_eig(mdl)
    H = bcs_bdg(mdl)
    F = eigen(H)
    return F.vectors, F.values, inv(F.vectors)
end

"""
    bcs_M0(mdl) -> M₀ = ⟨ΦΦ†⟩（热 Nambu 原始关联矩阵）

M₀ = V·diag(1/(1+e^{−βλ}))·V⁻¹ = I − f(𝓗)（整体式；⟨Φ†Φ⟩ = f(𝓗) 在本 Nambu
约定下逐元成立——单模解析与 4 维 Fock ED 双重验证；反常块符号已随之正确）。
"""
function bcs_M0(mdl)
    V, λ, Vinv = bcs_eig(mdl)
    return V * Diagonal(1 ./ (1 .+ exp.(-mdl.β .* λ))) * Vinv
end

# ------------------------- 严格参照 ------------------------- #

"推迟 Nambu 格林函数：𝒢^R(t) = −iθ(t)·S(t) 的杂质 2×2 块（态无关）。"
function exact_retarded_Nambu(mdl)
    V, λ, Vinv = bcs_eig(mdl)
    n = mdl.L + 1
    idx = [1, n + 1]                       # d↑（粒子）与 d̄↓†（空穴）
    Vimp, Vinvimp = V[idx, :], Vinv[:, idx]
    return t -> -im * (Vimp * Diagonal(exp.(-im .* λ .* t)) * Vinvimp)
end

"平衡态 greater：G₁>(t) = ⟨d↑(t)d↑†(0)⟩ = [S(t)M₀]_{d↑,d↑}（含初关联）。"
function exact_greater_eq(mdl)
    V, λ, Vinv = bcs_eig(mdl)
    M0 = bcs_M0(mdl)
    return t -> (V * (Diagonal(exp.(-im .* λ .* t)) * (Vinv * M0)))[1, 1]
end

"""
    exact_Gτ_ring(mdl, N) -> (𝒢[1..N], h)

Matsubara 环参照：𝒢(τ) = −⟨Φ(τ)Φ†(0)⟩ 的杂质 2×2 Nambu 块
（[1] = −⟨d↑(τ)d↑†⟩ greater 约定，同 secondorder/fermions；
 [2] = −⟨d↑(τ)d̄↓⟩ 反常 Matsubara；𝒢(τ) = −S(τ)M₀ 的 (1,1)/(1,n+1) 元）。
"""
function exact_Gτ_ring(mdl, N)
    V, λ, Vinv = bcs_eig(mdl)
    n = mdl.L + 1
    h = mdl.β / N
    M0 = bcs_M0(mdl)
    out = Vector{Vector{ComplexF64}}(undef, N)
    for m in 1:N
        A = V * (Diagonal(exp.(-λ .* m .* h)) * (Vinv * M0))    # S(τ_m)M₀
        out[m] = [-A[1, 1], -A[1, n + 1]]
    end
    return out, h
end

"误差阶数表。"
function order_table(io, name, Ns, errs)
    println(io, name)
    prev = 0.0
    for (N, e) in zip(Ns, errs)
        if prev == 0.0
            @printf(io, "   N=%4d   err=%.4e\n", N, e)
        else
            @printf(io, "   N=%4d   err=%.4e   ratio=%.2f   order=%.2f\n",
                    N, e, prev / e, log2(prev / e))
        end
        prev = e
    end
end
