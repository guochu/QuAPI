# =============================================================================
# secondorder/fermions/common.jl —— 费米 Toulouse（SIAM U=0）模型与严格解
#
# 模型：H = ϵ_d d†d + Σ_k ε_k c_k†c_k + Σ_k V_k (d†c_k + c_k†d)，离散费米浴。
# 严格 A(t) = [e^{-iHt}]₁₁（单粒子对角化）；严格 G(τ) 由热环公式给出；
# 严格平衡 G>(t) = −i[U(t)(I+γ₀ᵀ)]₁₁，γ₀ = 热占据对角（对角浴）。
# 冻结核：fermionic_Δt / fermionic_Δτ / fermionic_Δm（与 bosons 侧同骨架）。
# =============================================================================
using QuAPI, ImpurityModelBase, LinearAlgebra, Printf

"""
    fermi_toulouse(; β, ϵd, Γ, wlo, whi, Lmodes, μ) -> NamedTuple

费米 Toulouse 模型（U=0）：返回离散浴、单粒子矩阵与严格参照所需数据。
"""
function fermi_toulouse(; β=4.0, ϵd=1.0, Γ=0.15, wlo=-3.0, whi=3.0, Lmodes=48, μ=0.0)
    J = bounded(w -> Γ, wlo, whi)
    freqs = collect(range(wlo, whi; length=Lmodes + 1))
    bath = discretebath(Fermion, freqs, J; β=β, μ=μ)
    model = Toulouse(bath; ϵ_d=ϵd)
    ws, αs = frequencies(bath), spectrumvalues(bath)
    L = length(ws)
    hS = zeros(L + 1, L + 1); hS[1, 1] = ϵd
    hb = zeros(L + 1, L + 1)
    hi = zeros(ComplexF64, L + 1, L + 1)
    for k in 1:L
        hb[1+k, 1+k] = ws[k]
        hi[1, 1+k] = sqrt(αs[k]); hi[1+k, 1] = sqrt(αs[k])
    end
    return (; β, ϵd, Γ, μ, bath, model, ws, αs, L, hS, hb, hi, hfull=hS+hb+hi,
             Id=Matrix{Float64}(I, L+1, L+1))
end

"严格 A(t) = [e^{−iHt}]₁₁。"
function exact_return_amplitude(mdl)
    ev = eigen(Hermitian(mdl.hfull))
    return t -> sum(ev.vectors[1, k]^2 * exp(-im * ev.values[k] * t)
                    for k in 1:length(ev.values))
end

"严格 Matsubara G(τ_m)（自由环公式：G_m = [u^m (I−u^N)⁻¹]₁₁ 修正为直接闭式）。"
function exact_Gτ(mdl)
    return m -> begin
        h = mdl.β  # 占位：闭式由 exact_Gτ_ring 给出（见 imagtime.jl 用法）
        error("用 exact_Gτ_ring(mdl, N) 获得节点值序列")
    end
end

"""
    exact_Gτ_ring(mdl, N) -> (G_m for m=0..N−1 的 Vector, h)

严格热环（完整 H，无 Trotter）：G₀ = Tr[e^{−βH}]/… 的单粒子解
G(τ) 满足 G(τ−β) = −⟨d†(0)d(τ)⟩…；数值：构成 (L+1)×(L+1) 环
W = e^{−βh_full} 的块：G_m = [W^m·(I + W)^{-1}·…] —— 直接用
G_m = Σ_p [v_p v_p†]₁₁ e^{−λ_p τ}/(1+e^{−βλ_p})（谱求和，机器精度）。
"""
function exact_Gτ_ring(mdl, N)
    ev = eigen(Hermitian(mdl.hfull))
    h = mdl.β / N
    G = [sum(ev.vectors[1, p]^2 * exp(-ev.values[p] * m * h) /
             (1 + exp(-mdl.β * ev.values[p])) for p in 1:length(ev.values))
         for m in 0:N-1]
    return G, h
end

"热占据对角 γ₀（对角浴下杂质-浴单粒子热矩阵为对角分块）。"
function exact_gamma_diag(mdl)
    ev = eigen(Hermitian(mdl.hfull))
    nk = [1 / (exp(mdl.β * λ) + 1) for λ in ev.values]
    return ev.vectors * Diagonal(nk) * ev.vectors'
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
