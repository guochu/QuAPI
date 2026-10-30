# =============================================================================
# exactif/test/toulouse_greens.jl —— Toulouse（U=0 SIAM）严格核的
# 1/2/4 阶 Trotter Green 函数计算与误差标度（机器精度核数据版）
#
# 模型：H = Σε_k c†c + ΣV_k(d†c_k+h.c.) + ε_d d†d；ε_d ∈ kick（Ng–Reichman
#       分割），核 W(θ) 只含浴+耦合。浴：平台杂化 Γ₀（Toulouse），GL 节点。
# 核基计算（多项式，无指数）：IF 为高斯 ⟹ 双费米插入收缩 = 单粒子列递推
#   G>(t_N) = −i·r₁，r₀ = [1;n̄] = (I+γ₀ᵀ)e₀，r ← W(θ)·r 逐步作用。
# 核数据 = h₀ 本征分解闭式（机器精度，任意 M 稳定；极点/伴随矩阵法 M>8
#   病态已弃；Volterra/梯形法 8.5e-5 地板已弃）：
#   F(t) = Σⱼ v₁ⱼ² e^{−iλⱼt}
#   φ_t(E) = Σⱼ v₁ⱼ² (e^{−iλⱼt} − e^{−iEt}) / (i(E−λⱼ))
#   ψ_t(E_k,E_l) = Σⱼ v₁ⱼ² Ψ(λⱼ;E_k,E_l,t)   （闭式，对角 −iτ 结构）
# 阶：O1=(Φ(Δt)W(Δt))^N；O2=每胞 Φ(½)WΦ(½)；
#     O4=Suzuki-5 段 [p,p,θ−4p,p,p]，p=θ/(4−4^{1/3})（3 段只有 2 阶）；
#     O6=递归分形 S₆(θ)=S₄(p)S₄(p)S₄(θ−4p)S₄(p)S₄(p)，p=θ/(4−4^{1/5})
#     （25 个 S₂ 段/步；标量检验锁定 6 阶）。负系数段用 W(θ)=W(|θ|)†。
# 参照：节点浴含 ε_d 的本征分解（机器精度）。
# 输出：greens_error.dat；绘图 plot_greens.py → greens_error.png。
# =============================================================================
include(joinpath(@__DIR__, "..", "brute.jl"))
include(joinpath(@__DIR__, "..", "volterra.jl"))
using LinearAlgebra

function main()
    Γ0, Wb = 0.15, 3.0
    β, εd, T = 2.0, 1.0, 2.0
    M = 64
    xs, ws = gauss_legendre(M)
    Es = [Wb * x for x in xs]
    wΓ = [Wb * w * Γ0 for w in ws]
    Vsq = sqrt.(wΓ)
    n̄ = [1 / (1 + exp(β * e)) for e in Es]
    # ---- 核数据：h₀（ε=0）本征分解闭式 ----
    h0 = zeros(ComplexF64, M + 1, M + 1)
    for k in 1:M
        h0[1, 1+k] = Vsq[k]; h0[1+k, 1] = Vsq[k]; h0[1+k, 1+k] = Es[k]
    end
    ev = eigen(Hermitian(h0))
    λ = ev.values
    a2 = ev.vectors[1, :] .^ 2
    Fk(t) = sum(a2[j] * exp(-im * λ[j] * t) for j in 1:M+1)
    # I(E,λ,t) = ∫₀ᵗ e^{−iEτ}e^{−iλ(t−τ)}dτ
    Iint(E, λj, t) = (exp(-im * λj * t) - exp(-im * E * t)) / (im * (E - λj))
    φk(t, E) = sum(a2[j] * Iint(E, λ[j], t) for j in 1:M+1)
    # ψ 闭式：k≠l = [I(E_k)−I(E_l)]/(E_k−E_l)；k=l = −i e^{−iλt}∫₀ᵗτe^{−iaτ}dτ,
    #   ∫₀ᵗτe^{−iaτ}dτ = [(1+iat)e^{−iat}−1]/a²·(−1)… 以数值对矩阵校准（见自检）
    intτ(a, t) = (im == 0) ? t^2 / 2 : begin
        # ∫₀ᵗ τ e^{−iaτ} dτ = (1 − e^{−iat}(1 + i a t)) / a^2 · (−1/i²·…) → 用等价稳定式
        # 推导：∫τe^{cτ} = e^{cτ}(cτ−1)/c²，c=−ia ⟹ [(−iat−1)e^{−iat} + 1] / (−ia)²
        ((1 + im * a * t) * exp(-im * a * t) - 1) / (a^2 * im^2 * (-1))
    end
    function ψk(t, Eq, El)
        if Eq != El
            sum(a2[j] * (Iint(Eq, λ[j], t) - Iint(El, λ[j], t)) / (Eq - El)
                for j in 1:M+1)
        else
            sum(a2[j] * (-im) * exp(-im * λ[j] * t) *
                intτ(Eq - λ[j], t) for j in 1:M+1)
        end
    end
    # ---- 自检：闭式 vs 矩阵元（机器精度）----
    θc = 0.137
    Wc = exp(-im * θc * h0)
    φc = [φk(θc, e) for e in Es]
    ψc = [ψk(θc, Es[k], Es[l]) for k in 1:M, l in 1:M]
    chk1 = abs(Fk(θc) - Wc[1, 1])
    chk2 = maximum(abs.(-im * Vsq .* φc .- Wc[1, 2:end]))
    blk = [exp(-im*Es[k]*θc)*(k == l ? 1.0 : 0.0) - im*Vsq[k]*Vsq[l]*ψc[k,l]
           for k in 1:M, l in 1:M]
    chk3 = maximum(abs.(blk .- Wc[2:end, 2:end]))
    println("核数据自检（vs 矩阵）：F = $chk1  φ = $chk2  ψ = $chk3")
    # ---- W(θ) 列作用（正 θ 常规；负 θ 用 W(|θ|)†）----
    cache = Dict{Float64,Tuple{ComplexF64,Vector{ComplexF64},Matrix{ComplexF64}}}()
    Wdata(θ) = get!(cache, θ) do
        (Fk(θ), [φk(θ, e) for e in Es],
         [ψk(θ, Es[k], Es[l]) for k in 1:M, l in 1:M])
    end
    function apply_W!(r, θ)
        if θ < 0
            a = -θ
            Fθ, φv, ψm = Wdata(a)
            src = r[2:end] .* Vsq
            r0 = r[1]
            rn = zeros(ComplexF64, M + 1)
            rn[1] = conj(Fθ) * r0 + im * sum(src[k] * conj(φv[k]) for k in 1:M)
            for l in 1:M
                rn[l+1] = im * r0 * Vsq[l] * conj(φv[l]) +
                          r[l+1] * exp(im * Es[l] * a) +
                          im * Vsq[l] * sum(src[k] * conj(ψm[k, l]) for k in 1:M)
            end
            r .= rn
            return r
        end
        Fθ, φv, ψm = Wdata(θ)
        src = r[2:end] .* Vsq
        r0 = r[1]
        rn = zeros(ComplexF64, M + 1)
        rn[1] = r0 * Fθ - im * sum(src[k] * φv[k] for k in 1:M)
        for l in 1:M
            rn[l+1] = -im * r0 * Vsq[l] * φv[l] + r[l+1] * exp(-im * Es[l] * θ) -
                      im * Vsq[l] * sum(src[k] * ψm[l, k] for k in 1:M)
        end
        r .= rn
        return r
    end
    kick!(r, α) = (r[1] *= exp(-im * α * εd); r)
    # 递归展平：S₆ 的内层 S₄ 与中段都是完整阶梯（含自身 kick 结构），
    # kick 沿阶梯边界自然相消 ⟹ 最干净的组装 = 对每阶显式写递归：
    function apply_order!(r, order::Int, θ::Float64)
        if order == 1
            apply_W!(r, θ); kick!(r, θ)
        elseif order == 2
            kick!(r, θ / 2); apply_W!(r, θ); kick!(r, θ / 2)
        elseif order == 4
            p = θ / (4 - 4^(1 / 3))
            apply_order!(r, 2, p); apply_order!(r, 2, p)
            apply_order!(r, 2, θ - 4p)
            apply_order!(r, 2, p); apply_order!(r, 2, p)
        else
            p = θ / (4 - 4^(1 / 5))
            apply_order!(r, 4, p); apply_order!(r, 4, p)
            apply_order!(r, 4, θ - 4p)
            apply_order!(r, 4, p); apply_order!(r, 4, p)
        end
    end
    function G_trot(order::Int, Δt)
        N = round(Int, T / Δt)
        r = ComplexF64[1; n̄]
        for _ in 1:N
            apply_order!(r, order, Δt)
        end
        return -im * r[1]
    end
    # ---- 参照：节点浴含 ε_d 本征分解 ----
    hfull = copy(h0); hfull[1, 1] += εd
    evf = eigen(Hermitian(hfull))
    w0 = ComplexF64[1; n̄]
    c0 = evf.vectors' * w0
    gref = -im * sum(evf.vectors[1, j] * c0[j] * exp(-im * evf.values[j] * T)
                     for j in 1:M+1)
    # ---- 误差标度 ----
    open(joinpath(@__DIR__, "greens_error.dat"), "w") do io
        println(io, "# order  N  dt  abs_err")
        for order in [1, 2, 4, 6]
            prev = nothing
            for N in [5, 10, 20, 40, 80]
                Δt = T / N
                err = abs(G_trot(order, Δt) - gref)
                slope = prev === nothing ? "" : "  斜率 = $(round(log(prev/err)/log(2), digits=2))"
                println("order=$order N=$N Δt=$(round(Δt,digits=4)): |err| = ",
                        round(err, sigdigits=3), slope)
                println(io, "$order $N $Δt $err")
                prev = err
            end
        end
    end
    println("\n参照 G>_exact(T=$T) = ", round(gref, digits=6))
end
main()
