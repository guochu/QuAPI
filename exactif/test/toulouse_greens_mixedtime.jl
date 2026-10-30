# =============================================================================
# exactif/test/toulouse_greens_mixedtime.jl —— Toulouse（U=0）混合时间
# Green 函数：1/2/4/6 阶 Trotter + 误差标度
#
# L 轮廓两点的单粒子装配（本轮严格推导，多体逐位锚定）：
#   G(τᵢ, t_a) = e₁ᵀ · S(δτ)^m · F_P · K(δt)^a · e₁
# 其中 S = 虚时环单胞（Ng–Reichman 分割：Φ_τ 与 e^{−δτh₀} 按 Trotter 阶
# 组合），F_P = 环统计矩阵（S^{Nτ} 谱的等效频率 f=1/(1+e^{βω})），
# K = Keldysh 单胞（kick 半步夹 e^{−δt h₀}），m/a = τ/t 读数节点。
# 推导：Tr[d̂ X d̂† Y]/Z 中 ηₚX ηₚ† 的移位结构 ⟹ Tₚ = yₚ∏_{r≠p}(1+xᵣyᵣ)
# ⟹ mb = Σₚa2ₚe^{−λₚτᵢ}e^{−iλₚt_a}fₚ（ε_d=0 逐位 = 谱闭式）。
# 参照：完整 h₀+ε_d 谱的未 Trotter 化闭式。
# 输出：greens_mixedtime_error.dat；绘图 plot_greens_mixedtime.py。
# =============================================================================
include(joinpath(@__DIR__, "exactif_api.jl"))
using LinearAlgebra

function main_disabled()
    Γ0, Wb = 0.15, 3.0
    β, εd, t = 2.0, 0.3, 0.8
    M = 64
    Es, Vsq, _ = _quad_nodes(bounded(ω -> Γ0, -Wb, Wb), 0.0, M)
    h0 = zeros(ComplexF64, M + 1, M + 1)
    for k in 1:M
        h0[1, 1+k] = Vsq[k]; h0[1+k, 1] = conj(Vsq[k]); h0[1+k, 1+k] = Es[k]
    end
    hfull = copy(h0); hfull[1, 1] += εd
    evf = eigen(Hermitian(hfull))
    a2f = abs2.(evf.vectors[1, :])
    fF = [1 / (1 + exp(-β * evf.values[p])) for p in 1:M+1]
    # 参照 = 装配的连续极限（同一 L 轮廓两点，h_full 无 Trotter）：
    # G(τᵢ,t_a) = e₁ᵀe^{−τᵢh_f}n̄(h_f)e^{−t_a h_f}e₁ = Σₚa2fₚe^{−λₚτᵢ}e^{−iλₚt_a}fFₚ
    Gref(τ, tt) = sum(a2f[p] * exp(-evf.values[p] * τ) *
                      exp(-im * evf.values[p] * tt) * fF[p] for p in 1:M+1)
    Φm(α) = begin
        D = Matrix{ComplexF64}(I, M+1, M+1); D[1, 1] = exp(-α * εd); D
    end
    function build(order, θ)
        order == 1 && return Φm(θ) * exp(-θ * h0)
        order == 2 && return Φm(θ/2) * exp(-θ * h0) * Φm(θ/2)
        if order == 4
            p = θ / (4 - 4^(1/3))
            return build(2, p) * build(2, p) * build(2, θ - 4p) * build(2, p) * build(2, p)
        else
            p = θ / (4 - 4^(1/5))
            return build(4, p) * build(4, p) * build(4, θ - 4p) * build(4, p) * build(4, p)
        end
    end
    e1 = vcat([1.0], zeros(M))
    function G_trot(order::Int, Nτ::Int, Nt::Int, τ_frac::Float64, t_frac::Float64)
        δτ, δt = β / Nτ, t / Nt
        S = build(order, δτ)
        K = order == 1 ? Φm(δt) * exp(-δt * h0) : build(order, δt)
        m = round(Int, τ_frac * Nτ); a = round(Int, t_frac * Nt)
        PN = S^Nτ
        F = eigen(PN)
        ω = -real.(log.(F.values)) ./ β
        F_P = F.vectors * Diagonal([1 / (1 + exp(β * ω[j])) for j in 1:M+1]) / F.vectors
        return real((e1' * S^m * F_P * K^a * e1))
    end
    open(joinpath(@__DIR__, "greens_mixedtime_error.dat"), "w") do io
        println(io, "# order  Ntau  Nt  abs_err")
        for order in [1, 2, 4, 6]
            prev = nothing
            for N in [2, 4, 8, 16]
                Nτ = N * 5; Nt = N
                err = abs(G_trot(order, Nτ, Nt, 0.5, 0.5) - real(Gref(β/2, t/2)))
                slope = prev === nothing ? "" : "  斜率 = $(round(log(prev/err)/log(2), digits=2))"
                println("order=$order Nτ=$Nτ Nt=$Nt: |err| = ", round(err, sigdigits=3), slope)
                println(io, "$order $Nτ $Nt $err")
                prev = err
            end
        end
    end
    println("\n参照 G(β/2, t/2) = ", round(real(Gref(β/2, t/2)), digits=6))
end
# main()  # 自检模式禁用主程序

# ---- ε_d=0 自检（追加）----
function selfcheck()
    Γ0, Wb = 0.15, 3.0
    β, εd, t = 2.0, 0.0, 0.8
    M = 16
    Es, Vsq, _ = _quad_nodes(bounded(ω -> Γ0, -Wb, Wb), 0.0, M)
    h0 = zeros(ComplexF64, M+1, M+1)
    for k in 1:M; h0[1,1+k] = Vsq[k]; h0[1+k,1] = conj(Vsq[k]); h0[1+k,1+k] = Es[k]; end
    hfull = copy(h0)
    ev0 = eigen(Hermitian(h0)); evf = eigen(Hermitian(hfull))
    a2f = abs2.(evf.vectors[1, :])
    fF = [1/(1+exp(-β*evf.values[p])) for p in 1:M+1]
    Gref(τ, tt) = sum(a2f[p]*exp(-evf.values[p]*τ)*exp(-im*evf.values[p]*tt)*fF[p] for p in 1:M+1)
    Φm(α) = Matrix{ComplexF64}(I, M+1, M+1)   # ε_d=0 ⟹ 恒等
    build(order, θ) = order == 1 ? exp(-θ*h0) :
        order == 2 ? exp(-θ*h0) :
        (p = θ/(4 - 4^(1/3)); S2(θ) = exp(-θ*h0);
         order == 4 ? S2(p)*S2(p)*S2(θ-4p)*S2(p)*S2(p) :
         (p4 = θ/(4 - 4^(1/5)); v = [p, p, θ-4p, p, p];
          S4(q) = (pq = q/(4-4^(1/3)); S2(pq)*S2(pq)*S2(q-4pq)*S2(pq)*S2(pq));
          S4(p4)*S4(p4)*S4(θ-4p4)*S4(p4)*S4(p4)))
    e1 = vcat([1.0], zeros(M))
    function G_trot(order, Nτ, Nt, m, a)
        δτ, δt = β/Nτ, t/Nt
        S = build(order, δτ); K = build(order, δt)
        PN = S^Nτ
        F = eigen(PN)
        ω = -real.(log.(F.values)) ./ β
        F_P = F.vectors * Diagonal([1/(1+exp(β*ω[j])) for j in 1:M+1]) / F.vectors
        return real((e1' * S^m * F_P * K^a * e1))
    end
    println("ε_d=0 自检（各阶 vs 参照）：")
    for order in [1, 2, 4, 6]
        for N in [2, 4]
            Nτ = N*5; Nt = N
            gv = G_trot(order, Nτ, Nt, Nτ÷2, Nt÷2)
            gr = real(Gref(Nτ÷2 * β/Nτ, Nt÷2 * t/Nt))
            println("  order=$order Nτ=$Nτ Nt=$Nt: trot=$(round(gv,digits=6))  ref=$(round(gr,digits=6))  diff=$(round(abs(gv-gr),sigdigits=3))")
        end
    end
end
selfcheck()
