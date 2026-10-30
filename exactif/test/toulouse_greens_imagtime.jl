# =============================================================================
# exactif/test/toulouse_greens_imagtime.jl —— Toulouse（U=0）严格核的
# 虚时 Green 函数：1/2/4/6 阶 Trotter + 误差标度（实时版镜像）
#
# 虚时严格核（Ng–Reichman 分割的虚时对应）：
#   单胞 P(δτ) 由 kick Φ(α)=e^{−αε_d n̂} 与核 W(θ)=e^{−θh₀} 组装（h₀ = 浴+耦合，
#   ε_d 全在 kick）；Trotter 环 P^N ≈ e^{−βh_full}。
# Green（单粒子严格公式，dbg17/dbg18 双锚定）：
#   G(τ_m) = [P^m·(I − n̄_P)]₁₁，n̄_P = V·diag(f(ω_j))·V⁻¹，
#   ω_j = −log|Λ_j(P^N)|/β 的等效频率（f 由环谱给出——虚时统计进入 n̄，
#   与实时版的 (I+γ₀ᵀ) 初值地位对应）。
# 阶：O1/O2 的环读数误差逐点相同（多体直算证实：环周期闭合把 Lie 方案的
#   O(δτ) 边界 kick 误差均匀化——虚时环基线即 2 阶）；O4/O6 = 递归 Suzuki
#   （apply_order! 风格；负系数段用 e^{−θh₀} 直接支持 θ<0）。
# 参照：热环谱求和 G(τ)=Σ_p v₁ₚ²e^{−λₚτ}/(1+e^{−βλₚ})（机器精度，
#   与 performance/secondorder 同式）。
# 输出：greens_imagtime_error.dat；绘图 plot_greens_imagtime.py。
# =============================================================================
include(joinpath(@__DIR__, "..", "brute.jl"))
include(joinpath(@__DIR__, "..", "volterra.jl"))
using LinearAlgebra

function main()
    Γ0, Wb = 0.15, 3.0
    β, εd = 2.0, 1.0
    M = 64
    xs, ws = gauss_legendre(M)
    Es = [Wb * x for x in xs]
    wΓ = [Wb * w * Γ0 for w in ws]
    Vsq = sqrt.(wΓ)
    h0 = zeros(ComplexF64, M + 1, M + 1)
    for k in 1:M
        h0[1, 1+k] = Vsq[k]; h0[1+k, 1] = Vsq[k]; h0[1+k, 1+k] = Es[k]
    end
    hfull = copy(h0); hfull[1, 1] += εd
    evf = eigen(Hermitian(hfull))
    Gref(τ) = sum(evf.vectors[1, p]^2 * exp(-evf.values[p] * τ) /
                  (1 + exp(-β * evf.values[p])) for p in 1:M+1)
    Φm(α) = begin
        D = Matrix{ComplexF64}(I, M+1, M+1); D[1, 1] = exp(-α * εd); D
    end
    # ---- 环统计：单胞幂的等效频率 → n̄_P ----
    function ring_occup(order, N)
        δτ = β / N
        S = ring_cell(order, δτ)
        PN = S^N
        F = eigen(PN)
        ω = -real.(log.(F.values)) ./ β   # 保留复对数实部（保 kick 相位差，阶不虚高）
        ft = [1 / (1 + exp(β * ω[j])) for j in 1:M+1]
        return F.vectors * Diagonal(ft) / F.vectors
    end
    # ---- 单胞组装（递归 Suzuki；kick 与 W(θ) 交替，θ<0 直接支持）----
    function ring_cell(order, θ)
        if order == 1
            return Φm(θ) * exp(-θ * h0)
        elseif order == 2
            return Φm(θ / 2) * exp(-θ * h0) * Φm(θ / 2)
        elseif order == 4
            p = θ / (4 - 4^(1/3))
            return ring_cell(2, p) * ring_cell(2, p) * ring_cell(2, θ - 4p) *
                   ring_cell(2, p) * ring_cell(2, p)
        else
            p = θ / (4 - 4^(1/5))
            return ring_cell(4, p) * ring_cell(4, p) * ring_cell(4, θ - 4p) *
                   ring_cell(4, p) * ring_cell(4, p)
        end
    end
    function G_trot(order::Int, N::Int, m::Int)
        δτ = β / N
        P = ring_cell(order, δτ)
        nb = ring_occup(order, N)
        return real((P^m * (I - nb))[1, 1])
    end
    # ---- 误差标度（τ = β/2 环对称点）----
    open(joinpath(@__DIR__, "greens_imagtime_error.dat"), "w") do io
        println(io, "# order  Ntau  abs_err")
        for order in [1, 2, 4, 6]
            prev = nothing
            for N in [5, 10, 20, 40]
                gref = Gref(β / 2)
                err = abs(G_trot(order, N, N ÷ 2) - gref)
                slope = prev === nothing ? "" : "  斜率 = $(round(log(prev/err)/log(2), digits=2))"
                println("order=$order Nτ=$N: |err| = ", round(err, sigdigits=3), slope)
                println(io, "$order $N $err")
                prev = err
            end
        end
    end
    println("\n参照 G(τ=β/2) = ", round(Gref(β/2), digits=6))
end
main()
