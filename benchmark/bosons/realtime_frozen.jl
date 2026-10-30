# =============================================================================
# benchmark/bosons/realtime_frozen.jl —— 实时（实轴、直积初态）冻结核收敛阶
#
# 与 realtime.jl 的区别：这里初态是直积态 ρ_S(n_sys) ⊗ ρ_bath(β)，实轴演化，
# 严格解用 toulouse_neq_greater_lesser（非平衡 G>、G<），而不是混合轮廓的平衡解。
#
# 关键恒等式（直积初态）：G>(t) = -i(1+n_sys)A(t)，G<(t) = +i n_sys A(t)，
# A(t) = [e^{-iHt}]₁₁。故只测 A(t) 即可，冻结核求解器见 common.jl。
#
# 冻结核的杂交核直接取自 QuAPI 的 bosonic_Δt：K^R = −(Δ⁺⁺ + Δ⁺₋)（因果
# 下三角，含等时），是精确恒等式（热因子 (1+n)−n 相消，等时为 t>t′ 三角形）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. benchmark/bosons/realtime_frozen.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "performance", "common", "common.jl"))

function main(; Lmodes=48, nsys=0.3, T=2.0, Ns=[16, 32, 64, 128])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Aex = exact_return_amplitude(mdl)
    gl = toulouse_neq_greater_lesser(mdl.model; nsys=nsys)
    io = stdout
    println(io, "=" ^ 74)
    println(io, "实时（实轴、直积初态 n_sys=$nsys）冻结核收敛阶")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)

    # ---------------- 0) 直积初态恒等式自检 ----------------
    println(io, "\n[0] 直积初态恒等式自检：G> = -i(1+n_sys)A，G< = +i n_sys A")
    e = 0.0
    for t in (0.3, 0.7, 1.1, 1.5, 1.9)
        g, l = gl(t); a = Aex(t)
        e = max(e, abs(g + im*(1+nsys)*a), abs(l - im*nsys*a))
    end
    @printf(io, "    max|G>_exact + i(1+n_sys)A|, |G<_exact − i n_sys A| = %.2e\n", e)

    # ---------------- 1) 冻结核阶数 ----------------
    println(io, "\n[1] 冻结核（推迟 Dyson，cell 双积分核）对精确 A(t) 的误差")
    es = Float64[]
    for N in Ns
        vals, h = frozen_retarded_A(mdl, N, T)
        ms = round.(Int, [0.25, 0.5, 0.75] .* N)
        push!(es, maximum(abs.(vals[ms .+ 1] .- [Aex(m*h) for m in ms])) / maximum(abs.(Aex(m*h) for m in ms)))
    end
    order_table(io, "    frozen retarded", Ns, es)

    # ---------------- 2) 严格 cell 对照 ----------------
    println(io, "\n[2] 对照：严格 cell 实时链（同一 A(t)）")
    for order in (1, 2, 4)
        eo = Float64[]
        for N in Ns
            vals, h = chain_return_amplitude(mdl, order, N, T)
            ms = round.(Int, [0.25, 0.5, 0.75] .* N)
            push!(eo, maximum(abs.(vals[ms .+ 1] .- [Aex(m*h) for m in ms])) / maximum(abs.(Aex(m*h) for m in ms)))
        end
        order_table(io, "    strict cell, order=$order", Ns, eo)
    end

    # ---------------- 3) 核恒等式：−(Δ⁺⁺+Δ⁺₋) = 态无关推迟核（闭式） ----------------
    println(io, "\n[3] 核恒等式验证：−(Δ⁺⁺+Δ⁺₋) vs 逐模式闭式推迟核（N=8，下三角含等时）")
    N = 8; h = T / N
    corr = bosonic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+)
    X = branch(corr; b1=:+, b2=:-)
    closed(lag) = sum(α * (lag == 0 ? (1 - im*w*h - exp(-im*w*h)) / w^2 :
                                   exp(-im*w*lag*h) * 2*(1 - cos(w*h)) / w^2)
                      for (w, α) in zip(mdl.ws, mdl.αs))
    d = 0.0; sc = 0.0
    for lag in 0:N
        v = -(P[lag+1, 1] + X[lag+1, 1])
        d = max(d, abs(v - closed(lag))); sc = max(sc, abs(closed(lag)))
    end
    @printf(io, "    max|−(Δ⁺⁺+Δ⁺₋) − K_closedform| = %.2e   (核尺度 %.2e, 相对 %.1e)\n",
            d, sc, d / sc)
    println(io, "    Δ⁺⁺ = −(1+n)R, Δ⁺₋ = +nR ⟹ −(Δ⁺⁺+Δ⁺₋) = R：热因子严格相消，")
    println(io, "    与 retarded 核的态无关性一致；等时对角同为 t>t′ 三角形约定。")

    println(io, "\n结论：实轴直积初态、核直接取自 bosonic_Δt 的 −(Δ⁺⁺+Δ⁺₋)，冻结核为真 1 阶；")
    println(io, "      真 2/4 阶来自严格 cell 链；冻结核到 2/4 阶的组合路径（顶点移位装配 +")
    println(io, "      Richardson 外推）见 trotter_frozen.jl。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
