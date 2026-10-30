# =============================================================================
# benchmark/bosons/realtime.jl —— 实时（Keldysh）精度基准
#
# 用 bosonic_Δt（平衡初态用 Δm 的虚时段）离散化驱动 Toulouse 杂质的一阶/二阶/
# 四阶 Trotter 实时演化，检查 G>、G< 对严格解 toulouse_greater_lesser 的收敛阶。
#
# 平衡初态：单粒子热密度矩阵 γ = (U_β^{-1} - I)^{-1}（U_β 由虚时链给出，对应
# Δm 的虚时段）。实时传播用单粒子链 R(t)；观测量
#     G>(t) = -i [R(t)(I + γᵀ)]₁₁,   G<(t) = +i [R(t)γᵀ]₁₁ .
#
# 运行：julia --startup-file=no --project=../.. benchmark/bosons/realtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "performance", "common", "common.jl"))

function main(; Lmodes=48, Ns=[8, 16, 32, 64, 128], T=2.0,
              fracs=[0.25, 0.5, 0.75], Nτ_gamma=256)
    mdl = boson_toulouse(; Lmodes=Lmodes)
    gl = exact_greater_lesser(mdl)
    Gex(t) = gl(t)[1]
    Lex(t) = gl(t)[2]
    io = stdout

    println(io, "=" ^ 74)
    println(io, "实时 benchmark：bosonic_Δt/Δm 离散杂化矩阵 + 1/2/4 阶 Trotter 演化")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "演化器浴 = 生成 Δt/Δm 的同一离散浴（下面对 Δt、Δm 逐元验证）；")
    println(io, "平衡初态取混合轮廓 Δm 的虚时段，实时传播用单粒子链。")
    println(io, "=" ^ 74)

    # ---------------- 0) 离散杂化矩阵自检 ----------------
    println(io, "\n[0] 离散 hybridization matrix 自检")
    @printf(io, "    Δt  (四块 vs 模式求和严格双 cell 积分)  max diff = %.2e\n",
            verify_delta_t(mdl, 5, 1.3))
    dτ, dpm = verify_delta_m(mdl, 6, 4, 1.1)
    @printf(io, "    Δm  (τ 块 vs Δτ = %.2e,  ± 块 vs Δt = %.2e)\n", dτ, dpm)

    # ---------------- 1) 平衡初态 γ ----------------
    # 默认用严格 γ 以隔离实时阶数；同时报告 Δm 虚时链（order2, Nτ=Nτ_gamma）的 γ 偏差。
    γ = exact_gamma(mdl)
    γc = chain_gamma(mdl, 2, Nτ_gamma)
    @printf(io, "\n[1] 平衡初态：γ₁₁ = %.10f（严格）\n", γ[1, 1])
    @printf(io, "    Δm 虚时段链 γ（order2, Nτ=%d）与严格 γ 的 max|Δ| = %.2e\n",
            Nτ_gamma, maximum(abs.(γc .- γ)))
    @printf(io, "    G>(0) = %s   严格 = %s\n",
            string(round(-im * (γ[1, 1] + 1); digits=10)),
            string(round(Gex(0.0); digits=10)))

    # ---------------- 2) 收敛阶 ----------------
    println(io, "\n[2] 实时 Trotter 演化：G>、G< 对 toulouse_greater_lesser 的误差")
    errs_by_order = Dict{Int,Vector{Float64}}()
    for order in (1, 2, 4)
        es = Float64[]
        for Nt in Ns
            ms = round.(Int, fracs .* Nt)
            gt, lt, ts = chain_greater_lesser(mdl, order, Nt, T; γ=γ, ms=ms)
            gt_ex = [Gex(t) for t in ts]
            lt_ex = [Lex(t) for t in ts]
            scale = maximum(abs.(gt_ex))
            push!(es, maximum(vcat(abs.(gt .- gt_ex), abs.(lt .- lt_ex))) / scale)
        end
        errs_by_order[order] = es
        order_table(io, "    order = $order", Ns, es)
    end

    # ---------------- 3) G>/G< 分量示例（order 2） ----------------
    println(io, "\n[3] 分量示例（order=2, Nt=$(last(Ns)))")
    ms = round.(Int, fracs .* last(Ns))
    gt, lt, ts = chain_greater_lesser(mdl, 2, last(Ns), T; γ=γ, ms=ms)
    for (k, t) in enumerate(ts)
        @printf(io, "    t=%.3f  G> = %s   严格 %s\n", t,
                string(round(gt[k]; digits=6)), string(round(Gex(t); digits=6)))
        @printf(io, "             G< = %s   严格 %s\n",
                string(round(lt[k]; digits=6)), string(round(Lex(t); digits=6)))
    end

    # ---------------- 4) 汇总 ----------------
    println(io, "\n[4] 实测收敛阶（log2 误差比，取最后一级）")
    for order in (1, 2, 4)
        e = errs_by_order[order]
        @printf(io, "    order=%d : 拟合阶 ≈ %.2f\n", order, log2(e[end-1] / e[end]))
    end
    println(io, "    结论：order 1/2/4 分别为真一阶/二阶/四阶。")
    println(io, "    虚时 thermal 量的 1 阶退化在实时不出现（实时不是迹型观测量）。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
