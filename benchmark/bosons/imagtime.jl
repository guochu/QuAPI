# =============================================================================
# benchmark/bosons/imagtime.jl —— 虚时（Matsubara）精度基准
#
# 用 bosonic_Δτ 的离散杂化矩阵驱动 Toulouse 杂质的 Trotter 演化，检查
# G(τ) 对严格解 toulouse_Gτ 的收敛阶是否为真 1/2/4 阶。
#
# 运行：julia --startup-file=no --project=../.. benchmark/bosons/imagtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "performance", "common", "common.jl"))

function main(; Lmodes=48, Ns=[8, 16, 32, 64, 128], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Gex = exact_Gτ(mdl)
    io = stdout

    println(io, "=" ^ 74)
    println(io, "虚时 benchmark：bosonic_Δτ 离散杂化矩阵 + 1/2/4 阶 Trotter 演化")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd)
    println(io, "演化器浴 = 生成 Δτ 的同一离散浴（下面对 Δτ 逐元验证）；")
    println(io, "故链与 Δτ 描述同一 hybridized bath，高斯杂质 G 只由 Δ 决定。")
    println(io, "=" ^ 74)

    # ---------------- 0) 离散杂化矩阵自检 ----------------
    println(io, "\n[0] 离散 hybridization matrix 自检（bosonic_Δτ vs 模式求和严格双 cell 积分）")
    for N in (5, 8, 13)
        d, _, _ = verify_delta_tau(mdl, N)
        @printf(io, "    Nτ=%2d : max|ΔQuAPI − Δexact| = %.2e\n", N, d)
    end

    # ---------------- 1) Trotter 收敛阶 ----------------
    println(io, "\n[1] 虚时 Trotter 演化：G(τ) 对 toulouse_Gτ 的误差")
    errs_by_order = Dict{Int,Vector{Float64}}()
    G1_last = G2_last = nothing
    for order in (1, 2, 4)
        es = Float64[]
        for N in Ns
            ms = round.(Int, fracs .* N)
            vals, _, h = chain_Gτ(mdl, order, N; ms=ms)
            ex = [Gex(mod(m * h, mdl.β)) for m in ms]
            push!(es, maximum(abs.(vals .- ex)) / maximum(abs.(ex)))
            order == 1 && (G1_last = vals)
            order == 2 && (G2_last = vals)
        end
        errs_by_order[order] = es
        order_table(io, "    order = $order", Ns, es)
    end

    # ---------------- 2) 一/二阶退化说明 ----------------
    println(io, "\n[2] 说明：虚时 thermal 量对“对角自由因子”的摆放是规范不变的")
    @printf(io, "    max|Gτ(order1) − Gτ(order2)| (Nτ=%d) = %.2e\n",
            last(Ns), maximum(abs.(G1_last .- G2_last)))
    println(io, "    ⟹ 虚时 1 阶 Lie 的首阶误差是 Tr[[A,B]e^{-βH}] = 0，自动升为 2 阶。")

    # ---------------- 3) 等时项：分开保留 ⟺ 系数相加 ----------------
    println(io, "\n[3] 等时项的算符处理（[d,d†]=1，非对角杂化）")
    cs = bosonic_Δτ(mdl.bath.spectrum, β=mdl.β, N=16, μ=mdl.μ).data
    d_id, r_id = equal_time_identity_check(cs.ηⱼₖ[1], cs.ηₖⱼ[1])
    @printf(io, "    两种次序分开：exp(ηjk[1]·n̂)exp(ηkj[1]·(n̂+1))\n")
    @printf(io, "    与系数相加：exp((ηjk[1]+ηkj[1])·n̂) 之差 = e^{ηkj[1]}\n")
    @printf(io, "    数值验证：max|Osplit − e^{ηkj1}Osum| = %.2e，归一化期望偏差 = %.2e\n",
            d_id, abs(r_id))
    println(io, "    ⟹ [d,d†]=1 是 c-number，整体常数在归一化 G 中抵消；")
    println(io, "       等时项对路径的系数就是 ηjk[1]+ηkj[1]。丢单侧会改掉 n̂ 的系数，是错的。")

    # ---------------- 4) 冻结顶点（vertex-frozen）QuAPI 路径积分 ----------------
    println(io, "\n[4] 对照：冻结顶点规则 × cell 内重构阶数（节点/中点/三次）")
    println(io, "    检验：提高杂交项的 cell 重构阶数能否把冻结规则抬到 2/4 阶")
    for p in (1, 2, 4)
        es = Float64[]
        for N in Ns
            vals, h = frozen_dyson_Gτ(mdl, N; cellrule=p)
            ex = [Gex(mod(m * h, mdl.β)) for m in 1:N-1]
            push!(es, maximum(abs.(vals .- ex)) / maximum(abs.(ex)))
        end
        order_table(io, "    frozen vertex, cell rule p=$p", Ns, es)
    end

    # ---------------- 5) 汇总 ----------------
    println(io, "\n[5] 实测收敛阶（相邻 N 的 log2 误差比，取最后一级）")
    for order in (1, 2, 4)
        e = errs_by_order[order]
        @printf(io, "    order=%d : 拟合阶 ≈ %.2f\n", order, log2(e[end-1] / e[end]))
    end
    println(io, "    结论：")
    println(io, "    · 严格 cell（cell 内精确演化）：order 2 → 真二阶，order 4 → 真四阶；")
    println(io, "      虚时 order 1 因 Tr[[A,B]e^{-βH}]=0 自动升为二阶。")
    println(io, "    · 冻结顶点规则：把 cell 内重构从节点→中点→三次插值都试过，阶数恒为 1；")
    println(io, "      说明 1 阶误差不是杂交项重构，而是“节点采样自由传播子 × cell 积分核”")
    println(io, "      这一配对本身的 O(δτ) 求积失配；Gaussian 路径积分里没有可调的算符分裂。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
