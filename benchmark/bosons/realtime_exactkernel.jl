# =============================================================================
# benchmark/bosons/realtime_exactkernel.jl —— 精确核（exact-cell IF）实时求解器
#
# 参考：docs/NgReichman2023.pdf（PRB 107, 125103）——"对 Trotter 化动力学精确
# 求 IF"的思想：不近似 B = 浴+耦合 的 cell 演化，而是把 bath 自由度精确
# （Schur）消掉，得到 Trotter 链的**精确** IF 核。装配后对链机器精度，从而
# 直接继承组合分裂的阶数（Strang → 2 阶，Suzuki → 4 阶），无需 Richardson。
#
# 单粒子层面（Toulouse 玻色）：W = e^{-ihB} 分块（w₀; w; v; M），节点递推
#   a_m = u·w₀·a_{m-1} + Σ_{d≥2} u·(w M^{d-2} v)·a_{m-d},   u = e^{-iϵ_d h}
# （trotter_frozen.jl [0](a) 已验证对严格 Strang 链机器精度）。本文件把它
# 做成可复用的求解器：
#   · 精确核生成元 exactkernel_gens：(w₀, {κ(d)}) —— **bath-only**（不含 ϵ_d），
#     即 IF 的 H_S 无关性：一次计算，任意系统哈密顿量复用（[0](b) 实测）；
#   · 节点装配 exactkernel_rt：与 frozen_rt_shifted 同骨架、换核（N×N 三角
#     系统，O(N²)）；
#   · 块行进 exact_cell_march：W 的 (1+L) 分块直接行进，支持任意宽度布局
#     （含 Suzuki 的负 γ 步），精确 cell 按宽度缓存（O(L³)·宽度数 + O(L²)·步数）。
# 冻结核 = κ 的 V² 截断（= −u·G[d−1]，顶点移位装配；V⁴ 亏损的完整分析见
# trotter_frozen.jl）。费米子在包内已有同族实现：exactfermionic.jl 的
# exact_fermionic_Δτ/Δt（Schur：Δ_exact = (w₀−1)·S_kinetic − C·M⁻¹·C）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. benchmark/bosons/realtime_exactkernel.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "performance", "common", "trotter_frozen.jl"))   # 传递引入 common.jl；复用 frozen_rt_shifted 作对照

# ------------------------------------------------------------------ #
# 精确核生成元（bath-only）：κ(1)=w₀（次对角），κ(d)=w·M^{d-2}·v（d≥2）
# 迭代 urow ← urow·M，每 lag O(L²)；不含 ϵ_d（IF 的 H_S 无关性）
# ------------------------------------------------------------------ #
function exactkernel_gens(mdl, N, T)
    h = T / N
    W = exp(-im * h * (mdl.hb + mdl.hi))
    w₀ = W[1, 1]
    v = W[2:end, 1]
    M = W[2:end, 2:end]
    urow = reshape(W[1, 2:end], 1, :)
    κ = zeros(ComplexF64, N)                # κ[d]（d≥2 有值）
    for d in 2:N
        κ[d] = (urow * v)[1]
        urow = urow * M
    end
    return (; w₀, κ, h)
end

"""
    exactkernel_rt(gens, ϵd, N, T) -> (a[1..N], h)

精确核 + 顶点移位节点装配（与 frozen_rt_shifted 同骨架，核换成 w₀/κ）：

    a_m − u·w₀·a_{m−1} − Σ_{d≥2} u·κ(d)·a_{m−d} = 0,  a₀ = 1

对严格 Strang 链机器精度（[0](a)）；u 相位与源项由 ϵd 提供——
gens 本身与 ϵd 无关，可跨 H_S 复用（[0](b)）。
"""
function exactkernel_rt(gens, ϵd, N, T)
    h = T / N
    u = exp(-im * ϵd * h)
    A = Matrix{ComplexF64}(I, N, N)
    b = zeros(ComplexF64, N)
    for m in 1:N
        if m ≥ 2
            A[m, m-1] -= u * gens.w₀
        else
            b[m] += u * gens.w₀             # 耦合到已知初值 a₀=1
        end
        for d in 2:m
            if m - d ≥ 1
                A[m, m-d] -= u * gens.κ[d]
            else
                b[m] += u * gens.κ[d]       # d=m 项耦合 a₀=1
            end
        end
    end
    return A \ b, h
end

# ------------------------------------------------------------------ #
# 块行进：精确 cell W(w) 的 (1+L) 分块直接行进（半步自由相位夹逼）。
# 宽度可任取（含 Suzuki 负 γ 步）；W 按宽度缓存 ⟹ O(L³)×宽度数 + O(L²)×步数。
# ------------------------------------------------------------------ #
function exact_cell_march(mdl, widths)
    cache = Dict{Float64,Matrix{ComplexF64}}()
    W(w) = get!(cache, w) do
        exp(-im * w * (mdl.hb + mdl.hi))
    end
    a = 1.0 + 0.0im
    β = zeros(ComplexF64, mdl.L)
    for w in widths
        fp = exp(-im * mdl.ϵd * w / 2)
        a *= fp
        Ww = W(w)
        anew = Ww[1, 1] * a + sum(Ww[1, 2:end] .* β)   # 吸收：无共轭
        β = Ww[2:end, 1] * a + Ww[2:end, 2:end] * β
        a = anew * fp
    end
    return a
end

const SUZUKI_W5 = (SUZUKI_θ, SUZUKI_θ, SUZUKI_γ, SUZUKI_θ, SUZUKI_θ)   # 4θ+γ=1

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[8, 16, 32, 64], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Aex = exact_return_amplitude(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "精确核（exact-cell IF，NgReichman2023 思想）：实时 Toulouse 求解器")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)

    # ---------------- [0] 结构自检 ----------------
    println(io, "\n[0] 结构自检")
    # (a) 三方一致：核装配 vs 块行进 vs 严格 Strang 链（机器精度）
    for N in (8, 16)
        gens = exactkernel_gens(mdl, N, T)
        a, _ = exactkernel_rt(gens, mdl.ϵd, N, T)
        march = exact_cell_march(mdl, fill(T / N, N))
        st = real_step(mdl, 2, T / N)
        chain = [(st^m)[1, 1] for m in 1:N]
        @printf(io, "    (a) N=%2d: |装配−链|=%.1e  |行进−链|=%.1e  （机器精度）\n",
                N, maximum(abs.(a .- chain)), abs(march - chain[N]))
    end
    # (b) IF 复用（H_S 无关性）：同一批 gens 装配不同 ϵ_d，各对相应严格链
    let N = 16
        gens1 = exactkernel_gens(mdl, N, T)                    # 由 ϵ_d=1.0 的模型算出
        mdl2 = boson_toulouse(; Lmodes=Lmodes, ϵd=0.6)
        gens2 = exactkernel_gens(mdl2, N, T)
        @printf(io, "    (b) 核的 H_S 无关性：两个 ϵ_d 的 κ(4) 之差 = %.1e；",
                abs(gens1.κ[4] - gens2.κ[4]))
        a1, _ = exactkernel_rt(gens1, mdl2.ϵd, N, T)           # 复用 gens1 解 ϵ_d=0.6
        st2 = real_step(mdl2, 2, T / N)
        chain2 = [(st2^m)[1, 1] for m in 1:N]
        @printf(io, "复用解 vs 链：max|Δ| = %.1e\n", maximum(abs.(a1 .- chain2)))
    end
    # (c) 冻结核 = V² 截断：|κ(4) + u·G[3]| 的 Γ 减半 ÷4 标度（详见 trotter_frozen [0](b)）
    print(io, "    (c) 冻结核 = κ 的 V² 截断：")
    for Γ in (0.15, 0.075)
        m2 = boson_toulouse(; Lmodes=Lmodes, Γ=Γ)
        gens = exactkernel_gens(m2, 64, 64 * 0.125)
        corr = bosonic_Δt(m2.bath.spectrum; β=m2.β, N=64, t=64 * 0.125, μ=m2.μ)
        P = branch(corr; b1=:+, b2=:+); X = branch(corr; b1=:+, b2=:-)
        G3 = -(P[4, 1] + X[4, 1])
        u = exp(-im * m2.ϵd * 0.125)
        @printf(io, "Γ=%.3f → δκ=%.1e；", Γ, abs(u * gens.κ[4] + u * G3))
    end
    println(io, "（÷4 ⟹ V⁴）")

    # ---------------- [1] 阶数阶梯 ----------------
    println(io, "\n[1] 实时 A(t)：精确核装配直接继承分裂阶数（无需 Richardson）")
    es = Float64[]
    for N in Ns
        gens = exactkernel_gens(mdl, N, T)
        vals = exactkernel_rt(gens, mdl.ϵd, N, T)[1][ms_of(N)]
        push!(es, maximum(abs.(vals .- [Aex(m * T / N) for m in ms_of(N)])))
    end
    order_table(io, "    精确核 + Strang 节点装配", Ns, es)
    es = Float64[]
    for N in Ns
        st = real_step(mdl, 2, T / N)
        push!(es, maximum(abs.([(st^m)[1, 1] for m in ms_of(N)] .-
                               [Aex(m * T / N) for m in ms_of(N)])))
    end
    order_table(io, "    严格 cell 链 order=2（参照，应逐位同上）", Ns, es)
    esA = Float64[]; esB = Float64[]
    for N in Ns
        h = T / N
        push!(esA, abs(exact_cell_march(mdl, fill(h, N)) - Aex(T)))
        push!(esB, abs(exact_cell_march(mdl, repeat(collect(SUZUKI_W5), N) .* h) - Aex(T)))
    end
    order_table(io, "    精确核块行进（均匀，终点）", Ns, esA)
    order_table(io, "    精确核块行进 Suzuki(θ,θ,γ,θ,θ)（终点）", Ns, esB)
    es = Float64[]
    for N in Ns
        st = real_step(mdl, 4, T / N)
        push!(es, abs((st^N)[1, 1] - Aex(T)))
    end
    order_table(io, "    严格 cell 链 order=4（参照，应同上）", Ns, es)
    es = Float64[]
    for N in Ns
        vals = frozen_rt_shifted(mdl, N, T)[1][ms_of(N)]
        push!(es, maximum(abs.(vals .- [Aex(m * T / N) for m in ms_of(N)])))
    end
    order_table(io, "    冻结核装配（对照：V² 截断，1 阶）", Ns, es)

    # ---------------- [2] 成本与结论 ----------------
    println(io, "\n[2] 成本（N=2048，各重复取均值）与结论")
    let N = 2048, reps = 5
        tg = @elapsed (for _ in 1:reps; exactkernel_gens(mdl, N, T); end)
        gens = exactkernel_gens(mdl, N, T)
        ts = @elapsed (for _ in 1:reps; exactkernel_rt(gens, mdl.ϵd, N, T); end)
        tm = @elapsed (for _ in 1:reps; exact_cell_march(mdl, fill(T / N, N)); end)
        @printf(io, "    核生成 O(N·L²)=%.1fms + 装配求解 O(N²)=%.1fms；块行进 O(N·L²)=%.1fms\n",
                1e3tg / reps, 1e3ts / reps, 1e3tm / reps)
    end
    println(io, "\n结论（实测）：")
    println(io, "  · 精确核（w₀, κ(d)=wM^{d-2}v）为 bath-only：一次计算，跨 ϵ_d 复用（[0](b)），")
    println(io, "    即 NgReichman2023 的 exact-IF 对象在玻色单粒子层面的实现（费米子：exactfermionic.jl）。")
    println(io, "  · 装配/行进对严格 Strang 链机器精度（[0](a)）⟹ Strang → 真 2 阶；Suzuki(θ,θ,γ,θ,θ)")
    println(io, "    块行进 → 真 4 阶——分裂阶数被精确核完整继承，无需 Richardson。")
    println(io, "  · 冻结核 = κ 的 V² 截断，渐近 1 阶（对照行）；其 2/4 阶路径 = Richardson（trotter_frozen.jl）")
    println(io, "    或换用本文件的精确核。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
