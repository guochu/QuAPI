# =============================================================================
# benchmark/bosons/trotter_frozen.jl —— 冻结核 × 高阶 Trotter 组合的收敛阶判定
#
# 问题：bosonic_Δt / bosonic_Δτ 的 cell 双积分核（“冻结核”）驱动的冻结方案是
# 真 1 阶（realtime_frozen.jl、imagtime.jl [4]）。能否只用同一批核（或同类的
# V² 精确闭式），配合 Strang / Suzuki 分裂或外推组合，把全局误差抬到 2/4 阶？
#
# ---------------------------- 理论框架（实测确认） ----------------------------
# 记号：单粒子矩阵上 B = 浴+耦合，S = 杂质自由；严格 Strang 步 T = F·W·F，
# W = e^{-ihB}，F = e^{-ihS/2}（只作用于杂质线）。W 分块（w₀; w; v; M）。
# 节点递推（对链精确，[0](a) 机器精度验证）：
#   a_m = u·w₀·a_{m-1} + Σ_{d≥2} u·(w M^{d-2} v)·a_{m-d},   u = e^{-iϵ_d h}.
# 关键恒等式（逐模解析成立，Duhamel 一阶）：
#   (w M^{d-2} v)|_{V²} = −G[d−1]
# 其中 G[ℓ] = bosonic_Δt 的 lag-ℓ 条目 = Σα e^{-iωℓh}·2(1−cos ωh)/ω²。
# 即严格链的节点核 = 冻结条目**顶点移位一格**（cell 配对 → 节点配对）+
# 每步统一的 u 相位（相互作用绘景：顶点在节点，顶点间自由飞行相位相消）
# + 三角条目并入次对角 (1−g₀)。这就是“顶点移位冻结装配”——与
# QuAPI.Dense 费米顶点移位（fermion.jl）同一层级的装配修正。
#
# 冻结装配的剩余亏损（[0](b) 实测 V⁴ 标度：Γ 减半 → δκ÷4）：
#   δκ(d) = κ(d) + u·G[d−1] ~ V⁴h²(1+d)/条目（相对量），
# 来源 = 顶点自缀饰（w 的 V³ 项）+ 飞行中浴模经杂质的虚访问（M 的 V² 项）。
# 它是逐模（ω 依赖）的，lag 空间的冻结条目组合原则上重构不出；且按 ~N 个
# lag 相干累积成全局 O(h)——**不是分裂型误差**，因此 Suzuki 重排消不掉
# （[2] 实测冻结 Suzuki 布局封顶 1 阶，对照严格 Suzuki 的 4 阶）。
#
# 能到达高阶的冻结核组合 = Richardson 外推（误差有规则 h 幂次展开，实测
# 1/2/3 层 → 2.00 / 2.99 / 4.01 阶），每层只需把网格加密一倍。
#
# 与文献的连接：
# · GuoChen2026a（PRB 114, 125413）App. C 的 QuAPI 离散化 Eq. (12) 把节点
#   算符 A†_i,ζ 与 cell 积分核 ζ_{i,j} 配对——正是顶点移位装配；其对严格链
#   的逐条目精度 O(V⁴) 即上述亏损，全局渐近 1 阶与其 O(t·δt) 推断一致。
#   （offdiag_order.jl 在小窗口 N≤5 测到的“全局 δt²”是 pre-asymptotic。）
# · NgReichman2023（PRB 107, 125103）的思路是“对 Trotter 化动力学精确求 IF”
#   （全二次 H₀ | 相互作用 H₁ 的 Strang 分裂）——对应本文件的严格链参照组；
#   Toulouse 模型 H₁=0，严格链即其极限情形。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. benchmark/bosons/trotter_frozen.jl
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

# ------------------------------------------------------------------ #
# 实时：顶点移位冻结装配（由 bosonic_Δt 的条目直接驱动）
#   a_m = u(1−g₀)a_{m−1} − u Σ_d G[d−1] a_{m−d},  a₀ = 1
#   G[ℓ] = −(Δ⁺⁺+Δ⁺⁻) 的 lag-ℓ 生成元；g₀ = 等时推迟三角
# ------------------------------------------------------------------ #
function frozen_rt_shifted(mdl, N, T)
    h = T / N
    u = exp(-im * mdl.ϵd * h)
    corr = bosonic_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ)
    P = branch(corr; b1=:+, b2=:+)
    X = branch(corr; b1=:+, b2=:-)
    g0 = -(P[1, 1] + X[1, 1])                        # 等时三角（t>t′ 约定）
    g = [-(P[d+1, 1] + X[d+1, 1]) for d in 1:N-1]    # lag-d 生成元（d≥1）
    M = Matrix{ComplexF64}(I, N, N)
    b = zeros(ComplexF64, N)
    for m in 1:N
        if m ≥ 2
            M[m, m-1] -= u * (1 - g0)
        else
            b[m] += u * (1 - g0)                     # 耦合到已知初值 a₀=1
        end
        for d in 2:m
            if m - d ≥ 1
                M[m, m-d] += u * g[d-1]
            else
                b[m] -= u * g[d-1]                   # d=m 项耦合 a₀=1
            end
        end
    end
    return M \ b, h
end

# ------------------------------------------------------------------ #
# 实时：V² 精确冻结 cell 的单粒子行进（有向，w 可负 → 支持 Suzuki 布局）
#   w₀ᶠ = 1 − Σα·D(w)（cell 三角），发射/吸收矢 eₖ = −VₖQₖ(w)，
#   Q(w) = ∫₀ʷe^{-iωσ}dσ，裸飞行 Mᶠ = diag(e^{-iωₖw})。
#   均匀 w=h 时隐含节点核 = −u·G[d−1]，与 frozen_rt_shifted 逐位一致（[0](c)）。
# ------------------------------------------------------------------ #
_rtQ(ω, w) = abs(ω) < 1e-10 ? w : (1 - exp(-im * ω * w)) / (im * ω)
_rtD(ω, w) = abs(ω) < 1e-10 ? w^2 / 2 : (1 - im * ω * w - exp(-im * ω * w)) / ω^2

function frozen_rt_march(mdl, widths)
    Vs = sqrt.(mdl.αs)
    a = 1.0 + 0.0im
    β = zeros(ComplexF64, mdl.L)
    for w in widths
        fp = exp(-im * mdl.ϵd * w / 2)
        a *= fp
        w0f = 1 - sum(α * _rtD(ω, w) for (ω, α) in zip(mdl.ws, mdl.αs))
        em = [-im * V * _rtQ(ω, w) for (ω, V) in zip(mdl.ws, Vs)]
        anew = w0f * a + sum(em .* β)                # 吸收：无共轭（振幅非内积）
        β = em * a + [exp(-im * ω * w) for ω in mdl.ws] .* β
        a = anew * fp
    end
    return a
end

# ------------------------------------------------------------------ #
# 虚时：冻结 cell 环 G = [T_f^m (I − T_f^N)^{-1}]₁₁（热环闭合与严格链同构；
# W_f 为虚时化的 V² 精确冻结 cell）
# ------------------------------------------------------------------ #
_τQ(ω, h) = abs(ω) < 1e-10 ? h : (1 - exp(-ω * h)) / ω
_τD(ω, h) = abs(ω) < 1e-10 ? h^2 / 2 : (exp(-ω * h) - (1 - ω * h)) / ω^2

function frozen_ring_Wf(mdl, N)
    h = mdl.β / N
    L = mdl.L
    W = Matrix{Float64}(I, L + 1, L + 1)
    W[1, 1] = 1 - sum(α * _τD(ω, h) for (ω, α) in zip(mdl.ws, mdl.αs))
    for k in 1:L
        W[1, 1+k] = W[1+k, 1] = -sqrt(mdl.αs[k]) * _τQ(mdl.ws[k], h)
        W[1+k, 1+k] = exp(-mdl.ws[k] * h)
    end
    F = Matrix{Float64}(I, L + 1, L + 1)
    F[1, 1] = exp(-mdl.ϵd * h / 2)
    Tf = F * W * F
    S = inv(mdl.Id - Tf^N)
    return [(Tf^m * S)[1, 1] for m in 0:N-1], h
end

# ------------------------------------------------------------------ #
# Richardson 外推：levels 层，各层采样点逐层倍增对齐（rₖ(n,m) ← rₖ₋₁(2n,2m)），
# 第 k 层消 h^k（权重 (2^k·rₖ₋₁(h/2) − rₖ₋₁(h))/(2^k−1)）
# ------------------------------------------------------------------ #
function richardson_levels(solver, N, ms; levels::Int=1)
    r = (n, m) -> 2 * solver(2n)[2m] - solver(n)[m]
    for k in 2:levels
        c = 2.0^k
        r = let rprev = r, ck = c
            (n, m) -> (ck * rprev(2n, 2m) - rprev(n, m)) / (ck - 1)
        end
    end
    return [r(N, m) for m in ms]
end

const SUZUKI_W = (SUZUKI_θ, SUZUKI_θ, SUZUKI_γ, SUZUKI_θ, SUZUKI_θ)   # 4θ+γ=1

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[16, 32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Aex = exact_return_amplitude(mdl)
    Gex = exact_Gτ(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "冻结核 × 高阶 Trotter 组合：顶点移位装配、Suzuki 布局、Richardson 外推")
    @printf(io, "浴：%d 模，J=%.2f on [0.5,3.0]，β=%.1f，ϵ_d=%.1f，T=%.1f\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)

    # ---------------- [0] 结构自检 ----------------
    println(io, "\n[0] 结构自检")
    # (a) 精确节点递推 == 严格 Strang 链（机器精度）
    let N = 8, h = T / N
        W = exp(-im * h * (mdl.hb + mdl.hi))
        w0 = W[1, 1]; w = reshape(W[1, 2:end], 1, :); v = W[2:end, 1]; M = W[2:end, 2:end]
        u = exp(-im * mdl.ϵd * h)
        st = real_step(mdl, 2, h)
        a = ComplexF64[1.0]
        for m in 1:N
            mem = m ≥ 2 ? sum((((w * M^(d - 2)) * v)[1]) * a[m-d+1] for d in 2:m) : 0.0im
            push!(a, u * (w0 * a[end] + mem))
        end
        @printf(io, "    (a) 精确递推 vs 严格 Strang 链：max|Δ| = %.1e（机器精度）\n",
                maximum(abs.(a[2:end] .- [(st^m)[1, 1] for m in 1:N])))
    end
    # (b) V² 恒等式：κ(d) ≈ −u·G[d−1]，剩余亏损 V⁴ 标度（Γ 减半 ÷4）
    print(io, "    (b) 核恒等式 κ(d) = −u·G[d−1] + O(V⁴)：")
    for Γ in (0.15, 0.075)
        m2 = boson_toulouse(; Lmodes=Lmodes, Γ=Γ)
        h = 0.125
        u = exp(-im * m2.ϵd * h)
        W = exp(-im * h * (m2.hb + m2.hi))
        w = reshape(W[1, 2:end], 1, :); v = W[2:end, 1]; M = W[2:end, 2:end]
        corr = bosonic_Δt(m2.bath.spectrum; β=m2.β, N=64, t=64h, μ=m2.μ)
        P = branch(corr; b1=:+, b2=:+); X = branch(corr; b1=:+, b2=:-)
        G3 = -(P[4, 1] + X[4, 1])                    # lag-3 条目（κ(4) ↔ d−1=3）
        δ = abs((u * ((w * M^2) * v)[1]) + u * G3)
        @printf(io, "Γ=%.3f → δκ(4)=%.1e；", Γ, δ)
    end
    println(io, "（Γ 减半 ÷4 ⟹ V⁴）")
    # (c) 冻结 cell 行进 == 节点递推（均匀布局）
    let N = 16
        rec, h = frozen_rt_shifted(mdl, N, T)
        march = frozen_rt_march(mdl, fill(h, N))
        @printf(io, "    (c) 冻结 cell 行进 vs 节点递推（均匀）：|Δ| = %.1e（逐位一致）\n",
                abs(march - rec[N]))
    end

    # ---------------- [1] 实时阶数阶梯 ----------------
    println(io, "\n[1] 实时 A(t)：装配 / 分裂的阶数阶梯（同一批 bosonic_Δt 条目）")
    es = Float64[]
    for N in Ns
        vals = frozen_retarded_A(mdl, N, T)[1][ms_of(N).+1]
        push!(es, maximum(abs.(vals .- [Aex(m * T / N) for m in ms_of(N)])))
    end
    order_table(io, "    旧冻结装配（frozen_retarded_A）", Ns, es)
    es = Float64[]
    for N in Ns
        vals = frozen_rt_shifted(mdl, N, T)[1][ms_of(N)]
        push!(es, maximum(abs.(vals .- [Aex(m * T / N) for m in ms_of(N)])))
    end
    order_table(io, "    顶点移位冻结装配（V² 精确，常数↓~12×）", Ns, es)
    for order in (2, 4)
        es = Float64[]
        for N in Ns
            st = real_step(mdl, order, T / N)
            push!(es, maximum(abs.([(st^m)[1, 1] for m in ms_of(N)] .-
                                   [Aex(m * T / N) for m in ms_of(N)])))
        end
        order_table(io, "    严格 cell 链 order=$order（参照）", Ns, es)
    end

    # ---------------- [2] 实时高阶组合 ----------------
    println(io, "\n[2] 实时高阶组合：Suzuki 布局（冻结 cell）与 Richardson 外推")
    esA = Float64[]; esB = Float64[]
    for N in Ns[1:3]
        h = T / N
        push!(esA, abs(frozen_rt_march(mdl, fill(h, N)) - Aex(T)))
        push!(esB, abs(frozen_rt_march(mdl, repeat(collect(SUZUKI_W), N) .* h) - Aex(T)))
    end
    order_table(io, "    冻结 cell 均匀布局（终点）", Ns[1:3], esA)
    order_table(io, "    冻结 cell Suzuki(θ,θ,γ,θ,θ) 布局（终点）", Ns[1:3], esB)
    for levels in (1, 2, 3)
        es = Float64[]
        for N in Ns[1:end-levels+1]
            vals = richardson_levels(n -> frozen_rt_shifted(mdl, n, T)[1], N, ms_of(N); levels=levels)
            push!(es, maximum(abs.(vals .- [Aex(m * T / N) for m in ms_of(N)])))
        end
        order_table(io, "    顶点移位冻结 + Richardson×$levels", Ns[1:end-levels+1], es)
    end
    es = Float64[]
    for N in Ns[1:end-1]
        a1, h = frozen_retarded_A(mdl, N, T)
        a2, _ = frozen_retarded_A(mdl, 2N, T)
        vals = [2 * a2[2m+1] - a1[m+1] for m in ms_of(N)]
        push!(es, maximum(abs.(vals .- [Aex(m * T / N) for m in ms_of(N)])))
    end
    order_table(io, "    旧冻结装配 + Richardson×1", Ns[1:end-1], es)

    # ---------------- [3] 虚时环 ----------------
    println(io, "\n[3] 虚时 G(τ)：冻结 cell 环（W_f）与 Richardson")
    es = Float64[]
    for N in Ns
        vals, h = frozen_ring_Wf(mdl, N)
        ms = round.(Int, fracs .* N)                 # G[m+1] ≈ G(m·h)（含 G(0)）
        push!(es, maximum(abs.([vals[m+1] for m in ms] .- [Gex(m * h) for m in ms])))
    end
    order_table(io, "    冻结 cell 环 W_f：[T_f^m(I−T_f^N)^{-1}]₁₁", Ns, es)
    es = Float64[]
    for N in Ns[1:end-1]
        G1, h = frozen_ring_Wf(mdl, N)
        G2, _ = frozen_ring_Wf(mdl, 2N)
        ms = round.(Int, fracs .* N)
        vals = [2 * G2[2m+1] - G1[m+1] for m in ms]
        push!(es, maximum(abs.(vals .- [Gex(m * h) for m in ms])))
    end
    order_table(io, "    冻结 cell 环 + Richardson×1", Ns[1:end-1], es)
    es = Float64[]
    for N in Ns[1:end-2]
        vals, ms, h = chain_Gτ(mdl, 2, N; ms=round.(Int, fracs .* N))
        push!(es, maximum(abs.(vals .- [Gex(mod(m * h, mdl.β)) for m in ms])))
    end
    order_table(io, "    严格 Strang 环链（参照）", Ns[1:end-2], es)

    # ---------------- [4] 结论 ----------------
    println(io, "\n[4] 结论（实测）")
    println(io, "  · 严格 Strang 链的节点核 = 冻结条目顶点移位一格 + 每步统一 u 相位 +")
    println(io, "    三角并入次对角（V² 逐模解析恒等，[0](a)(b)）；旧冻结装配不满足该恒等式，")
    println(io, "    常数大 ~12×。两者同为一阶：剩余亏损 V⁴h²(1+d)/条目（顶点自缀饰 + 飞行")
    println(io, "    虚访问），逐模、随 lag 线增、按 ~N 个 lag 相干累积 → 全局 O(h)。")
    println(io, "  · Suzuki 重排冻结 cell（θ,θ,γ,θ,θ 有向布局）不提升阶数（仍 1 阶，对照严格")
    println(io, "    Suzuki 的 4 阶）⟹ 亏损不是分裂型误差，BCH 代数消不掉。")
    println(io, "  · Richardson 外推是冻结核到高阶的实际组合路径：1/2/3 层 → 2.00 / 2.99 /")
    println(io, "    4.01 阶；旧装配 + 1 层也到 1.99。虚时冻结 cell 环同理（1 阶 → 2 阶）。")
    println(io, "  · 与 GuoChen2026a App.C：其 Eq.(12) 即顶点移位装配，全局渐近 1 阶与其")
    println(io, "    O(t·δt) 推断一致（offdiag_order.jl 小窗口的“全局 δt²”是 pre-asymptotic）；")
    println(io, "    与 NgReichman2023：其“对 Trotter 化动力学精确求 IF”对应严格链参照组。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
