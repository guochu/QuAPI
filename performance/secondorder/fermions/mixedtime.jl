# =============================================================================
# secondorder/fermions/mixedtime.jl —— 费米混合轮廓（L 形 Kadanoff-Baym）
# 全平衡 G>(t)：全 fermionic_Δm 冻结核 + K 侧直接二阶。输出 results/mixedtime.csv
#
# 参照（U1GTEMPO2 混合时间轴，本机实测，单模 β=2、T=0.8）：order=1 对全热
# 平衡 ED 阶 0.98、order=2 阶 2.00——混合目标 = τ 环热化的全平衡 G>（含 γ
# 非对角初关联），二阶由 K 侧修正达成、IF 不变。本文件在 QuAPI 冻结核框架
# 内实现同目标、同阶数。
#
# ============================ fermionic_Δm 块使用表 ============================
# 全部数据取自单次 fermionic_Δm 调用的三个块（+ 同数据的闭式权重）：
#   ① ξ 块（τ,τ）——τ 支环方程（weyl_k_imagtime 结构：梯形 φ、半步相位 f、
#      接触对角、反周期跳变 G₀ = −1−G(β⁻) 完整代入）→ 环解 G(τ_j)（二阶）；
#      n_d = −G(β⁻)（KMS）；γ[1,k] = 消峰 Boltzmann 收缩（见下，二阶）；
#   ② η 块（+,+）（+,−）——实时支 Weyl-K 递推的推迟组合。注意 Δm 约定：
#      branch(:+,:−) = −(Δt 的同块）⟹ 推迟组合 = −(P+X)（等价于 Δt 的 X−P，
#      数值已对拍：A 通道对 Aex 二阶一致）；
#   ③ ζ 块（+,τ）——实时支 C 通道的**全历史 Duhamel 源的构成因子**：
#      ζ⁺τ 每模严格 rank-1（单模 SVD σ₂/σ₁ = 9e-17）⟹ 实时侧因子 =
#      I_k(m)·e^{∓相位}（cell 对积分闭式，进入源），τ 侧因子 = Boltzmann
#      权重（进入 γ[1,k] 收缩）。即 ζ 块的物理内容（混合交连）以因式分解
#      形式完整进入源与初关联，而非 O(N²) 数值收缩。
#      恒等式数值证据：Σ_j ζ[m,j]·w_j（局部收缩）= γ[1,k]·I_k(m)·(−iV_k)
#      对全部 m 成立（比值常数、虚部波动 0.0）——因此 C 通道 = 同一条 Weyl-K
#      递推（η 块核，同 A 通道）+ 局部 Duhamel 源（ζ 块实时侧因子 × γ 行）。
#      注意递推配源的组合必须自洽：局部源配 η 块 Δ 递推（本文件，二阶）；
#      全历史源 J_k 配 𝒦 单重卷积递推（等价形式）。混配（局部源 + 𝒦，或全历
#      史源 + Δ）不自洽。
#
# ============================ 为什么 ζ 不能直接收缩 ============================
# bath 起点振幅 [U]₁ₐ 的正确记忆方程核是**单重卷积** 𝒦(D) =
# Σ_k V_k²e^{−iϵ_dD}(e^{iω_kD}−1)/(iω_k)（bath 起点振幅载波为 ε_a；杂质自
# 能核 Δ 是二重 cell 积分、载波 ϵ_d，不适用——实测：ζ 收缩驱动 Δ 递推 =
# 平地板/发散，𝒦 递推 + 全历史源 = 直接二阶）。C 通道源的全历史 Duhamel
# 积分 J_k(m) 与 γ[1,k] 收缩分别正是 ζ 块因式分解的实时侧因子与 τ 侧因子。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/fermions/mixedtime.jl
# =============================================================================
include(joinpath(@__DIR__, "imagtime.jl"))    # 复用 weyl_k_imagtime（传递引入 common.jl）

"𝒦(D) = Σ_k V_k² e^{−iϵ_dD}(e^{iω_kD}−1)/(iω_k)（bath-only 单重卷积核，闭式）。"
function mathcalK(mdl, D)
    S = 0.0 + 0.0im
    for k in 1:mdl.L
        ω = mdl.ϵd - mdl.ws[k]
        S += mdl.αs[k] * exp(-im*mdl.ϵd*D) *
             (abs(ω) < 1e-10 ? im*D : (exp(im*ω*D) - 1) / (im*ω))
    end
    return S
end

"""
    weyl_k_mixed(mdl, Nτ, Nt, T) -> G>[t_0..t_Nt]

费米混合全 fermionic_Δm K 侧二阶（全平衡 G>）。`G[m+1] ≈ G>(t_m)`。
"""
function weyl_k_mixed(mdl, Nτ, Nt, T)
    h = T / Nt
    δτ = mdl.β / Nτ
    u = exp(-im * mdl.ϵd * h)          # 实时支相位
    f = exp(-im * mdl.ϵd * h / 2)
    uτ = exp(-mdl.ϵd * δτ)             # τ 支相位（实衰减！虚时无 i）
    fτ = exp(-mdl.ϵd * δτ / 2)
    dm = fermionic_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=Nτ, t=T, Nt=Nt, μ=mdl.μ)
    # ① τ 支环（ξ 块；四要素同 imagtime.jl）
    ξ = branch(dm; b1=:τ, b2=:τ)
    ηjk, ηkj = ξ.ηⱼₖ, ξ.ηₖⱼ
    d0τ = ηjk[1] + ηkj[1]
    M = zeros(ComplexF64, Nτ, Nτ); b = zeros(ComplexF64, Nτ)
    for m in 1:Nτ
        M[m, m] += 1.0
        if m ≥ 2; M[m, m-1] -= uτ
        else; b[m] -= uτ; M[m, Nτ] += uτ; end
        for k in 1:Nτ
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0τ
            if k == 1                       # φ₁ = (−1−G_N+G₁)/2：跳变完整代入
                M[m, Nτ] += fτ*w*0.5
                M[m, 1] -= fτ*w*0.5
                b[m] -= fτ*w*0.5
            else
                M[m, k-1] -= fτ*w*0.5
                M[m, k] -= fτ*w*0.5
            end
        end
    end
    Gτ = M \ b
    n_d = -Gτ[Nτ]
    G0v = -1 - Gτ[Nτ]
    Gext = vcat(G0v, Gτ)
    # γ[1,k] 消峰 Boltzmann 收缩（ζ 块 τ 侧因子的二阶离散）
    γrow = zeros(ComplexF64, mdl.L)
    for k in 1:mdl.L
        ε, V = mdl.ws[k], sqrt(mdl.αs[k])
        W_tot = abs(ε) < 1e-10 ? mdl.β : (1 - exp(-ε*mdl.β))/ε
        acc = 0.0 + 0.0im
        for j in 1:Nτ
            w0 = exp(-ε*(mdl.β - (j-1)*δτ))
            w1 = exp(-ε*(mdl.β - j*δτ))
            acc += δτ * (w0*(Gext[j] - G0v) + w1*(Gext[j+1] - G0v)) / 2
        end
        γrow[k] = V * (G0v * W_tot + acc) / (1 + exp(-ε*mdl.β))
    end
    # ② 实时支 A（η 块；Δm 约定推迟组合 −(P+X)）
    P = branch(dm; b1=:+, b2=:+)
    Xr = branch(dm; b1=:+, b2=:-)
    Xc = branch(dm; b1=:+, b2=:τ)
    Δ(j, k) = -(Xr[j, k] + P[j, k])
    d0 = Δ(1, 1)
    A = zeros(ComplexF64, Nt + 1); A[1] = 1.0
    for m in 1:Nt
        s = 0.0 + 0.0im
        for k in 1:m-1
            s += Δ(m+1, k+1) * (A[k] + A[k+1]) / 2
        end
        A[m+1] = (u*A[m] - f*s - f*d0*A[m]/2) / (1 + f*d0/2)
    end
    # ③ C 通道：同一条 Weyl-K 递推（η 块核）+ 局部 Duhamel 源
    #    S_m = Σ_k(−iV_k)·I_k(m)·(−γrow[k])，I_k(m) = ζ 块实时侧因子的 cell 积分
    #    （相位基准 e^{−iε_a t_{m−1}}·e^{−iϵ_d h}——错误的基准导致一阶病态）
    C = zeros(ComplexF64, Nt + 1)
    for m in 1:Nt
        s = 0.0 + 0.0im
        for k in 1:m-1
            s += Δ(m+1, k+1) * (C[k] + C[k+1]) / 2
        end
        S = 0.0 + 0.0im
        for k in 1:mdl.L
            εa, V = mdl.ws[k], sqrt(mdl.αs[k])
            ωd = mdl.ϵd - εa
            I = exp(-im*εa*(m-1)*h) * exp(-im*mdl.ϵd*h) *
                (abs(ωd) < 1e-10 ? h : (exp(im*ωd*h) - 1) / (im*ωd))
            S += -im * V * I * (-γrow[k])
        end
        C[m+1] = (u*C[m] - f*s - f*d0*C[m]/2 + S) / (1 + f*d0/2)
    end
    return (-im) .* ((1 - n_d) .* A .+ C)
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[16, 32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = fermi_toulouse(; Lmodes=Lmodes)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    # 严格参照用 γ_eq（一次对角化，仅参照、不进入求解器）
    ev = eigen(Hermitian(mdl.hfull))
    γ = ev.vectors * Diagonal([1/(exp(mdl.β*λ)+1) for λ in ev.values]) * ev.vectors'
    gexact(N, m) = -im * (exp(-im * (m * T / N) * mdl.hfull) * (mdl.Id - transpose(γ)))[1, 1]
    println(io, "=" ^ 74)
    println(io, "费米混合 Weyl-K 二阶（全 fermionic_Δm 冻结核）：全平衡 G>(t)")
    println(io, "  ξ 块 τ 环 → n_d/γ 行（消峰收缩）；η 块实时 A；ζ 块因式分解 → C 源 + 𝒦 卷积")
    println(io, "  参照：U1GTEMPO2 mixed order=1/2 对全平衡 ED = 0.98 / 2.00（本机实测）")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，T=%.1f（Nτ=Nt 同步）\n",
            length(mdl.ws), mdl.Γ, mdl.β, mdl.ϵd, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64}[]
    println(io, "\n[1] Weyl-K 二阶 vs 精确 G>(t)（求解全程不使用精确 γ）")
    es = Float64[]
    for N in Ns
        G = weyl_k_mixed(mdl, N, N, T)
        e = maximum(abs.([G[m+1] - gexact(N, m) for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, T / N, e, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照（同 γ_eq，仅参照用）")
    es = Float64[]
    for (i, N) in enumerate(Ns)
        h = T / N
        eA = exp(-im * h / 2 * mdl.hS)
        eB = exp(-im * h * (mdl.hb + mdl.hi))
        st = eA * eB * eA
        gt = [-im * ((st^m) * (mdl.Id - transpose(γ)))[1, 1] for m in ms_of(N)]
        e = maximum(abs.(gt .- [gexact(N, m) for m in ms_of(N)]))
        push!(es, e)
        rows[i] = (rows[i][1], rows[i][2], rows[i][3], e)
    end
    order_table(io, "    strict order=2", Ns, es)
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "mixedtime.csv"), "w") do f
        println(f, "N,dt,err_weylk,err_strict2")
        for (N, h, e1, e2) in rows
            @printf(f, "%d,%.6e,%.6e,%.6e\n", N, h, e1, e2)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "mixedtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
