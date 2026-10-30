# =============================================================================
# secondorder/bcs/imagtime.jl —— BCS 浴虚时 K 侧二阶（Weyl-K，Nambu 4 分量环）
# 输出画图数据 results/imagtime.csv
#
# 与 secondorder/fermions/imagtime.jl 同骨架，环未知量从标量升级为 Nambu 2×2：
#   𝒢(τ) = −⟨Φ(τ)Φ†(0)⟩ 的杂质块，Φ = (d↑, c↑… ; d̄↓†, c̄↓†…)（同 common.jl 基）。
#   四分量 𝒢[p,q]（p = τ 端算符行，q = 0 端列）由矩阵 Dyson 环方程闭合：
#   (𝒢₀⁻¹⊛𝒢)[p,j] = δ + Σ_s K[p,s]⊛𝒢[s,j]，核 = bcs_Δτ 四块（cn,nc,cc,nn）的
#   环用法（ηⱼₖ 前向/ηₖⱼ 反向/等时接触），梯形 φ、半步相位 f、跳变完整代入。
#   行相位：粒子行 (e^{−ϵ_d δτ}, e^{−ϵ_d δτ/2})，空穴行 (e^{+ϵ_d δτ}, +ϵ_d/2)
#   （𝒢₂₂ ∝ e^{+ϵτ} 增长，V=0 已验证）。跳变常数：𝒢₁₁/𝒢₂₂ 行 G₀ = −1−G(β⁻)
#   （反对易 δ 源），反常分量 G₀ = −G(β⁻)（不同模无 δ，连续跨界）。
#   提取符号：方案的反常列 𝒢[1][2]/𝒢[2][1] 相对参照差一个整体负号
#   （列腿序 Grassmann 约定），输出取 sgn(p,q) = (p==q ? +1 : −1)。
# 参照：𝒢(τ_m) = −[S(τ_m)M₀] 杂质块（谱分解，已与全 Fock ED/包机制对拍）；
#   严格参照 = 精确离散热环 ρ(τ_m) = T^m(I−T^N)⁻¹M₀，𝒢 = −ρ 杂质块。
#
# 复 Δ：环方程天然支持（双侧边值结构中共轭自洽；数值验证复 Δ = 0.2+0.45im
# 同样 2.00 阶，见 main() [4] 段）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/bcs/imagtime.jl
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

"""
    weyl_k_bcs_imagtime(mdl, N) -> (𝒢[1..N], h)

BCS 虚时 K 侧二阶（Nambu 4 分量矩阵环）。`𝒢[m] ≈ −⟨Φ(τ_m)Φ†(0)⟩ 杂质块`。
"""
function weyl_k_bcs_imagtime(mdl, N)
    h = mdl.β / N
    u1 = exp(-mdl.ϵd * h);  f1 = exp(-mdl.ϵd * h / 2)
    u2 = exp(+mdl.ϵd * h);  f2 = exp(+mdl.ϵd * h / 2)
    corr = bcs_Δτ(mdl.bath.spectrum; β=mdl.β, N=N, μ=mdl.μ, Δ=mdl.Δ)
    J = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    K = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    J[1,1], K[1,1] = corr.uu.data.ηⱼₖ, corr.uu.data.ηₖⱼ      # cn（正常 ↑↑）
    J[2,2], K[2,2] = corr.dd.data.ηⱼₖ, corr.dd.data.ηₖⱼ      # nc（正常 ↓↓）
    J[1,2], K[1,2] = corr.ud.data.ηⱼₖ, corr.ud.data.ηₖⱼ      # cc（反常 ā↑ā↓）
    J[2,1], K[2,1] = corr.du.data.ηⱼₖ, corr.du.data.ηₖⱼ      # nn（反常 a↓a↑）
    wpq(p, q, m, k) = k < m ? J[p, q][m-k+1] :
                      k > m ? K[p, q][k-m+1] : J[p, q][1] + K[p, q][1]
    M = zeros(ComplexF64, 4N, 4N)
    b = zeros(ComplexF64, 4N)
    δrow(p, j) = (p == j) ? -1.0 : 0.0        # 跳变常数（反对易 δ 源仅在对角行）
    for p in 1:2, j in 1:2, m in 1:N
        r = ((p-1)*2 + (j-1))*N + m
        up, fp = p == 1 ? (u1, f1) : (u2, f2)
        M[r, r] += 1.0
        if m ≥ 2
            M[r, ((p-1)*2+(j-1))*N + m-1] -= up
        else
            cp = δrow(p, j)
            cp != 0 && (b[r] += up * cp)
            M[r, ((p-1)*2+(j-1))*N + N] += up
        end
        for k in 1:N, s in 1:2
            wv = wpq(p, s, m, k)
            if k == 1
                cpk = δrow(s, j)
                M[r, ((s-1)*2+(j-1))*N + N] += fp * wv * 0.5
                M[r, ((s-1)*2+(j-1))*N + 1] -= fp * wv * 0.5
                cpk != 0 && (b[r] += fp * wv * 0.5 * cpk)
            else
                M[r, ((s-1)*2+(j-1))*N + k-1] -= fp * wv * 0.5
                M[r, ((s-1)*2+(j-1))*N + k]   -= fp * wv * 0.5
            end
        end
    end
    x = M \ b
    # 组装 2×2 Nambu 矩阵（反常分量带提取符号）
    return [ComplexF64[x[((p-1)*2+(q-1))*N + m] * (p == q ? 1.0 : -1.0)
                       for p in 1:2, q in 1:2] for m in 1:N], h
end

"""
    exact_Gτ_Nambu(mdl, N) -> (𝒢[1..N], h)

连续谱严格参照：𝒢(τ_m) = −[S(τ_m)M₀] 的杂质 2×2 Nambu 块。
"""
function exact_Gτ_Nambu(mdl, N)
    V, λ, Vinv = bcs_eig(mdl)
    n = mdl.L + 1
    h = mdl.β / N
    M0 = bcs_M0(mdl)
    out = Vector{Matrix{ComplexF64}}(undef, N)
    for m in 1:N
        A = V * (Diagonal(exp.(-λ .* m .* h)) * (Vinv * M0))
        out[m] = -A[[1, n+1], [1, n+1]]
    end
    return out, h
end

"""
    strict_ring_bcs(mdl, N) -> (𝒢[1..N], h)

严格 cell 链参照：Strang 分裂 cell st = eA·eB·eA（同 realtime.jl）的离散演化
M(τ_m) = st^m·M₀，𝒢 = −M(τ_m) 杂质块。
"""
function strict_ring_bcs(mdl, N)
    h = mdl.β / N
    𝓗 = bcs_bdg(mdl)
    n = mdl.L + 1
    hA = zeros(ComplexF64, 2n, 2n)
    hA[1, 1] = mdl.ϵd
    hA[n+1, n+1] = -mdl.ϵd
    st = exp(-h / 2 * hA) * exp(-h * (𝓗 - hA)) * exp(-h / 2 * hA)
    M0 = bcs_M0(mdl)
    idx = [1, n+1]
    out = [-(st^m * M0)[idx, idx] for m in 1:N]
    return out, h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Δ=0.5, Ns=[16, 32, 64, 128, 256], fracs=[0.25, 0.5, 0.75])
    mdl = bcs_toulouse(; Lmodes=Lmodes, Δ=Δ)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "BCS 浴虚时 Weyl-K 二阶（Nambu 4 分量矩阵环 + K 侧）：G(τ)")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，Δ=%.2f\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64}[]
    println(io, "\n[1] Weyl-K 二阶 vs 连续谱严格参照（分量 𝒢₁₁ = −⟨d↑(τ)d↑†⟩）")
    es = Float64[]
    for N in Ns
        G, h = weyl_k_bcs_imagtime(mdl, N)
        ex, _ = exact_Gτ_Nambu(mdl, N)
        e = maximum(abs.([G[m][1, 1] - ex[m][1, 1] for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, h, e, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照（Strang cell 热环 M = st^m·M₀）")
    es = Float64[]
    for (i, N) in enumerate(Ns)
        Gs, h = strict_ring_bcs(mdl, N)
        ex, _ = exact_Gτ_Nambu(mdl, N)
        e = maximum(abs.([Gs[m][1, 1] - ex[m][1, 1] for m in ms_of(N)]))
        push!(es, e)
        rows[i] = (rows[i][1], rows[i][2], rows[i][3], e)
    end
    order_table(io, "    strict（精确离散环）", Ns, es)
    # 反常分量（额外报告，不入 CSV）
    println(io, "\n[3] 反常分量 𝒢₁₂ = −⟨d↑(τ)d↓⟩ 收敛（报告）")
    for N in Ns
        G, h = weyl_k_bcs_imagtime(mdl, N)
        ex, _ = exact_Gτ_Nambu(mdl, N)
        e = maximum(abs.([G[m][1, 2] - ex[m][1, 2] for m in ms_of(N)]))
        @printf(io, "    N=%4d   err_F=%.4e\n", N, e)
    end
    println(io, "\n[4] 复数 Δ 验证（Δ = 0.2+0.45im；报告，不入 CSV）")
    mdlc = bcs_toulouse(; Lmodes=Lmodes, Δ=0.2 + 0.45im)
    es = Float64[]
    for N in Ns
        G, h = weyl_k_bcs_imagtime(mdlc, N)
        ex, _ = exact_Gτ_Nambu(mdlc, N)
        e = maximum(abs.([G[m][1, 1] - ex[m][1, 1] for m in ms_of(N)]))
        push!(es, e)
    end
    order_table(io, "    Weyl-K 二阶（复 Δ）", Ns, es)
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "imagtime.csv"), "w") do fh
        println(fh, "N,dtau,err_weylk,err_strict2")
        for (N, h, e1, e2) in rows
            @printf(fh, "%d,%.6e,%.6e,%.6e\n", N, h, e1, e2)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "imagtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
