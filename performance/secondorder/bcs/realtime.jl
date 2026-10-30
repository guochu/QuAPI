# =============================================================================
# secondorder/bcs/realtime.jl —— BCS 浴实时 K 侧二阶（Weyl-K，Nambu 矩阵化）
# 输出画图数据 results/realtime.csv
#
# 模型：H = Σ_σ ϵ_d d_σ†d_σ + Σ_kσ ε_k c_kσ†c_kσ + Σ_kσ V_k(d_σ†c_kσ + h.c.)
#         + Δ Σ_k (c_k↑†c_k↓† + c_k↓c_k↑)（μ=0；Δ 实数或复数均支持）。
#
# 与 secondorder/fermions/realtime.jl 同骨架，核从标量升级为 2×2 Nambu 矩阵：
#   Nambu 基 Φ = (d↑, c₁↑..c_L↑ ; d̄↓†, c̄₁↓†..c̄_L↓†)，杂质 h_𝓝 = diag(ϵ_d,−ϵ_d)，
#   推迟 Dyson 方程 (i∂_t − h_𝓝)𝒜 = Σ^R⊛𝒜 + δI 的隐式梯形离散：
#   粒子行相位 (u, f) = (e^{−iϵ_d h}, e^{−iϵ_d h/2})，空穴行 (e^{+iϵ_d h}, e^{+iϵ_d h/2})；
#   接触 = 2×2 等时矩阵 D₀（显式 −f D₀A[m]/2 + 隐式 (1+fD₀/2)⁻¹，按列解）。
# 推迟组合（Δt 约定；数值扫描 + Δ=0 退化 + 实/复 Δ 一致性三重验证）：
#   正常槽位（cn/nc 组合块）：X − P（同正常费米子；u²,v² 实 ⟹ 相位不敏感）；
#   反常槽位（cc/nn 组合块）：P − X = −_uvφ₊ + _uvφ₋（±E 双极点未共轭权重、
#   热因子 1）。注意：其正确性依赖 fermionic_Δt 反向扇区的负滞后路线（不共轭
#   谱权重）——旧版共轭存储在复 Δ 下会让组合块转置部分携带 conj(_uv)（实 Δ 无
#   感）；该核心问题已修复（见 src/realtime/fermionic.jl 注释）。
# 参照：𝒜(t) = 杂质 2×2 Nambu 块 of e^{−i𝓗_BdG t}（= −iθ·G^R_𝓝，与高斯态无关；
#   𝓗 已由 Fock 对易子逐行验证）。Δ=0 退化为两个解耦正常通道（已验证）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/bcs/realtime.jl
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

"""
    weyl_k_bcs_realtime(mdl, N, T) -> (𝒜[1..N], h)

BCS 实时 K 侧二阶（Nambu 矩阵隐式梯形）。`𝒜[m] ≈ e^{−i𝓗t_m} 的杂质 2×2 块`。
"""
function weyl_k_bcs_realtime(mdl, N, T)
    h = T / N
    u = exp(-im * mdl.ϵd * h)
    f = exp(-im * mdl.ϵd * h / 2)
    uI = exp(im * mdl.ϵd * h)
    fI = exp(im * mdl.ϵd * h / 2)
    corr = bcs_Δt(mdl.bath.spectrum; β=mdl.β, N=N, t=T, μ=mdl.μ, Δ=mdl.Δ)
    P1, X1 = branch(corr.uu, :+, :+), branch(corr.uu, :+, :-)
    P2, X2 = branch(corr.dd, :+, :+), branch(corr.dd, :+, :-)
    P3, X3 = branch(corr.ud, :+, :+), branch(corr.ud, :+, :-)
    P4, X4 = branch(corr.du, :+, :+), branch(corr.du, :+, :-)
    kcn(j, k) = -(X1[j, k] + P1[j, k])   # 正常 ↑↑（统一约定；≡ 旧 X−P）
    knc(j, k) = -(X2[j, k] + P2[j, k])   # 正常 ↓↓（统一约定；≡ 旧 X−P）
    kcc(j, k) = P3[j, k] + X3[j, k]      # 反常 𝒩₁₂（统一约定；≡ 旧 P−X）
    kdu(j, k) = P4[j, k] + X4[j, k]      # 反常 𝒩₂₁（统一约定；≡ 旧 P−X）
    𝒜 = [zeros(ComplexF64, 2, 2) for _ in 1:N+1]
    𝒜[1] = Matrix{ComplexF64}(I, 2, 2)
    D0 = [kcn(1, 1) kcc(1, 1); kdu(1, 1) knc(1, 1)]
    rowphase = [(u, f), (uI, fI)]
    for m in 1:N
        hist = zeros(ComplexF64, 2, 2)
        for k in 1:m-1
            K = [kcn(m+1, k+1) kcc(m+1, k+1); kdu(m+1, k+1) knc(m+1, k+1)]
            hist += K * (𝒜[k] + 𝒜[k+1]) / 2
        end
        RHS = zeros(ComplexF64, 2, 2)
        for p in 1:2
            up, fp = rowphase[p]
            RHS[p, :] .= up .* 𝒜[m][p, :] .- fp .* (hist[p, :] .+ (D0 * 𝒜[m] / 2)[p, :])
        end
        Msys = Matrix{ComplexF64}(I, 2, 2)
        for p in 1:2, r in 1:2
            Msys[p, r] += rowphase[p][2] * D0[p, r] / 2
        end
        for q in 1:2
            𝒜[m+1][:, q] = Msys \ RHS[:, q]
        end
    end
    return 𝒜[2:end], h
end

"""
    strict_chain_bcs(mdl, N, T) -> (𝒜[1..N], h)

严格 cell 链参照：BdG 单粒子问题 𝓗 = hA + hB（hA = 杂质 Nambu 块，hB = 其余）
的 Strang 分裂 eA·eB·eA；返回杂质 2×2 块（Lie ≡ Strang 节点值，费米扇区简并）。
"""
function strict_chain_bcs(mdl, N, T)
    h = T / N
    𝓗 = bcs_bdg(mdl)
    n = mdl.L + 1
    hA = zeros(ComplexF64, 2n, 2n)
    hA[1, 1] = mdl.ϵd
    hA[n+1, n+1] = -mdl.ϵd
    hB = 𝓗 - hA
    eA = exp(-im * h / 2 * hA)
    eB = exp(-im * h * hB)
    st = eA * eB * eA
    idx = [1, n+1]
    out = [(st^m)[idx, idx] for m in 1:N]
    return out, h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Δ=0.5, T=2.0, Ns=[16, 32, 64, 128, 256, 512],
              fracs=[0.25, 0.5, 0.75])
    mdl = bcs_toulouse(; Lmodes=Lmodes, Δ=Δ)
    ex = exact_retarded_Nambu(mdl)          # −i𝒢^R_𝓝 = e^{−i𝓗t} 杂质块
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "BCS 浴实时 Weyl-K 二阶（Nambu 矩阵核 + K 侧）vs 严格 cell 链：𝒜(t)")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，Δ=%.2f，T=%.1f\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64}[]
    err_of(𝒜s, N) = maximum(maximum(abs.([𝒜s[m][p, q] - (ex(m * T / N) ./ (-im))[p, q]
                                          for m in ms_of(N), p in 1:2, q in 1:2])))
    println(io, "\n[1] Weyl-K 二阶（Nambu 矩阵核）vs 精确 e^{−i𝓗t} 杂质块")
    es = Float64[]
    for N in Ns
        𝒜s, h = weyl_k_bcs_realtime(mdl, N, T)
        e = err_of(𝒜s, N)
        push!(es, e)
        push!(rows, (N, h, e, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照（BdG Strang 分裂）")
    es = Float64[]
    for (i, N) in enumerate(Ns)
        𝒜s, h = strict_chain_bcs(mdl, N, T)
        e = err_of(𝒜s, N)
        push!(es, e)
        rows[i] = (rows[i][1], rows[i][2], rows[i][3], e)
    end
    order_table(io, "    strict order=2", Ns, es)
    println(io, "\n[3] 复数 Δ 验证（Δ = 0.2+0.45im；报告，不入 CSV）")
    mdlc = bcs_toulouse(; Lmodes=Lmodes, Δ=0.2 + 0.45im)
    exc = exact_retarded_Nambu(mdlc)
    es = Float64[]
    for N in Ns
        𝒜s, h = weyl_k_bcs_realtime(mdlc, N, T)
        e = maximum(maximum(abs.([𝒜s[m][p, q] - (exc(m * T / N) ./ (-im))[p, q]
                                  for m in ms_of(N), p in 1:2, q in 1:2])))
        push!(es, e)
    end
    order_table(io, "    Weyl-K 二阶（复 Δ）", Ns, es)
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "realtime.csv"), "w") do fh
        println(fh, "N,dt,err_weylk,err_strict2")
        for (N, h, e1, e2) in rows
            @printf(fh, "%d,%.6e,%.6e,%.6e\n", N, h, e1, e2)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "realtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
