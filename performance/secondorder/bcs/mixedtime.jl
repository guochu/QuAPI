# =============================================================================
# secondorder/bcs/mixedtime.jl —— BCS 浴混合轮廓（L 形 Kadanoff-Baym）
# 全平衡 G₁>(t) = ⟨d↑(t)d↑†(0)⟩：单次 bcs_Δm 冻结核 + K 侧直接二阶。
# 输出 results/mixedtime.csv
#
# 与 secondorder/fermions/mixedtime.jl 同目标、同结构（τ 环热化的全平衡 G>，
# 含初关联；数据全部取自 Δm 轮廓对象）：
#   ① ξ 块（τ,τ）——τ 支 4 分量 Nambu 环（imagtime.jl 同构方程；branch(τ,τ)
#      场 ≡ bcs_Δτ 环数据，实/复 Δ 均逐元验证相等）→ 环解 𝒢[m]；
#      初值数据：⟨d↑d↑†⟩ = 1+𝒢₁₁[N]、反常 ⟨d̄↓†d↑†⟩ = 𝒢₂₁[N]（跳变反解）。
#   ② γ 初关联收缩（族 (1,1)，正常侧消峰公式的 Nambu 推广，2.00 对拍 M₀）：
#      未知 g=(⟨c_k↑(τ)d↑†⟩,⟨c_k↓†(τ)d↑†⟩) 满足
#        ∂g = [[−ε_k,−Δ],[−Δ*,+ε_k]]·g + V_k·(𝒢₁₁(τ), −𝒢₂₁(τ))
#      （源符号 (+,−)、无 j 因子；[AB,C] 恒等式 + BdG 行方程双推 + 真值函数
#      有限差分 1e-12 验证）。严格反周期 g(β⁻)=−g(0⁺)（不同模 {c,d†}=0 ⟹ 无
#      跳变）、τ=0⁺ 源值用环跳变代入、消峰分拆
#      X₀ = −(I+e^{Mβ})⁻¹[M⁻¹(e^{Mβ}−I)s₀ + CN-Σ]（Δ=0 退化为 fermions 侧
#      γrow = V[G₀W_tot+∫e^{−ε(β−s)}(G−G₀)]/(1+e^{−εβ})）。
#   ③ η 块（+,+）（+,−）——实时支推迟核（Δm 约定，全部取自组合块 branch 矩阵）：
#      正常槽位 −(P+X)（≡ Δt 约定的 X−P；fermion 侧同款翻译）；
#      反常槽位 P+X（≡ Δt 约定的 P−X，即 −_uvφ₊ + _uvφ₋，复 Δ 正确——依赖
#      MixedCorrelationFunction 的 ηⱼₖ₂/ηₖⱼ₂ 补全存储与 branch 的无共轭派生，
#      核心修复见 src/mixedtime/mixedcorr.jl）。
#   ④ C 通道：同一矩阵 Weyl-K 递推（零初值）+ 逐 cell Nambu Duhamel 源
#      Src_k(m) = −i·[expm(h·[[−iD(ϵ_d), 𝒱_k·e^{−i𝓗_bk t_{m−1}}],[0, −i𝓗_bk]])]₁₂·γ_k，
#      𝒱_k = diag(V_k,−V_k)（空穴耦合 −V，BdG −hᵀ 块）、𝓗_bk = [[ε_k,Δ],[Δ*,−ε_k]]
#      （厄米；复 Δ 下 [2,1] 取 conj(Δ)）；4×4 block-exp 是 cell 内双相位飞行
#      积分的精确闭式（正常侧 I_k(m) 的推广，即 ζ 交连块的因式分解形式——与
#      fermions 侧 ζ rank-1 因式分解同源，故求解器不直接收缩 ζ 块）。
#   ⑤ 组装：G₁>(t_m) = M₀₁₁·𝒜₁₁ + F₀·𝒜₁₂ + C₁（retarded 部分 O(h²)、C 通道
#      携带浴起点初关联）。
#
# 参照：exact_greater_eq = [S(t)M₀]₁₁（已与全 Fock ED 对拍）；严格 cell 链 =
# [(st)^m·M₀]₁₁（st = 实时 Strang cell，继承 2 阶；G₁> 有大虚部，勿 real 截断）。
# 复 Δ：全链条支持（γ/C 源矩阵按 BdG 厄米结构取 conj(Δ)；基准内置 [4] 验证段）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. secondorder/bcs/mixedtime.jl
# =============================================================================
include(joinpath(@__DIR__, "common.jl"))

"""
    _bcs_ring_from_dm(dm, Nτ, ϵd, β) -> (𝒢[1..Nτ], δτ)

τ 支 4 分量 Nambu 环（imagtime.jl 同构方程），数据取自单次 bcs_Δm 调用的
ξ 块（branch(τ,τ) 场，已验证 ≡ bcs_Δτ 环数据）。
"""
function _bcs_ring_from_dm(dm, Nτ, ϵd, β)
    δτ = β / Nτ
    J = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    K = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    for (p, q, blk) in ((1, 1, dm.uu), (1, 2, dm.ud), (2, 1, dm.du), (2, 2, dm.dd))
        ξ = branch(blk; b1=:τ, b2=:τ)
        J[p, q], K[p, q] = ξ.ηⱼₖ, ξ.ηₖⱼ
    end
    u1 = exp(-ϵd * δτ);  f1 = exp(-ϵd * δτ / 2)
    u2 = exp(+ϵd * δτ);  f2 = exp(+ϵd * δτ / 2)
    wpq(p, q, m, k) = k < m ? J[p, q][m-k+1] :
                      k > m ? K[p, q][k-m+1] : J[p, q][1] + K[p, q][1]
    M = zeros(ComplexF64, 4Nτ, 4Nτ)
    b = zeros(ComplexF64, 4Nτ)
    δrow(p, j) = (p == j) ? -1.0 : 0.0
    for p in 1:2, j in 1:2, m in 1:Nτ
        r = ((p-1)*2 + (j-1))*Nτ + m
        up, fp = p == 1 ? (u1, f1) : (u2, f2)
        M[r, r] += 1.0
        if m ≥ 2
            M[r, ((p-1)*2+(j-1))*Nτ + m-1] -= up
        else
            cp = δrow(p, j)
            cp != 0 && (b[r] += up * cp)
            M[r, ((p-1)*2+(j-1))*Nτ + Nτ] += up
        end
        for k in 1:Nτ, s in 1:2
            wv = wpq(p, s, m, k)
            if k == 1
                cpk = δrow(s, j)
                M[r, ((s-1)*2+(j-1))*Nτ + Nτ] += fp * wv * 0.5
                M[r, ((s-1)*2+(j-1))*Nτ + 1] -= fp * wv * 0.5
                cpk != 0 && (b[r] += fp * wv * 0.5 * cpk)
            else
                M[r, ((s-1)*2+(j-1))*Nτ + k-1] -= fp * wv * 0.5
                M[r, ((s-1)*2+(j-1))*Nτ + k]   -= fp * wv * 0.5
            end
        end
    end
    x = M \ b
    return [ComplexF64[x[((p-1)*2+(q-1))*Nτ + m] * (p == q ? 1.0 : -1.0)
                       for p in 1:2, q in 1:2] for m in 1:Nτ], δτ
end

"""
    gamma_nambu(G, mdl, N; j) -> γ (2×L)

族 (1,j) 初关联收缩：γ[1,k] = ⟨c_k↑(0)Φ_j†⟩、γ[2,k] = ⟨c_k↓†(0)Φ_j†⟩
（二阶收敛于 M₀ 对应元；消峰分拆 + 反周期闭合，见文件头 ②）。
"""
function gamma_nambu(G, mdl, N; j::Int)
    h = mdl.β / N
    L = mdl.L
    G0v = ComplexF64[-1-G[N][1,1]  -G[N][1,2]; -G[N][2,1]  -1-G[N][2,2]]
    γ = zeros(ComplexF64, 2, L)
    for k in 1:L
        ε, Δ, V = mdl.ws[k], mdl.Δ, mdl.Vs[k]
        Mk = ComplexF64[-ε -Δ; -conj(Δ) ε]
        P = exp(h * Mk)
        PN = exp(mdl.β * Mk)
        s0 = V .* [G0v[1, j], -G0v[2, j]]
        δs = Vector{Vector{ComplexF64}}(undef, N + 1)
        δs[1] = [0.0, 0.0]
        for i in 1:N
            δs[i+1] = V .* [G[i][1, j], -G[i][2, j]] .- s0
        end
        acc = zeros(ComplexF64, 2)
        for i in 0:N-1
            acc = P * acc + (h / 2) .* (P * δs[i+1] + δs[i+2])
        end
        A0 = Mk \ (PN - Matrix{ComplexF64}(I, 2, 2))
        γ[:, k] .= -((Matrix{ComplexF64}(I, 2, 2) + PN) \ (A0 * s0 + acc))
    end
    return γ
end

"矩阵 Weyl-K 递推（retarded 核组合 + 行相位 + 可选源与初值）。"
function _bcs_march(kcn, knc, kcc, kdu, phases, N; src = nothing,
                    init = Matrix{ComplexF64}(I, 2, 2))
    As = [zeros(ComplexF64, 2, 2) for _ in 1:N+1]
    As[1] = init
    D0 = [kcn(1, 1) kcc(1, 1); kdu(1, 1) knc(1, 1)]
    for m in 1:N
        hist = zeros(ComplexF64, 2, 2)
        for k in 1:m-1
            hist += [kcn(m+1, k+1) kcc(m+1, k+1); kdu(m+1, k+1) knc(m+1, k+1)] *
                    (As[k] + As[k+1]) / 2
        end
        RHS = zeros(ComplexF64, 2, 2)
        for p in 1:2
            up, fp = phases[p]
            RHS[p, :] .= up .* As[m][p, :] .- fp .* (hist[p, :] .+ (D0 * As[m] / 2)[p, :])
        end
        src !== nothing && (RHS .+= src[m])
        Msys = Matrix{ComplexF64}(I, 2, 2)
        for p in 1:2, r in 1:2
            Msys[p, r] += phases[p][2] * D0[p, r] / 2
        end
        for q in 1:2
            As[m+1][:, q] = Msys \ RHS[:, q]
        end
    end
    return As
end

"""
    weyl_k_bcs_mixed(mdl, Nτ, Nt, T) -> G₁>[t_0..t_Nt]

BCS 混合轮廓 K 侧二阶（全平衡 G₁>(t) = ⟨d↑(t)d↑†(0)⟩，不含 −i）。`G[m+1] ≈ G₁>(t_m)`。
数据：单次 bcs_Δm 调用（ξ 环 + cn/nc 实时核）+ 反常槽位的 uv/vu 原始加权 Δm 块。
"""
function weyl_k_bcs_mixed(mdl, Nτ, Nt, T)
    h = T / Nt
    dm = bcs_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=Nτ, t=T, Nt=Nt, μ=mdl.μ, Δ=mdl.Δ)
    # ① τ 环 + 初值数据
    G, _ = _bcs_ring_from_dm(dm, Nτ, mdl.ϵd, mdl.β)
    M11 = 1 + G[Nτ][1, 1]           # ⟨d↑d↑†⟩ = 1−n_d
    F0 = G[Nτ][2, 1]                # 反常 ⟨d̄↓†d↑†⟩
    # ② γ 收缩（族 (1,1)）
    γ = gamma_nambu(G, mdl, Nτ; j = 1)
    # ③ 实时核（Δm 约定）：全部取自组合块 branch 矩阵——
    #    正常槽位 −(P+X)（≡ Δt 的 X−P；fermion 侧同款翻译）；
    #    反常槽位 P+X（≡ Δt 的 P−X，即 −_uvφ₊ + _uvφ₋）。依赖 branch(:+,:−) 的
    #    无共轭派生（MixedCorrelationFunction 的 ηⱼₖ₂ 补全向量；复 Δ 正确）。
    P1 = branch(dm.uu; b1=:+, b2=:+); X1 = branch(dm.uu; b1=:+, b2=:−)
    P2 = branch(dm.dd; b1=:+, b2=:+); X2 = branch(dm.dd; b1=:+, b2=:−)
    P3 = branch(dm.ud; b1=:+, b2=:+); X3 = branch(dm.ud; b1=:+, b2=:−)
    P4 = branch(dm.du; b1=:+, b2=:+); X4 = branch(dm.du; b1=:+, b2=:−)
    kcn(j, k) = -(X1[j, k] + P1[j, k])
    knc(j, k) = -(X2[j, k] + P2[j, k])
    kcc(j, k) = P3[j, k] + X3[j, k]
    kdu(j, k) = P4[j, k] + X4[j, k]
    u = exp(-im * mdl.ϵd * h);  f = exp(-im * mdl.ϵd * h / 2)
    uI = exp(im * mdl.ϵd * h);  fI = exp(im * mdl.ϵd * h / 2)
    phases = [(u, f), (uI, fI)]
    𝒜 = _bcs_march(kcn, knc, kcc, kdu, phases, Nt)
    # ④ C 通道源：Src_k(m) = −i·[expm(h·[[−iD, 𝒱e^{−i𝓗_bk t_{m−1}}],[0, −i𝓗_bk]])]₁₂·γ_k
    D4 = ComplexF64[-im*mdl.ϵd 0; 0 im*mdl.ϵd]
    src = [zeros(ComplexF64, 2, 2) for _ in 1:Nt]
    for k in 1:mdl.L
        ε, Δ, V = mdl.ws[k], mdl.Δ, mdl.Vs[k]
        Hbk = ComplexF64[ε Δ; conj(Δ) -ε]
        B4 = -im * Hbk
        𝒱 = ComplexF64[V 0; 0 -V]
        Eb = Matrix{ComplexF64}(I, 2, 2)      # e^{−i𝓗_bk t_{m−1}}
        for m in 1:Nt
            G4 = [D4 𝒱*Eb; zeros(ComplexF64, 2, 2) B4]
            src[m] .+= (-im * exp(h * G4)[1:2, 3:4]) * γ[:, k]
            Eb *= exp(-im * h * Hbk)
        end
    end
    C = _bcs_march(kcn, knc, kcc, kdu, phases, Nt; src = src,
                   init = zeros(ComplexF64, 2, 2))
    # ⑤ 组装
    return [M11 * 𝒜[m+1][1, 1] + F0 * 𝒜[m+1][1, 2] + C[m+1][1] for m in 0:Nt]
end

"""
    strict_mixed_bcs(mdl, N, T) -> G₁>[t_1..t_N]

严格 cell 链参照：G₁>(t_m) = [(st)^m·M₀]₁₁，st = 实时 Strang cell（2 阶）。
"""
function strict_mixed_bcs(mdl, N, T)
    h = T / N
    𝓗 = bcs_bdg(mdl)
    n = mdl.L + 1
    hA = zeros(ComplexF64, 2n, 2n)
    hA[1, 1] = mdl.ϵd
    hA[n+1, n+1] = -mdl.ϵd
    st = exp(-im * h / 2 * hA) * exp(-im * h * (𝓗 - hA)) * exp(-im * h / 2 * hA)
    M0 = bcs_M0(mdl)
    return [(st^m * M0)[1, 1] for m in 1:N], h
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Δ=0.5, T=2.0, Ns=[16, 32, 64, 128, 256],
              fracs=[0.25, 0.5, 0.75])
    mdl = bcs_toulouse(; Lmodes=Lmodes, Δ=Δ)
    gex = exact_greater_eq(mdl)
    io = stdout
    ms_of(N) = round.(Int, fracs .* N)
    println(io, "=" ^ 74)
    println(io, "BCS 浴混合 Weyl-K 二阶（单次 bcs_Δm：ξ 环 + γ 收缩 + Duhamel C 通道）：全平衡 G₁>(t)")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，Δ=%.2f，T=%.1f（Nτ=Nt 同步）\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64}[]
    println(io, "\n[1] Weyl-K 二阶 vs 精确 [S(t)M₀]₁₁（求解全程不使用精确 M₀/γ）")
    es = Float64[]
    for N in Ns
        G1 = weyl_k_bcs_mixed(mdl, N, N, T)
        e = maximum(abs.([G1[m+1] - gex(m * T / N) for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, T / N, e, NaN))
    end
    order_table(io, "    Weyl-K 二阶", Ns, es)
    println(io, "\n[2] 严格 cell 链参照（[(st)^m·M₀]₁₁，M₀ 仅参照用）")
    es = Float64[]
    for (i, N) in enumerate(Ns)
        Gs, h = strict_mixed_bcs(mdl, N, T)
        e = maximum(abs.([Gs[m] - gex(m * T / N) for m in ms_of(N)]))
        push!(es, e)
        rows[i] = (rows[i][1], rows[i][2], rows[i][3], e)
    end
    order_table(io, "    strict order=2", Ns, es)
    # γ 收缩质量（额外报告，不入 CSV）
    println(io, "\n[3] γ 收缩 vs M₀（报告；族 (1,1)+(1,2)）")
    n = mdl.L + 1
    M0 = bcs_M0(mdl)
    for N in Ns
        dm = bcs_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=N, t=T, Nt=N, μ=mdl.μ, Δ=mdl.Δ)
        G, _ = _bcs_ring_from_dm(dm, N, mdl.ϵd, mdl.β)
        e = 0.0
        for j in 1:2
            slot = j == 1 ? 1 : n + 1
            γ = gamma_nambu(G, mdl, N; j = j)
            for k in 1:mdl.L
                e = max(e, abs(γ[1, k] - M0[1+k, slot]), abs(γ[2, k] - M0[n+1+k, slot]))
            end
        end
        @printf(io, "    N=%4d   err_γ=%.4e\n", N, e)
    end
    println(io, "\n[4] 复数 Δ 验证（Δ = 0.2+0.45im；报告，不入 CSV）")
    mdlc = bcs_toulouse(; Lmodes=Lmodes, Δ=0.2 + 0.45im)
    gexc = exact_greater_eq(mdlc)
    es = Float64[]
    for N in Ns
        G1 = weyl_k_bcs_mixed(mdlc, N, N, T)
        e = maximum(abs.([G1[m+1] - gexc(m * T / N) for m in ms_of(N)]))
        push!(es, e)
    end
    order_table(io, "    Weyl-K 二阶（复 Δ）", Ns, es)
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "mixedtime.csv"), "w") do fh
        println(fh, "N,dt,err_weylk,err_strict2")
        for (N, h, e1, e2) in rows
            @printf(fh, "%d,%.6e,%.6e,%.6e\n", N, h, e1, e2)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "mixedtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
