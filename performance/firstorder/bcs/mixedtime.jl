# =============================================================================
# firstorder/bcs/mixedtime.jl —— BCS 浴混合轮廓：一阶冻结核（naive）vs 严格
# 1/2 阶 Trotter。对照：严格解（全平衡 G₁>(t) = [S(t)M₀]₁₁）。
# 输出 results/mixedtime.csv
#
# 与 secondorder/bcs/mixedtime.jl 同骨架、同核数据（单次 bcs_Δm 调用 + 原始
# uv/vu 加权 Δm 块），但全部 K 侧要素退回 naive：
#   ① τ 支 4 分量 Nambu 环：ξ 块（branch(τ,τ)）朴素矩形权重（无梯形 0.5、
#      无跳变常数代入）；跳变源只在 Nambu 对角分量（imagtime.jl 同构）；
#   ② γ 初关联收缩：梯形 + 跳变延拓端点，无消峰分拆（一阶病态来源）；
#   ③ 实时 march：face-value 递推（无半步 f、无接触对角隐式、无梯形）：
#      X_{m+1} = P·(X_m − Σ_{k≤m−1} 𝒦(m+1,k+1)·X_k + [C] h·S_m)，
#      正常槽位 −(P+X)（Δm 约定）、反常槽位 Ŵ_uv − conj(Ŵ_vu)（原始块）；
#   ④ C 通道源：矩形 Duhamel 卷积 S_k(m) = −i·h·Σ_{j<m} E_D(m−j)𝒱_k E_b(j)·γ_k
#      （E_D = diag(e^{∓iϵ_d··})，E_b = e^{−i𝓗_bk t}，𝒱_k = diag(V_k,−V_k)）；
#   ⑤ 组装：G₁>(t_m) = M11·𝒜₁₁ + F0·𝒜₁₂ + C₁（M11 = 1+𝒢₁₁[Nτ]、F0 = 𝒢₂₁[Nτ]，
#      跳变反解）。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bcs/mixedtime.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "bcs", "common.jl"))

"""
    bcs_mixed_naive(mdl, Nτ, Nt, T) -> G₁>[t_0..t_Nt]

BCS 混合一阶冻结核（naive 装配）。`G[m+1] ≈ G₁>(t_m)`（不含 −i）。
"""
function bcs_mixed_naive(mdl, Nτ, Nt, T)
    h = T / Nt
    δτ = mdl.β / Nτ
    dm = bcs_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=Nτ, t=T, Nt=Nt, μ=mdl.μ, Δ=mdl.Δ)
    # ① τ 环（4 分量 Nambu，朴素矩形；数据取自单次 bcs_Δm 的 ξ 块）
    J = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    K = Matrix{Vector{ComplexF64}}(undef, 2, 2)
    for (p, q, blk) in ((1, 1, dm.uu), (1, 2, dm.ud), (2, 1, dm.du), (2, 2, dm.dd))
        ξ = branch(blk; b1=:τ, b2=:τ)
        J[p, q], K[p, q] = ξ.ηⱼₖ, ξ.ηₖⱼ
    end
    uτ = [exp(-mdl.ϵd * δτ), exp(mdl.ϵd * δτ)]
    fτ = [exp(-mdl.ϵd * δτ / 2), exp(mdl.ϵd * δτ / 2)]
    M = zeros(ComplexF64, 4Nτ, 4Nτ)
    b = zeros(ComplexF64, 4Nτ)
    idx(i, j) = 2*(i-1) + j
    for m in 1:Nτ, i in 1:2, j in 1:2
        s = 4*(m-1) + idx(i, j)
        M[s, s] += 1.0
        if m ≥ 2
            M[s, 4*(m-2)+idx(i, j)] -= uτ[i]
        else
            M[s, 4*(Nτ-1)+idx(i, j)] += uτ[i]
            i == j && (b[s] = -uτ[i])
        end
        for k in 1:Nτ
            for q in 1:2                              # 𝚺̂ 反常槽位带 −1
                w = k < m ? J[i, q][m-k+1] : k > m ? K[i, q][k-m+1] :
                    J[i, q][1] + K[i, q][1]
                sgn = i == q ? 1.0 : -1.0
                M[s, 4*(k-1)+idx(q, j)] -= fτ[i] * sgn * w
            end
        end
    end
    x = M \ b
    𝒢 = [ComplexF64[x[4*(m-1)+idx(p, q)] for p in 1:2, q in 1:2] for m in 1:Nτ]
    M11 = 1 + 𝒢[Nτ][1, 1]
    F0 = 𝒢[Nτ][2, 1]
    # ② γ 初关联收缩（族 (1,1)；梯形 + 跳变端点延拓，无消峰）
    G0v = ComplexF64[-1-𝒢[Nτ][1, 1]  -𝒢[Nτ][1, 2];
                     -𝒢[Nτ][2, 1]    -1-𝒢[Nτ][2, 2]]
    γ = zeros(ComplexF64, 2, mdl.L)
    for k in 1:mdl.L
        ε, Δ, V = mdl.ws[k], mdl.Δ, mdl.Vs[k]
        Mk = ComplexF64[-ε -Δ; -conj(Δ) ε]
        PN = exp(mdl.β * Mk)
        s0 = V .* [G0v[1, 1], -G0v[2, 1]]
        sj = [V .* [𝒢[m][1, 1], -𝒢[m][2, 1]] for m in 1:Nτ]
        acc = zeros(ComplexF64, 2)
        Q0 = PN
        for m in 1:Nτ
            Q1 = exp(Mk * (mdl.β - m * δτ))
            acc += (δτ / 2) .* (Q0 * (m == 1 ? s0 : sj[m-1]) + Q1 * sj[m])
            Q0 = Q1
        end
        γ[:, k] .= -((Matrix{ComplexF64}(I, 2, 2) + PN) \ acc)
    end
    # ③ 实时核（Δm 约定）：正常 −(P+X)；反常用原始 uv/vu 加权 Δm 块
    #    Ŵ = X+P、k₁₂ = Ŵ_uv − conj(Ŵ_vu)（混合轮廓必需——组合块 cc/nn 的
    #    naive march 配对不收敛，见 secondorder/bcs/mixedtime.jl 机制说明）
    P1 = branch(dm.uu; b1=:+, b2=:+); X1 = branch(dm.uu; b1=:+, b2=:-)
    P2 = branch(dm.dd; b1=:+, b2=:+); X2 = branch(dm.dd; b1=:+, b2=:-)
    fu2, fuv, fvu, fv2 = QuAPI.get_all_fs(mdl.bath.spectrum, mdl.Δ)
    disp = ϵ -> sqrt(ϵ^2 + abs2(mdl.Δ))
    uv_m = fermionic_Δm(fuv, disp; β=mdl.β, Nτ=Nτ, t=T, Nt=Nt, μ=mdl.μ)
    vu_m = fermionic_Δm(fvu, disp; β=mdl.β, Nτ=Nτ, t=T, Nt=Nt, μ=mdl.μ)
    Puv = branch(uv_m; b1=:+, b2=:+); Xuv = branch(uv_m; b1=:+, b2=:-)
    Pvu = branch(vu_m; b1=:+, b2=:+); Xvu = branch(vu_m; b1=:+, b2=:-)
    K22(j, k) = [-(X1[j, k] + P1[j, k])   (Xuv[j, k] + Puv[j, k]) - conj(Xvu[j, k] + Pvu[j, k]);
                 (Xvu[j, k] + Pvu[j, k]) - conj(Xuv[j, k] + Puv[j, k])   -(X2[j, k] + P2[j, k])]
    Pph = Diagonal([exp(-im * mdl.ϵd * h), exp(im * mdl.ϵd * h)])
    # ④ A / C 通道（face-value 递推；C 带逐 cell 矩形 Duhamel 源）
    #    cellSrc(m) = −i·h·Σ_k 𝒱_k·e^{−i𝓗_bk t_m}·γ_k（bath 相位冻结在 cell 左端点，
    #    −i 因子同 secondorder Duhamel 约定；E_D 传播由 Pph 承担）
    A = [zeros(ComplexF64, 2, 2) for _ in 1:Nt+1]
    C = [zeros(ComplexF64, 2) for _ in 1:Nt+1]
    A[1] = Matrix{ComplexF64}(I, 2, 2)
    Ebs = [Matrix{ComplexF64}(I, 2, 2) for _ in 1:mdl.L]    # e^{−i𝓗_bk t_m}
    for m in 1:Nt
        sA = zeros(ComplexF64, 2, 2)
        sC = zeros(ComplexF64, 2)
        for k in 1:m-1
            Kk = K22(m+1, k+1)
            sA += Kk * A[k]
            sC += Kk * C[k]
        end
        cellSrc = zeros(ComplexF64, 2)
        for k in 1:mdl.L
            ε, Δ, V = mdl.ws[k], mdl.Δ, mdl.Vs[k]
            Hbk = ComplexF64[ε Δ; conj(Δ) -ε]
            𝒱 = ComplexF64[V 0; 0 -V]
            cellSrc += -im * h * (𝒱 * (Ebs[k] * γ[:, k]))
            Ebs[k] = exp(-im * h * Hbk) * Ebs[k]           # e^{−i𝓗_bk t_m}
        end
        A[m+1] = Pph * (A[m] - sA)
        C[m+1] = Pph * (C[m] - sC + cellSrc)
    end
    # ⑤ 组装
    return [M11 * A[m+1][1, 1] + F0 * A[m+1][1, 2] + C[m+1][1] for m in 0:Nt]
end

"""
    strict_mixed_bcs(mdl, N, T, order) -> G₁>[t_1..t_N]

严格 cell 链参照：G₁>(t_m) = [(st)^m·M₀]₁₁，st = 实时 Lie / Strang cell。
"""
function strict_mixed_bcs(mdl, N, T, order::Int)
    h = T / N
    𝓗 = bcs_bdg(mdl)
    n = mdl.L + 1
    hA = zeros(ComplexF64, 2n, 2n)
    hA[1, 1] = mdl.ϵd
    hA[n+1, n+1] = -mdl.ϵd
    eA = exp(-im * h / 2 * hA)
    eB = exp(-im * h * (𝓗 - hA))
    st = order == 1 ? exp(-im * h * hA) * eB : eA * eB * eA
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
    println(io, "BCS 浴混合：一阶冻结核（naive：矩形环+无消峰γ+face-value march）")
    println(io, "  vs 严格 1/2 阶 Trotter（全平衡 G₁>(t)，求解不用精确 M₀/γ）")
    @printf(io, "浴：%d 模，J=%.2f on [−3,3]，β=%.1f，ϵ_d=%.1f，Δ=%.2f，T=%.1f（Nτ=Nt 同步）\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ, T)
    println(io, "=" ^ 74)
    rows = Tuple{Int,Float64,Float64,Float64,Float64}[]
    println(io, "\n[1] 一阶冻结核 vs 精确 [S(t)M₀]₁₁")
    es = Float64[]
    for N in Ns
        G = bcs_mixed_naive(mdl, N, N, T)
        e = maximum(abs.([G[m+1] - gex(m * T / N) for m in ms_of(N)]))
        push!(es, e)
        push!(rows, (N, T / N, e, NaN, NaN))
    end
    order_table(io, "    一阶冻结核（naive 装配）", Ns, es)
    println(io, "\n[2] 严格 Trotter 参照（M₀ 仅参照用）")
    for order in (1, 2)
        es = Float64[]
        for (i, N) in enumerate(Ns)
            Gs, _ = strict_mixed_bcs(mdl, N, T, order)
            e = maximum(abs.([Gs[m] - gex(m * T / N) for m in ms_of(N)]))
            push!(es, e)
            rows[i] = (rows[i][1], rows[i][2], rows[i][3],
                       order == 1 ? e : rows[i][4], order == 2 ? e : rows[i][5])
        end
        order_table(io, "    strict order=$order", Ns, es)
    end
    println(io, "\n[3] Δ=0 退化（报告）：反常链应关闭，退化为正常费米 naive 混合")
    mdl0 = bcs_toulouse(; Lmodes=Lmodes, Δ=0.0)
    gex0 = exact_greater_eq(mdl0)
    for N in (16, 32)
        G = bcs_mixed_naive(mdl0, N, N, T)
        e = maximum(abs.([G[m+1] - gex0(m * T / N) for m in ms_of(N)]))
        @printf(io, "    N=%3d   err=%.4e\n", N, e)
    end
    out = joinpath(@__DIR__, "results")
    mkpath(out)
    open(joinpath(out, "mixedtime.csv"), "w") do fh
        println(fh, "N,dt,err_frozen1,err_strict1,err_strict2")
        for (N, h, a, c, d) in rows
            @printf(fh, "%d,%.6e,%.6e,%.6e,%.6e\n", N, h, a, c, d)
        end
    end
    println(io, "\n数据已写入 $(joinpath(out, "mixedtime.csv"))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
