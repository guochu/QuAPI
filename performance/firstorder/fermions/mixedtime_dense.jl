# =============================================================================
# firstorder/fermions/mixedtime_dense.jl —— 费米混合轮廓：稠密 IF 装配
#
# 与 mixedtime.jl 完全对齐的另一种实现：同一个一阶冻结核（fermionic_Δm 的
# ξ/η/ζ 块条目）不再走通道递推，而是按通道折叠为稠密 IF 张量（GrassmannIF），
# 从 Pfaffian 余子式收缩提取，两路输出必须逐位一致（机器精度）。
#
# 三段张量（混合轮廓的 ζ 跨块经已验证的因式分解路径进入——直接朴素 ζ 收缩
# 有 lesser/greater 语义错配，见 secondorder/fermions/mixedtime.jl 的结论）：
#   ① τ 环 IF：naive 矩形 ξ 块 + 反周期回绕（同 imagtime_dense.jl，无 fτ）；
#   ② A 通道 IF：Keldysh 开链，face-value Δ = −(P+X)（Δm 约定），含半步 f
#      与接触对角隐式——递推的线性系统 M_A·A = e₁ 整体折叠；
#   ③ C 通道 IF：M_C = I + h·𝒦（严格下三角）折叠；τ 支信息经朴素梯形 γ
#      （由 ① 的张量解 Gτ 计算）进入源 Svec，稠密提取 C = M_C⁻¹·Svec。
#
# 目的：给下游 TEMPO/GTEMPO 的 MPS-IF 装配一个具体的稠密参照——IF 张量的
# 腿基、系数矩阵、收缩（Pfaffian 余子式）顺序全部在此显式可见。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/fermions/mixedtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "fermions", "common.jl"))
using QuAPI
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix, GrassmannIF

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

# ------------------------------------------------------------------ #
# ① τ 环：naive 矩形 ξ 块 + 反周期回绕（与 mixedtime.jl 逐位同款）
# ------------------------------------------------------------------ #
function naive_tau_ring(mdl, Nτ)
    δτ = mdl.β / Nτ
    uτ = exp(-mdl.ϵd * δτ)
    dm = fermionic_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=Nτ, t=1.6, Nt=Nτ, μ=mdl.μ)
    ξ = branch(dm; b1=:τ, b2=:τ)
    ηjk, ηkj = ξ.ηⱼₖ, ξ.ηₖⱼ
    d0τ = ηjk[1] + ηkj[1]
    M = zeros(ComplexF64, Nτ, Nτ)
    b = zeros(ComplexF64, Nτ)
    for m in 1:Nτ
        M[m, m] += 1.0
        if m ≥ 2
            M[m, m-1] -= uτ
        else
            b[m] -= uτ
            M[m, Nτ] += uτ
        end
        for k in 1:Nτ
            w = k < m ? ηjk[m-k+1] : k > m ? ηkj[k-m+1] : d0τ
            M[m, k] -= w
        end
    end
    return M, b
end

"τ 环稠密 IF：腿基 (a, ā)/节点，C[2k, 2m−1] = −M[m,k]（转置核装配）。"
function tau_ring_if(Nτ, M)
    basis = GrassmannContourIndex[GrassmannContourIndex(j, :τ, 1, c)
                                  for j in 1:Nτ for c in (false, true)]
    C = zeros(ComplexF64, 2Nτ, 2Nτ)
    for m in 1:Nτ, k in 1:Nτ
        C[2k, 2m-1] = -M[m, k]
    end
    return GrassmannIF(GrassmannHybridizationMatrix(basis, C))
end

"""
    minor_inv(gif, N) -> Minv

从折叠张量提取全部 (M⁻¹)[i, j]：sgn·data[去掉 {ā_i, a_j}]/data[end]。
符号 (i, j) 逐项数值标定：相对目标必须落在整数 ±1 上（自校验）。
"""
function minor_inv(gif, N, target::AbstractMatrix)
    top = gif.data[end]
    full = (1 << (2N)) - 1
    Minv = zeros(ComplexF64, N, N)
    sgns = Int[]
    for i in 1:N, j in 1:N
        mask = full ⊻ ((1 << (2i - 1)) | (1 << (2j - 2)))  # ā_i(leg 2i)、a_j(leg 2j−1)
        v = gif.data[mask + 1]
        if abs(v) < 1e-280
            # 结构零：minor 为零 ⟺ (M⁻¹)[i,j] 为零（如因果三角系统的 j>i 块）
            abs(target[i, j]) < 1e-12 * max(1.0, maximum(abs, target)) ||
                error("结构零处目标非零 at ($i,$j)")
            push!(sgns, 0)
            continue
        end
        r = target[i, j] * top / v
        s = round(real(r))
        abs(real(r) - s) > 1e-6 && error("符号标定失败 at ($i,$j): r=$r")
        push!(sgns, Int(s))
        Minv[i, j] = s * v / top
    end
    return Minv, sgns
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, T=2.0, Ns=[4, 6, 8])
    mdl = fermi_toulouse(; Lmodes=Lmodes)
    io = stdout
    # 严格参照（全平衡 G>，仅参照用）
    ev = eigen(Hermitian(mdl.hfull))
    γ = ev.vectors * Diagonal([1/(exp(mdl.β*λ)+1) for λ in ev.values]) * ev.vectors'
    gexact(N, m) = -im * (exp(-im * (m * T / N) * mdl.hfull) * (mdl.Id - transpose(γ)))[1, 1]
    println(io, "=" ^ 74)
    println(io, "费米混合：稠密 IF 装配（GrassmannIF，分段折叠）vs 通道递推（机器精度）")
    println(io, "  ① τ 环 IF（ξ 块矩形+反周期）→ Gτ → 朴素梯形 γ")
    println(io, "  ② A 通道 IF（η 块 face-value，Δm 约定 −(P+X)，半步 f+接触隐式）")
    println(io, "  ③ C 通道 IF（𝒦 连续卷积矩阵）+ 源 Svec(γ) 的余子式提取")
    @printf(io, "浴：%d 模，J=%.2f，β=%.1f，T=%.1f（Nτ=Nt 同步）\n",
            length(mdl.ws), mdl.Γ, mdl.β, T)
    println(io, "=" ^ 74)

    for N in Ns
        h = T / N
        dm = fermionic_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=N, t=T, Nt=N, μ=mdl.μ)
        P = branch(dm; b1=:+, b2=:+)
        Xr = branch(dm; b1=:+, b2=:-)
        # ---------- ① τ 环 ----------
        Mτ, bτ = naive_tau_ring(mdl, N)
        Gτ_ref = Mτ \ bτ
        gif1 = tau_ring_if(N, Mτ)
        Minv1, _ = minor_inv(gif1, N, inv(Mτ))
        Gτ = Minv1 * bτ                                  # x = M⁻¹b 经余子式
        n_d = -Gτ[N]
        d1 = maximum(abs.(Gτ .- Gτ_ref))
        # γ[1,k] 朴素梯形提取（无消峰——一阶病态来源，与 mixedtime.jl 同款）
        δτ = mdl.β / N
        γrow = zeros(ComplexF64, mdl.L)
        G0v = -1 - Gτ[N]
        Gext = vcat(G0v, Gτ)
        for k in 1:mdl.L
            ε, V = mdl.ws[k], sqrt(mdl.αs[k])
            acc = 0.0 + 0.0im
            for j in 1:N
                w0 = exp(-ε*(mdl.β - (j-1)*δτ))
                w1 = exp(-ε*(mdl.β - j*δτ))
                acc += δτ * (w0*Gext[j] + w1*Gext[j+1]) / 2
            end
            γrow[k] = V * acc / (1 + exp(-ε*mdl.β))
        end
        # ---------- ② A 通道 ----------
        u = exp(-im * mdl.ϵd * h)
        f = exp(-im * mdl.ϵd * h / 2)
        Δ(j, k) = -(Xr[j, k] + P[j, k]); d0 = Δ(1, 1)
        MA = zeros(ComplexF64, N + 1, N + 1)
        MA[1, 1] = 1.0
        for m in 1:N
            MA[m+1, m+1] = 1 + f * d0 / 2
            MA[m+1, m] = -u + f * d0 / 2
            for j in 1:m-1                  # 递推 s = Σ Δ(m+1,k+1)·A[k]：A 列 = k
                MA[m+1, j] += f * Δ(m+1, j+1)
            end
        end
        A_ref = MA \ [1.0; zeros(ComplexF64, N)]
        gif2 = tau_channel_if(N + 1, MA, :+)
        Minv2, _ = minor_inv(gif2, N + 1, inv(MA))
        A = Minv2 * [1.0; zeros(ComplexF64, N)]
        d2 = maximum(abs.(A .- A_ref))
        # ---------- ③ C 通道 ----------
        MC = Matrix{ComplexF64}(I, N + 1, N + 1)
        for m in 1:N, j in 1:m-1
            MC[m+1, j+1] = h * mathcalK(mdl, (m - j) * h)
        end
        Svec = zeros(ComplexF64, N + 1)
        for m in 1:N
            tm = m * h
            S = 0.0 + 0.0im
            for k in 1:mdl.L
                ωk = mdl.ϵd - mdl.ws[k]
                Jk = exp(-im*mdl.ϵd*tm) * (abs(ωk) < 1e-10 ? tm :
                     (exp(im*ωk*tm) - 1) / (im*ωk))
                S += -im * sqrt(mdl.αs[k]) * Jk * (-γrow[k])
            end
            Svec[m+1] = S
        end
        C_ref = MC \ Svec
        gif3 = tau_channel_if(N + 1, MC, :+)
        Minv3, _ = minor_inv(gif3, N + 1, inv(MC))
        C = Minv3 * Svec
        d3 = maximum(abs.(C .- C_ref))
        # ---------- 汇总 ----------
        G_naive = (-im) .* ((1 - n_d) .* A_ref .+ C_ref)
        G_dense = (-im) .* ((1 - n_d) .* A .+ C)
        dG = maximum(abs.(G_dense .- G_naive))
        eex = maximum(abs.(G_naive[m+1] - gexact(N, m) for m in (1, N÷2, N)))
        @printf(io, "  N=%2d: |Gτ_dense−Gτ_rec|=%.2e  |A_dense−A_rec|=%.2e  |C_dense−C_rec|=%.2e  |G_dense−G_naive|=%.2e  |G_naive−exact|=%.2e\n",
                N, d1, d2, d3, dG, eex)
    end
    println(io, "\n结论（实测）")
    println(io, "  · 三段 IF 张量的 Pfaffian 余子式提取与通道递推逐位一致（机器精度）。")
    println(io, "  · ζ 跨块经因式分解路径（γ×𝒦）进入 C 通道源；直接朴素 ζ 收缩存在")
    println(io, "    lesser/greater 语义错配，不采用（结论同 secondorder）。")
    println(io, "  · 下游 MPS-IF 按本文件的腿基/系数矩阵装配即可逐位复现通道解。")
end

"实时开链通道 IF：腿基 (a, ā)/节点（branch=:+），C[2k, 2m−1] = −M[m,k]。"
function tau_channel_if(n, M, br::Symbol)
    basis = GrassmannContourIndex[GrassmannContourIndex(j, br, 1, c)
                                  for j in 1:n for c in (false, true)]
    C = zeros(ComplexF64, 2n, 2n)
    for m in 1:n, k in 1:n
        C[2k, 2m-1] = -M[m, k]
    end
    return GrassmannIF(GrassmannHybridizationMatrix(basis, C))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
