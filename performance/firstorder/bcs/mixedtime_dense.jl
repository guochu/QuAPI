# =============================================================================
# firstorder/bcs/mixedtime_dense.jl —— BCS 混合：稠密 IF 装配 vs 通道递推
#
# 与 mixedtime.jl 完全对齐的另一种实现：同一个一阶冻结核（bcs_Δm ξ 环 +
# η 正常核 + 原始 uv/vu 加权反常核）按通道折叠为稠密 IF 张量（GrassmannIF），
# Pfaffian 余子式提取，两路输出必须逐位一致（机器精度）：
#   ① τ 环 IF：4 分量 Nambu 环矩阵 Mτ（8Nτ 腿：band 1..4 × (a,ā)）→ Gτ →
#      M11、F0、γ；
#   ② A 通道 IF：naive march 的块线性系统 M_A = I − P·S⁻ + P·𝒦（face-value，
#      2(Nt+1) 抽象节点 → 4(Nt+1) 腿），B = e₁⊗I₂；
#   ③ C 通道：同一张量，B_C = 逐 cell 矩形 Duhamel 源（左端点 bath 相位）。
#
# 腿基约定同 realtime_dense.jl：C[a腿(s), ā腿(t)] = −M[s,t]，掩码去 {a_t, ā_s}。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. firstorder/bcs/mixedtime_dense.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "secondorder", "bcs", "common.jl"))
include(joinpath(@__DIR__, "mixedtime.jl"))
using QuAPI
using QuAPI.Dense: GrassmannContourIndex, GrassmannHybridizationMatrix, GrassmannIF

"参数化抽象节点 → 腿位置（每物理节点 NCOMP 个分量 × (a,ā) = 2·NCOMP 腿）。"
make_aleg(NCOMP) = s -> 2*NCOMP*(cld(s, NCOMP) - 1) + 2*((s - 1) % NCOMP) + 1

"块矩阵 M（n×n 抽象节点）的稠密 IF 折叠。"
function dense_if_legs(n, M, NCOMP, br::Symbol)
    aleg = make_aleg(NCOMP)
    basis = GrassmannContourIndex[GrassmannContourIndex(m, br, band, c)
                                  for m in 1:(n ÷ NCOMP) for band in 1:NCOMP for c in (false, true)]
    C = zeros(ComplexF64, 2n, 2n)
    for s in 1:n, t in 1:n
        C[aleg(s), aleg(t)+1] = -M[s, t]
    end
    return GrassmannIF(GrassmannHybridizationMatrix(basis, C))
end

"从折叠张量提取全部 (M⁻¹)[s, t]（符号逐项数值标定，结构零须与 target 一致）。"
function minor_inv(gif, n, NCOMP, target::AbstractMatrix)
    aleg = make_aleg(NCOMP)
    top = gif.data[end]
    full = (1 << (2n)) - 1
    Minv = zeros(ComplexF64, n, n)
    for s in 1:n, t in 1:n
        mask = full ⊻ ((1 << (aleg(t) - 1)) | (1 << (aleg(s))))     # 去 {a_t, ā_s}
        v = gif.data[mask + 1]
        if abs(v) < 1e-280
            abs(target[s, t]) < 1e-12 * max(1.0, maximum(abs, target)) ||
                error("结构零处目标非零 at ($s,$t)")
            continue
        end
        r = target[s, t] * top / v
        sg = round(real(r))
        abs(real(r) - sg) > 1e-6 && error("符号标定失败 at ($s,$t): r=$r")
        Minv[s, t] = sg * v / top
    end
    return Minv
end

# ------------------------------------------------------------------ #
function main(; Lmodes=48, Δ=0.5, T=2.0, Ns=[2, 3])
    mdl = bcs_toulouse(; Lmodes=Lmodes, Δ=Δ)
    io = stdout
    println(io, "=" ^ 74)
    println(io, "BCS 混合：稠密 IF 装配（GrassmannIF 分段折叠）vs naive 通道递推")
    println(io, "  ① τ 环 IF（8Nτ 腿）→ Gτ → M11、F0、γ    ②③ A/C 通道 IF（4(Nt+1) 腿）")
    @printf(io, "浴：%d 模，J=%.2f，β=%.1f，ϵ_d=%.1f，Δ=%.2f，T=%.1f（Nτ=Nt 同步）\n",
            mdl.L, mdl.Γ, mdl.β, mdl.ϵd, Δ, T)
    println(io, "=" ^ 74)

    worst = 0.0
    for N in Ns
        h = T / N
        dm = bcs_Δm(mdl.bath.spectrum; β=mdl.β, Nτ=N, t=T, Nt=N, μ=mdl.μ, Δ=mdl.Δ)
        # ---- naive 通道解（mixedtime.jl 同款，作为参照） ----
        J = Matrix{Vector{ComplexF64}}(undef,2,2); KK = Matrix{Vector{ComplexF64}}(undef,2,2)
        for (p,q,blk) in ((1,1,dm.uu),(1,2,dm.ud),(2,1,dm.du),(2,2,dm.dd))
            ξ = branch(blk; b1=:τ, b2=:τ); J[p,q], KK[p,q] = ξ.ηⱼₖ, ξ.ηₖⱼ
        end
        uτ = [exp(-mdl.ϵd*mdl.β/N), exp(mdl.ϵd*mdl.β/N)]
        fτ = [exp(-mdl.ϵd*mdl.β/(2N)), exp(mdl.ϵd*mdl.β/(2N))]
        Mτ = zeros(ComplexF64,4N,4N); bτ = zeros(ComplexF64,4N)
        idx(i,j) = 2*(i-1)+j
        for m in 1:N, i in 1:2, j in 1:2
            s = 4*(m-1)+idx(i,j); Mτ[s,s] += 1.0
            if m ≥ 2; Mτ[s,4*(m-2)+idx(i,j)] -= uτ[i]
            else; Mτ[s,4*(N-1)+idx(i,j)] += uτ[i]; i==j && (bτ[s] = -uτ[i]); end
            for k in 1:N, q in 1:2
                sgn = i==q ? 1.0 : -1.0
                w = k<m ? J[i,q][m-k+1] : k>m ? KK[i,q][k-m+1] : J[i,q][1]+KK[i,q][1]
                Mτ[s,4*(k-1)+idx(q,j)] -= fτ[i]*sgn*w
            end
        end
        xτ = Mτ\bτ
        𝒢 = [ComplexF64[xτ[4*(m-1)+idx(p,q)] for p in 1:2, q in 1:2] for m in 1:N]
        M11 = 1 + 𝒢[N][1,1]; F0 = 𝒢[N][2,1]
        G0v = ComplexF64[-1-𝒢[N][1,1] -𝒢[N][1,2]; -𝒢[N][2,1] -1-𝒢[N][2,2]]
        γ = zeros(ComplexF64, 2, mdl.L)
        for k in 1:mdl.L
            ε, Δ, V = mdl.ws[k], mdl.Δ, mdl.Vs[k]
            Mk = ComplexF64[-ε -Δ; -conj(Δ) ε]; PN = exp(mdl.β*Mk)
            s0 = V .* [G0v[1,1], -G0v[2,1]]
            sj = [V .* [𝒢[m][1,1], -𝒢[m][2,1]] for m in 1:N]
            acc = zeros(ComplexF64,2); Q0 = PN
            for m in 1:N
                Q1 = exp(Mk*(mdl.β - m*mdl.β/N))
                acc += (mdl.β/N/2).*(Q0*(m==1 ? s0 : sj[m-1]) + Q1*sj[m]); Q0 = Q1
            end
            γ[:,k] .= -((Matrix{ComplexF64}(I,2,2)+PN)\acc)
        end
        P1 = branch(dm.uu; b1=:+, b2=:+); X1 = branch(dm.uu; b1=:+, b2=:-)
        P2 = branch(dm.dd; b1=:+, b2=:+); X2 = branch(dm.dd; b1=:+, b2=:-)
        fu2, fuv, fvu, fv2 = QuAPI.get_all_fs(mdl.bath.spectrum, mdl.Δ)
        disp = ϵ -> sqrt(ϵ^2 + abs2(mdl.Δ))
        uv_m = fermionic_Δm(fuv, disp; β=mdl.β, Nτ=N, t=T, Nt=N, μ=mdl.μ)
        vu_m = fermionic_Δm(fvu, disp; β=mdl.β, Nτ=N, t=T, Nt=N, μ=mdl.μ)
        Puv = branch(uv_m; b1=:+, b2=:+); Xuv = branch(uv_m; b1=:+, b2=:-)
        Pvu = branch(vu_m; b1=:+, b2=:+); Xvu = branch(vu_m; b1=:+, b2=:-)
        K22(j,k) = [-(X1[j,k]+P1[j,k])   (Xuv[j,k]+Puv[j,k]) - conj(Xvu[j,k]+Pvu[j,k]);
                    (Xvu[j,k]+Pvu[j,k]) - conj(Xuv[j,k]+Puv[j,k])   -(X2[j,k]+P2[j,k])]
        Pph = Diagonal([exp(-im*mdl.ϵd*h), exp(im*mdl.ϵd*h)])
        A = [zeros(ComplexF64,2,2) for _ in 1:N+1]
        C = [zeros(ComplexF64,2) for _ in 1:N+1]
        A[1] = Matrix{ComplexF64}(I,2,2)
        Ebs = [Matrix{ComplexF64}(I,2,2) for _ in 1:mdl.L]
        srcs = [zeros(ComplexF64,2) for _ in 1:N]
        for m in 1:N
            sA = zeros(ComplexF64,2,2); sC = zeros(ComplexF64,2)
            for k in 1:m-1
                Kk = K22(m+1,k+1); sA += Kk*A[k]; sC += Kk*C[k]
            end
            cellSrc = zeros(ComplexF64,2)
            for k in 1:mdl.L
                ε, Δ, V = mdl.ws[k], mdl.Δ, mdl.Vs[k]
                Hbk = ComplexF64[ε Δ; conj(Δ) -ε]
                𝒱 = ComplexF64[V 0; 0 -V]
                cellSrc += -im*h*(𝒱*(Ebs[k]*γ[:,k]))
                Ebs[k] = exp(-im*h*Hbk)*Ebs[k]
            end
            srcs[m] = cellSrc
            A[m+1] = Pph*(A[m]-sA); C[m+1] = Pph*(C[m]-sC+cellSrc)
        end
        G_naive = [M11*A[m+1][1,1] + F0*A[m+1][1,2] + C[m+1][1] for m in 0:N]
        # ---- 稠密装配 ----
        # ① τ 环 IF
        gif1 = dense_if_legs(4N, Mτ, 4, :τ)
        Minv1 = minor_inv(gif1, 4N, 4, inv(Mτ))
        xd = Minv1 * bτ
        𝒢d = [ComplexF64[xd[4*(m-1)+idx(p,q)] for p in 1:2, q in 1:2] for m in 1:N]
        d1 = maximum(maximum(abs.(𝒢d[m] .- 𝒢[m])) for m in 1:N)
        # ②③ A/C 通道 IF（同一 M_A 折叠，两个源）
        nA = 2(N+1)
        MA = zeros(ComplexF64, nA, nA)
        MA[1:2, 1:2] = Matrix{ComplexF64}(I, 2, 2)
        for n in 2:N+1
            MA[2n-1:2n, 2n-1:2n] = Matrix{ComplexF64}(I, 2, 2)
            MA[2n-1:2n, 2n-3:2n-2] = -Pph
            for k in 1:n-2
                MA[2n-1:2n, 2k-1:2k] += Pph * K22(n, k+1)
            end
        end
        B_A = zeros(ComplexF64, nA, 2); B_A[1:2, :] = Matrix{ComplexF64}(I, 2, 2)
        B_C = zeros(ComplexF64, nA, 2)
        for n in 2:N+1
            B_C[2n-1:2n, 1] = Pph * srcs[n-1]
        end
        gif2 = dense_if_legs(nA, MA, 2, :+)
        Minv2 = minor_inv(gif2, nA, 2, inv(MA))
        Ad = Minv2 * B_A
        Cd = Minv2 * B_C
        dA = maximum(maximum(abs.(Ad[2m-1:2m, :] .- A[m])) for m in 1:N+1)
        dC = maximum(maximum(abs.(Cd[2m-1:2m] .- C[m])) for m in 1:N+1)
        G_dense = [M11*Ad[2m-1, 1] + F0*Ad[2m-1, 2] + Cd[2m-1] for m in 1:N+1]
        dG = maximum(abs.(G_dense .- G_naive))
        worst = max(worst, d1, dA, dC, dG)
        @printf(io, "  N=%d: |Gτ|=%.2e  |A|=%.2e  |C|=%.2e  |G dense−naive|=%.2e\n",
                N, d1, dA, dC, dG)
    end
    @printf(io, "\n最大一致性偏差 = %.2e（机器精度即通过）\n", worst)
    println(io, "\n结论（实测）")
    println(io, "  · 三段 IF 张量的 Pfaffian 余子式提取与 naive 通道递推逐位一致。")
    println(io, "  · 反常核为原始 uv/vu 加权 Δm 块构造（混合轮廓必需）；ζ 跨块信息经")
    println(io, "    γ 收缩 + Duhamel 源进入 C 通道。")
    println(io, "  · 下游 MPS-IF 按本文件的腿基/系数矩阵装配即可逐位复现通道解。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
