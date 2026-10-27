# 现行费米核（fermionic_Δτ/Δt 冻结源）的 G 标度检验 —— 目标：G⁴ 平台
# （与玻色 non-diagonal PT 同级：O(G²) 系数精确 + 残差按 G⁴ 标度，G 减半比值 ≈16）
#
# 参照：exact_fermionic_Δτ/Δt 提供**严格方案**的机器精确值（同骨架、同腿基），
# 冻结核与其之差 = 纯核 G 截断误差（方案误差完全抵消）。
#      |Z_frozen − Z_exact| ~ G⁴ ⟹ 现行核 O(g²) 内容与严格 IF 一致（达标）
#      ~ G²（比值≈4）          ⟹ O(g²) 内容有缺口（对角项/常数），需修
# 另打印对物理链 (B·A)^N 的误差作背景（方案误差，两边同阶）。
using QuAPI, QuAPI.Dense, ImpurityModelBase, LinearAlgebra, Printf

const _I2 = Matrix{Float64}(I, 2, 2)
const _Z2 = [1.0 0.0; 0.0 -1.0]
const _SP2 = ComplexF64[0 0; 1 0]
_fock_cdag(K, m) = kron([i < m ? _Z2 : (i == m ? _SP2 : _I2) for i in 1:K+1]...)
function fock_chain_ops(εs, Vs, ϵ_d)
    K = length(εs)
    ddag = _fock_cdag(K, 1); dd = ddag'
    H_B = zeros(ComplexF64, 2^(K+1), 2^(K+1))
    for k in 1:K
        cdag = _fock_cdag(K, k+1); cc = cdag'
        H_B .+= εs[k] * (cdag * cc) .+ Vs[k] * (ddag * cc .+ cdag * dd)
    end
    return H_B, ϵ_d * (ddag * dd)
end

function tau_skeleton!(C, N, x)
    for i in 1:N; C[2i, 2i-1] += -1; end
    for i in 2:N; C[2i, 2(i-1)-1] += x; end
    C[2, 2N-1] += -x
    return C
end
Ztau(C, basis, N, pref) =
    (-1.0)^N * Zvalue(GrassmannIF(GrassmannHybridizationMatrix(basis, C))) * pref

prt(j0, br, ā) = 4 * j0 + (br === :+ ? 0 : 2) + (ā ? 2 : 1)
rtbasis(N) = [GrassmannContourIndex(j0 + 1, br, 1, c)
              for j0 in 0:N for br in (:+, :-) for c in (false, true)]
function mono_from_seq(seq)
    s = 0
    for i in 1:length(seq), j in (i+1):length(seq); seq[i] > seq[j] && (s += 1); end
    mask = 0
    for p in seq; mask |= 1 << (p - 1); end
    mask, isodd(s) ? -1 : 1
end
function merge_sign(A, B)
    s = 0
    while B > 0
        lb = B & -B; B -= lb
        s += count_ones(A & ~(2 * lb - 1))
    end
    isodd(s) ? -1 : 1
end
sf(terms) = [isempty(seq) ? (0, ComplexF64(c)) :
             (m -> (m[1], ComplexF64(c) * m[2]))(mono_from_seq(seq)) for (seq, c) in terms]
function gds(P, S)
    out = zeros(ComplexF64, length(P))
    for (m, c) in S, m2 in 0:length(P)-1
        v = P[m2+1]; iszero(v) && continue
        (m & m2) == 0 || continue
        out[(m | m2) + 1] += c * v * merge_sign(m2, m)
    end
    out
end
function rt_skeleton!(C, N, x)
    for j0 in 0:N, br in (:+, :-)
        C[prt(j0, br, true), prt(j0, br, false)] += -1
    end
    for k in 1:N
        C[prt(k, :+, true), prt(k-1, :+, false)] += x
        C[prt(k-1, :-, true), prt(k, :-, false)] += conj(x)
    end
    return C
end
function rho_from_C(C, pref, N, ρ0)
    gif = GrassmannIF(GrassmannHybridizationMatrix(rtbasis(N), ComplexF64.(C)))
    P = ComplexF64.(gif.data)
    one_ = 1.0 + 0.0im
    P = gds(P, sf([(Int[], ρ0[1,1]), ([prt(0,:+,true)], -ρ0[2,1]),
                   ([prt(0,:-,false)], -ρ0[1,2]),
                   ([prt(0,:+,true), prt(0,:-,false)], ρ0[2,2])]))
    ρ = zeros(ComplexF64, 2, 2)
    for a in 0:1, b in 0:1
        seq = Int[]
        a == 1 && push!(seq, prt(N,:+,false))
        b == 1 && push!(seq, prt(N,:-,true))
        Q = gds(P, sf([(seq, (-one_)^(a+b))]))
        ρ[a+1,b+1] = (a == 1 && b == 0 ? -1 : 1) * Q[end]
    end
    return pref .* ρ
end

const EPSD = 0.6
const BETA = 1.2
const RHO0 = ComplexF64[0.62 0.2-0.15im; 0.2+0.15im 0.38]

println("=" ^ 76)
println("现行费米核（冻结源）vs 严格方案（exact 参照）：G 标度（目标 G⁴ ⟹ 比值≈16）")
println("=" ^ 76)

println("-- 虚时 τ（单模 DiracDelta ω=0.9）--")
for N in (2, 4)
    println("  N = $N:")
    prev = nothing
    for G in (0.4, 0.2, 0.1)
        spec = DiracDelta(0.9; α = G^2)
        x = exp(-EPSD * BETA / N)
        hfr = GrassmannHybridizationMatrix(fermionic_Δτ(spec; β = BETA, N = N))
        Zfr = Ztau(tau_skeleton!(Matrix{ComplexF64}(hfr.Δ), N, x), hfr.basis, N, 1.0)
        ex = exact_fermionic_Δτ(spec; β = BETA, N = N)
        Zex = Ztau(tau_skeleton!(copy(ex.hmat.Δ), N, x), ex.hmat.basis, N, ex.prefactor)
        H_B, H_A = fock_chain_ops([0.9], [G], EPSD)
        Zch = tr((exp(-(BETA/N) .* H_B) * exp(-(BETA/N) .* H_A))^N) / (1 + exp(-BETA * 0.9))
        d = abs(Zfr - Zex)
        @printf("    G=%.2f: |Zfr−Zex| = %.3e%s   [物理链: |Zfr−Zch|=%.2e |Zex−Zch|=%.2e]\n",
                G, d, prev === nothing ? "          " : "  比值=$(round(prev/d, digits=1))",
                abs(Zfr - Zch), abs(Zex - Zch))
        prev = d
    end
end

println("-- 实时 RT（单模 DiracDelta ω=0.9, N=2, T=0.9）--")
for N in (2, 3)
    T = 0.9
    println("  N = $N, T = $T:")
    prev = nothing
    for G in (0.4, 0.2, 0.1)
        spec = DiracDelta(0.9; α = G^2)
        x = exp(-im * EPSD * T / N)
        hfr = GrassmannHybridizationMatrix(fermionic_Δt(spec; β = BETA, t = T, N = N))
        ρfr = rho_from_C(rt_skeleton!(Matrix{ComplexF64}(hfr.Δ), N, x), 1.0, N, RHO0)
        ex = exact_fermionic_Δt(spec; β = BETA, t = T, N = N)
        ρex = rho_from_C(rt_skeleton!(copy(ex.hmat.Δ), N, x), ex.prefactor, N, RHO0)
        H_B, H_A = fock_chain_ops([0.9], [G], EPSD)
        U = (exp(-im*(T/N) .* H_B) * exp(-im*(T/N) .* H_A))^N
        n = 1/(exp(BETA*0.9)+1)
        R = U * kron(RHO0, ComplexF64[1-n 0; 0 n]) * U'
        T4 = reshape(R, 2, 2, 2, 2)
        ρch = [sum(T4[i,a,i,b] for i in 1:2) for a in 1:2, b in 1:2]
        d = norm(ρfr .- ρex)
        @printf("    G=%.2f: ‖ρfr−ρex‖ = %.3e%s   [物理链: ‖ρfr−ρch‖=%.2e ‖ρex−ρch‖=%.2e]\n",
                G, d, prev === nothing ? "          " : "  比值=$(round(prev/d, digits=1))",
                norm(ρfr .- ρch), norm(ρex .- ρch))
        prev = d
    end
end

println("=" ^ 76)
println("修正后标准装配（GrassmannHybridizationMatrix 顶点移位）vs 严格方案：")
println("  (a) 原始 Z：残差 = 缺 O(g²) 常数（比值≈4）   (b) 除以 exact prefactor：应 G⁴（≈16）")
println("=" ^ 76)

println("-- 虚时 τ --")
for N in (2, 4)
    println("  N = $N:")
    prevb = nothing
    for G in (0.4, 0.2, 0.1)
        spec = DiracDelta(0.9; α = G^2)
        x = exp(-EPSD * BETA / N)
        hfr = GrassmannHybridizationMatrix(fermionic_Δτ(spec; β = BETA, N = N))
        Zfr = Ztau(tau_skeleton!(Matrix{ComplexF64}(hfr.Δ), N, x), hfr.basis, N, 1.0)
        ex = exact_fermionic_Δτ(spec; β = BETA, N = N)
        Zex = Ztau(tau_skeleton!(copy(ex.hmat.Δ), N, x), ex.hmat.basis, N, ex.prefactor)
        d = abs(Zfr * ex.prefactor - Zex)     # (b)：补上缺失的 O(g²) 常数后
        @printf("    G=%.2f: (b) |Zfr·P − Zex| = %.3e%s\n", G, d,
                prevb === nothing ? "          " : "  比值=$(round(prevb/d, digits=1))")
        prevb = d
    end
end

println("-- 实时 RT（T=0.9，迹归一化 ρ——常数自动消去）--")
for N in (2, 3)
    T = 0.9
    println("  N = $N, T = $T:")
    prev = nothing
    for G in (0.4, 0.2, 0.1)
        spec = DiracDelta(0.9; α = G^2)
        x = exp(-im * EPSD * T / N)
        hfr = GrassmannHybridizationMatrix(fermionic_Δt(spec; β = BETA, t = T, N = N))
        ρfr = rho_from_C(rt_skeleton!(Matrix{ComplexF64}(hfr.Δ), N, x), 1.0, N, RHO0)
        ρfr ./= tr(ρfr)
        ex = exact_fermionic_Δt(spec; β = BETA, t = T, N = N)
        ρex = rho_from_C(rt_skeleton!(copy(ex.hmat.Δ), N, x), ex.prefactor, N, RHO0)
        ρex ./= tr(ρex)
        d = norm(ρfr .- ρex)
        @printf("    G=%.2f: ‖ρ̂fr − ρ̂ex‖ = %.3e%s\n", G, d,
                prev === nothing ? "          " : "  比值=$(round(prev/d, digits=1))")
        prev = d
    end
end

println("-- 混合 Δm（全迹 Z，补 prefactor 后应 G⁴）--")
for (Nτ, Nt) in ((2, 2), (3, 1))
    T = 0.8
    println("  Nτ = $Nτ, Nt = $Nt:")
    prev = nothing
    for G in (0.4, 0.2, 0.1)
        spec = DiracDelta(0.9; α = G^2)
        n = Nτ + 2 * (Nt + 1)
        ipos(c, ā) = 2 * c - (ā ? 0 : 1)
        cτ(j) = j; cp(jj) = Nτ + 2jj - 1; cm(jj) = Nτ + 2jj
        function mxC0()
            C0 = zeros(ComplexF64, 2n, 2n)
            for c in 1:n; C0[ipos(c, true), ipos(c, false)] -= 1; end
            C0[ipos(cp(1), true), ipos(cτ(Nτ), false)] += 1
            C0[ipos(cm(Nt + 1), true), ipos(cp(Nt + 1), false)] += -1
            C0[ipos(cτ(1), true), ipos(cm(1), false)] += 1
            for k in 2:Nτ;  C0[ipos(cτ(k), true), ipos(cτ(k - 1), false)] += 1; end
            for jj in 1:Nt
                C0[ipos(cp(jj + 1), true), ipos(cp(jj), false)] += 1
                C0[ipos(cm(jj), true), ipos(cm(jj + 1), false)] += 1
            end
            return C0
        end
        basis = vcat([GrassmannContourIndex(j, :τ, 1, c) for j in 1:Nτ for c in (false, true)],
                     [GrassmannContourIndex(jj, br, 1, c)
                      for jj in 1:Nt+1 for br in (:+, :-) for c in (false, true)])
        hfr = GrassmannHybridizationMatrix(fermionic_Δm(spec; β = BETA, Nτ = Nτ, t = T, Nt = Nt))
        Zfr = Zvalue(GrassmannIF(GrassmannHybridizationMatrix(
            basis, mxC0() + Matrix{ComplexF64}(hfr.Δ))))
        ex = exact_fermionic_Δm(spec; β = BETA, Nτ = Nτ, t = T, Nt = Nt)
        Zex = Zvalue(GrassmannIF(GrassmannHybridizationMatrix(
            basis, mxC0() + copy(ex.hmat.Δ)))) * ex.prefactor
        d = abs(Zfr * ex.prefactor - Zex)
        @printf("    G=%.2f: |Zfr·P − Zex| = %.3e%s\n", G, d,
                prev === nothing ? "          " : "  比值=$(round(prev/d, digits=1))")
        prev = d
    end
end
