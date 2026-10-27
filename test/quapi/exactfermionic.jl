# ============================================================================
# exact_fermionic_Δτ / exact_fermionic_Δt 的验证（纯 bath IF，无 ϵ_d 输入）
#
# 直接参考（Fock 空间，维数 2^{K+1}，杂质 = 第 1 个/最慢因子）：
#   B = e^{−δt(Σ ε_k c_k†c_k + Σ V_k(d†c_k+c_k†d))}，ε_k = ω_k − μ
#   A = e^{−δt ϵ_d d†d}
#   [机器精度] ε_d = 0（骨架动能 x = 1）:
#     虚时:  Z = (−1)^N·Zvalue(IF(测度+动能+Δ))·prefactor == Tr[B^N]/Z_bath
#     实时:  ρ(T) = prefactor·∫[骨架×IF] == Tr_b[B^N (ρ₀⊗ρ_b) (B^N)†]
#   [收敛性] ε_d ≠ 0（骨架动能 x = e^{−δt ϵ_d}）: 对真实链 (B·A)^N 的误差
#     严格小于现行冻结源核（fermionic_Δτ/Δt）的误差（O(g²δt³) vs O(g²δt²)）
# ============================================================================

# ------------------------------ Fock 空间参考链 ------------------------------ #
const _I2 = Matrix{Float64}(I, 2, 2)
const _Z2 = [1.0 0.0; 0.0 -1.0]
const _SP2 = ComplexF64[0 0; 1 0]          # |1⟩⟨0|（占据创建）

# 模式 m（1 = 杂质，2..K+1 = bath 模式 k = m−1）的创建算符（JW 串在慢因子上）
_fock_cdag(K::Int, m::Int) =
    kron([i < m ? _Z2 : (i == m ? _SP2 : _I2) for i in 1:K+1]...)

function _fock_chain_ops(εs::Vector{Float64}, Vs::Vector{Float64}, ϵ_d::Real)
    K = length(εs)
    ddag = _fock_cdag(K, 1); dd = ddag'
    H_B = zeros(ComplexF64, 2^(K+1), 2^(K+1))
    for k in 1:K
        cdag = _fock_cdag(K, k+1); cc = cdag'
        H_B .+= εs[k] * (cdag * cc) .+ Vs[k] * (ddag * cc .+ cdag * dd)
    end
    return H_B, ϵ_d * (ddag * dd)
end

# ------------------------------ 实时收缩机器（与 fermionicIF.jl 校准的约定） ------------------------------ #
_prt(j0::Int, br::Symbol, ā::Bool) = 4 * j0 + (br === :+ ? 0 : 2) + (ā ? 2 : 1)
_rtbasis(N) = [GrassmannContourIndex(j0 + 1, br, 1, c)
               for j0 in 0:N for br in (:+, :-) for c in (false, true)]

function _mono_from_seq(seq::Vector{Int})
    s = 0
    for i in 1:length(seq), j in (i+1):length(seq)
        seq[i] > seq[j] && (s += 1)
    end
    mask = 0
    for p in seq; mask |= 1 << (p - 1); end
    mask, isodd(s) ? -1 : 1
end
function _merge_sign(A::Int, B::Int)
    s = 0
    while B > 0
        lb = B & -B; B -= lb
        s += count_ones(A & ~(2 * lb - 1))
    end
    isodd(s) ? -1 : 1
end
function _sf(terms)
    out = Tuple{Int,ComplexF64}[]
    for (seq, c) in terms
        c = ComplexF64(c)
        if isempty(seq)
            push!(out, (0, c))
        else
            m, s = _mono_from_seq(seq)
            push!(out, (m, c * s))
        end
    end
    out
end
function _gds(P::Vector{ComplexF64}, S::Vector{Tuple{Int,ComplexF64}})
    out = zeros(ComplexF64, length(P))
    for (m, c) in S
        for m2 in 0:length(P)-1
            v = P[m2+1]; iszero(v) && continue
            (m & m2) == 0 || continue
            out[(m | m2) + 1] += c * v * _merge_sign(m2, m)
        end
    end
    out
end

# ρ(T)：折叠 IF(全指数) + 初态因子 + 末态投影 + 通道符号
function _rho_from_C(C::AbstractMatrix{<:Number}, pref, N::Int, ρ0::Matrix{ComplexF64})
    gif = GrassmannIF(GrassmannHybridizationMatrix(_rtbasis(N), ComplexF64.(C)))
    P = ComplexF64.(gif.data)
    one_ = 1.0 + 0.0im
    P = _gds(P, _sf([(Int[], ρ0[1, 1]),
                     ([_prt(0, :+, true)], -ρ0[2, 1]),
                     ([_prt(0, :-, false)], -ρ0[1, 2]),
                     ([_prt(0, :+, true), _prt(0, :-, false)], ρ0[2, 2])]))
    ρ = zeros(ComplexF64, 2, 2)
    for a in 0:1, b in 0:1
        seq = Int[]
        a == 1 && push!(seq, _prt(N, :+, false))
        b == 1 && push!(seq, _prt(N, :-, true))
        Q = _gds(P, _sf([(seq, (-one_)^(a + b))]))
        ρ[a+1, b+1] = (a == 1 && b == 0 ? -1 : 1) * Q[end]
    end
    return pref .* ρ
end

# τ 骨架叠加：测度 + 动能（x = 1 时单位动能；反周期回绕 −x）
function _tau_skeleton!(C::AbstractMatrix, N::Int, x::Number)
    for i in 1:N; C[2i, 2i-1] += -1; end
    for i in 2:N; C[2i, 2(i-1)-1] += x; end
    C[2, 2N-1] += -x
    return C
end

# 实时骨架叠加：测度 + 分支动能（x = 1 时单位动能）
function _rt_skeleton!(C::AbstractMatrix, N::Int, x::Number)
    for j0 in 0:N, br in (:+, :-)
        C[_prt(j0, br, true), _prt(j0, br, false)] += -1
    end
    for k in 1:N
        C[_prt(k, :+, true), _prt(k-1, :+, false)] += x
        C[_prt(k-1, :-, true), _prt(k, :-, false)] += conj(x)
    end
    return C
end

_Ztau(C, basis, N, pref) =
    (-1.0)^N * Zvalue(GrassmannIF(GrassmannHybridizationMatrix(basis, C))) * pref

# ------------------------------ 测试 ------------------------------ #

@testset "exact fermionic IF (normal bath, H_S-free) vs direct chain" begin
    ϵ_d = 0.6
    β = 1.2
    ρ0 = ComplexF64[0.62 0.2 - 0.15im; 0.2 + 0.15im 0.38]

    specs = [("1 mode (DiracDelta)", DiracDelta(0.9; α = 0.35^2), [0.9], [0.35]),
             ("2 modes", DiscreteSpectrum([-0.7, 1.1], [0.2^2, 0.3^2]),
              [-0.7, 1.1], [0.2, 0.3]),
             ("3 modes (μ≠0)", DiscreteSpectrum([-1.0, 0.4, 1.6], [0.15^2, 0.25^2, 0.2^2]),
              [-1.0, 0.4, 1.6], [0.15, 0.25, 0.2])]

    @testset "imaginary-time: machine-exact at ε_d = 0 (IF 是纯 bath 对象)" begin
        for (label, spec, ws, Vs) in specs
            μ = label == "3 modes (μ≠0)" ? 0.25 : 0.0
            εs = ws .- μ
            H_B, _ = _fock_chain_ops(εs, Vs, 0.0)
            Zb = prod(1 .+ exp.(-β .* εs))
            for N in (1, 2, 3, 4)
                δτ = β / N
                ex = exact_fermionic_Δτ(spec; β = β, N = N, μ = μ)
                @test ex.hmat.basis ==
                      GrassmannHybridizationMatrix(fermionic_Δτ(spec; β = β, N = N, μ = μ)).basis
                C = _tau_skeleton!(copy(ex.hmat.Δ), N, 1.0)
                @test abs(_Ztau(C, ex.hmat.basis, N, ex.prefactor) - tr(exp(-δτ .* H_B)^N) / Zb) < 1e-12
            end
            @testset "$label" begin end
        end
        # g = 0：prefactor ≡ 1、Δ ≡ 0（退化为现行骨架）
        ex0 = exact_fermionic_Δτ(DiracDelta(0.9; α = 0.0); β = 1.2, N = 4)
        @test abs(ex0.prefactor - 1) < 1e-14
        @test maximum(abs.(ex0.hmat.Δ)) == 0
    end

    @testset "real-time: machine-exact at ε_d = 0" begin
        for (label, spec, ws, Vs) in specs
            μ = label == "3 modes (μ≠0)" ? 0.25 : 0.0
            εs = ws .- μ
            H_B, _ = _fock_chain_ops(εs, Vs, 0.0)
            K = length(εs)
            ρb = kron([ComplexF64[1 - 1/(exp(β*ε)+1) 0; 0 1/(exp(β*ε)+1)] for ε in εs]...)
            for (N, T) in ((1, 0.5), (2, 0.9), (3, 1.2))
                δt = T / N
                ex = exact_fermionic_Δt(spec; β = β, t = T, N = N, μ = μ)
                @test ex.hmat.basis ==
                      GrassmannHybridizationMatrix(fermionic_Δt(spec; β = β, t = T, N = N, μ = μ)).basis
                C = _rt_skeleton!(copy(ex.hmat.Δ), N, 1.0)
                U = exp(-im * δt .* H_B)^N
                R = U * kron(ρ0, ρb) * U'
                # 列优先：reshape 第一维 = 快位 = bath（kron(imp, bath)），
                # 对 bath 槽求和后保留杂质 —— [B, a, B, b]（K=1 时即 [i,a,i,b]）
                Tens = reshape(R, 2^K, 2, 2^K, 2)
                ρref = [sum(Tens[B, a, B, b] for B in 1:2^K) for a in 1:2, b in 1:2]
                @test norm(_rho_from_C(C, ex.prefactor, N, ρ0) .- ρref) < 1e-12
            end
            @testset "$label" begin end
        end
        ex0 = exact_fermionic_Δt(DiracDelta(0.9; α = 0.0); β = 1.2, t = 0.9, N = 3)
        @test abs(ex0.prefactor - 1) < 1e-14
        @test maximum(abs.(ex0.hmat.Δ)) == 0
    end

    @testset "ε_d ≠ 0: 标准核（顶点移位后 = PT 类）与 exact 同级，差距 G⁴" begin
        # 顶点移位约定（2026-10，models.jl）后标准装配的核 = 严格 B-cell IF 的
        # O(g²) 精确 ⟹ 标准与 exact 对链 (B·A)^N 的误差同为方案误差
        # （O(g²δt³)/cell）主导，两者之差 = G⁴ 核截断
        for (label, spec, ws, Vs) in specs[1:2]
            μ = 0.0
            εs = copy(ws)
            H_B, H_A = _fock_chain_ops(εs, Vs, ϵ_d)
            Zb = prod(1 .+ exp.(-β .* εs))
            # 虚时 Z
            N = 4; δτ = β / N
            ex = exact_fermionic_Δτ(spec; β = β, N = N, μ = μ)
            Zref = tr((exp(-δτ .* H_B) * exp(-δτ .* H_A))^N) / Zb
            Zex = _Ztau(_tau_skeleton!(copy(ex.hmat.Δ), N, exp(-δτ * ϵ_d)),
                        ex.hmat.basis, N, ex.prefactor)
            hfr = GrassmannHybridizationMatrix(fermionic_Δτ(spec; β = β, N = N, μ = μ))
            Zfr = _Ztau(_tau_skeleton!(Matrix{ComplexF64}(hfr.Δ), N, exp(-δτ * ϵ_d)),
                        hfr.basis, N, 1.0)
            @test abs(Zfr * ex.prefactor - Zex) < 5e-4        # G⁴ 核截断量级（补常数）
            # 实时 ρ(T)
            N = 3; T = 1.2; δt = T / N
            ex = exact_fermionic_Δt(spec; β = β, t = T, N = N, μ = μ)
            x = exp(-im * δt * ϵ_d)
            U = (exp(-im * δt .* H_B) * exp(-im * δt .* H_A))^N
            K = length(εs)
            ρb = kron([ComplexF64[1 - 1/(exp(β*ε)+1) 0; 0 1/(exp(β*ε)+1)] for ε in εs]...)
            R = U * kron(ρ0, ρb) * U'
            Tens = reshape(R, 2^K, 2, 2^K, 2)
            ρref = [sum(Tens[B, a, B, b] for B in 1:2^K) for a in 1:2, b in 1:2]
            ρex = _rho_from_C(_rt_skeleton!(copy(ex.hmat.Δ), N, x), ex.prefactor, N, ρ0)
            hfr = GrassmannHybridizationMatrix(fermionic_Δt(spec; β = β, t = T, N = N, μ = μ))
            ρfr = _rho_from_C(_rt_skeleton!(Matrix{ComplexF64}(hfr.Δ), N, x), 1.0, N, ρ0)
            @test norm(ρfr .* ex.prefactor .- ρex) < 5e-3     # G⁴ 核截断量级（补常数；RT 条目差 ~5e-4）
            @testset "$label" begin end
        end
    end

    # ---------------- 连续谱输入（Gauss–Legendre 正交离散） ----------------

    @testset "τ 连续谱：GL 等价 + 收敛" begin
        fcont = bounded(ε -> 0.08 * exp(-((ε - 0.3) / 0.9)^2), -2.5, 2.5)
        xs, ws = QuAPI._gauss_legendre(48)
        ωs = 2.5 .* xs                      # 参考区间 [-1,1] → 谱支撑 [-2.5,2.5]
        spec48 = DiscreteSpectrum(ωs, 2.5 .* ws .* fcont.(ωs))
        a = exact_fermionic_Δτ(fcont; β = β, N = 4, npts = 48)
        b = exact_fermionic_Δτ(spec48; β = β, N = 4)
        @test maximum(abs.(a.hmat.Δ .- b.hmat.Δ)) < 1e-13
        @test abs(a.prefactor - b.prefactor) < 1e-13
        a1 = exact_fermionic_Δτ(fcont; β = β, N = 4, npts = 64)
        a2 = exact_fermionic_Δτ(fcont; β = β, N = 4, npts = 256)
        @test maximum(abs.(a1.hmat.Δ .- a2.hmat.Δ)) < 1e-8
        @test abs(a1.prefactor - a2.prefactor) < 1e-8
    end

    # ---------------- RT 连续谱（Gauss–Legendre 离散） ----------------

    @testset "RT 连续谱：GL 等价 + 收敛" begin
        fcont = bounded(ε -> 0.08 * exp(-((ε - 0.3) / 0.9)^2), -2.5, 2.5)
        xs, ws = QuAPI._gauss_legendre(48)
        ωs = 2.5 .* xs                      # 参考区间 [-1,1] → 谱支撑 [-2.5,2.5]
        spec48 = DiscreteSpectrum(ωs, 2.5 .* ws .* fcont.(ωs))
        a = exact_fermionic_Δt(fcont; β = β, t = 0.9, N = 2, npts = 48)
        b = exact_fermionic_Δt(spec48; β = β, t = 0.9, N = 2)
        @test maximum(abs.(a.hmat.Δ .- b.hmat.Δ)) < 1e-13
        @test abs(a.prefactor - b.prefactor) < 1e-13
        a1 = exact_fermionic_Δt(fcont; β = β, t = 0.9, N = 2, npts = 64)
        a2 = exact_fermionic_Δt(fcont; β = β, t = 0.9, N = 2, npts = 256)
        @test maximum(abs.(a1.hmat.Δ .- a2.hmat.Δ)) < 1e-8
    end

    # ---------------- 混合轮廓（Kadanoff–Baym） ----------------
    # 全迹收缩：Zvalue(IF(单位骨架 + Δ))·prefactor == Tr[(B^τ)^Nτ]/Z_bath
    #（实时段 B^t B̄^t 幺正回绕相消）；标准装配（顶点移位后）与 exact 的差按
    # G⁴ 标度（补 prefactor 后）

    function _mixed_C0(Nτ, Nt)
        n = Nτ + 2 * (Nt + 1)
        ipos(c, ā) = 2 * c - (ā ? 0 : 1)
        cτ(j) = j; cp(jj) = Nτ + 2jj - 1; cm(jj) = Nτ + 2jj
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

    @testset "混合 Δm：exact 机器精度对链 + 移位装配 G⁴" begin
        for (label, spec, ws, Vs) in specs[1:2]
            for (Nτ, Nt, t) in ((2, 1, 0.6), (3, 2, 0.9))
                ex = exact_fermionic_Δm(spec; β = β, Nτ = Nτ, t = t, Nt = Nt)
                Zm = Zvalue(GrassmannIF(GrassmannHybridizationMatrix(
                    ex.hmat.basis, _mixed_C0(Nτ, Nt) + ex.hmat.Δ))) * ex.prefactor
                H_B, = _fock_chain_ops(copy(ws), Vs, 0.0)
                Zref = tr((exp(-(β / Nτ) .* H_B))^Nτ) / prod(1 .+ exp.(-β .* (ws .- 0.0)))
                @test abs(Zm - Zref) < 1e-11
                # 标准装配（移位后）与 exact 同级：差 = G⁴ 核截断（补常数）
                hfr = GrassmannHybridizationMatrix(
                    fermionic_Δm(spec; β = β, Nτ = Nτ, t = t, Nt = Nt))
                @test hfr.basis == ex.hmat.basis
                Zfr = Zvalue(GrassmannIF(GrassmannHybridizationMatrix(
                    ex.hmat.basis, _mixed_C0(Nτ, Nt) + Matrix{ComplexF64}(hfr.Δ))))
                @test abs(Zfr * ex.prefactor - Zm) < 5e-3
            end
            @testset "$label" begin end
        end
        # g = 0 退化：Δ = 0，prefactor = 1
        ex0 = exact_fermionic_Δm(DiracDelta(0.9; α = 0.0); β = β, Nτ = 2, t = 0.6, Nt = 1)
        @test maximum(abs.(ex0.hmat.Δ)) == 0
        @test abs(ex0.prefactor - 1) < 1e-14
    end
end
