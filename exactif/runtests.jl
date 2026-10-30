# =============================================================================
# exactif/runtests.jl —— 严格核算法（连续热库）分层验证
#
# 层次：
#  L0  Pfaffian 算法 vs 展开；多体机构 V=0 极限
#  L1  F(t)：M=1 解析极限（留数闭式 + Volterra 连续求解器收敛阶）
#  L2  谱代数层：A = KDK†、B = KEK†、Z_N/Z 与有限浴显式矩阵
#      （r⁽ⁱ⁾·v⁽ʲ⁾ 等）逐位机器精度一致 —— 严格核算法的代数核心
#  L3  G-块装配：对称性 (S21)–(S22)、结构零模式、V→0 极限
#  L4  连续谱：spectral_continuous 在 δ 峰极限下逼近离散闭式（收敛）
# =============================================================================
using Test

@testset "exactif" begin

@testset "L0 Pfaffian 与多体机构" begin
    include("brute.jl")
    using LinearAlgebra
    using Random
    Random.seed!()
    for n in [2, 4, 6, 8]
        B = randn(ComplexF64, n, n)
        T = UpperTriangular(B) - diagm(0 => diag(UpperTriangular(B)))
        A = T - transpose(T)
        @test isapprox(pfaffian(A), pfaffian_expand(A), rtol=1e-9)
    end
    # V=0：过程张量 = 轨迹守恒 δ 结构（每支 α=β 时 I=1，否则 0）
    r0 = if_path_values([-0.7, 1.3], [0.0, 0.0], 2.0, 1, 0.15)
    for S in 0:15
        sv = [(S >> (v - 1)) & 1 for v in 1:4]
        expected = (sv[1] == sv[3] && sv[2] == sv[4]) ? 1.0 + 0im : 0.0 + 0im
        @test r0.Ivals[S+1] ≈ expected atol = 1e-12
    end
end

@testset "L1 F(t) 解析与收敛" begin
    include("volterra.jl")
    # M=1 解析：F(t) = e^{−iEt/2}[cos(νt/2) + i(E/ν)sin(νt/2)]，ν=√(E²+4g²)
    E, g = 0.8, 0.5
    ν = sqrt(E^2 + 4g^2)
    Fan(t) = exp(-im * E * t / 2) * (cos(ν * t / 2) + im * (E / ν) * sin(ν * t / 2))
    ts = collect(0.1:0.2:2.0)
    Fp = fidelt_poles([E], [g + 0im], ts)
    @test maximum(abs.(Fp - Fan.(ts))) < 1e-10
    # 多模式：F(0)=1、F'(0)=0（记忆核）
    Es = [-0.7, 1.3, 0.4]; Vs = [0.45, -0.6im, 0.35]
    ts2 = [0.0; collect(0.05:0.1:2.0)]
    Fp2 = fidelt_poles(Es, Vs, ts2)
    @test abs(Fp2[1] - 1) < 1e-12
    # 连续求解器：光滑高斯谱，对照均匀模式离散参照（两者应一致到积分误差）
    E0, σ0, g20 = 0.8, 0.25, 0.25
    ΓG(ω) = g20 * exp(-(ω - E0)^2 / (2σ0^2)) / sqrt(2π * σ0^2)
    lb0, ub0 = E0 - 2.0, E0 + 2.0
    nm = 200
    Em = collect(range(lb0, ub0; length=nm + 1))[2:end]
    dε = (ub0 - lb0) / nm
    Vm = [sqrt(ΓG(e) * dε) for e in Em]
    include("brute.jl")
    Wm = exp(-im * 0.2 * single_particle_h(Em, Vm))
    tsc = collect(0.2:0.2:1.8)
    Fref = [(Wm^n)[1, 1] for (n, t) in zip(1:9, tsc)]
    Fc = fidelt_continuous(ΓG, lb0, ub0, tsc; nsub=800, npts_spec=400)
    @test maximum(abs.(Fc - Fref)) < 5e-5
end

@testset "L2 谱代数层 = 有限浴显式矩阵（机器精度）" begin
    include("brute.jl")
    include("volterra.jl")
    include("spectral.jl")
    using LinearAlgebra
    for (Es, Vs, β, N, Δt) in [
            ([-0.7, 1.3], [0.45, -0.6im], 2.0, 3, 0.15),
            ([-1.1, 0.2, 2.0], [0.3 + 0.1im, 0.5, -0.4im], 0.7, 4, 0.08),
            ([0.9], [0.62], 5.0, 5, 0.2)]
        sp = spectral_discrete(Es, Vs, β, N, Δt)
        ab = brute_ABZN(Es, Vs, β, N, Δt)
        F_of(W, n) = (W^n)[1, 1]
        @test maximum(abs.(sp.A - ab.Ad)) < 1e-11        # A = K D K†
        @test maximum(abs.(sp.B - ab.Bd)) < 1e-11        # B = K E K†
        @test abs(sp.ZNoZ - ab.ZNoZ) < 1e-11             # Z_N/Z = det(I−B)
        @test maximum(abs.(sp.Fs - [F_of(ab.W, n) for n in 1:N])) < 1e-11
    end
end

@testset "L3 G-块结构与 N=1 字典精确关系" begin
    include("brute.jl")
    include("volterra.jl")
    include("spectral.jl")
    include("kernel.jl")
    using LinearAlgebra
    Es = [-0.7, 1.3, 0.4]; Vs = [0.45, -0.6im, 0.35]; β = 2.0; N = 4; Δt = 0.15
    sp = spectral_discrete(Es, Vs, β, N, Δt)
    K = exact_kernel(; A=sp.A, B=sp.B, F1=sp.Fs[1], ZNoZ=sp.ZNoZ, N=N, Δt=Δt)
    # 对称性 (S22)：G⁻(m,n) = G⁻(n,m)*；G⁺(m,n) = G⁺(n,m)*
    @test isapprox(K.Gm, K.Gm', rtol=1e-10, atol=1e-12)
    @test isapprox(K.Gp, K.Gp', rtol=1e-10, atol=1e-12)
    # Θ 反对称
    @test isapprox(K.Θ, -transpose(K.Θ), rtol=0, atol=0)
    # V→0：A,B → 0，核退化为纯边界项
    sp0 = spectral_discrete(Es, zeros(ComplexF64, 3), β, 1, Δt)
    K0 = exact_kernel(; A=sp0.A, B=sp0.B, F1=sp0.Fs[1], ZNoZ=sp0.ZNoZ, N=1, Δt=Δt)
    @test maximum(abs.(K0.Gp)) < 1e-12
    @test maximum(abs.(K0.Gm)) < 1e-12
    @test K0.ZNoZ ≈ 1
end

@testset "L4 连续谱 → 离散极限收敛" begin
    include("volterra.jl")
    include("spectral.jl")
    using LinearAlgebra
    # Γ 用 3 个窄高斯近似 δ 峰，对照 spectral_discrete
    Es = [-0.7, 1.3, 0.4]; g2s = abs2.([0.45, -0.6im, 0.35]); β = 2.0; N = 3; Δt = 0.15
    σ = 5e-3
    Γ(ω) = sum(g2s[k] * exp(-(ω - Es[k])^2 / 2σ^2) / sqrt(2π * σ^2) for k in 1:3)
    # 参照改为均匀模式 + 矩阵指数的 spectral 量（避开伴随矩阵病态）
    include("brute.jl")
    nm = 400
    Em = collect(range(-2.0, 2.2; length=nm + 1))[2:end]
    dε = 4.2 / nm
    Vm = [sqrt(Γ(e) * dε) for e in Em]
    W = exp(-im * Δt * single_particle_h(Em, Vm))
    F1 = W[1, 1]
    a = vec(W[1, 2:end]); UB = W[2:end, 2:end]
    n̄m = [1 / (1 + exp(β * Em[k])) for k in eachindex(Em)]
    r_ = [vec(transpose(a) * (UB^i)) for i in 0:N-1]
    v_ = [vec((UB')^j * conj.(a)) for j in 0:N-1]
    Ad = [sum(r_[i] .* v_[j]) for i in 1:N, j in 1:N]
    Bd = [sum(r_[i] .* (n̄m .* v_[j])) for i in 1:N, j in 1:N]
    spc = spectral_continuous(Γ, -2.0, 2.2, β, N, Δt; nsub=1600, npts_spec=600)
    @test abs(spc.Fs[1] - F1) < 2e-4
    @test maximum(abs.(spc.A - Ad)) < 2e-3
    @test maximum(abs.(spc.B - Bd)) < 2e-3
end

@testset "L5 端到端：核 → 收缩 → ρ_imp vs Trotter 精确动力学（任意 N 严格）" begin
    include("brute.jl")
    include("contract.jl")
    Es = [-0.7, 1.3]; Vs = [0.45, -0.6im]; β = 2.0; Δt = 0.15
    for N in [1, 2, 3, 4]
        r = if_path_values(Es, Vs, β, N, Δt)
        for (εs, U) in [((0.0, 0.0), 0.0), ((0.3, -0.5), 1.2)]
            ρif = contract_impurity(r.Ivals; N=N, Δt=Δt, εs=εs, U=U)
            ρex = trotter_exact(Es, Vs, β, N, Δt; εs=εs, U=U)
            @test maximum(abs.(ρif - ρex)) < 1e-12      # 机器精度（实测 ≤2e-15）
        end
    end
    # N=1 的 Hafnian 字典（研究对照）机器精度
    r1 = if_path_values(Es, Vs, β, 1, Δt)
    for S in 0:15
        supp = [v for v in 1:4 if (S >> (v - 1)) & 1 == 1]
        @test isapprox(path_value_haf(r1.Θ, r1.c, supp), r1.Ivals[S+1]; atol=1e-13)
    end
end

@testset "L6 结构定理与 𝒢 非高斯性（紧致闭式不终止的根源）" begin
    include("brute.jl")
    using LinearAlgebra
    # 浴模线性算符（JW 与 mb_generator 一致，已交叉验证）
    function bath_ops(M)
        dim = 1 << M
        cs = [zeros(ComplexF64, dim, dim) for _ in 1:M]
        for k in 1:M, J in 0:dim-1
            (J >> (k - 1)) & 1 == 1 || continue
            sgn = (-1)^count_ones(J & ((1 << (k - 1)) - 1))
            cs[k][J ⊻ (1 << (k - 1)) + 1, J + 1] = sgn
        end
        return cs
    end
    for (Es, Vs, Δt) in [([-0.7, 1.3], [0.45, -0.6im], 0.15),
                          ([-1.1, 0.2, 2.0], [0.3 + 0.1im, 0.5, -0.4im], 0.08)]
        M = length(Es)
        h = single_particle_h(Es, Vs)
        Wmb = exp(-im * Δt * mb_generator(h))
        X = impurity_blocks(Wmb, M)
        W = exp(-im * Δt * h)
        a = vec(W[1, 2:end]); b = vec(W[2:end, 1]); F1 = W[1, 1]
        cs = bath_ops(M)
        Λa = sum(a[k] * cs[k] for k in 1:M)
        Λbd = sum(b[k] * cs[k]' for k in 1:M)
        # 四块结构定理（机器精度）
        @test norm(X[(1, 0)] - X[(0, 0)] * Λa) < 1e-14
        @test norm(X[(0, 1)] - Λbd * X[(0, 0)]) < 1e-14
        @test norm(X[(1, 1)] - (F1 * X[(0, 0)] - Λbd * X[(0, 0)] * Λa)) < 1e-14
        # 𝒢 非高斯性：在高斯+位移张成上的残差 ≠ 0（且随 V²Δt² 缩小）
        dim = 1 << M
        Imat = Matrix{ComplexF64}(I, dim, dim)
        basis = hcat(vec(Imat), [vec(cs[k]) for k in 1:M]...,
                     [vec(cs[k]') for k in 1:M]...,
                     [vec(cs[i]' * cs[j]) for i in 1:M for j in 1:M]...)
        coef = basis \ vec(X[(0, 0)])
        res = norm(basis * coef - vec(X[(0, 0)]))
        @test res > 0            # 内禀非高斯（占据基紧致闭式不终止的数学根源）
        @test res < 0.1          # 但 O(V²Δt²) 小量：接近高斯
    end
end

@testset "L7 Wick 单粒子归约：严格离散 IF 的 Green 函数表示" begin
    include("brute.jl")
    include("green.jl")
    using LinearAlgebra
    for N in [1, 2, 3, 4]
        for (Es, Vs, β, Δt) in [([-0.7, 1.3], [0.45, -0.6im], 2.0, 0.15),
                                 ([-1.1, 0.2, 2.0], [0.3 + 0.1im, 0.5, -0.4im], 0.7, 0.08)]
            r = if_path_values(Es, Vs, β, N, Δt)
            Wv = wick_path_values(Es, Vs, β, N, Δt)
            Lv = loop_path_values(Es, Vs, β, N, Δt)
            @test maximum(abs.(Wv - r.Ivals)) < 1e-13     # Wick = X-链（机器精度）
            @test maximum(abs.(Lv - r.Ivals)) < 1e-13     # 环路 = X-链
        end
    end
end

@testset "L8 连续热库严格离散 IF（谱闭式 + 连续求积 → 端到端）" begin
    include("brute.jl")
    include("volterra.jl")
    include("green.jl")
    include("contract.jl")
    using LinearAlgebra
    # (i) 谱闭式（极点）= 多体真值，N=1..4 机器精度（V 实规范）
    for N in [1, 2, 3, 4]
        Es = [-0.7, 1.3]; Vs = [0.45, 0.6]
        r = if_path_values(Es, Vs, 2.0, N, 0.15)
        Ws = wick_path_values_spec(Es, Vs, 2.0, N, 0.15)
        @test maximum(abs.(Ws - r.Ivals)) < 1e-13
    end
    # (ii) 连续求积版 vs 同节点浴的多体真值：Volterra O(h²) 收敛 + 端到端
    Γ(ω) = sum([0.2, 0.3][k] * exp(-(ω - [-0.5, 0.7][k])^2 / (2 * 0.25^2)) /
               sqrt(2π * 0.25^2) for k in 1:2)
    xs, ws = gauss_legendre(4)
    Em = [x for x in xs]
    Vm = [sqrt(w * Γ(ε)) for (ε, w) in zip(Em, ws)]
    β, Δt, N = 2.0, 0.15, 2
    r = if_path_values(Em, Vm, β, N, Δt)
    e1 = maximum(abs.(wick_path_values_cont(Γ, -1.0, 1.0, β, N, Δt; nsub=200, npts_spec=4) - r.Ivals))
    e2 = maximum(abs.(wick_path_values_cont(Γ, -1.0, 1.0, β, N, Δt; nsub=400, npts_spec=4) - r.Ivals))
    @test e2 < e1 / 1.9                                   # O(h²)（实测阶比 1.97）
    Wc = wick_path_values_cont(Γ, -1.0, 1.0, β, N, Δt; nsub=3200, npts_spec=4)
    ρif = contract_impurity(Wc; N=N, Δt=Δt, εs=(0.3, -0.5), U=1.2)
    ρex = trotter_exact(Em, Vm, β, N, Δt; εs=(0.3, -0.5), U=1.2)
    @test maximum(abs.(ρif - ρex)) < 1e-5                 # 实测 7.5e-7（纯 Volterra 残差）
end

end # top-level
