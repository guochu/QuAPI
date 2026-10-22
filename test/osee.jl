println("------------------------------------")
println("|   GaussianStates extension: OSEE  |")
println("------------------------------------")

@testset "GaussianStates extension" begin
    # ---------------- contour Δ 组装（经 bdg 重载间接验证） ----------------
    corrτ = fermionic_Δτ(semicircular(), β=5.0, N=8, μ=0.3)
    Δτ_ref = [index(corrτ, i, j; b1=:τ, b2=:τ) for i in 1:8, j in 1:8]

    corrt = fermionic_Δt(semicircular(), β=2.0, t=1.0, N=3, μ=0.2)
    Δt_ref = Matrix{ComplexF64}(undef, 8, 8)
    for (i, b1) in enumerate((:+, :-)), (j, b2) in enumerate((:+, :-))
        Δt_ref[(i-1)*4+1:i*4, (j-1)*4+1:j*4] =
            [index(corrt, i2, j2; b1=b1, b2=b2) for i2 in 1:4, j2 in 1:4]
    end

    corrm = fermionic_Δm(semicircular(), β=2.0, Nτ=4, t=1.0, Nt=3, μ=0.2)
    Δm_ref = Matrix{ComplexF64}(undef, 12, 12)
    for (i, b1) in enumerate((:+, :-, :τ)), (j, b2) in enumerate((:+, :-, :τ))
        ni = (b1 == :τ) ? 4 : 4
        nj = (b2 == :τ) ? 4 : 4
        oi = (i - 1) * 4
        oj = (j - 1) * 4
        Δm_ref[oi+1:oi+ni, oj+1:oj+nj] =
            [index(corrm, i2, j2; b1=b1, b2=b2) for i2 in 1:ni, j2 in 1:nj]
    end

    # ---------------- 费米：fermionic_bdg_matrix(corr) → 系数态 OSEE ----------
    for (corr, Δref) in ((corrτ, Δτ_ref), (corrt, Δt_ref), (corrm, Δm_ref))
        W = fermionic_bdg_matrix(corr)
        n = size(Δref, 1)
        @test size(W) == (2n, 2n)
        # 反对称块结构：W[ā_p, a_q] = Δ_pq ⟹ 左下 = Δ、右上 = −Δᵀ
        @test W == [zeros(n, n) -transpose(Δref); Δref zeros(n, n)]
        @test W ≈ -transpose(W)
        s = FermionicGaussianState(W)
        @test !s.isbdg
        Γ = GaussianStates.covariance(s)                # 未导出，限定名（惰性触发）
        @test Γ isa Matrix{Float64}                # 任意向量 Majorana 协方差恒实
        @test maximum(abs, Γ^2 + I) < 1e-10        # 纯态
        @test osee(s, collect(1:2n)) ≈ 0.0 atol = 1e-10   # 无 bipartition
        # 互补对称（纯态）：前缀割 = 余下割
        k = div(n, 2)
        @test osee(s, vcat(1:k, n+1:n+k)) ≈
              osee(s, vcat(k+1:n, n+k+1:2n)) atol = 1e-10
        @test_throws BoundsError osee(s, [0])
        @test_throws ArgumentError osee(s, [1, 1])
    end

    # 跨 site 单配对解析（矩阵级，验证 W → 态 → osee 的约定链）：
    # Δ[ā₂, a₁] = d 只配对 ā₂–a₁ ⟹ site 割 1|2 劈开该配对，两割同熵
    d = 0.7 + 0.2im
    Δ = [0 0; d 0]
    n = 2
    W = zeros(ComplexF64, 2n, 2n); W[n+1:2n, 1:n] .= Δ
    s2 = FermionicGaussianState(W - transpose(W))
    ε = abs(d)^2 / (1 + abs(d)^2)
    h2(p) = -p * log(p) - (1 - p) * log(1 - p)
    @test osee(s2, [1, 3]) ≈ h2(ε) atol = 1e-12
    @test osee(s2, [2, 4]) ≈ h2(ε) atol = 1e-12
    @test osee(s2, [1, 3]; α=2) ≈ -log(ε^2 + (1 - ε)^2) atol = 1e-12

    # ---------------- 玻色：bosonic_bdg_matrix(corr) → 系数态 OSEE ----------
    # 玻色 bath 谱须全正（BE 分布要求 ε > μ），semicircular 含负能量不适用
    f0b = spectrum(ε -> exp(-ε), lb=0.5, ub=8.0)
    corr_b = bosonic_Δτ(f0b, β=2.0, N=4)
    Wb = bosonic_bdg_matrix(corr_b)
    Δb = real.([index(corr_b, i, j; b1=:τ, b2=:τ) for i in 1:4, j in 1:4])
    @test Wb == [zeros(4, 4) exp(-Δb); transpose(exp(-Δb)) zeros(4, 4)]
    @test Wb == transpose(Wb) && Wb isa Matrix{Float64}
    @test maximum(abs, eigvals(Symmetric(Wb))) < 1        # 谱 ⊂ (−1,1)
    sb = BosonicGaussianState(Wb)
    @test !sb.isbdg
    @test_throws ArgumentError entropy(sb)                # 非密度矩阵
    Γb = GaussianStates.covariance(sb)
    Ω4 = GaussianStates._omega(8)                         # 8 模（16×16）的辛形式
    ν = abs.(eigvals(Ω4 * Γb))
    @test maximum(abs.(ν .- 0.5)) < 1e-10                 # 纯态：辛本征值全 ½
    @test osee(sb, collect(1:8)) ≈ 0.0 atol = 1e-10
    # site p ↔ 模式对 [p, n+p]；纯态互补对称
    @test osee(sb, vcat([1, 2], 4 .+ [1, 2])) ≈
          osee(sb, vcat(3:4, 7:8)) atol = 1e-10
    # 复 Δ（实时 contour）报错
    @test_throws ArgumentError bosonic_bdg_matrix(corrt)

    # 玻色 TMSV 解析（矩阵级）：M = m·[[0,1],[1,0]] ⟹ 配对 (c₁,d₂)、(c₂,d₁)
    # 跨 site ⟹ site 割（模式对 [p, n+p]）= 单模热熵（tanh r = m）
    m = 0.45
    M = m .* [0.0 1.0; 1.0 0.0]
    Wm = [zeros(2, 2) M; transpose(M) zeros(2, 2)]
    sm = BosonicGaussianState(Wm)
    n̄(x) = x^2 / (1 - x^2)
    SvN(x) = (n̄(x) + 1) * log(n̄(x) + 1) - n̄(x) * log(n̄(x))
    # site 割 {c_p,d_p} 同时劈开两条配对（c₁–d₂ 与 c₂–d₁）⟹ 2 × 单模热熵
    @test osee(sm, [1, 3]) ≈ 2SvN(m) atol = 1e-12        # site 1 = 模式 {1, 3}
    @test osee(sm, [2, 4]) ≈ 2SvN(m) atol = 1e-12        # 互补割同熵
    @test osee(sm, [1, 2, 3, 4]) ≈ 0.0 atol = 1e-12
    @test osee(sm, [1, 3]; α=2) ≈ 2log(1 + 2n̄(m)) atol = 1e-12
end
