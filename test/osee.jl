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

# ---------------- 稠密 IF 张量 vs 快速 OSEE（GTEMPO 默认 ordering） ----------------
@testset "稠密 IF 张量 vs 快速 OSEE" begin
    # 前缀割 SVD 熵（plain SVD：连续前缀割下与 Z₂/协方差约定一致）
    svent(M) = (sv = svdvals(M); p = abs2.(sv); p ./= sum(p);
                -sum(x * log(x) for x in p if x > 1e-15))
    cases = (
        ("IM", fermionic_Δτ(semicircular(), β=5.0, N=4, μ=0.3)),
        ("RT", fermionic_Δt(semicircular(), β=2.0, t=1.0, N=2, μ=0.2)),
        ("MX", fermionic_Δm(semicircular(), β=2.0, Nτ=2, t=0.5, Nt=1, μ=0.2)),
    )
    for (name, corr) in cases
        T = dense_fermionicIF(corr)
        n = size(QuAPI._contour_delta(corr), 1)
        @test T isa DenseFermionicIF
        @test length(T.basis) == 2n && length(T.data) == 2^(2n)
        @test norm(T.data) > 0
        # basis：链序节点 × (a, ā)，band=1（BCS 多 band 场景待扩展）
        b1, j1 = QuAPI._if_legorder(corr)[1]
        @test T.basis[1:2] == [FermionicLeg(j1, b1, 1, false),
                               FermionicLeg(j1, b1, 1, true)]
        @test all(l -> l.band == 1, T.basis)
        s = FermionicGaussianState(fermionic_bdg_matrix(corr))
        chain = QuAPI._if_legorder(corr)             # GTEMPO ordering 链序节点
        nodes = QuAPI._contour_nodes(corr)           # Δ 行序节点（fast 模式指标）
        v = T.data
        for k in 1:n-1                                # 链前缀 k 个节点块
            modes = [findfirst(==(bi), nodes) for bi in chain[1:k]]
            S_dense = svent(reshape(v, 2^(2k), 2^(2(n - k))))
            @test osee(s, vcat(modes, n .+ modes)) ≈ S_dense atol = 1e-10
        end
        # 整割 = 0（纯系数向量）
        @test svent(reshape(v, 2^(2n), 1)) ≈ 0.0 atol = 1e-12
    end
    # GTEMPO 默认 ordering 链序（三种 contour）
    @test QuAPI._if_legorder(fermionic_Δτ(semicircular(), β=5.0, N=3, μ=0.0)) ==
          [(:τ, 3), (:τ, 2), (:τ, 1)]                 # 虚时降序
    @test QuAPI._if_legorder(fermionic_Δt(semicircular(), β=2.0, t=1.0, N=2, μ=0.0)) ==
          [(:+, 3), (:-, 3), (:+, 2), (:-, 2), (:+, 1), (:-, 1)]  # 时间降序，+ 先 − 后
    @test QuAPI._if_legorder(fermionic_Δm(semicircular(), β=2.0, Nτ=2, t=0.5, Nt=1, μ=0.0)) ==
          [(:τ, 2), (:τ, 1), (:+, 1), (:-, 1), (:+, 2), (:-, 2)]  # 虚时降序 ++ 实时升序
    # 阶数上限：n = 16 ⟹ L = 32 > 30 提前报错（不进入 2^32 系数构造）
    big = fermionic_Δτ(semicircular(), β=5.0, N=16, μ=0.0)
    @test_throws ArgumentError dense_fermionicIF(big)
end

# ---------------- 稠密玻色 IF（路径张量，hyb·(a†+a) 耦合） ----------------
@testset "稠密玻色 IF 路径张量" begin
    λsbm = eigvals(Symmetric([1.0 0.0; 0.0 -1.0]))   # σ_z 本征值（自旋-玻色）
    # 逐元素对照公式 I[v] = exp(Σ_pq Δ_pq v_p v_q)，v 取 λ（链序）
    # 张量维度 p（链序）的指标 I[p] 选路径值 λ[I[p]]；Δ 按行序节点索引
    function ref_IF(corr, λ)
        Δ = QuAPI._contour_delta(corr)
        nodes = QuAPI._contour_nodes(corr)
        chain = QuAPI._if_legorder(corr)
        d, n = length(λ), length(chain)
        rowof = [findfirst(==(bi), nodes) for bi in chain]   # 链位置 → Δ 行号
        T = ones(ComplexF64, ntuple(_ -> d, n))
        for I in CartesianIndices(T)
            v = zeros(n)
            for p in 1:n
                v[rowof[p]] = λ[I[p]]
            end
            T[I] = exp(dot(v, Δ, v))
        end
        return T
    end
    # 自旋-玻色（σ_z 本征值）：三种 contour 全体元素对照
    f0b = spectrum(ε -> exp(-ε), lb=0.5, ub=8.0)
    for (name, corr) in (("IM", bosonic_Δτ(f0b, β=2.0, N=4)),
                         ("RT", bosonic_Δt(f0b, β=2.0, t=1.0, N=2)),
                         ("MX", bosonic_Δm(f0b, β=2.0, Nτ=2, t=0.5, Nt=1)))
        T = dense_bosonicIF(corr, λsbm)
        @test T isa DenseBosonicIF
        @test T.d == 2 && length(T.data) == 2^length(T.basis)
        # basis：链序节点，一条路径腿
        b1, j1 = QuAPI._if_legorder(corr)[1]
        @test T.basis[1] == BosonicLeg(j1, b1)
        @test T.data ≈ vec(ref_IF(corr, λsbm)) rtol = 1e-12
    end
    # 一般本征值矢量（非 ±1；含零）
    λg = [0.3, 1.7, -0.9]
    T4 = dense_bosonicIF(bosonic_Δτ(f0b, β=2.0, N=3), λg)
    @test T4.data ≈ vec(ref_IF(bosonic_Δτ(f0b, β=2.0, N=3), λg)) rtol = 1e-12
    @test T4.d == 3 && length(T4.data) == 27 && length(T4.basis) == 3
    # 负检验：链序（GTEMPO ordering）确实影响结果——升序链的参考与实现不同，
    # 保证上面的逐元素对照不是平凡通过。用实时 corr（虚时 Δ 有 KMS 时间反转
    # 对称 Δ = JΔJ，升降序恰好重合，不能用于此检验）
    corr3 = bosonic_Δt(f0b, β=2.0, t=1.0, N=2)
    Tdown = dense_bosonicIF(corr3, λsbm)
    function ref_IF_wrong(corr, λ)   # 故意用升序链
        Δ = QuAPI._contour_delta(corr)
        nodes = QuAPI._contour_nodes(corr)
        chain = reverse(QuAPI._if_legorder(corr))
        d, n = length(λ), length(chain)
        rowof = [findfirst(==(bi), nodes) for bi in chain]
        T = ones(ComplexF64, ntuple(_ -> d, n))
        for I in CartesianIndices(T)
            v = zeros(n)
            for p in 1:n
                v[rowof[p]] = λ[I[p]]
            end
            T[I] = exp(dot(v, Δ, v))
        end
        return T
    end
    @test !(Tdown.data ≈ vec(ref_IF_wrong(corr3, λsbm)))
    # 空 λ 报错
    @test_throws ArgumentError dense_bosonicIF(corr3, Float64[])
    # 元素数上限：d=2, n=31 ⟹ 2^31 > 2^30
    big = bosonic_Δτ(f0b, β=2.0, N=31)
    @test_throws ArgumentError dense_bosonicIF(big, λsbm)
    # d=4, n=16 ⟹ 4^16 = 2^32 > 2^30
    big16 = bosonic_Δτ(f0b, β=2.0, N=16)
    @test_throws ArgumentError dense_bosonicIF(big16, [1.0, 2.0, 3.0, 4.0])
end
