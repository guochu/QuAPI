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
        T = GrassmannIF(GrassmannHybridizationMatrix(corr))
        n = size(QuAPI.Dense._contour_delta(corr), 1)
        chain = QuAPI.Dense._if_legorder(corr)       # GTEMPO ordering 链序节点
        nodes = QuAPI.Dense._contour_nodes(corr)     # Δ 行序节点（fast 模式指标）
        @test T isa GrassmannADT
        @test length(T.basis) == 2n && length(T.data) == 2^(2n)
        @test norm(T.data) > 0
        # 元素类型由 correlation function 决定（不写死）：虚时为实、实时含虚部
        @test eltype(T.data) == eltype(QuAPI.Dense._contour_delta(corr))
        @test name == "IM" ? eltype(T.data) <: Real : eltype(T.data) <: Complex
        # basis 逐位置固化：非混合为链序节点 × (a, ā)；混合的实时支按
        # band、branch（− 先 + 后）、(a, ā)（与 GTEMPO 默认 ordering 一致）
        expected = GrassmannContourIndex[]
        if name == "MX"
            for j in isize(corr):-1:1, c in (false, true)
                push!(expected, GrassmannContourIndex(j, :τ, 1, c))
            end
            for j in 1:rsize(corr), b in (:-, :+), c in (false, true)
                push!(expected, GrassmannContourIndex(j, b, 1, c))
            end
        else
            append!(expected, [GrassmannContourIndex(j, b, 1, c) for (b, j) in chain
                                                         for c in (false, true)])
        end
        @test T.basis == expected
        s = FermionicGaussianState(fermionic_bdg_matrix(corr))
        v = T.data
        for k in 1:n-1                                # 链前缀 k 个节点块
            pref = unique([(bi.branch, bi.j) for bi in T.basis[1:2k]])  # 前 k 个节点
            @test length(pref) == k
            modes = [findfirst(==(bi), nodes) for bi in pref]
            S_dense = svent(reshape(v, 2^(2k), 2^(2(n - k))))
            @test osee(s, vcat(modes, n .+ modes)) ≈ S_dense atol = 1e-10
        end
        # 整割 = 0（纯系数向量）
        @test svent(reshape(v, 2^(2n), 1)) ≈ 0.0 atol = 1e-12
    end
    # GTEMPO 默认 ordering 链序（三种 contour）
    @test QuAPI.Dense._if_legorder(fermionic_Δτ(semicircular(), β=5.0, N=3, μ=0.0)) ==
          [(:τ, 3), (:τ, 2), (:τ, 1)]                 # 虚时降序
    @test QuAPI.Dense._if_legorder(fermionic_Δt(semicircular(), β=2.0, t=1.0, N=2, μ=0.0)) ==
          [(:+, 3), (:-, 3), (:+, 2), (:-, 2), (:+, 1), (:-, 1)]  # 时间降序，+ 先 − 后
    @test QuAPI.Dense._if_legorder(fermionic_Δm(semicircular(), β=2.0, Nτ=2, t=0.5, Nt=1, μ=0.0)) ==
          [(:τ, 2), (:τ, 1), (:+, 1), (:-, 1), (:+, 2), (:-, 2)]  # 虚时降序 ++ 实时升序
    # 阶数上限：n = 16 ⟹ L = 32 > 30 提前报错（不进入 2^32 系数构造）
    big = fermionic_Δτ(semicircular(), β=5.0, N=16, μ=0.0)
    @test_throws ArgumentError GrassmannIF(GrassmannHybridizationMatrix(big))
end

# ---------------- 稠密费米 IF：BCS bath（2 band） ----------------
@testset "稠密费米 IF（BCS：2 band）" begin
    β, N = 2.0, 3
    spec = DiscreteSpectrum([1.0], [0.25])
    corr = fermionic_Δτ(spec, β=β, N=N, μ=0.0)
    z = corr - corr                                  # 零块
    bcs = BCSCorrelationFunction(; uu=corr, ud=z, du=z, dd=z)
    B = GrassmannIF(GrassmannHybridizationMatrix(bcs))
    n = N
    @test B isa GrassmannADT
    # 每节点 4 条腿：band 1/2，每个 band 内 (a, ā) 相邻，链序 τ 降序
    @test length(B.basis) == 4n
    @test length(B.data) == 2^(4n)
    @test B.basis == GrassmannContourIndex[
        GrassmannContourIndex(j, b, band, c)
        for (b, j) in QuAPI.Dense._if_legorder(corr) for band in (1, 2)
        for c in (false, true)]
    # 只有 ↑ 的正常块非零：↓ 腿全偶的系数应与单 flavor 张量逐一相等
    A = GrassmannIF(GrassmannHybridizationMatrix(corr))   # 单 flavor
    @test length(A.basis) == 2n
    for i in 0:(2^(2n) - 1)
        mask_b = 0
        for (k, l) in enumerate(A.basis)
            isodd(i >> (k - 1)) || continue
            lb = GrassmannContourIndex(l.j, l.branch, 1, l.conj)
            mask_b |= 1 << (findfirst(==(lb), B.basis) - 1)
        end
        @test B.data[mask_b + 1] ≈ A.data[i + 1] atol = 1e-14
    end
    # 规模上限：BCS 每节点 4 条腿 ⟹ n = 8 即 L = 32 > 30
    bigc = fermionic_Δτ(spec, β=β, N=8, μ=0.0)
    bigbcs = BCSCorrelationFunction(; uu=bigc, ud=bigc - bigc, du=bigc - bigc,
                                    dd=bigc - bigc)
    @test_throws ArgumentError GrassmannIF(GrassmannHybridizationMatrix(bigbcs))
end

# ---------------- AbstractContourIndex：腿指标上的 correlation function 访问 ----------------
@testset "AbstractContourIndex 访问" begin
    corrτ = fermionic_Δτ(semicircular(), β=5.0, N=4, μ=0.3)
    i = GrassmannContourIndex(2, :τ, 1, false)
    j = GrassmannContourIndex(3, :τ, 1, true)
    @test branch(i) == :τ && branch(j) == :τ
    @test i isa AbstractContourIndex && BosonicContourIndex(1, :+) isa AbstractContourIndex
    # a-ā 取 Δ，ā-a 取 -Δ，a-a / ā-ā 为零
    @test corrτ[i, j] == index(corrτ, 2, 3; b1=:τ, b2=:τ)
    @test corrτ[j, i] == -index(corrτ, 3, 2; b1=:τ, b2=:τ)
    @test iszero(corrτ[i, GrassmannContourIndex(3, :τ, 1, false)])
    @test iszero(corrτ[j, GrassmannContourIndex(3, :τ, 1, true)])
    # 普通费米 bath：不再要求两个 band 相等
    @test corrτ[i, GrassmannContourIndex(3, :τ, 2, true)] ==
          index(corrτ, 2, 3; b1=:τ, b2=:τ)
    # 实时 contour：分支由指标携带
    corrt = fermionic_Δt(semicircular(), β=2.0, t=1.0, N=2, μ=0.2)
    ip = GrassmannContourIndex(2, :+, 1, false)
    jm = GrassmannContourIndex(3, :-, 1, true)
    @test corrt[ip, jm] == index(corrt, 2, 3; b1=:+, b2=:-)
    @test corrt[jm, ip] == -index(corrt, 3, 2; b1=:-, b2=:+)
    # 玻色（无 band）
    corrb = bosonic_Δτ(spectrum(ε -> exp(-ε), lb=0.5, ub=8.0), β=2.0, N=3)
    bi = BosonicContourIndex(1, :τ)
    bj = BosonicContourIndex(2, :τ)
    @test corrb[bi, bj] == index(corrb, 1, 2; b1=:τ, b2=:τ)
    @test corrb[bj, bi] == index(corrb, 2, 1; b1=:τ, b2=:τ)
    # BCS：band 奇数 = 上自旋，偶数 = 下自旋（允许带间块）
    bcs = bcs_Δτ(DiscreteSpectrum([1.0], [0.25]), β=2.0, N=3, μ=0.0, Δ=0.5)
    kup = GrassmannContourIndex(1, :τ, 1, false)
    kdn = GrassmannContourIndex(2, :τ, 2, true)
    @test bcs[kup, kdn] == index(bcs[:↑, :↓], 1, 2; b1=:τ, b2=:τ)
    @test bcs[kdn, kup] == index(bcs[:↓, :↑], 2, 1; b1=:τ, b2=:τ)
    @test bcs[kup, GrassmannContourIndex(3, :τ, 3, true)] ==
          index(bcs[:↑, :↑], 1, 3; b1=:τ, b2=:τ)     # band 3 = 上自旋
end

# ---------------- 稠密 IF 乘法（两个 basis 可以不同） ----------------
@testset "稠密 IF 乘法" begin
    # Grassmann：两腿因子 exp(c ξᵢξⱼ) = 1 + c ξᵢξⱼ 乘进大的 IF
    corr = fermionic_Δτ(semicircular(), β=2.0, N=2)
    gif = GrassmannIF(GrassmannHybridizationMatrix(corr))
    oneif = GrassmannADT(GrassmannContourIndex[], [1.0])
    @test (gif * oneif).data ≈ gif.data
    @test (oneif * gif).data ≈ gif.data
    for (i, j) in ((1, 4), (4, 1), (2, 3))
        c = 0.37
        term = GrassmannADT([gif.basis[i], gif.basis[j]], [1.0, 0.0, 0.0, c])
        prod = gif * term
        @test prod.basis == gif.basis
        # 独立参考：把 c 加进指数矩阵再折叠
        C = QuAPI.Dense._if_coefficients(corr)
        C[i, j] += c
        @test prod.data ≈ QuAPI.Dense._fold_grassmann(C) rtol = 1e-12
    end
    # 玻色 ADT：共享节点上的逐点乘积
    bcorr = bosonic_Δτ(spectrum(ε -> exp(-ε), lb=0.5, ub=8.0), β=2.0, N=3)
    A = BosonicIF(BosonicHybridizationMatrix(bcorr), [0.5, -0.5])
    B = BosonicADT([A.basis[1]], 2, [1.0, 2.0])
    C = A * B
    @test C.basis == A.basis
    @test reshape(C.data, 2, 2, 2) ≈
          reshape(A.data, 2, 2, 2) .* reshape([1.0, 2.0], 2, 1, 1)
    @test Zvalue(A) ≈ sum(A.data)
    @test Zvalue(oneif) ≈ 1.0
    # 玻色 PT：按矩阵乘积合成（行 = 各节点 ket，列 = 各节点 bra）
    P1 = BosonicIF(BosonicHybridizationMatrix(bcorr; split=true), Matrix(Diagonal([0.5, -0.5])))
    P2 = BosonicIF(BosonicHybridizationMatrix(bcorr; split=true), [0.2 0.1; -0.1 0.3])
    n = length(P1.basis)
    perm = vcat(collect(1:2:2n), collect(2:2:2n))
    Mp(Q) = reshape(permutedims(reshape(Q.data, ntuple(_ -> 2, 2n)), perm), 2^n, 2^n)
    Cpt = P1 * P2
    @test Cpt.basis == P1.basis
    @test Mp(Cpt) ≈ Mp(P1) * Mp(P2)
end

# ---------------- 稠密 IF 索引与边缘化 ----------------
@testset "稠密 IF 索引与 margin" begin
    # Grassmann：x[legs...] 是这些腿组成的单项式系数，换位带反对易符号
    corr = fermionic_Δτ(semicircular(), β=2.0, N=2)
    gif = GrassmannIF(GrassmannHybridizationMatrix(corr))
    b = gif.basis
    L = length(b)
    @test gif[] ≈ gif.data[1]                     # 常数项
    @test gif[b[1]] ≈ gif.data[2]                 # ξ₁
    @test gif[b[1], b[2]] ≈ gif.data[4]           # ξ₁ξ₂
    @test gif[b[2], b[1]] ≈ -gif.data[4]          # 反对易
    @test gif[b[1], b[1]] == 0                    # ξ² = 0
    # Grassmann margin：对未列出的腿的占据状态求和
    @test margin(gif).data[1] ≈ sum(gif.data)
    M = margin(gif, b[1], b[3])
    @test M.basis == [b[1], b[3]]
    @test sum(M.data) ≈ sum(gif.data)
    @test M.data[1] ≈ sum(gif.data[m + 1] for m in 0:(2^L - 1)
                           if (m & 1) == 0 && (m & 4) == 0)
    @test M.data[end] ≈ sum(gif.data[m + 1] for m in 0:(2^L - 1)
                             if (m & 1) == 1 && (m & 4) == 4)

    # 玻色 ADT：未列出的节点求和
    bcorr = bosonic_Δτ(spectrum(ε -> exp(-ε), lb=0.5, ub=8.0), β=2.0, N=3)
    A = BosonicIF(BosonicHybridizationMatrix(bcorr), [0.5, -0.5])
    mA = margin(A, A.basis[2])
    @test mA.basis == [A.basis[2]]
    @test mA.data ≈ vec(sum(reshape(A.data, 2, 2, 2), dims=(1, 3)))
    @test margin(A).data ≈ [Zvalue(A)]
end

# ---------------- 稠密玻色 ADT（对角耦合：路径权重） ----------------
@testset "稠密玻色 ADT 路径权重" begin
    λsbm = eigvals(Symmetric([1.0 0.0; 0.0 -1.0]))   # σ_z 本征值（自旋-玻色）
    # 逐元素对照公式 I[v] = exp(Σ_pq Δ_pq v_p v_q)，v 取 λ（链序）
    # 张量维度 p（链序）的指标 I[p] 选路径值 λ[I[p]]；Δ 按行序节点索引
    function ref_IF(corr, λ)
        Δ = QuAPI.Dense._contour_delta(corr)
        nodes = QuAPI.Dense._contour_nodes(corr)
        chain = QuAPI.Dense._if_legorder(corr)
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
        A = BosonicIF(BosonicHybridizationMatrix(corr), λsbm)
        n = length(QuAPI.Dense._if_legorder(corr))
        chain = QuAPI.Dense._if_legorder(corr)
        @test A isa BosonicADT
        @test A.d == 2 && length(A.data) == 2^n
        @test eltype(A.data) == eltype(QuAPI.Dense._contour_delta(corr))
        # basis 逐位置固化
        @test A.basis == BosonicContourIndex[BosonicContourIndex(j, b) for (b, j) in chain]
        @test A.data ≈ vec(ref_IF(corr, λsbm)) rtol = 1e-12
    end
    # 一般本征值矢量（非 ±1；含零）
    λg = [0.3, 1.7, -0.9]
    A4 = BosonicIF(BosonicHybridizationMatrix(bosonic_Δτ(f0b, β=2.0, N=3)), λg)
    @test A4.data ≈ vec(ref_IF(bosonic_Δτ(f0b, β=2.0, N=3), λg)) rtol = 1e-12
    @test A4.d == 3 && length(A4.data) == 27 && length(A4.basis) == 3
    # 负检验：链序（GTEMPO ordering）确实影响结果——升序链的参考与实现不同，
    # 保证上面的逐元素对照不是平凡通过。用实时 corr（虚时 Δ 有 KMS 时间反转
    # 对称 Δ = JΔJ，升降序恰好重合，不能用于此检验）
    corr3 = bosonic_Δt(f0b, β=2.0, t=1.0, N=2)
    Tdown = BosonicIF(BosonicHybridizationMatrix(corr3), λsbm)
    function ref_IF_wrong(corr, λ)   # 故意用升序链
        Δ = QuAPI.Dense._contour_delta(corr)
        nodes = QuAPI.Dense._contour_nodes(corr)
        chain = reverse(QuAPI.Dense._if_legorder(corr))
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
    @test_throws ArgumentError BosonicIF(BosonicHybridizationMatrix(corr3), Float64[])
    # 元素数上限：d=2, n=31 ⟹ 2^31 > 2^30
    big = bosonic_Δτ(f0b, β=2.0, N=31)
    @test_throws ArgumentError BosonicIF(BosonicHybridizationMatrix(big), λsbm)
    # d=4, n=16 ⟹ 4^16 = 2^32 > 2^30
    big16 = bosonic_Δτ(f0b, β=2.0, N=16)
    @test_throws ArgumentError BosonicIF(BosonicHybridizationMatrix(big16), [1.0, 2.0, 3.0, 4.0])
end

# ---------------- 稠密玻色 PT（非对角耦合：IF = exp(H)） ----------------
@testset "稠密玻色 PT 矩阵元" begin
    f0b = spectrum(ε -> exp(-ε), lb=0.5, ub=8.0)
    λ = [1.0, -1.0]
    # 对角 hyb 极限：PT 退化为 ADT 权重的对角嵌入（exp(H) 在路径基下对角）
    for corr in (bosonic_Δτ(f0b, β=2.0, N=4),
                 bosonic_Δt(f0b, β=2.0, t=1.0, N=2),
                 bosonic_Δm(f0b, β=2.0, Nτ=2, t=0.5, Nt=1))
        n = length(QuAPI.Dense._if_legorder(corr))
        d = 2
        A = BosonicIF(BosonicHybridizationMatrix(corr), λ)
        P = BosonicIF(BosonicHybridizationMatrix(corr; split=true), Matrix(Diagonal(λ)))
        @test P isa BosonicPT
        @test P.d == d && length(P.data) == d^(2n) == 4^n
        @test P.basis == A.basis                       # 同一链序节点
        @test eltype(P.data) == eltype(A.data)
        T = reshape(P.data, ntuple(_ -> d, 2n))
        # 对角：(α₁,α₁,α₂,α₂,…) = ADT 权重
        W = Array{eltype(P.data)}(undef, ntuple(_ -> d, n))
        for I in CartesianIndices(W)
            W[I] = T[Tuple(vcat([[I[p], I[p]] for p in 1:n]...))...]
        end
        @test vec(W) ≈ A.data rtol = 1e-12
        # 非对角元为 0（纯对角嵌入）
        offdiag = maximum(abs, T[I] for I in CartesianIndices(T)
                          if any(I[2p - 1] != I[2p] for p in 1:n))
        @test offdiag < 1e-12
    end
    # 非对角（实对称）hyb：PT = 在本征基 ADT 上逐节点套正交变换 U
    # hyb = U·diag(λ)·Uᵀ，exp(H) = (⊗U)·exp(H_diag)·(⊗Uᵀ)
    corr = bosonic_Δτ(f0b, β=2.0, N=4)
    n = length(QuAPI.Dense._if_legorder(corr))
    d = 2
    θ = 0.37
    U = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    λnd = [0.7, -0.3]
    hybU = U * Diagonal(λnd) * U'
    P = BosonicIF(BosonicHybridizationMatrix(corr; split=true), hybU)
    w = reshape(BosonicIF(BosonicHybridizationMatrix(corr), λnd).data, ntuple(_ -> d, n))
    T = reshape(P.data, ntuple(_ -> d, 2n))
    Texp = Array{eltype(P.data)}(undef, ntuple(_ -> d, 2n))
    for I in CartesianIndices(Texp)
        α = [I[2p - 1] for p in 1:n]
        β = [I[2p] for p in 1:n]
        s = zero(eltype(Texp))
        for g in CartesianIndices(w)
            s += prod(U[α[p], g[p]] * U[β[p], g[p]] for p in 1:n) * w[g]
        end
        Texp[I] = s
    end
    @test T ≈ Texp atol = 1e-11
    # 复厄米 hyb：H 厄米 ⟹ exp(H) 厄米 ⟹ T[α;β] = conj(T[β;α])
    hybH = [0.5 0.4im; -0.4im -0.2]
    Pc = BosonicIF(BosonicHybridizationMatrix(corr; split=true), hybH)
    Tc = reshape(Pc.data, ntuple(_ -> d, 2n))
    swap = vcat([[2p, 2p - 1] for p in 1:n]...)
    @test Tc ≈ conj.(permutedims(Tc, swap)) atol = 1e-11
    # 非方阵 hyb / 规模上限
    @test_throws ArgumentError BosonicIF(BosonicHybridizationMatrix(corr; split=true), [1.0 2.0 3.0])
    big = bosonic_Δτ(f0b, β=2.0, N=16)          # d=2, n=16 ⟹ 2^32
    @test_throws ArgumentError BosonicIF(BosonicHybridizationMatrix(big; split=true), Matrix(Diagonal(λ)))
end
