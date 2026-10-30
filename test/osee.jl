println("------------------------------------")
println("|   GaussianStates extension: OSEE  |")
println("------------------------------------")

# 测试内小工具（只依赖导出接口，避免调用未导出内部函数）：
# 2n×2n 辛形式 Ω = ⊕ₖ [[0,1],[−1,0]]
function symplectic_form(n::Int)
    Ω = zeros(2n, 2n)
    for k in 1:n
        Ω[2k-1, 2k] = 1.0
        Ω[2k, 2k-1] = -1.0
    end
    return Ω
end

@testset "GaussianStates extension" begin
    # ---------------- Δ 提取约定（经构造器间接验证） ----------------
    # 新约定（GaussianStates v7：态直接由 corr 构造）：扩展读的 Δ 与稠密 IF
    # （GrassmannHybridizationMatrix，顶点移位）同一 —— 独立性体现在
    # "不调用扩展自身的 Δ 提取"；约定本身的正确性由下方解析锚、
    # 稠密 SVD 对照与 performance/fermionicG4.jl 的 V⁴ 平台检验背书。
    function delta_ref(corr)
        C = Matrix(GrassmannHybridizationMatrix(corr).Δ)
        n = size(C, 1) ÷ 2
        return [C[2i, 2j-1] for i in 1:n, j in 1:n]
    end

    corrτ = fermionic_Δτ(semicircular(), β=5.0, N=8, μ=0.3)
    Δτ_ref = delta_ref(corrτ)

    corrt = fermionic_Δt(semicircular(), β=2.0, t=1.0, N=3, μ=0.2)
    Δt_ref = delta_ref(corrt)

    corrm = fermionic_Δm(semicircular(), β=2.0, Nτ=4, t=1.0, Nt=3, μ=0.2)
    Δm_ref = delta_ref(corrm)

    # ---------------- 费米：FermionicGaussianState(corr) → 系数态 OSEE ----------
    for (corr, Δref) in ((corrτ, Δτ_ref), (corrt, Δt_ref), (corrm, Δm_ref))
        s = FermionicGaussianState(corr)
        n = size(Δref, 1)
        @test nmodes(s) == 2n                        # doubled：模 1..n = ket、n+1..2n = bra
        # 约定：corr 入口 = 腿矩阵 [ā_p, a_q] 位放 Δ 的矩阵级入口（与旧
        # W 块测试同一约定声明：W[ā_p, a_q] = Δ_pq ⟺ C[2p, 2q−1] = Δ_pq）
        C = zeros(eltype(Δref), 2n, 2n)
        for i in 1:n, j in 1:n
            C[2i, 2j-1] = Δref[i, j]
        end
        sref = FermionicGaussianState(GrassmannHybridizationMatrix(basis(corr), C))
        @test covariance(sref) ≈ covariance(s) atol = 1e-13
        Γ = covariance(s)
        @test Γ isa Matrix{Float64}                # 任意向量 Majorana 协方差恒实
        @test Γ ≈ -transpose(Γ)
        @test maximum(abs, Γ^2 + I) < 1e-10        # 纯态
        @test ispure(s)
        @test entropy(s, collect(1:2n)) ≈ 0.0 atol = 1e-10   # 无 bipartition
        # 互补对称（纯态）：前缀割 = 余下割
        k = div(n, 2)
        @test entropy(s, vcat(1:k, n+1:n+k)) ≈
              entropy(s, vcat(k+1:n, n+k+1:2n)) atol = 1e-10
        @test_throws BoundsError entropy(s, [0])
        @test_throws ArgumentError entropy(s, [1, 1])
    end

    # 跨 site 单配对解析（矩阵级，验证 Δ → 态 → entropy 的约定链）：
    # Δ[ā₂, a₁] = d 只配对 ā₂–a₁ ⟹ site 割 1|2 劈开该配对，两割同熵
    d = 0.7 + 0.2im
    b2 = [GrassmannContourIndex(j, :τ, 1, c) for j in 1:2 for c in (false, true)]
    Cd = zeros(ComplexF64, 4, 4)
    Cd[4, 1] = d                                     # C[ā₂, a₁]
    s2 = FermionicGaussianState(GrassmannHybridizationMatrix(b2, Cd))
    ε = abs(d)^2 / (1 + abs(d)^2)
    h2(p) = -p * log(p) - (1 - p) * log(1 - p)
    @test entropy(s2, [1, 3]) ≈ h2(ε) atol = 1e-12
    @test entropy(s2, [2, 4]) ≈ h2(ε) atol = 1e-12
    @test entropy(s2, [1, 3]; α=2) ≈ -log(ε^2 + (1 - ε)^2) atol = 1e-12
    # 稠密 IF 侧同值：I = 1 + d·ξ₄ξ₁，节点前缀割的 plain SVD 熵同 h2(ε)
    #（两条路径（协方差 vs SVD）在同一解析值上会合）
    gif2 = GrassmannIF(GrassmannHybridizationMatrix(b2, Cd))
    sv = svdvals(reshape(gif2.data, 4, 4))
    p = abs2.(sv); p ./= sum(p)
    @test -sum(x * log(x) for x in p if x > 1e-15) ≈ h2(ε) atol = 1e-12

    # ---------------- 玻色：BosonicGaussianState(corr) → 系数态 OSEE ----------
    # 玻色 bath 谱须全正（BE 分布要求 ε > μ），semicircular 含负能量不适用
    f0b = spectrum(ε -> exp(-ε), lb=0.5, ub=8.0)
    corr_b = bosonic_Δτ(f0b, β=2.0, N=4)
    sb = BosonicGaussianState(corr_b)
    # 约定：Δ 与 BosonicHybridizationMatrix（contour_ordering 行序）同源
    Δb = real.([branch(corr_b; b1=:τ, b2=:τ)[i, j] for i in 1:4, j in 1:4])
    hb = BosonicHybridizationMatrix([BosonicContourIndex(j, :τ) for j in 1:4], Δb)
    @test covariance(BosonicGaussianState(hb)) ≈ covariance(sb) atol = 1e-13
    # 显式 BosonicHybridizationMatrix 入口：split=false 与 correlation 入口一致；
    # split=true 的 BosonicSplitDelta(M1, M2) 是等时次序拆分，不是配对振幅，
    # 必须报错而不是静默拆成两个矩阵
    @test covariance(BosonicGaussianState(BosonicHybridizationMatrix(corr_b))) ≈
          covariance(sb) atol = 1e-13
    hsplit = BosonicHybridizationMatrix(corr_b; split=true)
    @test !(hsplit.Δ isa Tuple)                    # split 结果为专用类型，非裸元组
    @test_throws ArgumentError BosonicGaussianState(hsplit)
    @test_throws MethodError BosonicGaussianState(hsplit.Δ...)   # 不可 splat
    Γb = covariance(sb)
    @test Γb isa Matrix{Float64} && Γb ≈ transpose(Γb)
    Ω4 = symplectic_form(8)                         # 8 模（16×16）的辛形式
    ν = abs.(eigvals(Ω4 * Γb))
    @test maximum(abs.(ν .- 0.5)) < 1e-10           # 纯态：辛本征值全 ½
    @test ispure(sb)
    @test entropy(sb) ≈ 0.0 atol = 1e-10            # 纯态整态熵为 0（v7：entropy
    @test entropy(sb, collect(1:8)) ≈ 0.0 atol = 1e-10  # 对任意 GaussianState 有定义）
    # site p ↔ 模式对 [p, n+p]；纯态互补对称
    @test entropy(sb, vcat([1, 2], 4 .+ [1, 2])) ≈
          entropy(sb, vcat(3:4, 7:8)) atol = 1e-10
    @test_throws BoundsError entropy(sb, [0])
    @test_throws ArgumentError entropy(sb, [1, 1])
    # 复 Δ（实时 contour）报错
    @test_throws ArgumentError BosonicGaussianState(corrt)

    # 玻色解析（矩阵级）：Δ = diag(λ₁,λ₂) ⟹ M = diag(t₁,t₂)（t_k = e^{−λ_k}），
    # 系数态 = 两个独立双模压缩对 (c₁d₁)、(c₂d₂)；单模割劈开一对 ⟹ 单模热熵
    #（site 割 [p, n+p] 此时含完整配对对，熵为 0——交叉配对下的 site 割
    # 解析由费米单配对锚与下方 d→∞ PT 夹逼覆盖）
    λ1, λ2 = 1.1, 0.75
    hd = BosonicHybridizationMatrix([BosonicContourIndex(j, :τ) for j in 1:2],
                                    [λ1 0.0; 0.0 λ2])
    sd = BosonicGaussianState(hd)
    t1, t2 = exp(-λ1), exp(-λ2)
    n̄(x) = x^2 / (1 - x^2)
    SvN(x) = (n̄(x) + 1) * log(n̄(x) + 1) - n̄(x) * log(n̄(x))
    @test entropy(sd, [1]) ≈ SvN(t1) atol = 1e-12
    @test entropy(sd, [2]) ≈ SvN(t2) atol = 1e-12
    @test entropy(sd, [1, 2]) ≈ SvN(t1) + SvN(t2) atol = 1e-12   # ket 半空间
    @test entropy(sd, collect(1:4)) ≈ 0.0 atol = 1e-12
    @test entropy(sd, [1]; α=2) ≈ log((1 + t1^2) / (1 - t1^2)) atol = 1e-12
    # 旋转不变性：对称 Δ' = VΛVᵀ ⟹ M' = V diag(t) Vᵀ 奇异值不变，
    # ket|bra 割的熵只依赖 M 的奇异值（每个奇异值一对跨割双模压缩）
    θ = 0.37
    V = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    hr = BosonicHybridizationMatrix([BosonicContourIndex(j, :τ) for j in 1:2],
                                    V * [λ1 0.0; 0.0 λ2] * transpose(V))
    sr = BosonicGaussianState(hr)
    @test entropy(sr, [1, 2]) ≈ SvN(t1) + SvN(t2) atol = 1e-12
    @test entropy(sr, collect(1:4)) ≈ 0.0 atol = 1e-12
    @test ispure(sr)
end

# ---------------- 玻色 OSEE ↔ BosonicPT IF 的纠缠熵（d→∞ 收敛） ----------------
@testset "玻色 OSEE ↔ BosonicPT IF 熵（d→∞）" begin
    # 单模玻色杂质 + ohmic bath（TEMPO/osee 报告 §5.2、§6 的模型）：
    # N=2、δτ=β/2=0.1 ⟹ Δ 正定、λ_min ≈ 2.1e-5。高斯协方差公式给的是
    # **d→∞ 极限**；有限 d 的 BosonicPT 熵与它是同一条曲线的两端（§6），
    # 因此这里检验 d 扫描穿过高斯值。
    omic(w) = 2 * 0.1 * w * exp(-w / 2)
    specb = spectrum(w -> omic(w), lb=0.0, ub=2.0)
    corr = bosonic_Δτ(specb; β=0.2, N=2)

    # 全部走导出接口：BosonicHybridizationMatrix 的公共字段 basis / Δ
    h = BosonicHybridizationMatrix(corr)               # split=false ⟹ Δ 为矩阵
    Δ = h.Δ
    @test minimum(eigvals(Symmetric(Matrix(Δ)))) > 0   # 实 Δ、正定
    n = length(h.basis)                                # contour 节点数
    site = h.basis[1].j                                # 虚时 Δ 行序 (:τ,1..n)，链首为 n

    sb = BosonicGaussianState(h)
    Γb = covariance(sb)
    Ω = symplectic_form(2n)                            # Γ 为 4n×4n（2n 个 doubled 模）
    @test maximum(abs.(abs.(eigvals(Ω * Γb)) .- 0.5)) < 1e-6     # 纯态（λ_min 小 ⟹ 近奇异）
    S_fast = entropy(sb, [site, n + site])             # d→∞ 高斯 OSEE

    # BosonicPT（PT 耦合 a†/a）：hyb = 截断湮灭算符；data 指标序
    # (k₁,b₁,k₂,b₂)。OSEE 取 vectorized operator 的 ket|bra 二分，
    # 重排为 (k₁,k₂ | b₁,b₂) 后做 SVD。
    aop(nbas) = [k == j - 1 ? sqrt(k) : 0.0 for k in 1:nbas+1, j in 1:nbas+1]
    function pt_osee(d)
        P = BosonicIF(BosonicHybridizationMatrix(corr; split=true), aop(d - 1))
        T = reshape(P.data, d, d, d, d)
        sv = svdvals(reshape(permutedims(T, (1, 3, 2, 4)), d^2, d^2))
        p = abs2.(sv); p ./= sum(p)
        return -sum(x * log(x) for x in p if x > 1e-15)
    end
    @test BosonicIF(BosonicHybridizationMatrix(corr; split=true), aop(2)) isa BosonicPT

    # 有限 d 熵随 d 单调上升（§4.2：小 d 近似 ∝ d²），d≈25 处穿过高斯极限
    @test pt_osee(14) < pt_osee(18) < pt_osee(22) < pt_osee(26)
    @test pt_osee(24) < S_fast < pt_osee(26)               # 极限被相邻 d 夹住
    @test pt_osee(25) ≈ S_fast atol = 0.05                 # 交叉处 ~0.3%
end

# ---------------- 稠密 IF 张量 vs 快速 OSEE（统一 contour ordering） ----------------
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
        n = length(basis(corr)) ÷ 2                    # 单 flavor：每节点 2 条腿
        chain = contour_ordering(corr)                 # 统一节点序（= Δ 行列序）
        @test T isa GrassmannADT
        @test length(T.basis) == 2n && length(T.data) == 2^(2n)
        @test norm(T.data) > 0
        # 元素类型由 correlation function 决定（不写死）：虚时为实、实时含虚部
        @test eltype(T.data) == eltype(GrassmannHybridizationMatrix(corr).Δ)
        @test name == "IM" ? eltype(T.data) <: Real : eltype(T.data) <: Complex
        # basis 逐位置固化：节点按 contour_ordering，每节点 (a, ā)（band 内）
        expected = GrassmannContourIndex[]
        append!(expected, [GrassmannContourIndex(j, b, 1, c) for (b, j) in chain
                                                         for c in (false, true)])
        @test T.basis == expected
        s = FermionicGaussianState(corr)
        v = T.data
        for k in 1:n-1                                # 链前缀 k 个节点块
            pref = unique([(bi.branch, bi.j) for bi in T.basis[1:2k]])  # 前 k 个节点
            @test length(pref) == k
            modes = [findfirst(==(bi), chain) for bi in pref]
            S_dense = svent(reshape(v, 2^(2k), 2^(2(n - k))))
            @test entropy(s, vcat(modes, n .+ modes)) ≈ S_dense atol = 1e-10
        end
        # 整割 = 0（纯系数向量）
        @test svent(reshape(v, 2^(2n), 1)) ≈ 0.0 atol = 1e-12
    end
    # 统一节点序（三种 contour）
    @test contour_ordering(fermionic_Δτ(semicircular(), β=5.0, N=3, μ=0.0)) ==
          [(:τ, 1), (:τ, 2), (:τ, 3)]                 # 虚时升序
    @test contour_ordering(fermionic_Δt(semicircular(), β=2.0, t=1.0, N=2, μ=0.0)) ==
          [(:+, 1), (:-, 1), (:+, 2), (:-, 2), (:+, 3), (:-, 3)]  # 时间升序，+ 先 − 后
    @test contour_ordering(fermionic_Δm(semicircular(), β=2.0, Nτ=2, t=0.5, Nt=1, μ=0.0)) ==
          [(:τ, 1), (:τ, 2), (:+, 1), (:-, 1), (:+, 2), (:-, 2)]  # 虚时升序接实时升序
    # 阶数上限：n = 16 ⟹ L = 32 > 30 提前报错（不进入 2^32 系数构造）
    big = fermionic_Δτ(semicircular(), β=5.0, N=16, μ=0.0)
    @test_throws ArgumentError GrassmannIF(GrassmannHybridizationMatrix(big))
end

# ---------------- 稠密 IF 张量 vs 快速 OSEE（BCS：一般核） ----------------
@testset "稠密 IF 张量 vs 快速 OSEE（BCS 一般核）" begin
    # 一般配对核（含 ξξ/ξ̄ξ̄ 反常块）的系数态：湮灭子投影子一般路线 vs
    # 张量前缀割 SVD 熵（ground truth）。BCS bath：Δ ≠ 0，核四块全非零。
    svent(M) = (sv = svdvals(M); p = abs2.(sv); p ./= sum(p);
                -sum(x * log(x) for x in p if x > 1e-15; init=0.0))
    spec = DiscreteSpectrum([1.0], [0.25])
    corr = bcs_Δτ(spec, β=2.0, N=3, μ=0.0, Δ=0.5)
    T = GrassmannIF(GrassmannHybridizationMatrix(corr))
    L = length(T.basis)                            # 4n = 12 腿（2 band）
    @test T isa GrassmannADT && length(T.data) == 2^L
    s = FermionicGaussianState(corr)
    @test nmodes(s) == L
    # 模式映射：a 腿（c=false）在前半区（基序），ā 腿在后半区
    alp = findall(l -> l.conj == false, T.basis)
    pos_alp = Dict(p => i for (i, p) in enumerate(alp))
    pos_bar = Dict(p => i for (i, p) in enumerate(setdiff(1:L, alp)))
    modeof(p) = haskey(pos_alp, p) ? pos_alp[p] : length(alp) + pos_bar[p]
    for k in 1:2                                   # 前缀割：前 k 个节点（4k 条腿）
        S_dense = svent(reshape(T.data, 2^(4k), 2^(L - 4k)))
        cutmodes = [modeof(p) for p in 1:4k]
        @test entropy(s, cutmodes) ≈ S_dense atol = 1e-10
    end
    @test entropy(s, collect(1:L)) ≈ 0.0 atol = 1e-10   # 全割（纯系数态）
    # 复 Δ 的 BCS：同锚
    corrc = bcs_Δτ(spec, β=2.0, N=3, μ=0.0, Δ=0.2+0.45im)
    Tc = GrassmannIF(GrassmannHybridizationMatrix(corrc))
    sc = FermionicGaussianState(corrc)
    for k in 1:2
        S_dense = svent(reshape(Tc.data, 2^(4k), 2^(L - 4k)))
        cutmodes = [modeof(p) for p in 1:4k]
        @test entropy(sc, cutmodes) ≈ S_dense atol = 1e-10
    end
    # 纯度
    @test ispure(s) && ispure(sc)
end

# ---------------- 稠密 IF 张量 vs 快速 OSEE（BCS：一般核） ----------------
# 待实现：BCS 一般核（含 ξξ/ξ̄ξ̄ 反常块）的矢量化态是一般配对高斯态，
# 需要复 skew-Schur（Youla）分解下的一般配对协方差。锚已定义：BCS 稠密
# IF 张量（GrassmannIF 2-band，N ≤ 3 因 2^{4n} 上限）的前缀割 SVD 熵；
# 实现后按单 flavor 的同款模式对照（张量 SVD vs entropy）。

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
    # 每节点 4 条腿：band 1/2，每个 band 内 (a, ā) 相邻；节点按 contour_ordering（τ 升序）
    @test length(B.basis) == 4n
    @test length(B.data) == 2^(4n)
    @test B.basis == GrassmannContourIndex[
        GrassmannContourIndex(j, b, band, c)
        for (b, j) in contour_ordering(corr) for band in (1, 2)
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
    # 2026-10-30 访问器 τ 环回退 cell（节点值）读取（顶点重组使 g(τ) 误差
    # 常数恶化 ~30×，GTEMPO Toulouse 基准；判决见 contourindex.jl）：
    # i-ā / j-cross 交叉读 = branch[i.j, j.j] 与 −branch[i.j, j.j]（conj 符号）
    @test corrτ[i, j] == branch(corrτ; b1=:τ, b2=:τ)[2, 3]
    @test corrτ[j, i] == -branch(corrτ; b1=:τ, b2=:τ)[3, 2]
    @test iszero(corrτ[i, GrassmannContourIndex(3, :τ, 1, false)])
    @test iszero(corrτ[j, GrassmannContourIndex(3, :τ, 1, true)])
    # 普通费米 bath：不再要求两个 band 相等
    @test corrτ[i, GrassmannContourIndex(3, :τ, 2, true)] ==
          branch(corrτ; b1=:τ, b2=:τ)[2, 3]
    # 实时 contour：分支由指标携带（+− 块值不变；(−,+) 的 (3,2) 越窗为零）
    corrt = fermionic_Δt(semicircular(), β=2.0, t=1.0, N=2, μ=0.2)
    ip = GrassmannContourIndex(2, :+, 1, false)
    jm = GrassmannContourIndex(3, :-, 1, true)
    @test corrt[ip, jm] == branch(corrt; b1=:+, b2=:-)[2, 3]
    @test iszero(corrt[jm, ip])
    # 玻色（无 band）
    corrb = bosonic_Δτ(spectrum(ε -> exp(-ε), lb=0.5, ub=8.0), β=2.0, N=3)
    bi = BosonicContourIndex(1, :τ)
    bj = BosonicContourIndex(2, :τ)
    @test corrb[bi, bj] == branch(corrb; b1=:τ, b2=:τ)[1, 2]
    @test corrb[bj, bi] == branch(corrb; b1=:τ, b2=:τ)[2, 1]
    # BCS：band 奇数 = 上自旋，偶数 = 下自旋（允许带间块；τ 环同移位规则）
    bcs = bcs_Δτ(DiscreteSpectrum([1.0], [0.25]), β=2.0, N=3, μ=0.0, Δ=0.5)
    kup = GrassmannContourIndex(1, :τ, 1, false)
    kdn = GrassmannContourIndex(2, :τ, 2, true)
    @test bcs[kup, kdn] == branch(bcs[:↑, :↓]; b1=:τ, b2=:τ)[1, 2]
    @test bcs[kdn, kup] == -branch(bcs[:↓, :↑]; b1=:τ, b2=:τ)[2, 1]
    @test bcs[kup, GrassmannContourIndex(3, :τ, 3, true)] ==
          branch(bcs[:↑, :↑]; b1=:τ, b2=:τ)[1, 3]      # band 3 = 上自旋；k = 1 环回

    # contour 次序：分支 :+ < :- < :τ；:+ / :τ 内 j 升序，:- 内 j 降序
    @test GrassmannContourIndex(1, :+, 1, false) < GrassmannContourIndex(2, :+, 1, false)
    @test GrassmannContourIndex(2, :+, 1, false) < GrassmannContourIndex(1, :-, 1, false)
    @test GrassmannContourIndex(2, :-, 1, false) < GrassmannContourIndex(1, :-, 1, false)
    @test GrassmannContourIndex(1, :-, 1, false) < GrassmannContourIndex(1, :τ, 1, false)
    @test GrassmannContourIndex(1, :τ, 1, false) < GrassmannContourIndex(2, :τ, 1, false)
    # 同一节点：band 升序，ā 在 a 之前
    @test GrassmannContourIndex(1, :τ, 1, false) < GrassmannContourIndex(1, :τ, 2, false)
    @test GrassmannContourIndex(1, :τ, 1, true) < GrassmannContourIndex(1, :τ, 1, false)
    @test GrassmannContourIndex(1, :τ, 1, false) <= GrassmannContourIndex(1, :τ, 1, false)
    @test GrassmannContourIndex(1, :τ, 1, false) >= GrassmannContourIndex(1, :τ, 1, false)
    @test GrassmannContourIndex(2, :τ, 1, false) > GrassmannContourIndex(1, :τ, 1, false)
    @test GrassmannContourIndex(1, :τ, 1, false) == GrassmannContourIndex(1, :τ, 1, false)
    @test GrassmannContourIndex(1, :τ, 1, false) != GrassmannContourIndex(1, :τ, 2, false)

    @test BosonicContourIndex(1, :+) < BosonicContourIndex(2, :+)
    @test BosonicContourIndex(2, :+) < BosonicContourIndex(1, :-)
    @test BosonicContourIndex(2, :-) < BosonicContourIndex(1, :-)
    @test BosonicContourIndex(1, :-) < BosonicContourIndex(1, :τ)
    @test BosonicContourIndex(1, :τ) < BosonicContourIndex(2, :τ)
    @test BosonicContourIndex(1, :τ) <= BosonicContourIndex(1, :τ)
    @test BosonicContourIndex(1, :τ) == BosonicContourIndex(1, :τ)
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
        # 独立参考：把 c 加进指数矩阵再经公开构造函数折叠
        C = copy(GrassmannHybridizationMatrix(corr).Δ)
        C[i, j] += c
        @test prod.data ≈ GrassmannIF(GrassmannHybridizationMatrix(gif.basis, C)).data rtol = 1e-12
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
    # 玻色腿的补集索引
    @test complementindex(A.basis, A.basis[2]) == [A.basis[1], A.basis[3]]
end

# ---------------- 稠密玻色 ADT（对角耦合：路径权重） ----------------
@testset "稠密玻色 ADT 路径权重" begin
    λsbm = eigvals(Symmetric([1.0 0.0; 0.0 -1.0]))   # σ_z 本征值（自旋-玻色）
    # 逐元素对照公式 I[v] = exp(Σ_pq Δ_pq v_p v_q)，v 取 λ（统一节点序）
    # 张量维度 p 的指标 I[p] 选路径值 λ[I[p]]；Δ 的行列序与 contour_ordering 相同
    function ref_IF(corr, λ)
        Δ = BosonicHybridizationMatrix(corr).Δ
        chain = contour_ordering(corr)
        d, n = length(λ), length(chain)
        T = ones(ComplexF64, ntuple(_ -> d, n))
        for I in CartesianIndices(T)
            v = [λ[I[p]] for p in 1:n]
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
        chain = contour_ordering(corr)
        n = length(chain)
        @test A isa BosonicADT
        @test A.d == 2 && length(A.data) == 2^n
        @test eltype(A.data) == eltype(BosonicHybridizationMatrix(corr).Δ)
        # basis 逐位置固化
        @test A.basis == BosonicContourIndex[BosonicContourIndex(j, b) for (b, j) in chain]
        @test A.data ≈ vec(ref_IF(corr, λsbm)) rtol = 1e-12
    end
    # 一般本征值矢量（非 ±1；含零）
    λg = [0.3, 1.7, -0.9]
    A4 = BosonicIF(BosonicHybridizationMatrix(bosonic_Δτ(f0b, β=2.0, N=3)), λg)
    @test A4.data ≈ vec(ref_IF(bosonic_Δτ(f0b, β=2.0, N=3), λg)) rtol = 1e-12
    @test A4.d == 3 && length(A4.data) == 27 && length(A4.basis) == 3
    # 负检验：节点序确实影响结果——故意用逆序链的参考与实现不同，
    # 保证上面的逐元素对照不是平凡通过。用实时 corr（虚时 Δ 有 KMS 时间反转
    # 对称 Δ = JΔJ，升降序恰好重合，不能用于此检验）
    corr3 = bosonic_Δt(f0b, β=2.0, t=1.0, N=2)
    Tdown = BosonicIF(BosonicHybridizationMatrix(corr3), λsbm)
    function ref_IF_wrong(corr, λ)   # 故意用逆序链
        Δ = BosonicHybridizationMatrix(corr).Δ
        nodes = contour_ordering(corr)
        chain = reverse(nodes)
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
        n = length(contour_ordering(corr))
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
    n = length(contour_ordering(corr))
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
    # 一般非对角耦合（op1 ≠ op2，含复非厄米）：与「精确 cell 演化算符 + 偏迹
    # 环境」的严格单模 bath IF 对照——对角/厄米 hyb 机器精度，一般 hyb 残差
    # O(V⁴)（非高斯累积量；V=0.1 时 ~1e-6）
    begin
        ω, βs, Gc, nb, h = 1.3, 1.7, 0.1, 24, 0.3
        bth = diagm(1 => sqrt.(1.0:(nb - 1)))
        hb = ω * Diagonal(0.0:(nb - 1))
        pth = exp.(-βs * ω * (0.0:(nb - 1))); pth ./= sum(pth)
        spec = DiracDelta(ω, α = Gc^2)
        for hyb in ([0 1; 0 0],                                      # σ⁺（非厄米）
                    [0.2 + 0.3im 0.7; 0.4 - 0.1im 0.15 + 0.25im],    # 随机复非厄米
                    [0.3 + 0.4im 0; 0 -0.2 + 0.1im])                 # 复对角（应机器精度）
            op1, op2 = hyb, hyb'
            # 严格 1-cell IF：V = expm(-ih(hb + G(op1⊗b + op2⊗b†)))，
            # F[m,a,mm,b] = tr(⟨a|V|m⟩ ρ_b ⟨b|V|mm⟩†)，交错布局 (m,a,mm,b)
            V = exp(-im * h * (kron(Matrix(I, 2, 2), Matrix(hb)) +
                               Gc * (kron(op1, Matrix(bth)) + kron(op2, Matrix(bth')))))
            Vb = permutedims(reshape(V, nb, 2, nb, 2), (2, 1, 4, 3))
            Fex = [tr(Vb[a, :, m, :] * Diagonal(pth) * Vb[b, :, mm, :]')
                   for m in 1:2, a in 1:2, mm in 1:2, b in 1:2]
            corrq = bosonic_Δt(spec; β = βs, N = 1, t = h)
            hq = BosonicHybridizationMatrix(corrq; split = true)
            keep = trues(4); keep[1] = keep[2] = false       # 丢 t=0 窗口节点对
            idx = findall(keep)
            basis = [BosonicContourIndex(2, s) for s in (:+, :-)]
            h2 = BosonicHybridizationMatrix(basis,
                  QuAPI.Dense.BosonicSplitDelta(hq.Δ.M1[idx, idx], hq.Δ.M2[idx, idx]))
            P = reshape(BosonicIF(h2, hyb).data, 2, 2, 2, 2)
            if hyb == [0.3 + 0.4im 0; 0 -0.2 + 0.1im]         # 对角：严格 Gaussian
                @test maximum(abs.(P .- Fex)) < 1e-11
            else                                               # O(V⁴) ≈ 1e-6 @ V=0.1
                @test maximum(abs.(P .- Fex)) < 1e-5
            end
        end
    end
    # 非方阵 hyb / 规模上限
    @test_throws ArgumentError BosonicIF(BosonicHybridizationMatrix(corr; split=true), [1.0 2.0 3.0])
    big = bosonic_Δτ(f0b, β=2.0, N=16)          # d=2, n=16 ⟹ 2^32
    @test_throws ArgumentError BosonicIF(BosonicHybridizationMatrix(big; split=true), Matrix(Diagonal(λ)))
end
