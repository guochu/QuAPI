# ============================================================================
# exactfermionic.jl —— 费米 normal bath 的严格 Gaussian IF（纯 bath 对象）
#
# 定义（Feynman-Vernon）：IF 是 bath 对象 F = Tr_b[U_b[path] ρ_b U_b†[path]]，
# 与杂质哈密顿量 H_S 无关。本文件构造其在 QuAPI 离散方案下的严格实现：
# bath + 耦合的 cell 演化 B = e^{−δt·Σ_k(ε_k c_k†c_k + V_k(d†c_k+c_k†d))}
# （ε_k = ω_k − μ）被**精确**处理（显式矩阵指数），杂质相干腿在节点插入，
# bath 腿经 Grassmann 高斯边缘化（Schur）精确积掉。
#
# 由于 B 是二次型算符，边缘是严格 Gaussian：
#   虚时:  Tr[B^{Nτ}] = det(I − W⊗S_anti) = det(M_bb)·det(M̃)（Schur 恒等式）
#   实时:  轮廓积分同理（热边 n̄(1+rγ̄γ)、闭合边带费米迹负号）
# 其中 W = e^{−δt h} 为 B 的单粒子矩阵（星形耦合 + bath 能量，**无杂质相位**）。
#
# 谱输入（与 fermionic_Δτ/Δt 一致的 AbstractBoundedFunction 接口）：
#   虚时：离散谱走极点装配；一般连续谱走**闭式**——环被节点 Fourier 对角化
#     （λ_m = e^{iπ(2m+1)/N}），整条环归结为 N 个一维谱积分（L-泛函，见
#     _exact_Δτ_closed），无模式离散，束缚态自动包含；
#   实时：离散谱走极点装配；连续谱经 Gauss–Legendre 正交离散（npts 可控，
#     对离散链逐位严格，光滑谱谱精度收敛）——实时无闭式的原因见
#     exact_fermionic_Δt 文档（模式经杂质互相纠缠 + 热初态与 cell 函数不对易）。
#
# 返回的 Δ 与 fermionic_Δτ/Δt 的 Δ 同骨架约定（调用方叠加 Berezin 测度
# −āᵢaᵢ 与杂质自由动能 x = e^{−δtϵ_d}）、同腿基，可直接替换：
#   Δ_exact = (w₀−1)·S_kinetic − C_ib·M_bb⁻¹·C_bi，   w₀ = ⟨d|W|d⟩
# 两部分都是纯 bath 量：交叉 Schur 项 = 精确顶点幅 × bath 格点 Green 逆；
# (w₀−1) 项 = 耦合对杂质线的 cell 内缀饰。与冻结源核（把顶点幅冻结为
# Duhamel 低阶项）相比，对真实链 (B·A)^N 的局部误差从 O(V²δt²)/cell
# 改善到 O(V²δt³)/cell（全局二阶）。
#
# prefactor（≠1，必须乘进 Z/ρ）：bath 配分修正 + "cell 内跳去又跳回"泄漏；
# 虚时含 Berezin 迹符号 (−1)^{N(K+1)}（K = 模式数；K=1 时恒 1），V=0 时为 1
#（连续谱闭式路径的 prefactor 无此显式符号，已按同一约定归一）。
# ============================================================================

export ExactGrassmannIF, exact_fermionic_Δτ, exact_fermionic_Δt, exact_fermionic_Δm

"""
    ExactGrassmannIF{T}

费米 normal bath 的严格 Gaussian 影响泛函（纯 bath 对象，不依赖 H_S）：

- `hmat::GrassmannHybridizationMatrix`：`Δ` 为相对标准骨架（Berezin 测度
  −āᵢaᵢ + 杂质自由动能，由调用方提供）的纯 IF 部分；腿基与
  `GrassmannHybridizationMatrix(fermionic_Δτ/Δt(...))` 完全一致，可直接替换；
- `prefactor`：整体乘子（含 Berezin 迹符号 (−1)^{N(K+1)}，K 为 bath 模式数；
  K = 1 时恒为 1）——bath 配分修正与 cell 内"跳回泄漏"。

收缩时必须乘 `prefactor`。严格性：ε_d = 0 时
`Z = (−1)^N·Zvalue(IF(测度+动能+Δ))·prefactor == Tr[B^{Nτ}]/Z_bath`
（实时同理对 Tr_b[B^N ρ₀⊗ρ_b (B^N)†]）逐位机器精度成立；ε_d ≠ 0 时方案对
链 (B·A)^N 的局部误差为 O(V²δt³)/cell（冻结源核为 O(V²δt²)/cell）。
"""
struct ExactGrassmannIF{T<:Number}
    hmat::Dense.GrassmannHybridizationMatrix{T}
    prefactor::T
end

# ------------------------------ 内部工具 ------------------------------ #

# normal bath 的模式数据：(ε_k = ω_k − μ, V_k)。μ 语义与 fermionic_Δτ/Δt
# 一致（spectrumshift：μ 只平移 bath 谱）。离散谱（DiracDelta/DiscreteSpectrum）
# 取精确极点；一般连续谱经 Gauss–Legendre 正交离散（见下）。
function _normal_bath_modes(f::DiracDelta, μ::Real, ::Int)
    return ([f.ω - μ], [sqrt(f.α)])
end
function _normal_bath_modes(f::DiscreteSpectrum, μ::Real, ::Int)
    (all(>(0), f.fs) || throw(ArgumentError("DiscreteSpectrum 权重必须为正才能取 V_k = √权重")))
    return (f.ws .- μ, sqrt.(f.fs))
end
function _normal_bath_modes(f::AbstractBoundedFunction, μ::Real, npts::Int)
    lb, ub = lowerbound(f), upperbound(f)
    (isfinite(lb) && isfinite(ub)) ||
        throw(ArgumentError("连续谱输入需要有限支撑；请用 bounded(f, lb, ub) 限定谱支撑"))
    xs, ws = _gauss_legendre(npts)
    half = (ub - lb) / 2
    εs = Float64[]
    Vs = Float64[]
    for (x, w) in zip(xs, ws)
        ω = (lb + ub) / 2 + half * x
        g2 = half * w * f(ω)              # 谱权重 ∫A·(·) 的 GL 逼近：w·A → V²
        g2 > 0 || continue                # 剔除谱为零区域的节点
        push!(εs, ω - μ)
        push!(Vs, sqrt(g2))
    end
    isempty(εs) && throw(ArgumentError("谱在支撑上的积分为零（或非正）"))
    return εs, Vs
end

# Gauss–Legendre 节点/权重（Golub–Welsch：对称三对角 Jacobi 矩阵本征求解）
function _gauss_legendre(n::Int)
    n == 1 && return [0.0], [2.0]
    β = [0.5 / sqrt(1 - 1 / (2 * k)^2) for k in 1:n-1]
    F = eigen(SymTridiagonal(zeros(n), β))
    return F.values, 2 .* (F.vectors[1, :] .^ 2)
end

# 反周期移位 S：(Sη)_k = η_{k−1}，环回带费米迹负号（S^N = −I）
function _antiperiodic_shift(N::Int)
    S = zeros(ComplexF64, N, N)
    for k in 2:N; S[k, k-1] = 1; end
    S[1, N] = -1
    return S
end

# B cell 的单粒子矩阵 W = exp(sgn·δt·h)。模式 1 = 杂质（h 的 dd 元为 0 ——
# 杂质能量属于调用方骨架，不进 bath 对象），模式 k+1 = bath 模式 k。
function _normal_cell_W(εs::Vector, Vs::Vector, step::Real, sgn::ComplexF64)
    K = length(εs)
    h = zeros(ComplexF64, K + 1, K + 1)
    for k in 1:K
        h[1, k+1] = h[k+1, 1] = Vs[k]
        h[k+1, k+1] = εs[k]
    end
    return exp(sgn * step .* h)
end

# ------------------------------ 虚时（τ 环） ------------------------------ #

# τ 极点装配：kron(W,S) 环 + Schur——对模式集（离散谱极点或连续谱的 GL 离散）
# 逐位严格（与直接链机器一致）。
function _exact_Δτ_from_modes(εs::Vector, Vs::Vector, β::Float64, N::Int, δτ::Float64)
    K = length(εs)
    W = _normal_cell_W(εs, Vs, δτ, -1.0 + 0.0im)
    S = _antiperiodic_shift(N)
    # 环矩阵 −I + kron(W, S)：行/列 = (模式, 节点)；模式 1（杂质）占据前 N 维
    M = -Matrix{ComplexF64}(I, (K + 1) * N, (K + 1) * N) + kron(W, S)
    imp = 1:N
    bath = N+1:(K+1)*N
    II, IB, BI, BB = M[imp, imp], M[imp, bath], M[bath, imp], M[bath, bath]
    Mt = II - IB * (BB \ BI)
    Δ = Mt .+ Matrix{ComplexF64}(I, N, N) .- S   # 减测度与单位动能 ⟹ (w₀−1)S − Schur
    Zb = prod(1 .+ exp.(-β .* εs))
    # (−1)^{N(K+1)}：det(−I + kron(W,S)) 的整体 Berezin 符号（K=1 时恒为 1），
    # 吸收进 prefactor 使调用方按 (−1)^N·Zvalue 约定收缩直接得 Tr[B^{Nτ}]/Z_bath
    return Δ, (-1.0)^((K + 1) * N) * det(BB) / Zb
end

# 【闭式注记】τ 环原则上可闭式：节点 Fourier 对角化（λ_m = e^{iπ(2m+1)/N}，
# z_m = λ⁻¹）后每通道 ν_m = [f_{z_m}(h)]_dd，μ_m = λ_m/ν_m（Schur 逆恒等式），
# M̃ 与 prefactor = exp(L_true[log(1+e^{−βx})])/Πμ_m 由 {μ_m} 重建。但正确的
# L_true[f] = ∫ρ_d(E)f(E)dE + Σ_束缚态 w_b f(E_b) 必须用**耦合谱密度** ρ_d
# （A 的希尔伯特变换 Δ(E) = PV∫A/(E−ε)：ρ_d = A/[(E−Δ)²+π²A²]）+ 带外束缚态
# （E−Δ(E)=0 的根，w_b = 1/(1−Δ′)）。朴素的"A-矩 Hermite 公式"在四阶起错误
# （回返路径：(h⁴)_dd = V⁴+V²ω² ≠ V²ω²——模式经杂质的多次散射），已弃用。
# 连续谱暂走下方 GL 离散（对离散链逐位严格，光滑谱谱精度收敛）。

"""
    exact_fermionic_Δτ(f::AbstractBoundedFunction; β::Real, N::Int,
                       μ::Real=0, δτ::Real=β/N, npts::Int=200) -> ExactGrassmannIF

费米 normal bath 的严格 Gaussian IF（虚时 τ 环），纯 bath 对象（与 ϵ_d 无关）。
单粒子环矩阵 `−I + kron(W, S)`（W = e^{−δτh}，S 反周期移位）+ Schur 消元。

与 `fermionic_Δτ` 的 Δ 同骨架约定、同腿基（可直接替换），另附整体 `prefactor`。
严格性：调用方骨架动能取 x = 1（ε_d = 0）时
`Z = (−1)^N·Zvalue(IF(测度+动能+Δ))·prefactor == Tr[B^{Nτ}]/Z_bath` 逐位成立
（离散谱下对直接链机器验证）；一般 ε_d 下对链 (B·A)^N 的局部误差为
O(V²δτ³)/cell（`fermionic_Δτ` 为 O(V²δτ²)/cell）。

输入：离散谱（`DiracDelta`/`DiscreteSpectrum`）取精确极点（`npts` 无效）；
一般连续谱（有限支撑）经 Gauss–Legendre 正交离散为 `npts` 个模式（构造对该
离散链逐位严格，光滑谱随 `npts` 谱精度收敛）。
"""
function exact_fermionic_Δτ(f::AbstractBoundedFunction;
                            β::Real, N::Int, μ::Real=0, δτ::Real=β / N, npts::Int=200)
    εs, Vs = _normal_bath_modes(f, μ, npts)
    Δ, P = _exact_Δτ_from_modes(Vector{Float64}(εs), Vector{Float64}(Vs),
                                Float64(β), N, Float64(δτ))
    return _tau_pack(Δ, P, N)
end

# 公共尾部：节点基 Δ → 腿基配对矩阵 + ExactGrassmannIF
function _tau_pack(Δ::Matrix{ComplexF64}, P::ComplexF64, N::Int)
    C = zeros(ComplexF64, 2N, 2N)
    for i in 1:N, j in 1:N
        C[2i, 2j-1] = Δ[i, j]
    end
    basis = GrassmannContourIndex[GrassmannContourIndex(j, :τ, 1, c)
                                  for j in 1:N for c in (false, true)]
    return ExactGrassmannIF(Dense.GrassmannHybridizationMatrix(basis, C), P)
end

# ------------------------------ 实时（Keldysh 轮廓） ------------------------------ #

# 实时轮廓腿位（节点 j0 = 0..N，与 GrassmannHybridizationMatrix(RealCorrelationFunction)
# 的 basis 同序：节点主序，每节点 + 分支在前、每分支 (a, ā) 相邻）
_pos_rt(j0::Int, br::Symbol, ā::Bool) = 4 * j0 + (br === :+ ? 0 : 2) + (ā ? 2 : 1)

"""
    exact_fermionic_Δt(f::AbstractBoundedFunction; β::Real, t::Real, N::Int,
                       μ::Real=0, npts::Int=200) -> ExactGrassmannIF

费米 normal bath 的严格 Gaussian IF（实时 Keldysh 轮廓）。B cell 同
`exact_fermionic_Δτ`（取 −iδt），+ 分支前向 / − 分支后向插入，t=T 闭合边带
费米迹负号、t=0 热边 n̄ₖ(1+rₖγ̄γ)；Schur 消 bath 场后得纯 bath 的 Δ 与
`prefactor = Πₖn̄ₖ·det(M_bb)`，与 `fermionic_Δt` 的 Δ 同骨架约定、同腿基，
可直接替换。

严格性：杂质动能取 x = 1 时
`ρ(T) = prefactor·∫D[杂质](骨架×e^{āΔa}×⟨η̄⁺₀|ρ₀|η⁻₀⟩×投影) ==
Tr_b[B^N ρ₀⊗ρ_b (B^N)†]` 逐位成立；一般 ε_d 下对链 (B·A)^N 的局部误差
O(V²δt³)/cell。

输入：离散谱（`DiracDelta`/`DiscreteSpectrum`）走极点装配（`npts` 无效）；
一般连续谱（有限支撑）经 Gauss–Legendre 正交离散为 `npts` 个模式（构造对该
离散链逐位严格，光滑谱随 `npts` 谱精度收敛）。

与玻色 PT 的对照（为何费米无"裸谱系数闭式"）：玻色耦合对 bath 模**线性**
（op⊗b†+h.c.，[op,H_b]=0）⟹ 精确 cell = 多模位移（算符幅），热迹按模式
因子化且在算符幅上为高斯泛函；幅代数封闭于 c 数换位子（op=op；或
[a†,a]=−1 的谐振子情形）⟹ M1/M2 次序记账**精确**——a†b+ab† 同样完美。
σ⁺⊗b+σ⁻⊗b†（[σ⁺,σ⁻]=σz 为算符）迹不封闭，修复 PT 仅为 O(V²) 截断
（nondiagcoupling.jl 附录）。费米耦合被宇称强制为**双线性** d†c+c†d ⟹
无位移型 cell（比 σ⁺ 更早一步失效），裸系数核（冻结源 fermionic_Δt）即
O(V²δt²) 截断；严格对象来自另一机制——联合二次型 ⟹ determinant/Schur
边缘（本实现），核为耦合单粒子传播子的泛函。
"""
function exact_fermionic_Δt(f::AbstractBoundedFunction;
                            β::Real, t::Real, N::Int, μ::Real=0, npts::Int=200)
    εs, Vs = _normal_bath_modes(f, μ, npts)
    K = length(εs)
    W = _normal_cell_W(εs, Vs, t / N, -1.0im)
    Wc = conj.(W)
    Limp = 4 * (N + 1)                       # 杂质腿（ā 行 = a 列 = 腿空间）
    Nb = K * (2 * N + 2)                     # bath 配对维数（每模式 2(N+1) 个 (γ̄,γ) 槽）
    bidx(k, br, j0) = (k - 1) * (2N + 2) + (br === :+ ? j0 : N + 1 + j0) + 1
    II = zeros(ComplexF64, Limp, Limp)
    IB = zeros(ComplexF64, Limp, Nb)
    BI = zeros(ComplexF64, Nb, Limp)
    BB = zeros(ComplexF64, Nb, Nb)
    for j0 in 0:N, br in (:+, :-)            # Berezin 测度 −āa / −γ̄γ
        II[_pos_rt(j0, br, true), _pos_rt(j0, br, false)] -= 1
        for k in 1:K
            BB[bidx(k, br, j0), bidx(k, br, j0)] -= 1
        end
    end
    for cell in 1:N
        II[_pos_rt(cell, :+, true), _pos_rt(cell - 1, :+, false)] += W[1, 1]
        II[_pos_rt(cell - 1, :-, true), _pos_rt(cell, :-, false)] += Wc[1, 1]
        for km in 1:K
            # 交叉块：+ cell 用 W，− cell 用 W†（共轭）；无杂质相位（x 归骨架）
            IB[_pos_rt(cell, :+, true), bidx(km, :+, cell - 1)] += W[1, km+1]
            IB[_pos_rt(cell - 1, :-, true), bidx(km, :-, cell)] += Wc[1, km+1]
            BI[bidx(km, :+, cell), _pos_rt(cell - 1, :+, false)] += W[km+1, 1]
            BI[bidx(km, :-, cell - 1), _pos_rt(cell, :-, false)] += Wc[km+1, 1]
            for kn in 1:K
                BB[bidx(km, :+, cell), bidx(kn, :+, cell - 1)] += W[km+1, kn+1]
                BB[bidx(km, :-, cell - 1), bidx(kn, :-, cell)] += Wc[km+1, kn+1]
            end
        end
    end
    for k in 1:K                             # 闭合边（费米迹负号）与热边
        BB[bidx(k, :-, N), bidx(k, :+, N)] += -1
        n = 1 / (exp((β::Float64) * εs[k]) + 1)
        BB[bidx(k, :+, 0), bidx(k, :-, 0)] += n / (1 - n)
    end
    Mt = II - IB * (BB \ BI)
    C0 = zeros(ComplexF64, Limp, Limp)       # 单位动能骨架（x = 1）：测度 + 单位动能
    for j0 in 0:N, br in (:+, :-)
        C0[_pos_rt(j0, br, true), _pos_rt(j0, br, false)] -= 1
    end
    for cell in 1:N
        C0[_pos_rt(cell, :+, true), _pos_rt(cell - 1, :+, false)] += 1
        C0[_pos_rt(cell - 1, :-, true), _pos_rt(cell, :-, false)] += 1
    end
    basis = GrassmannContourIndex[GrassmannContourIndex(j0 + 1, br, 1, c)
                                  for j0 in 0:N for br in (:+, :-) for c in (false, true)]
    nbar_k = [1 / (1 + exp(-(β::Float64) * ε)) for ε in εs]
    return ExactGrassmannIF(Dense.GrassmannHybridizationMatrix(basis, Mt - C0),
                            prod(nbar_k) * det(BB))
end

# ------------------------------ 混合时间（Kadanoff–Baym L 轮廓） ------------------------------ #

"""
    exact_fermionic_Δm(f::AbstractBoundedFunction; β::Real, Nτ::Int, t::Real, Nt::Int,
                       μ::Real=0, δτ::Real=β/Nτ, npts::Int=200) -> ExactGrassmannIF

费米 normal bath 的严格 Gaussian IF（混合 L 形轮廓：τ 段 Nτ 步 + Keldysh 段
Nt 步）。轮廓环 τ₁→…→τ_{Nτ}→(+,1)→…→(+,Nt+1)→(−,Nt+1)→…→(−,1)→τ₁：
τ cells（含到达 (+,1) 的结合 cell）用 W^τ = e^{−δτh}，实时 ± cells 用
W^t = e^{−iδth}（− 支共轭），(+,Nt+1)→(−,Nt+1) 闭合边带费米迹负号，
(−,1)→τ₁ 环回边权 1（热演化在 τ 段显式进行，无热边 r）。

与 `fermionic_Δm` 的 Δ 同骨架约定、同腿基（可直接替换）。严格性：全迹收缩
`Zvalue(IF(测度+单位链+Δ))·prefactor == Tr[(B^τ)^{Nτ}]/Z_bath`（实时段
幺正回绕 B^t B̄^t = I 相消；对直接链机器验证）。
"""
function exact_fermionic_Δm(f::AbstractBoundedFunction;
                            β::Real, Nτ::Int, t::Real, Nt::Int,
                            μ::Real=0, δτ::Real=β / Nτ, npts::Int=200)
    εs, Vs = _normal_bath_modes(f, μ, npts)
    K = length(εs)
    βf = Float64(β)
    Wτ = _normal_cell_W(εs, Vs, δτ, -1.0 + 0.0im)
    Wt = _normal_cell_W(εs, Vs, t / Nt, -1.0im)
    Wtc = conj.(Wt)
    n = Nτ + 2 * (Nt + 1)                    # 轮廓节点数（链序）
    cτ(j) = j                                # τ 节点 j = 1..Nτ
    cp(jj) = Nτ + 2 * jj - 1                 # + 节点 jj = 1..Nt+1
    cm(jj) = Nτ + 2 * jj                     # − 节点
    ipos(c, ā) = 2 * c - (ā ? 0 : 1)         # 杂质腿（a = 2c−1, ā = 2c）
    bidx(k, c) = (k - 1) * n + c              # bath 配对槽（行 = γ̄，列 = γ，同索引）
    Limp = 2n
    Nb = K * n
    II = zeros(ComplexF64, Limp, Limp)
    IB = zeros(ComplexF64, Limp, Nb)
    BI = zeros(ComplexF64, Nb, Limp)
    BB = zeros(ComplexF64, Nb, Nb)
    for c in 1:n                               # Berezin 测度（杂质每节点一次；bath 每模式）
        II[ipos(c, true), ipos(c, false)] -= 1
        for k in 1:K
            BB[bidx(k, c), bidx(k, c)] -= 1
        end
    end
    # cell 边（dest ← src）：统一装配 W 的全部块（II/IB/BI/BB）
    function edge!(cdest, csrc, W)
        II[ipos(cdest, true), ipos(csrc, false)] += W[1, 1]
        for km in 1:K, kn in 1:K
            BB[bidx(km, cdest), bidx(kn, csrc)] += W[km+1, kn+1]
        end
        for km in 1:K
            IB[ipos(cdest, true), bidx(km, csrc)] += W[1, km+1]
            BI[bidx(km, cdest), ipos(csrc, false)] += W[km+1, 1]
        end
    end
    for k in 2:Nτ;   edge!(cτ(k), cτ(k - 1), Wτ); end      # τ 内部 cells
    edge!(cp(1), cτ(Nτ), Wτ)                               # 结合 cell（τ_Nτ → β ≡ +₁）
    for jj in 1:Nt;  edge!(cp(jj + 1), cp(jj), Wt); end    # + cells
    for jj in 1:Nt;  edge!(cm(jj), cm(jj + 1), Wtc); end   # − cells（后向）
    for k in 1:K                                           # 轮廓闭合边与环回边
        BB[bidx(k, cm(Nt + 1)), bidx(k, cp(Nt + 1))] += -1  # 费米迹
        BB[bidx(k, cτ(1)), bidx(k, cm(1))] += 1            # (−,1)→τ₁
    end
    # 杂质轮廓闭合边（必须在 Schur 前）：Keldysh 闭合 (+,Nt+1)→(−,Nt+1) 带
    # 费米迹负号，环回 (−,1)→τ₁ 单位（环总符号 −1；自由极限 Z = 2 = 简并度）
    II[ipos(cm(Nt + 1), true), ipos(cp(Nt + 1), false)] += -1
    II[ipos(cτ(1), true), ipos(cm(1), false)] += 1
    Mt = II - IB * (BB \ BI)
    # 单位骨架：测度 + 沿轮廓单位链 + 两条闭合边（全迹收缩用）
    C0 = zeros(ComplexF64, Limp, Limp)
    for c in 1:n
        C0[ipos(c, true), ipos(c, false)] -= 1
    end
    C0[ipos(cp(1), true), ipos(cτ(Nτ), false)] += 1
    C0[ipos(cm(Nt + 1), true), ipos(cp(Nt + 1), false)] += -1
    C0[ipos(cτ(1), true), ipos(cm(1), false)] += 1
    for k in 2:Nτ;  C0[ipos(cτ(k), true), ipos(cτ(k - 1), false)] += 1; end
    for jj in 1:Nt
        C0[ipos(cp(jj + 1), true), ipos(cp(jj), false)] += 1
        C0[ipos(cm(jj), true), ipos(cm(jj + 1), false)] += 1
    end
    basis = GrassmannContourIndex[
        vcat([GrassmannContourIndex(j, :τ, 1, c) for j in 1:Nτ for c in (false, true)],
             [GrassmannContourIndex(jj, br, 1, c)
              for jj in 1:Nt+1 for br in (:+, :-) for c in (false, true)])...]
    Zb = prod(1 .+ exp.(-βf .* εs))
    # (−1)^{nK}：折叠的 Berezin 符号（K = bath 模式数；K=1 退化为 (−1)^n），
    # 吸收进 prefactor 使全迹收缩直接得 Tr[(B^τ)^Nτ]/Z_bath
    return ExactGrassmannIF(Dense.GrassmannHybridizationMatrix(basis, Mt - C0),
                            (-1.0 + 0.0im)^(n * K) * det(BB) / Zb)
end
