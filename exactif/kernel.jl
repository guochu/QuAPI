# =============================================================================
# exactif/kernel.jl —— 严格核装配（Ng–Reichman 论文式 (S35)–(S38)，SM-II）
#
# 输入：spectral.jl 的 (F₁, A, B)（离散闭式或连续谱积分）。
# 输出：ExactKernel —— 论文意义下的严格离散影响泛函核：
#   I({η}) = (Z_N/Z)·exp[ 边界项 + Σ G-块配对 ]（式 (S14)），
#   Grassmann 变量按 (S39) 排序：η⁺_{n+}↔4n+1, η⁻_{n+}↔4n+2,
#                            η⁺_{(n+1)−}↔4n+3, η⁻_{(n+1)−}↔4n+4。
#
# 四个 G-块（式 (S35)–(S38)，按 SM 版式逐条实现；记
#   X ≡ (I−B_N)⁻¹(I−A_N)，A_n = A[1:n,1:n]）：
#   G⁺(m,n)  = −δ_{mn} + [(I−B_N)⁻¹]_{nm}                        (S37)
#   G⁺⁺₊(m,n) = −F₁ Σ_{k≤n} [I − X]_{mk} [(I−A_n)⁻¹]_{k,n−1}      (S35)
#   G⁺⁺₋(m,n) = F₁δ_{n,m−1} − F₁ Σ_{k≤m} [X]_{nk} [(I−A_m)⁻¹]_{k,m−1}  (S36)
#   G⁻(m,n)  = −|F₁|²[ −δ_{mn}[(I−A_m)⁻¹]_{m−1,m−1}
#                + Σ_{k<m}Σ_{ℓ<n} [(I−A_m)⁻¹]_{m−1,k}
#                  ·[(I−A_N)X(I−A_N)]_{kℓ}·[(I−A_n)⁻¹]_{ℓ,n−1} ]  (S38)
# 指标约定：G⁺⁺₊ 中 m=1..N−1, n=1..m；G⁺⁺₋ 中 m=1..N, n=0..m−1；
#           G⁺ 中 m,n=0..N−1；G⁻ 中 m,n=1..N（论文求和范围）。
#
# 对称性 (S21)–(S22)：G⁺⁺₊→G⁻⁻₊ 取复共轭等（runtests 检验 Hermitian 结构）。
#
# ⚠ 字典注记（本研究已数值确认，见 README §5）：论文 Grassmann 变量层的
# G-块与占据基过程张量的对角扇区配对值 I[{u,v}]/I[∅] 之间隔一层逐顶点
# 字典因子（非常数），因此 G-块的最终消费方是 MPS-IF 装配（论文主线），
# 而非直接与占据路径值画等号。
# =============================================================================
using LinearAlgebra: I

"""
    ExactKernel

Ng–Reichman 严格离散影响泛函核（单自旋扇区）：
- `N`、`Δt`、`F1`、`A`、`B`、`ZNoZ`（= Z_N/Z）
- `Gp`（G⁺，N×N，索引 m,n = 0..N−1 → 数组 [m+1,n+1]）
- `Gg`（G⁺⁺₊，(N−1)×N 稠密存放，范围 m=1..N−1, n=1..m）
- `Gl`（G⁺⁺₋，N×N，m=1..N, n=0..m−1）
- `Gm`（G⁻，N×N，m,n=1..N）
- `Θ`：按 (S14)+(S39) 装配的 4N×4N 反对称 Grassmann 核
"""
struct ExactKernel
    N::Int
    Δt::Float64
    F1::ComplexF64
    A::Matrix{ComplexF64}
    B::Matrix{ComplexF64}
    ZNoZ::ComplexF64
    Gp::Matrix{ComplexF64}
    Gg::Matrix{ComplexF64}
    Gl::Matrix{ComplexF64}
    Gm::Matrix{ComplexF64}
    Θ::Matrix{ComplexF64}
end

"""
    exact_kernel(; A, B, F1, ZNoZ, N, Δt) -> ExactKernel

由谱代数层 (F₁, A, B, Z_N/Z) 装配严格核。
"""
function exact_kernel(; A::AbstractMatrix, B::AbstractMatrix, F1::Number,
                      ZNoZ::Number, N::Int, Δt::Real)
    iB = (I - B) \ Matrix{ComplexF64}(I, N, N)
    Xn = iB * (I - A)
    # (S37)
    Gp = [(-((m == n) ? 1 : 0) + iB[n+1, m+1])::ComplexF64 for m in 0:N-1, n in 0:N-1]
    # (S35)：G⁺⁺₊(m,n)，m=1..N−1, n=1..m（论文 0 基矩阵指标 → 1 基：行 m+1，列 k）
    # G⁺⁺₊ = −F₁ Σ_{k=0}^{n−1}[I−X]_{mk}[iA_n]_{k,n−1}；k ≤ n−1 < m ⟹ δ 不触发
    Gg = zeros(ComplexF64, N, N)
    for m in 1:N-1, n in 1:m
        An = A[1:n, 1:n]
        iAn = (I - An) \ Matrix{ComplexF64}(I, n, n)
        Gg[m, n] = F1 * sum(Xn[m+1, k] * iAn[k, n] for k in 1:n)
    end
    # (S36)：G⁺⁺₋(m,n)，m=1..N, n=0..m−1
    Gl = zeros(ComplexF64, N, N)
    for m in 1:N, n in 0:m-1
        Am = A[1:m, 1:m]
        iAm = (I - Am) \ Matrix{ComplexF64}(I, m, m)
        s = zero(ComplexF64)
        for k in 1:m
            s += Xn[n+1, k] * iAm[k, m]
        end
        Gl[m, n+1] = F1 * ((n == m - 1) ? 1 : 0) - F1 * s
    end
    # (S38)：G⁻(m,n)，m,n=1..N
    Gm = zeros(ComplexF64, N, N)
    Yn = (I - A) * iB * (I - A)
    for m in 1:N, n in 1:N
        Am = A[1:m, 1:m]; An = A[1:n, 1:n]
        iAm = (I - Am) \ Matrix{ComplexF64}(I, m, m)
        iAn = (I - An) \ Matrix{ComplexF64}(I, n, n)
        s = zero(ComplexF64)
        for k in 1:m, ℓ in 1:n
            s += iAm[m, k] * Yn[k, ℓ] * iAn[ℓ, n]
        end
        Gm[m, n] = -abs2(F1) * (-((m == n) ? 1 : 0) * iAm[m, m] + s)
    end
    Θ = assemble_theta(Gp, Gg, Gl, Gm, F1, N)
    return ExactKernel(N, Float64(Δt), ComplexF64(F1), A, B, ComplexF64(ZNoZ),
                       Gp, Gg, Gl, Gm, Θ)
end

δ₀(m, k, N) = (m == k) ? 1.0 : 0.0

"""
    assemble_theta(Gp, Gg, Gl, Gm, F1, N) -> 4N×4N

按 (S14)+(S39) 装配 Grassmann 核 Θ（变量位置：
4n+1: η⁺_{n+}，4n+2: η⁻_{n+}，4n+3: η⁺_{(n+1)−}，4n+4: η⁻_{(n+1)−}）。

配对结构（(S14) 各求和项，1-基顶点编号）：
- 边界：Θ[β_{m}⁺(4m−1), α_{m−1}⁺(4m−3)] += F₁；共轭于反向支
- G⁺⁺₊：η⁺_{m+1−}η⁺_{n−1+} → Θ[4n−3, 4m+3] += G⁺⁺₊(m,n)
- G⁺⁺₋：η⁺_{m−1+}η⁺_{n+1−} → Θ[4m−3, 4n+3] += G⁺⁺₋(m,n)
- 反向支取复共轭（G⁻⁻₊ = (G⁺⁺₊)*, G⁻⁻₋ = (G⁺⁺₋)*）
- G⁻：η⁻_{m−1+}η⁺_{n−1+} → Θ[4n−3, 4m−1] += G⁻(m,n)
- G⁺：η⁻_{m+1−}η⁺_{n+1−} → Θ[4n+3, 4m+4] += G⁺(m,n)
（首指标 < 次指标按 (S39) 位置序；符号遵循 (S14) 的书写变量序，
与 MPS-IF 装配约定一致。）
"""
function assemble_theta(Gp, Gg, Gl, Gm, F1, N)
    nvert = 4N
    Θ = zeros(ComplexF64, nvert, nvert)
    add(u, v, x) = begin
        Θ[u, v] += x
        Θ[v, u] -= x
    end
    # 边界项：β_m⁺ × α_{m−1}⁺（位置 4m−1, 4m−3；4m−3 < 4m−1）
    for m in 1:N
        add(4m - 3, 4m - 1, F1)
        add(4m - 2, 4m, -conj(F1))
    end
    # G⁺⁺₊(m,n)：η⁺_{m+1−}(4m+3) × η⁺_{n−1+}(4n−3)，m=1..N−1, n=1..m
    for m in 1:N-1, n in 1:m
        add(4n - 3, 4m + 3, Gg[m, n])
        add(4n - 2, 4m + 4, conj(Gg[m, n]))
    end
    # G⁺⁺₋(m,n)：η⁺_{m−1+}(4m−3) × η⁺_{n+1−}(4n+3)，m=1..N, n=0..m−1
    for m in 1:N, n in 0:m-1
        add(4m - 3, 4n + 3, Gl[m, n+1])
        add(4m - 2, 4n + 4, conj(Gl[m, n+1]))
    end
    # G⁻(m,n)：η⁻_{m−1+}(4m−2) × η⁺_{n−1+}(4n−3)，m,n=1..N
    for m in 1:N, n in 1:N
        add(4n - 3, 4m - 2, -Gm[m, n])
    end
    # G⁺(m,n)：η⁻_{m+1−}(4m+4) × η⁺_{n+1−}(4n+3)，m,n=0..N−1
    for m in 0:N-1, n in 0:N-1
        add(4n + 3, 4m + 4, -Gp[m+1, n+1])
    end
    return Θ
end
