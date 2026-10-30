# =============================================================================
# benchmark/bosons/cellrule_galerkin.jl
#
# 诊断：把冻结核/杂交项做成“高阶 cell 规则”能否把 Δ-only 方案抬到 2/4 阶？
#
# 三种尝试（都用 bosonic_Δτ 的同一离散杂化）：
#   (A) 场重构：ẑ_i = Σ_j T[i,j] z_j，T = 节点 / 中点 / 三次插值 cell 平均
#       → 见 common.jl 的 cellrule_T；阶数恒为 1（imagtime.jl [4] 节已打印）。
#   (B) 一致质量（Galerkin）：K_ab = ∫∫ φ_a Δ φ_b，线性/三次 Hermite 基，
#       子格点 Gauss 求积（等时 cell 用三角形拆分），时间点随阶数增加
#       （线性 N 个自由度；Hermite 三次 2N 个自由度）。
#   (C) 精确自由传播子 + 一致质量 K（把 (B) 的自由项换成节点精确的 e^{-hϵ} 跳跃）。
#
# 结论（本文件运行输出）：
#   · 纯 Galerkin（F − K）在 ∂_τ 上用线性/三次基会产生 O(1) 的自由 G 误差
#     （一阶算符的 Galerkin 离散无耗散、正对称，无法表示 G0 在 τ=0 的跳变，
#     产生振荡）——这条路线在相干态路径积分里不可用。
#   · 换成“精确自由传播子 + 一致质量 K”，阶数仍为 1（与集中质量的冻结核相同）。
#   · 因此 Δ-only 的冻结核收缩最高只能到 1 阶；2/4 阶需要 cell 内精确演化
#     （严格 cell 链），它同样只由 Δ 决定。
#
# 运行：JULIA_PROBE_LIBSTDCXX=0 julia --startup-file=no --project=. benchmark/bosons/cellrule_galerkin.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "performance", "common", "common.jl"))

# Gauss-Legendre 节点/权重直接复用 common.jl 的 gauss01（已修正的 Golub–Welsch）
const GX, GW = gauss01(6)

phi(k, u, h, ::Val{:linear}) = k == 1 ? 1 - u : u
dphi(k, u, h, ::Val{:linear}) = k == 1 ? -1 / h : 1 / h
phi(k, u, h, ::Val{:hermite}) =
    k == 1 ? 1 - 3u^2 + 2u^3 : k == 2 ? h * (u - 2u^2 + u^3) :
    k == 3 ? 3u^2 - 2u^3 : h * (-u^2 + u^3)
dphi(k, u, h, ::Val{:hermite}) =
    k == 1 ? (-6u + 6u^2) / h : k == 2 ? (1 - 4u + 3u^2) :
    k == 3 ? (6u - 6u^2) / h : (-2u + 3u^2)
nbasis(::Val{:linear}) = 2
nbasis(::Val{:hermite}) = 4
ndof(V, N) = nbasis(V) == 2 ? N : 2N
function localdof(I, k, N, V)
    nbasis(V) == 2 && return mod(I - 1 + (k - 1), N) + 1
    node = mod(I - 1 + (k >= 3 ? 1 : 0), N) + 1
    return 2 * (node - 1) + (isodd(k) ? 1 : 2)
end
valdof(i, V) = nbasis(V) == 2 ? i : 2 * (i - 1) + 1

"T-ordered 点核 Δ(τ−τ')（模式求和）"
function deltapoint(mdl, s)
    r = 0.0
    for (w, α) in zip(mdl.ws, mdl.αs)
        n = 1 / (exp(mdl.β * w) - 1)
        r += α * (s >= 0 ? (1 + n) * exp(-w * s) : n * exp(-w * s))
    end
    r
end

"Galerkin 装配：返回 (F, K)；F 为 ∫φ_a(∂+ϵ)φ_b，K 为 ∫∫φ_aΔφ_b。"
function assemble(mdl, N, V)
    h = mdl.β / N; nb = nbasis(V); nd = ndof(V, N)
    F = zeros(nd, nd); K = zeros(nd, nd)
    for I in 1:N, a in 1:nb, b in 1:nb
        s = 0.0
        for (x, w) in zip(GX, GW)
            s += w * phi(a, x, h, V) * (dphi(b, x, h, V) + mdl.ϵd * phi(b, x, h, V))
        end
        F[localdof(I, a, N, V), localdof(I, b, N, V)] += s * h
    end
    for I in 1:N, J in 1:N, a in 1:nb, b in 1:nb
        gi = localdof(I, a, N, V); gj = localdof(J, b, N, V)
        s = 0.0
        if I == J
            for (x, wx) in zip(GX, GW), (y, wy) in zip(GX, GW)
                u = h * x; u2 = h * x * y
                s += wx * wy * h^2 * x * phi(a, x, h, V) * deltapoint(mdl, u - u2) * phi(b, x * y, h, V)
                v2 = h * x; v1 = h * x * y
                s += wx * wy * h^2 * x * phi(a, x * y, h, V) * deltapoint(mdl, v1 - v2) * phi(b, x, h, V)
            end
        else
            for (x, wx) in zip(GX, GW), (y, wy) in zip(GX, GW)
                s += wx * wy * h^2 * phi(a, x, h, V) *
                     deltapoint(mdl, (I - J) * h + h * x - h * y) * phi(b, y, h, V)
            end
        end
        K[gi, gj] += s
    end
    F, K
end

function errs(G, Gex, mdl, N, V)
    h = mdl.β / N; e = 0.0; sc = 0.0
    for i in 1:N, j in 1:N
        i == j && continue
        τ = mod((i - j) * h, mdl.β)
        e = max(e, abs(G[valdof(i, V), valdof(j, V)] - Gex(τ)))
        sc = max(sc, abs(Gex(τ)))
    end
    e / sc
end

function main(; Lmodes=48, Ns=[8, 16, 32, 64])
    mdl = boson_toulouse(; Lmodes=Lmodes)
    Gex = exact_Gτ(mdl)
    io = stdout
    println(io, "=" ^ 74)
    println(io, "高阶 cell 规则（一致质量 / 子格点）能否让 Δ-only 方案升阶？")
    println(io, "=" ^ 74)

    # ---- (B) 纯 Galerkin：先看自由 G（K=0）是否可用 ----
    println(io, "\n[B] 纯 Galerkin 自由项自检（K=0，线性/三次基）")
    n0 = 1 / (exp(mdl.β * mdl.ϵd) - 1)
    for V in (Val(:linear), Val(:hermite))
        N = 16; F, _ = assemble(mdl, N, V); G0 = inv(F)
        e0 = maximum(abs(G0[valdof(i,V), valdof(j,V)] - (1+n0)*exp(-mdl.ϵd*mod((i-j)*mdl.β/N, mdl.β)))
                     for i in 1:N, j in 1:N) / (1 + n0)
        @printf(io, "    %-8s  free-G 相对误差 = %.2e（应为 ~0；>1 表示自由项离散失效）\n", V, e0)
    end

    # ---- (C) 精确自由传播子 + 一致质量 K ----
    println(io, "\n[C] 精确自由传播子 + 一致质量 K（线性基）")
    es = Float64[]
    for N in Ns
        h = mdl.β / N
        A = zeros(N, N)
        for j in 1:N
            A[mod(j, N) + 1, j] = 1.0
        end
        Fex = Matrix{Float64}(I, N, N) - exp(-h * mdl.ϵd) * A
        _, K = assemble(mdl, N, Val(:linear))
        push!(es, errs(inv(Fex - K), Gex, mdl, N, Val(:linear)))
    end
    order_table(io, "    exact-free + consistent mass", Ns, es)

    # ---- (A) 场重构阶梯（节点/中点/三次）----
    println(io, "\n[A] 场重构阶梯（集中质量冻结核，T = 节点 / 中点 / 三次）")
    for p in (1, 2, 4)
        ep = Float64[]
        for N in Ns
            vals, h = frozen_dyson_Gτ(mdl, N; cellrule=p)
            ex = [Gex(mod(m * h, mdl.β)) for m in 1:N-1]
            push!(ep, maximum(abs.(vals .- ex)) / maximum(abs.(ex)))
        end
        order_table(io, "    cell rule p=$p", Ns, ep)
    end

    println(io, "\n结论：三种路线都没有超过 1 阶。Δ-only 冻结核收缩的 O(δτ) 来自")
    println(io, "      “节点精确自由传播子 × cell 积分核”这一配对本身；场重构、")
    println(io, "      一致质量、子格点求积都无法消除。2/4 阶需要 cell 内精确演化")
    println(io, "      （严格 cell 链），它同样只由 Δ 决定。")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
