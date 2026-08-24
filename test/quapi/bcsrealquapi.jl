println("------------------------------------")
println("|  Real-time QuAPI (BCS bath)      |")
println("------------------------------------")

using Test, LinearAlgebra
include("../../src/includes.jl")

# =====================================================================
# Real-time QuAPI for the Toulouse model with a BCS bath
#
#     H = Σ_σ ϵ_d d_σ†d_σ + Σ_kσ ω_k c_kσ†c_kσ + Σ_kσ V_k (d_σ†c_kσ + c_kσ†d_σ)
#         + Δ Σ_k (c_k↑†c_k↓† + c_k↓c_k↑),
#
# starting from the product state ρ(0) = ρ_imp ⊗ ρ_bath(β).  The discrete
# Keldysh weight is exp(½ξᵀBξ) × (initial polynomial) × (final polynomial)
# over the 4(N+1) Grassmann variables per flavor, and all Berezin integrals
# are evaluated exactly with Pfaffians (Wick's theorem), so no path
# enumeration is needed.
#
# Conventions (validated against exact diagonalization):
#   variables ξ = (ā⁺₀,a⁺₀,…,ā⁺_N,a⁺_N, ā⁻₀,a⁻₀,…,ā⁻_N,a⁻_N) per flavor,
#                 flavor ↓ offset by 4(N+1)
#   measure:     -āᵢaᵢ on every node (both branches)
#   forward:     e^{-im·δt·ϵ_d}·ā_{i+1}a_i        (i = 0..N-1)
#   backward:    e^{+im·δt·ϵ_d}·ā_i a_{i+1}       (i = 0..N-1)
#   influence:   exp[Σ_{uv} Σ_{ij} σ_u σ_v āᵢᵘ Δᵢⱼᵘᵛ aⱼᵛ] with σ₊=+1, σ₋=-1,
#                i.e. the cross-branch blocks Δ⁺⁻, Δ⁻⁺ enter with -1
#   BCS normal blocks: BOTH spin flavors use `corr.uu` (Δ_uu - Δ_vvᵀ).
#                G_↓↓ = G_↑↑ by spin symmetry of the singlet-paired bath,
#                whereas `corr.dd` = -Δ_uuᵀ at Δ=0 is the hole (Nambu 22)
#                partner and does NOT decouple to the single-flavor result.
#   BCS anomalous blocks: ā↑ Δ↑↓ ā↓ and a↓ Δ↓↑ a↑ from `corr.ud`, `corr.du`.
#   initial:     ⟨a⁺₀|ρ_imp|a⁻₀⟩,  final: ⟨α|a⁺_N⟩⟨ā⁻_N|β⟩,  monomials in
#                Keldysh contour order 0⁺ → N⁺ → N⁻ → 0⁻
#
# The BCS ED reference is the strict diagonalization of the quadratic
# (BdG) Hamiltonian from ImpurityModelBase.exactdiagonalizations.
# =====================================================================

# =====================================================================
# Pfaffian of an antisymmetric matrix (Parlett–Reid, O(n³))
# =====================================================================
function pfaffian(A::AbstractMatrix)
    n = size(A, 1)
    (iseven(n)) || throw(ArgumentError("Pfaffian requires even dimension"))
    n == 0 && return one(eltype(A))
    B = Matrix{ComplexF64}(A)
    p = 1.0 + 0.0im
    for k in 1:2:n-1
        if iszero(B[k, k+1])
            found = false
            for i in k+2:n
                if !iszero(B[k, i])
                    # swap rows/columns i <-> k+1 (each swap flips Pf by -1)
                    B[:, [k+1, i]] .= B[:, [i, k+1]]
                    B[[k+1, i], :] .= B[[i, k+1], :]
                    p = -p
                    found = true
                    break
                end
            end
            found || return zero(ComplexF64)
        end
        p *= B[k, k+1]
        for i in k+2:n, j in i:n
            B[i, j] = (B[i, j] * B[k, k+1] - B[i, k] * B[j, k+1] + B[i, k+1] * B[j, k]) / B[k, k+1]
            B[j, i] = -B[i, j]
        end
    end
    return p
end

# Gaussian Grassmann integral with monomial insertions:
#   ∫ exp(½ξᵀBξ) ξ_{i₁}···ξ_{i_{2k}} = Pf(B) · Pf(C),  C_{rs} = (B⁻¹)_{iᵣ,iₛ}
function gauss_integral(B::AbstractMatrix, monomial::Vector{Int})
    m = length(monomial)
    iseven(m) || return 0.0 + 0.0im
    m == 0 && return pfaffian(B)
    Binv = inv(B)
    C = Binv[monomial, monomial]
    return pfaffian(B) * pfaffian(C)
end

# =====================================================================
# Variable indexing
# =====================================================================
# single flavor: M = 4(N+1) variables.
#   ind(ā, u, i) = base(u) + 2i + 1 ; ind(a, u, i) = base(u) + 2i + 2
#   base(+) = 0, base(-) = 2(N+1)
# BCS: two flavors, flavor ↓ offset by 4(N+1)
function quark_single(N)
    (u, i) -> (u === :+ ? 0 : 2 * (N + 1)) + 2i + 1   # ā
end
function aindex_single(N)
    (u, i) -> (u === :+ ? 0 : 2 * (N + 1)) + 2i + 2   # a
end
function quark_bcs(N, σ)
    off = (σ == 2 ? 4 * (N + 1) : 0)
    (u, i) -> off + (u === :+ ? 0 : 2 * (N + 1)) + 2i + 1   # ā_σ
end
function aindex_bcs(N, σ)
    off = (σ == 2 ? 4 * (N + 1) : 0)
    (u, i) -> off + (u === :+ ? 0 : 2 * (N + 1)) + 2i + 2   # a_σ
end

# =====================================================================
# Build the antisymmetric matrix B for a single flavor
# =====================================================================
function build_B_single(N, β, ϵ_d, t, spec, μ)
    δt = t / N
    M = 4 * (N + 1)
    B = zeros(ComplexF64, M, M)
    ā = quark_single(N)
    a = aindex_single(N)

    # helper: add term c·ξ_p ξ_q  (½ξᵀBξ = Σ_{i<j} B[i,j] ξ_i ξ_j)
    addterm!(p, q, c) = (B[p, q] += c; B[q, p] -= c)

    # measures: -āᵢaᵢ
    for u in (:+, :-), i in 0:N
        addterm!(ā(u, i), a(u, i), -1.0)
    end

    # forward / backward hoppings
    x = exp(-im * δt * ϵ_d)
    for i in 0:N-1
        addterm!(ā(:+, i + 1), a(:+, i), x)        # x·ā⁺_{i+1}a⁺_i
        addterm!(ā(:-, i), a(:-, i + 1), conj(x))  # conj(x)·ā⁻_i a⁻_{i+1}
    end

    # influence functional: cross branches enter with -1 (Keldysh relative sign)
    corr = fermionic_Δt(spec, β=β, N=N, t=t, μ=μ)
    for u in (:+, :-), v in (:+, :-)
        G = Matrix(branch(corr, u, v))
        s = (u === v) ? 1.0 : -1.0
        for i in 0:N, j in 0:N
            c = s * G[i + 1, j + 1]
            iszero(c) && continue
            addterm!(ā(u, i), a(v, j), c)
        end
    end
    return B
end

# =====================================================================
# Build the antisymmetric matrix B for the BCS model (two flavors)
# =====================================================================
function build_B_bcs(N, β, ϵ_d, t, spec, Δ, μ)
    δt = t / N
    M = 8 * (N + 1)
    B = zeros(ComplexF64, M, M)
    ā1, a1 = quark_bcs(N, 1), aindex_bcs(N, 1)   # ↑
    ā2, a2 = quark_bcs(N, 2), aindex_bcs(N, 2)   # ↓

    addterm!(p, q, c) = (B[p, q] += c; B[q, p] -= c)

    # measures
    for u in (:+, :-), i in 0:N
        addterm!(ā1(u, i), a1(u, i), -1.0)
        addterm!(ā2(u, i), a2(u, i), -1.0)
    end

    # hoppings (both flavors, same ϵ_d)
    x = exp(-im * δt * ϵ_d)
    for i in 0:N-1
        for (ā, a) in ((ā1, a1), (ā2, a2))
            addterm!(ā(:+, i + 1), a(:+, i), x)
            addterm!(ā(:-, i), a(:-, i + 1), conj(x))
        end
    end

    # BCS influence functional (normal + anomalous blocks)
    corr = bcs_Δt(spec, β=β, N=N, t=t, Δ=Δ, μ=μ)
    for u in (:+, :-), v in (:+, :-)
        s = (u === v) ? 1.0 : -1.0
        Gnn_up = Matrix(branch(corr.uu, u, v))   # ā↑ Δ↑↑ a↑
        Gnn_dn = Matrix(branch(corr.uu, u, v))   # ā↓ Δ↓↓ a↓ (spin symmetry: G↓↓ = G↑↑)
        Gan_ud  = Matrix(branch(corr.ud, u, v))  # ā↑ Δ↑↓ ā↓
        Gan_du  = Matrix(branch(corr.du, u, v))  # a↓ Δ↓↑ a↑
        for i in 0:N, j in 0:N
            for (c, p1, p2) in ((s * Gnn_up[i+1, j+1], ā1(u, i), a1(v, j)),
                                (s * Gnn_dn[i+1, j+1], ā2(u, i), a2(v, j)),
                                (s * Gan_ud[i+1, j+1], ā1(u, i), ā2(v, j)),
                                (s * Gan_du[i+1, j+1], a2(u, i), a1(v, j)))
                iszero(c) && continue
                addterm!(p1, p2, c)
            end
        end
    end
    return B
end

# =====================================================================
# Reduced density matrix (single flavor) at time t
# =====================================================================
# ρ_imp given as 2×2 matrix in the Fock basis {|0⟩, |1⟩}.  Boundary monomials
# are written in Keldysh contour order 0⁺ → N⁺ → N⁻ → 0⁻:
#   ⟨a⁺₀|s1⟩ (init 0⁺) · ⟨α|a⁺_N⟩ (final N⁺) · ⟨a⁻_N|β⟩ (final N⁻) · ⟨s2|a⁻₀⟩ (init 0⁻)
function quapi_rt_single(N, β, ϵ_d, t, spec, ρ_imp, μ=0.0)
    B = build_B_single(N, β, ϵ_d, t, spec, μ)
    ā = quark_single(N)
    a = aindex_single(N)

    init_terms = ComplexF64[]
    init_plus = Vector{Vector{Int}}()
    init_minus = Vector{Vector{Int}}()
    push!(init_terms, ρ_imp[1, 1]); push!(init_plus, Int[]); push!(init_minus, Int[])
    push!(init_terms, ρ_imp[2, 1]); push!(init_plus, [ā(:+, 0)]); push!(init_minus, Int[])
    push!(init_terms, ρ_imp[1, 2]); push!(init_plus, Int[]); push!(init_minus, [a(:-, 0)])
    push!(init_terms, ρ_imp[2, 2]); push!(init_plus, [ā(:+, 0)]); push!(init_minus, [a(:-, 0)])

    final_plus = Dict(1 => Int[], 2 => [a(:+, N)])
    final_minus = Dict(1 => Int[], 2 => [ā(:-, N)])

    ρ = zeros(ComplexF64, 2, 2)
    for α in 1:2, β in 1:2
        val = 0.0 + 0.0im
        for (c, ip, im_) in zip(init_terms, init_plus, init_minus)
            mono = vcat(ip, final_plus[α], final_minus[β], im_)
            val += c * gauss_integral(B, mono)
        end
        ρ[α, β] = val
    end
    Z = ρ[1, 1] + ρ[2, 2]
    return ρ ./ Z, Z
end

# =====================================================================
# Reduced density matrix (BCS, two flavors) at time t
# =====================================================================
# ρ_imp given as 4×4 matrix in the Fock basis {|0⟩,|↑⟩,|↓⟩,|↑↓⟩} with
# |↑↓⟩ = ā↓†ā↑†|0⟩, so that ⟨a|↑↓⟩ = ā↓ā↑ and ⟨↑↓|a⟩ = a↑a↓.
function quapi_rt_bcs(N, β, ϵ_d, t, spec, Δ, ρ_imp, μ=0.0)
    B = build_B_bcs(N, β, ϵ_d, t, spec, Δ, μ)
    ā1, a1 = quark_bcs(N, 1), aindex_bcs(N, 1)
    ā2, a2 = quark_bcs(N, 2), aindex_bcs(N, 2)

    # ⟨a⁺₀|s1⟩ (init_plus) and ⟨s2|a⁻₀⟩ (init_minus) for the 4 Fock states
    bra_state_plus = Dict(1 => Int[],                              # |0⟩ → 1
                          2 => [ā1(:+, 0)],                        # |↑⟩ → ā↑
                          3 => [ā2(:+, 0)],                        # |↓⟩ → ā↓
                          4 => [ā2(:+, 0), ā1(:+, 0)])             # |↑↓⟩ → ā↓ā↑
    ket_state_minus = Dict(1 => Int[],                             # ⟨0| → 1
                           2 => [a1(:-, 0)],                       # ⟨↑| → a↑
                           3 => [a2(:-, 0)],                       # ⟨↓| → a↓
                           4 => [a1(:-, 0), a2(:-, 0)])            # ⟨↑↓| → a↑a↓

    # final: ⟨α|a⁺_N⟩ (final_plus) and ⟨a⁻_N|β⟩ (final_minus)
    bra_final_plus = Dict(1 => Int[],
                          2 => [a1(:+, N)],
                          3 => [a2(:+, N)],
                          4 => [a1(:+, N), a2(:+, N)])
    ket_final_minus = Dict(1 => Int[],
                           2 => [ā1(:-, N)],
                           3 => [ā2(:-, N)],
                           4 => [ā2(:-, N), ā1(:-, N)])

    ρ = zeros(ComplexF64, 4, 4)
    for α in 1:4, β in 1:4
        val = 0.0 + 0.0im
        for s1 in 1:4, s2 in 1:4
            c = ρ_imp[s1, s2]
            iszero(c) && continue
            mono = vcat(bra_state_plus[s1], bra_final_plus[α], ket_final_minus[β], ket_state_minus[s2])
            val += c * gauss_integral(B, mono)
        end
        ρ[α, β] = val
    end
    Z = ρ[1, 1] + ρ[2, 2] + ρ[3, 3] + ρ[4, 4]
    return ρ ./ Z, Z
end

# =====================================================================
# Exact ED reference: BCS Toulouse model (quadratic ED, ImpurityModelBase)
# =====================================================================
# Without interactions the Toulouse model coupled to a BCS bath is a
# quadratic Hamiltonian, so its real-time dynamics can be computed exactly
# with the strict diagonalization routines of
# ImpurityModelBase.exactdiagonalizations (the cdm obeys dρ/dt = -i[hᵀ, ρ],
# so the eigencache is built from the transposed BdG matrix):
#   m    = Toulouse(DiscreteBCSBath(ωs, αs, Δ; β, μ=0); ϵ_d)
#   h    = cmatrix(m)             # 2L×2L BdG coefficient matrix
#   ρ₀   = separablecdm(m, nsys)  # impurity occupation nsys ⊗ thermal bath
#   ρ(t) = timeevo(ρ₀, hᵀ, +im·t) # exp(+im·t·hᵀ)·ρ₀·exp(-im·t·hᵀ)
# with ⟨n_σ(t)⟩ = ρ(t)[σ, σ] (cdm entries encode ⟨c_j†c_i⟩).  Note the
# DiscreteBCSBath constructor takes the spectral density αs (so the coupling
# in the BdG matrix is √αs = Vs, matching the QuAPI spectral function).
function bcs_ed_imb(β, ϵ_d, ωs, Vs, Δ, nsys, ts; μ=0.0)
    (μ == 0) || throw(ArgumentError("BCS bath should have μ=0"))
    m = Toulouse(DiscreteBCSBath(ωs, Vs .^ 2, Δ; β=β, μ=μ); ϵ_d=ϵ_d)
    h = cmatrix(m)
    cache = eigencache(transpose(h))
    ρ₀ = separablecdm(m, nsys)
    out = Tuple{Float64, Float64}[]
    for t in ts
        ρt = timeevo(ρ₀, transpose(h), im * t, cache)
        push!(out, (real(ρt[1, 1]), real(ρt[2, 2])))   # (⟨n_↑⟩, ⟨n_↓⟩)
    end
    return out
end

# =====================================================================
# Tests
# =====================================================================
@testset "real-time QuAPI (BCS bath)" begin
    β, ϵ_d, μ = 2.0, 0.5, 0.0
    ρ0 = Matrix{ComplexF64}(I, 4, 4)
    ρ0[1, 1] = 1.0; ρ0[2, 2] = 0.0; ρ0[3, 3] = 0.0; ρ0[4, 4] = 0.0   # |0⟩⟨0|

    # Free-impurity limit (λ = 0, Δ = 0): the two spin flavors decouple and
    # evolve exactly as ρ(t) = U(t) ρ_imp U(t)† with U = diag(1, e^{-iϵt},
    # e^{-iϵt}, e^{-2iϵt}), pinning down the boundary-monomial convention.
    @testset "free impurity (λ=0, Δ=0)" begin
        spec0 = DiracDelta(ω=1.0, α=0.0)
        for (t, N) in ((0.5, 2), (1.0, 3))
            ρb, Zb = quapi_rt_bcs(N, β, ϵ_d, t, spec0, 0.0, ρ0, μ)
            @test Zb ≈ 1.0 atol=1.0e-8   # unitary evolution preserves tr ρ
            @test ρb ≈ ρ0 atol=1.0e-8    # the vacuum is stationary
            # coherent (all-equal) initial state: tests the off-diagonal /
            # cross-branch structure of the 4×4 reduced density matrix
            ρc = fill(0.25 + 0.0im, 4, 4)
            ρb, _ = quapi_rt_bcs(N, β, ϵ_d, t, spec0, 0.0, ρc, μ)
            U = Diagonal([1.0, exp(-im * ϵ_d * t), exp(-im * ϵ_d * t), exp(-2im * ϵ_d * t)])
            @test ρb ≈ U * ρc * U' atol=1.0e-8
        end
    end

    # Δ = 0: the two spin channels must decouple to two identical
    # single-flavor Toulouse problems.  This pins down the ↓-channel IF block
    # (both normal blocks must equal `corr.uu`; `corr.dd` is the hole partner
    # and does NOT reproduce the single-flavor result at Δ=0).
    @testset "Δ=0 decoupling vs single flavor" begin
        V = 0.25
        spec = DiracDelta(ω=1.0, α=V^2)
        t, N = 1.5, 8
        ρb, Zb = quapi_rt_bcs(N, β, ϵ_d, t, spec, 0.0, ρ0, μ)
        ρs, Zs = quapi_rt_single(N, β, ϵ_d, t, spec, [1.0 0.0; 0.0 0.0], μ)
        @test real(ρb[2, 2] + ρb[4, 4]) ≈ real(ρs[2, 2]) atol=1.0e-8   # ⟨n_↑⟩
        @test real(ρb[3, 3] + ρb[4, 4]) ≈ real(ρs[2, 2]) atol=1.0e-8   # ⟨n_↓⟩
        @test Zb ≈ Zs^2 rtol=1.0e-8                                     # independent flavors
    end

    # BCS with pairing gap vs exact diagonalization (K=1 bath mode).  The
    # normal block of both channels is `corr.uu`; the anomalous blocks
    # ā↑Δ↑↓ā↓ + a↓Δ↓↑a↑ enter through `corr.ud` / `corr.du`.  The ED is the
    # quadratic (BdG) diagonalization of ImpurityModelBase, with the initial
    # state |0⟩⟨0| ⊗ ρ_bath(β) corresponding to the impurity cdm occupation
    # nsys = 0.
    @testset "BCS vs ED" begin
        for (V, Δ) in ((0.25, 0.6), (0.25, 0.2 + 0.45im))
            spec = DiracDelta(ω=1.0, α=V^2)
            t, N = 1.5, 8
            ed = bcs_ed_imb(β, ϵ_d, [1.0], [V], Δ, 0.0, [t]; μ=μ)[1]
            ρb, Zb = quapi_rt_bcs(N, β, ϵ_d, t, spec, Δ, ρ0, μ)
            @test real(ρb[2, 2] + ρb[4, 4]) ≈ ed[1] atol=0.02   # ⟨n_↑⟩
            @test real(ρb[3, 3] + ρb[4, 4]) ≈ ed[2] atol=0.02   # ⟨n_↓⟩
            # convergence: finer N reduces the discretization error
            errs = Float64[]
            for Nc in (4, 8, 16)
                ρc, _ = quapi_rt_bcs(Nc, β, ϵ_d, t, spec, Δ, ρ0, μ)
                push!(errs, abs(real(ρc[2, 2] + ρc[4, 4]) - ed[1]))
            end
            @test errs[end] < errs[1]
        end
    end
end
