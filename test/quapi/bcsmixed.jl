println("------------------------------------")
println("|   Mixed-time QuAPI (BCS bath)    |")
println("------------------------------------")

using Test, LinearAlgebra

# =====================================================================
# Mixed-time (Kadanoff / L-shaped) QuAPI for the Toulouse model
#
#     H = Σ_σ ϵ_d d_σ†d_σ + Σ_kσ ω_k c_kσ†c_kσ
#         + Σ_kσ V_k (d_σ†c_kσ + c_kσ†d_σ)
#         + Δ Σ_k (c_k↑†c_k↓† + c_k↓c_k↑) ,
#
# starting from the thermal state.  The contour is the L-shaped Kadanoff
# contour: imaginary branch τ = 0 → β glued to the closed Keldysh part
# t = 0 → t → 0.  Following the validated spin-boson mixed-time path sum
# (test/quapi/mixed.jl):
#
#   - real + / − nodes j = 0..Nt (t_j = jδt), imaginary nodes
#     τ_i = (i−1)δτ, with τ_1 identified with the −0 node (τ=0);
#   - the imaginary chain runs from −0 through τ_2..τ_Nτ and closes onto
#     +0 with the step e^{−δτϵ_d};
#   - the Keldysh turn-around at t closes as −ā⁻(t)a⁺(t).
#
# The bath influence functional is the exact Grassmann bilinear of
# `bcs_Δm` (all nine branch blocks; `GrassmannHybridizationMatrix`).
# The same Gaussian integral can be evaluated in two ways:
#   - dense: fold the coefficient matrix into a `GrassmannADT` and read
#     `Zvalue` / `complementindex` (feasible for ~20 legs);
#   - Pfaffian: O(L³), needed for the larger contours below.
# The thermal Green's function is read off as
#
#     G(τ_m) = −⟨d(τ_m)d†(0)⟩ = A_{ā(−0), a(τ_m)} / Z ,
#
# m = 1..Nτ−1 (the τ=0 entry is the contact term).
# =====================================================================

# Mixed-time path integral for a `bands`-flavor correlation function:
# returns (Z, G(τ_m), δτ), τ_m = mδτ, m = 1..Nτ−1.
function _mixed_quapi(corr, β, ϵ_d; Nt, Nτ, t, bands, dense::Bool=false)
    δt = t / Nt
    δτ = β / Nτ
    x = exp(-im * δt * ϵ_d)
    xc = conj(x)
    xτ = exp(-δτ * ϵ_d)

    h = GrassmannHybridizationMatrix(corr)
    legs = h.basis
    li = Dict(l => k for (k, l) in enumerate(legs))
    C0 = h.Δ
    L0 = length(legs)

    # τ_1 (τ=0) is the same contour point as −0
    id = collect(1:L0)
    for band in 1:bands, c in (false, true)
        id[li[GrassmannContourIndex(1, :τ, band, c)]] =
            id[li[GrassmannContourIndex(1, :-, band, c)]]
    end
    uniq = sort(unique(id))
    ren = Dict(u => k for (k, u) in enumerate(uniq))
    L = length(uniq)
    # representatives of the merged legs, used by the dense ADT construction
    merged_legs = GrassmannContourIndex[]
    if dense
        merged_legs = Vector{GrassmannContourIndex}(undef, L)
        for l in legs
            u = ren[id[li[l]]]
            isassigned(merged_legs, u) || (merged_legs[u] = l)
        end
    end
    C = zeros(ComplexF64, L, L)
    for i in 1:L0, j in 1:L0
        c = C0[i, j]
        iszero(c) && continue
        p, q = ren[id[i]], ren[id[j]]
        p == q && continue          # ξ² = 0 for coincident legs
        C[p, q] += c
    end
    leg(branch, j, band, c) = ren[id[li[GrassmannContourIndex(j, branch, band, c)]]]

    for band in 1:bands
        # Berezin measures −āa
        for j in 1:Nt+1
            C[leg(:+, j, band, true), leg(:+, j, band, false)] += -1.0
            C[leg(:-, j, band, true), leg(:-, j, band, false)] += -1.0
        end
        for j in 2:Nτ
            C[leg(:τ, j, band, true), leg(:τ, j, band, false)] += -1.0
        end
        # real-time kinetics x ā_{j+1}a_j (forward) and x̄ ā_j a_{j+1} (backward)
        for j in 1:Nt
            C[leg(:+, j + 1, band, true), leg(:+, j, band, false)] += x
            C[leg(:-, j, band, true), leg(:-, j + 1, band, false)] += xc
        end
        C[leg(:-, Nt + 1, band, true), leg(:+, Nt + 1, band, false)] += -1.0
        # imaginary chain −0 → τ_2 → … → τ_Nτ → +0
        for j in 1:Nτ-1
            C[leg(:τ, j + 1, band, true), leg(:τ, j, band, false)] += xτ
        end
        C[leg(:+, 1, band, true), leg(:τ, Nτ, band, false)] += xτ
    end

    if dense
        gif = GrassmannIF(GrassmannHybridizationMatrix(merged_legs, C))
        Z = Zvalue(gif)
    else
        M = zeros(ComplexF64, L, L)
        for i in 1:L, j in 1:L
            i == j && continue
            M[i, j] = C[i, j] - C[j, i]
        end
        Z = pfaffian(M)
        W = inv(M)
    end
    sgn(i, j) = isodd(i + j) ? -1.0 : 1.0
    G = ComplexF64[]
    for m in 1:Nτ-1
        i = leg(:-, 1, 1, true)     # ā at τ=0
        j = leg(:τ, m + 1, 1, false) # a at τ=mδτ (band ↑)
        if dense
            push!(G, gif[complementindex(gif.basis, merged_legs[i], merged_legs[j])...] / Z)
        else
            push!(G, -sgn(i, j) * W[i, j])
        end
    end
    return Z, G, δτ
end

"""
    quapi_bcs_mixed(β, ϵ_d, spectrum; Nt, Nτ, t, Δ, μ=0.0, dense=false)

Mixed-time QuAPI for the Toulouse model with a BCS bath.  Returns
`(Z, G, δτ)` with `G[m] = −⟨d(τ_m)d†(0)⟩` at `τ_m = mδτ` (`m = 1..Nτ−1`).
`dense=true` constructs the `GrassmannADT` directly (small contours only).
"""
function quapi_bcs_mixed(β, ϵ_d, spectrum; Nt, Nτ, t, Δ, μ=0.0, dense::Bool=false)
    corr = bcs_Δm(spectrum, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ, Δ=Δ)
    return _mixed_quapi(corr, β, ϵ_d; Nt=Nt, Nτ=Nτ, t=t, bands=2, dense=dense)
end

# single-flavor version, used for the Δ=0 decoupling check
function quapi_normal_mixed(β, ϵ_d, spectrum; Nt, Nτ, t, μ=0.0, dense::Bool=false)
    corr = fermionic_Δm(spectrum, β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ)
    return _mixed_quapi(corr, β, ϵ_d; Nt=Nt, Nτ=Nτ, t=t, bands=1, dense=dense)
end

@testset "mixed-time QuAPI (BCS bath, Toulouse)" begin
    β, ϵ_d, μ = 2.0, 0.5, 0.0

    # ---- dense ADT construction (Nτ=2, Nt=1: 20 legs after the τ=0 merge) ----
    @testset "dense ADT construction" begin
        spec = DiracDelta(ω=1.0, α=0.25^2)
        Nt, Nτ, t = 1, 2, 0.1
        for Δ in (0.0, 0.6)
            Zd, Gd, _ = quapi_bcs_mixed(β, ϵ_d, spec; Nt=Nt, Nτ=Nτ, t=t,
                                        Δ=Δ, μ=μ, dense=true)
            Zp, Gp, _ = quapi_bcs_mixed(β, ϵ_d, spec; Nt=Nt, Nτ=Nτ, t=t,
                                        Δ=Δ, μ=μ)
            @test Zd ≈ Zp atol=1.0e-9
            @test Gd ≈ Gp atol=1.0e-9
        end
        # free impurity: exact
        Z, G, δτ = quapi_bcs_mixed(β, ϵ_d, DiracDelta(ω=1.0, α=0.0);
                                   Nt=Nt, Nτ=Nτ, t=t, Δ=0.0, μ=μ, dense=true)
        @test Z ≈ (1 + exp(-β * ϵ_d))^2 atol=1.0e-9
        τs = [m * δτ for m in 1:Nτ-1]
        @test G ≈ [exp(-τ * ϵ_d) / (1 + exp(-β * ϵ_d)) for τ in τs] atol=1.0e-9
    end

    # ---- free impurity (V = 0, Δ = 0): exact ----
    @testset "free impurity (λ=0, Δ=0)" begin
        for (Nt, Nτ, t) in ((1, 4, 1.0e-4), (2, 8, 1.0e-3))
            Z, G, δτ = quapi_bcs_mixed(β, ϵ_d, DiracDelta(ω=1.0, α=0.0);
                                       Nt=Nt, Nτ=Nτ, t=t, Δ=0.0, μ=μ)
            @test Z ≈ (1 + exp(-β * ϵ_d))^2 atol=1.0e-9
            τs = [m * δτ for m in 1:Nτ-1]
            Gan = [exp(-τ * ϵ_d) / (1 + exp(-β * ϵ_d)) for τ in τs]
            @test G ≈ Gan atol=1.0e-9
        end
    end

    # ---- Δ = 0: the two spin flavors decouple ----
    @testset "Δ=0 decoupling" begin
        spec = DiracDelta(ω=1.0, α=0.25^2)
        Nt, Nτ = 2, 8
        Zn, Gn, _ = quapi_normal_mixed(β, ϵ_d, spec; Nt=Nt, Nτ=Nτ, t=1.0e-3, μ=μ)
        Zb, Gb, _ = quapi_bcs_mixed(β, ϵ_d, spec; Nt=Nt, Nτ=Nτ, t=1.0e-3, Δ=0.0, μ=μ)
        @test Zb ≈ Zn^2 atol=1.0e-9          # independent flavors
        @test Gb ≈ Gn atol=1.0e-9            # ↑ and ↓ Green's functions coincide
    end

    # ---- coupled BCS bath vs exact diagonalization ----
    @testset "coupled vs ED" begin
        V = 0.25
        spec = DiracDelta(ω=1.0, α=V^2)
        t, Nt, Nτ = 0.1, 4, 16

        for Δ in (0.6, 0.2 + 0.45im)
            Z, G, δτ = quapi_bcs_mixed(β, ϵ_d, spec; Nt=Nt, Nτ=Nτ, t=t, Δ=Δ, μ=μ)
            # Toulouse exact diagonalization (loose: the QuAPI hybridization
            # carries the usual O(δτ²) discretization error, cf. bcsimag.jl)
            bath = DiscreteBCSBath([1.0], [V], Δ; β=β, μ=μ)
            m = Toulouse(bath; ϵ_d=ϵ_d)
            τs = [m * δτ for m in 1:Nτ-1]
            Gan = toulouse_Gτ(m, τs)
            @test maximum(abs.(G .- Gan)) < 0.06
        end
    end
end
