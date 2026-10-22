println("------------------------------------")
println("|  Imaginary-time QuAPI (BCS bath) |")
println("------------------------------------")

using Test, LinearAlgebra

# The influence functional — both flavors and the anomalous pairing blocks —
# comes from the shared `QuAPI.DenseFermionicIF` machinery
# (`QuAPI._if_coefficients` in the GTEMPO 2-flavor leg order); only the
# impurity kinetic terms and the Berezin measures are added here.

"""
    quapi_bcs_gt(β, ϵ_d, spectrum::AbstractBoundedFunction; N, Δ, μ=0.0)

Imaginary-time QuAPI for the Toulouse (resonant-level) model with a BCS bath,

    H = Σ_σ ϵ_d d_σ†d_σ + Σ_kσ ω_k c_kσ†c_kσ + Σ_kσ V_k (d_σ†c_kσ + c_kσ†d_σ)
        + Δ Σ_k (c_k↑†c_k↓† + c_k↓c_k↑),

whose impurity (two spin flavors ↑, ↓) couples to the bath only through the
normal part while the pairing gap Δ enters the bath propagators.  Following
src/bcs/bcscorr.jl the discretized path weight is

    exp[ā↑(M_imp+Δ↑↑)a↑ + ā↓(M_imp+Δ↓↓)a↓ + ā↑Δ↑↓ā↓ + a↓Δ↓↑a↑]

with M_imp the anti-periodic kinetic matrix (x = e^{-δτϵ_d}) and Δ↑↑/Δ↓↓ the
normal, Δ↑↓/Δ↓↑ the anomalous influence-functional blocks from `bcs_Δτ`.  The
Δ part is assembled by `QuAPI._if_coefficients` in the GTEMPO 2-flavor leg
order (per node: a↑, ā↑, a↓, ā↓ with band = 1, 2), the impurity terms are
added to the same exponent matrix, and the tensor is folded once with
`_fold_grassmann`.  Z and the Matsubara Green's function in the standard
convention G(τ_m) = -⟨d↑(τ_m)d↑†(0)⟩ are the Berezin integrals read off the
coefficients

    Z      = (-1)^(2N) · (all-odd top coefficient)
    G(τ_m) = (coefficient with legs ā↑₀, a↑_m even) / (top coefficient).

Returns (Z, G, δτ) with τ_m = m·δτ.  ImpurityModelBase implements no analytic
solution for BCS baths, so the reference is the strict exact-diagonalization
Green's function `toulouse_Gτ` of the `Toulouse` class, which uses the
raw-expectation convention G = +⟨d(τ)d†(0)⟩, i.e. the negative of the
convention returned here.
"""
function quapi_bcs_gt(β, ϵ_d, spectrum::AbstractBoundedFunction; N, Δ, μ=0.0)
    δτ = β / N
    x = exp(-δτ * ϵ_d)
    corr = bcs_Δτ(spectrum, β=β, N=N, μ=μ, Δ=Δ)

    # influence functional (normal + anomalous blocks, both flavors)
    C = QuAPI._if_coefficients(corr)
    basis = QuAPI._if_legs(corr)
    li = Dict(l => k for (k, l) in enumerate(basis))
    leg(j, band, c) = li[GrassmannContourIndex(j, :τ, band, c)]
    for band in (1, 2)
        for i in 1:N                              # Berezin measure −āᵢaᵢ
            C[leg(i, band, true), leg(i, band, false)] += -1.0
        end
        for i in 2:N                              # forward kinetic x āᵢaᵢ₋₁
            C[leg(i, band, true), leg(i - 1, band, false)] += x
        end
        C[leg(1, band, true), leg(N, band, false)] += -x    # anti-periodic wrap
    end

    data = QuAPI._fold_grassmann(C)
    L = length(basis)
    ctop = data[1 << L]                           # all-odd top coefficient
    # Z = +ctop in this ordering (the (-1)^(2N) measure sign and the
    # τ-descending leg reversal cancel); the G numerator monomial
    # {ā↑₁, a↑_{m+1}} changes order under the same relabelling, hence the
    # extra minus sign.
    c = [QuAPI._grassmann_coefficient(data, leg(1, 1, true), leg(m + 1, 1, false))
         for m in 0:N-1]
    return ctop, -c ./ ctop, δτ
end

@testset "imaginary-time QuAPI (BCS bath, Grassmann tensor) vs Toulouse exact diagonalization" begin
    β, ϵ_d, μ = 2.0, 0.5, 0.0

    # Free-impurity limit (λ = 0, Δ = 0): the two spin flavors decouple and the
    # Grassmann-tensor ADT must reproduce the exact results to machine precision,
    # pinning down the extraction convention (including the (-1)^(2N) factor).
    @testset "free impurity (λ=0, Δ=0)" begin
        for N in (2, 3)   # odd and even N; N ≤ 3 keeps the test fast (4N legs)
            Z, G, δτ = quapi_bcs_gt(β, ϵ_d, DiracDelta(ω=1.0, α=0.0); N=N, Δ=0.0, μ=μ)
            @test Z ≈ (1 + exp(-β * ϵ_d))^2 atol=1.0e-10
            τs = [(m) * δτ for m in 0:N-1]
            Gan = [-exp(-τ * ϵ_d) / (1 + exp(-β * ϵ_d)) for τ in τs]
            @test G ≈ Gan atol=1.0e-10
        end
    end

    # Each case is (label, ωs, fs, Δ): bath spectrum Σ_k fs[k] δ(ε-ωs[k]) with
    # pairing gap Δ.  ImpurityModelBase has no analytic solution for BCS baths,
    # so the reference is the strict exact-diagonalization Green's function of
    # the `Toulouse` class (free-fermion BdG diagonalization of the 2L-mode
    # Hamiltonian, `toulouse_Gτ`), evaluated at the same imaginary times.
    cases = (("1 mode, real Δ", [1.0], [0.25], 0.6),
             ("1 mode, complex Δ", [1.0], [0.25], 0.2 + 0.45im),
             ("2 modes, real Δ", [1.0, 2.0], [0.25, 0.09], 0.5))
    for (label, ωs, fs, Δ) in cases
        @testset "$label" begin
            Ncoarse, Nfine = 2, 3   # N ≤ 3 keeps the test fast (4N legs)
            spec = DiscreteSpectrum(ωs, fs)
            bath = DiscreteBCSBath(ωs, fs, Δ; β=β, μ=0)
            m = Toulouse(bath; ϵ_d=ϵ_d)

            Z, G, δτ = quapi_bcs_gt(β, ϵ_d, spec; N=Nfine, Δ=Δ, μ=μ)
            τs = [(k) * δτ for k in 0:Nfine-1]
            # toulouse_Gτ uses the raw-expectation convention G = +⟨d(τ)d†(0)⟩,
            # the negative of the standard convention returned by quapi_bcs_gt.
            Gan = -real.(toulouse_Gτ(m, τs))
            # N = 3 (δτ = β/3) keeps the discretization error below ~0.065 for
            # all three cases, so the bound is set to 0.08.
            @test maximum(abs.(real.(G) .- Gan)) < 0.08

            # convergence: finer discretization (smaller δτ) ⇒ smaller error
            _, Gc, _ = quapi_bcs_gt(β, ϵ_d, spec; N=Ncoarse, Δ=Δ, μ=μ)
            τsc = [(k) * β / Ncoarse for k in 0:Ncoarse-1]
            Ganc = -real.(toulouse_Gτ(m, τsc))
            @test maximum(abs.(real.(G) .- Gan)) < maximum(abs.(real.(Gc) .- Ganc))
        end
    end
end
