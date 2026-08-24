println("------------------------------------")
println("|  Imaginary-time QuAPI (BCS bath) |")
println("------------------------------------")

using Test, LinearAlgebra

# Dense Grassmann-tensor helpers shared with imagquapi.jl (dense coefficient
# array over 2^L monomials instead of the GrassmannTensors `hadamardproduct`
# fusion tree, which was the dominant test-suite cost).  Guarded so that this
# file can also be included after imagquapi.jl without redefining functions.
if !isdefined(@__MODULE__, :_gt_from_coeffs)
    function _gt_from_coeffs(coeffs::AbstractMatrix, L::Int)
        T = zeros(ComplexF64, 1 << L)
        T[1] = 1.0
        for i in 1:L, j in 1:L
            c = coeffs[i, j]
            iszero(c) && continue
            a, b = minmax(i, j)
            cc = i < j ? c : -c
            bita, bitb = 1 << (a - 1), 1 << (b - 1)
            # in-place fold of the factor 1 + cc·ξₐξ_b; descending target order
            # keeps the source T[mask] unmodified within this factor (m < m2)
            for m2 in (1 << L - 1):-1:0
                (m2 & bita) == 0 || continue
                (m2 & bitb) == 0 || continue
                m = m2 | bita | bitb
                # fermionic exchange sign: move ξₐξ_b past the odd legs below a, b
                sgn = isodd(count_ones(m2 & (bita - 1)) + count_ones(m2 & (bitb - 1))) ? -1.0 : 1.0
                T[m + 1] += cc * sgn * T[m2 + 1]
            end
        end
        return T
    end

    function _gt_coef(T, L, tup::Tuple)
        mask = (1 << L) - 1
        for k in tup
            mask &= ~(1 << (k - 1))
        end
        return T[mask + 1]
    end
end

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
normal, Δ↑↓/Δ↓↑ the anomalous influence-functional blocks from `bcs_Δτ`.  It is
built as a Grassmann tensor over 4N legs — flavor ↑ on legs 2i-1 (ā), 2i (a),
flavor ↓ on legs 2N+2i-1 (ā), 2N+2i (a) — by multiplying every pairwise factor
1 + c·ξᵢξⱼ with the `hadamardproduct` (@hadamard) interface.  Z and the Matsubara
Green's function in the standard convention G(τ_m) = -⟨d↑(τ_m)d↑†(0)⟩ are the
Berezin integrals read off the Grassmann tensor coefficients

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
    Δup = Matrix(corr.uu.data)   # normal ↑↑ (ā↑ Δ a↑)
    Δdn = Matrix(corr.dd.data)   # normal ↓↓ (a↓ Δ ā↓)
    Δud = Matrix(corr.ud.data)   # anomalous ↑↓ (ā↑ Δ ā↓)
    Δdu = Matrix(corr.du.data)   # anomalous ↓↑ (a↓ Δ a↑)

    # coefficient matrix of the bilinear exponent (see docstring): legs
    # ā↑ᵢ = 2i-1, a↑ᵢ = 2i, ā↓ᵢ = 2N+2i-1, a↓ᵢ = 2N+2i.  Coincident pairs
    # accumulate with +=, matching the multiplicative merge of the individual
    # factors in the original tensor formulation.
    L = 4N
    M = zeros(ComplexF64, L, L)
    # ↑ block: kinetic M_imp + normal Δ↑↑
    for i in 1:N                                  # diagonal: (M_imp)ᵢᵢ + Δ↑↑ᵢᵢ
        M[2i - 1, 2i] += Δup[i, i] - 1.0
    end
    for i in 2:N                                  # forward kinetic: x ā↑ᵢ a↑ᵢ₋₁
        M[2i - 1, 2i - 2] += x
    end
    M[1, 2N] += -x                                # anti-periodic wrap: -x ā↑₁ a↑ₙ
    for i in 1:N, j in 1:N                        # normal influence off-diagonal
        i != j && (M[2i - 1, 2j] += Δup[i, j])
    end
    # ↓ block: kinetic M_imp + normal Δ↓↓ (legs shifted by 2N)
    for i in 1:N
        M[2N + 2i - 1, 2N + 2i] += Δdn[i, i] - 1.0
    end
    for i in 2:N
        M[2N + 2i - 1, 2N + 2i - 2] += x
    end
    M[2N + 1, 4N] += -x
    for i in 1:N, j in 1:N
        i != j && (M[2N + 2i - 1, 2N + 2j] += Δdn[i, j])
    end
    # anomalous blocks: ā↑_i Δ↑↓_ij ā↓_j and a↓_i Δ↓↑_ij a↑_j
    for i in 1:N, j in 1:N
        M[2i - 1, 2N + 2j - 1] += Δud[i, j]
    end
    for i in 1:N, j in 1:N
        M[2N + 2i, 2j] += Δdu[i, j]
    end

    T = _gt_from_coeffs(M, L)
    ctop = T[1 << L]                              # all-odd top coefficient
    c = zeros(ComplexF64, N)
    for m in 0:N-1
        c[m+1] = _gt_coef(T, L, (1, 2m + 2))
    end
    return (-1)^(2N) * ctop, c ./ ctop, δτ
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
