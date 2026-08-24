println("------------------------------------")
println("|      Imaginary-time QuAPI        |")
println("------------------------------------")

using Test, LinearAlgebra

"""
    quapi_adt(β, ε, Δ, ωs, λs; N)

Imaginary-time QuAPI for the spin-boson model

    H = ε/2 σz + Δ/2 σx + Σ_k ω_k a_k†a_k + (σz/2) Σ_k λ_k (a_k† + a_k)

The Feynman-Vernon influence functional is discretized on N time slices
(δτ = β/N) and the full sum over the 2^N spin paths is assembled into an
augmented density tensor (ADT) of rank N and physical dimension 2.  Returns
the partition function Z, the Matsubara Green's function C[m+1] = ⟨σz(τ_m)σz(0)⟩
at τ_m = m δτ, and δτ.

Sign convention: following src/QuAPI.jl, the minus sign in the influence
functional exponent I[s] = exp[-∫∫ s Δ s] is absorbed into the hybridization
matrix η, so the discretized path weight is exp[+sᵀ η s].
"""
function quapi_adt(β, ε, Δ, ωs, λs; N)
    σx = [0.0 1.0; 1.0 0.0]
    σz = [1.0 0.0; 0.0 -1.0]
    δτ = β / N

    # discretized bath correlation η[j+1, k+1] = Δ_{jk}
    corr = bosonic_Δτ(DiscreteSpectrum(ωs, λs .^ 2), β=β, N=N, μ=0)
    η = corr.data

    # short-time propagator of the isolated spin
    Hs = ε/2 * σz + Δ/2 * σx
    W = exp(-δτ * Hs)

    # enumerate the 2^N spin paths s_j ∈ {±1/2}
    masks = UInt32.(0:2^N-1)
    s = Matrix{Float64}(undef, N, 2^N)
    for j in 0:N-1
        bit = (masks .>> j) .& 1
        s[j+1, :] .= ifelse.(bit .== 0, 0.5, -0.5)
    end

    # influence functional I[s] = exp[+sᵀηs]
    IF = exp.(sum(s .* (η * s), dims=1))

    # product of short-time spin propagators along the path (periodic s_N = s_0)
    prodW = ones(2^N)
    for j in 0:N-1
        i1 = (masks .>> j) .& 1 .+ 1
        i2 = (masks .>> ((j + 1) % N)) .& 1 .+ 1
        prodW .*= [W[i1[m], i2[m]] for m in 1:2^N]
    end

    # augmented density tensor: rank-N, physical dimension 2 per site
    ADT = reshape(IF[:] .* prodW, ntuple(_ -> 2, N))
    Z = sum(ADT)

    # ⟨σz(τ_m)σz(0)⟩ obtained by contracting the ADT with σz at sites m and 0
    σ0 = ifelse.((masks .& 1) .== 0, 1.0, -1.0)     # σz(s_0)
    C = zeros(N)
    for m in 0:N-1
        σm = ifelse.((masks .>> m) .& 1 .== 0, 1.0, -1.0)
        C[m+1] = sum(ADT .* reshape(σm .* σ0, size(ADT))) / Z
    end
    return Z, C, δτ
end

"""
    ed_hamiltonian(ε, Δ, ωs, λs; d)

Full spin-boson Hamiltonian on the truncated Fock space (d boson levels per
mode), to be diagonalized exactly as a reference.
"""
function ed_hamiltonian(ε, Δ, ωs, λs; d)
    σx = [0.0 1.0; 1.0 0.0]
    σz = [1.0 0.0; 0.0 -1.0]
    I2 = Matrix{Float64}(I, 2, 2)
    K = length(ωs)

    a  = [bosonaoperator(K, k; d=d) for k in 1:K]
    ad = [bosonadagoperator(K, k; d=d) for k in 1:K]
    n  = [bosondensityoperator(K, k; d=d) for k in 1:K]
    Ib = Matrix{Float64}(I, d^K, d^K)

    Hs = ε/2 * σz + Δ/2 * σx
    Hb = sum(ωs[k] * n[k] for k in 1:K)
    Hc = kron(σz/2, sum(λs[k] * (ad[k] + a[k]) for k in 1:K))
    H = kron(Hs, Ib) + kron(I2, Hb) + Hc
    return (H + H') / 2
end

"""
    quapi_ibm(β, ε, ωs, λs; N)

Independent-boson limit of the spin-boson model (tunneling Δ = 0) solved with
the discretized QuAPI influence functional.  Since σz is conserved when Δ = 0,
the only non-trivial Matsubara correlator is the spin-flip (fermionic) one

    G(τ_m) = ⟨d(τ_m) d†(0)⟩,   d = σ_-,  d† = σ_+.

It is given by the single-flip path integral: the spin takes +1/2 on [0, τ_m)
and -1/2 on [τ_m, β) (one domain wall), with weight

    W[1,1]^m W[2,2]^(N-m) × exp[+sᵀηs] / Z,

where W = exp(-δτ (ε/2) σz) is the diagonal short-time spin propagator and Z
receives contributions only from the two constant spin paths (since W is
diagonal).  Returns Z, G and δτ.

Note on normalization: exp[+sᵀηs] is built with the free-bath partition
function factored out (Ση = βΣλ²/ω exactly), so the returned Z is the true
partition function divided by Z_bath; the correlators (ratios) are unaffected.
"""
function quapi_ibm(β, ε, ωs, λs; N)
    δτ = β / N
    σz = [1.0 0.0; 0.0 -1.0]
    corr = bosonic_Δτ(DiscreteSpectrum(ωs, λs .^ 2), β=β, N=N, μ=0)
    η = corr.data
    W = exp(-δτ * (ε/2) * σz)                    # diagonal (Δ = 0)
    # partition function: only the constant spin paths contribute
    Z = (W[1,1]^N + W[2,2]^N) * exp(sum(η) / 4)
    # single-flip path: s_j = +1/2 for j < m, -1/2 for j ≥ m
    G = zeros(N)
    for m in 0:N-1
        sp = [(j < m ? 0.5 : -0.5) for j in 0:N-1]
        G[m+1] = W[1,1]^m * W[2,2]^(N-m) * exp(dot(sp, η, sp)) / Z
    end
    return Z, G, δτ
end

# ---------------------------------------------------------------------------
# Fermionic QuAPI (Toulouse / resonant-level model) with a Grassmann-tensor ADT
# ---------------------------------------------------------------------------
# Dense rank-L Grassmann tensor: T[mask+1] = coefficient of the monomial whose
# odd legs are given by the L-bit word `mask` (bit k set ⇔ leg k+1 odd).  This
# avoids the GrassmannTensors `hadamardproduct` fusion-tree machinery, which is
# O(3^L) per pairwise factor and was the dominant test-suite cost.
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

"""
    quapi_rlm_gt(β, ϵ_d, spectrum::AbstractBoundedFunction; N, μ=0.0)

Imaginary-time QuAPI for the Toulouse (resonant-level) model

    H = ϵ_d d†d + Σ_k ω_k c_k†c_k + Σ_k V_k (d†c_k + c_k†d)

with the discretized path weight exp[Σᵢⱼ āᵢ (M_imp + Δ)ᵢⱼ aⱼ] (M_imp the
anti-periodic kinetic matrix with x = e^{-δτϵ_d}, Δ the fermionic influence-
functional matrix from `fermionic_Δτ`) built as a Grassmann tensor: each
pairwise factor 1 + c·āᵢaⱼ is multiplied with the `hadamardproduct` (@hadamard)
interface of the GrassmannTensors package.  Z and the Matsubara Green's function
in the standard convention G(τ_m) = -⟨d(τ_m)d†(0)⟩ are the Berezin integrals
read off the Grassmann tensor coefficients

    Z      = (-1)^N · (all-odd top coefficient)
    G(τ_m) = (coefficient with legs ā₀, a_m even) / (top coefficient).

Returns (Z, G, δτ) with τ_m = m·δτ.  Note that `toulouse_Gτ` from
ImpurityModelBase uses the raw-expectation convention G = +⟨d(τ)d†(0)⟩, i.e.
the negative of the convention returned here.
"""
function quapi_rlm_gt(β, ϵ_d, spectrum::AbstractBoundedFunction; N, μ=0.0)
    δτ = β / N
    x = exp(-δτ * ϵ_d)
    corr = fermionic_Δτ(spectrum, β=β, N=N, μ=μ)
    Δ = Matrix(corr.data)

    # coefficient matrix of the bilinear exponent: e^{Σ āᵢ Mᵢⱼ aⱼ}
    # = Πᵢⱼ (1 + Mᵢⱼ āᵢaⱼ) by Grassmann nilpotency; legs āᵢ = 2i-1, aⱼ = 2j.
    # Coincident pairs (e.g. ā₁aₙ from the anti-periodic wrap and from the
    # influence functional) accumulate with +=, matching the multiplicative
    # merge of the individual factors in the original tensor formulation.
    L = 2N
    M = zeros(ComplexF64, L, L)
    for i in 1:N                                  # diagonal: (M_imp)ᵢᵢ + Δᵢᵢ
        M[2i - 1, 2i] += Δ[i, i] - 1.0
    end
    for i in 2:N                                  # forward kinetic: x āᵢ aᵢ₋₁
        M[2i - 1, 2i - 2] += x
    end
    M[1, 2N] += -x                                # anti-periodic wrap: -x ā₁ aₙ
    for i in 1:N, j in 1:N                        # influence functional off-diagonal
        i != j && (M[2i - 1, 2j] += Δ[i, j])
    end

    T = _gt_from_coeffs(M, L)
    ctop = T[1 << L]                              # all-odd top coefficient
    c = zeros(Float64, N)
    for m in 0:N-1
        c[m+1] = _gt_coef(T, L, (1, 2m + 2))
    end
    return (-1)^N * ctop, c ./ ctop, δτ
end

# Convenience wrapper for a discrete bath spectrum Σ_k λ_k² δ(ε-ω_k).
quapi_rlm_gt(β, ϵ_d, ωs::AbstractVector{<:Real}, λs::AbstractVector{<:Real}; N, μ=0.0) =
    quapi_rlm_gt(β, ϵ_d, DiscreteSpectrum(collect(ωs), collect(λs .^ 2)); N=N, μ=μ)

@testset "imaginary-time QuAPI (spin-boson) vs exact diagonalization" begin
    β, ε, Δ = 2.0, 0.5, 0.5
    d = 12   # boson Fock cutoff for ED
    cases = (("1 mode", [1.0], [0.5]),
             ("2 modes", [1.0, 2.0], [0.5, 0.3]))

    for (label, ωs, λs) in cases
        @testset "$label" begin
            K = length(ωs)
            H = ed_hamiltonian(ε, Δ, ωs, λs; d=d)
            cache = eigencache(H)
            σz_full = kron([1.0 0.0; 0.0 -1.0], Matrix{Float64}(I, d^K, d^K))

            # ED reference at τ = β/2 (used for the convergence check)
            ged_mid = correlation_2op_1τ(H, σz_full, σz_full, [β/2], cache; β=β, reverse=false)[1]

            # QuAPI on the fine grid: full curve vs ED
            N = 12
            Z, C, δτ = quapi_adt(β, ε, Δ, ωs, λs; N=N)
            @test Z > 0
            @test C[1] ≈ 1.0 atol=1.0e-10    # ⟨σz(0)σz(0)⟩ = 1 by construction

            τs = [(m) * δτ for m in 0:N-1]
            ged = correlation_2op_1τ(H, σz_full, σz_full, τs, cache; β=β, reverse=false)
            @test C ≈ ged atol=2.0e-3

            # convergence: finer discretization (smaller δτ) ⇒ smaller error
            Ncoarse = 6
            _, Ccoarse, _ = quapi_adt(β, ε, Δ, ωs, λs; N=Ncoarse)
            err_coarse = abs(Ccoarse[Ncoarse÷2+1] - ged_mid)
            err_fine   = abs(C[N÷2+1] - ged_mid)
            @test err_fine < err_coarse
            @test err_fine < 1.0e-3
        end
    end
end

@testset "independent-boson limit (Δ=0) vs ImpurityModelBase analytic solution" begin
    β, ε = 2.0, 0.5
    d = 12   # boson Fock cutoff for ED
    cases = (("1 mode", [1.0], [0.5]),
             ("2 modes", [1.0, 2.0], [0.5, 0.3]))

    for (label, ωs, λs) in cases
        @testset "$label" begin
            N = 16
            Z, G, δτ = quapi_ibm(β, ε, ωs, λs; N=N)
            τs = [(m) * δτ for m in 0:N-1]

            # ImpurityModelBase analytic solution for the independent-boson model.
            # Our spin-boson Hamiltonian maps onto a fermionic impurity (n=(σz+1)/2)
            # with level ϵ_d = ε + Δ, where Δ = Σλ²/ω is the polaron shift absorbed
            # inside the analytic solution, so that its renormalized level ϵ_d - Δ = ε
            # equals the bare level of the spin-boson Hamiltonian.
            Δ = sum(λs .^ 2 ./ ωs)
            spec = DiscreteSpectrum(ωs, λs .^ 2)
            Gan = [independentbosons_Gτ(spec, τ, β=β, ϵ_d=ε + Δ, bands=1) for τ in τs]
            @test G ≈ Gan atol=1.0e-8

            # cross-check against exact diagonalization of the same spin-boson
            # Hamiltonian (Δ_tun = 0), i.e. ⟨d(τ)d†(0)⟩
            H = ed_hamiltonian(ε, 0.0, ωs, λs; d=d)
            cache = eigencache(H)
            K = length(ωs)
            σm = [0.0 0.0; 1.0 0.0]
            σp = [0.0 1.0; 0.0 0.0]
            d_full  = kron(σm, Matrix{Float64}(I, d^K, d^K))
            dp_full = kron(σp, Matrix{Float64}(I, d^K, d^K))
            Ged = correlation_2op_1τ(H, d_full, dp_full, τs, cache; β=β, reverse=false)
            @test G ≈ Ged atol=1.0e-8

            # partition function with the free-bath Z factored out (see quapi_ibm)
            @test Z ≈ 2 * cosh(β*ε/2) * exp(β * sum(λs .^ 2 ./ (4 .* ωs))) atol=1.0e-8
        end
    end
end

@testset "imaginary-time QuAPI (Toulouse, Grassmann tensor) vs ImpurityModelBase analytic solution" begin
    β, ϵ_d, μ = 2.0, 0.5, 0.0

    # Free-impurity limit (λ = 0): the Grassmann-tensor ADT must reproduce the
    # exact results to machine precision, pinning down the extraction convention
    # (including the (-1)^N factor of the Berezin integral in Z).
    @testset "free impurity (λ=0)" begin
        for N in (3, 4)   # odd and even N check the (-1)^N sign; N ≤ 4 keeps the test fast (2N legs)
            Z, G, δτ = quapi_rlm_gt(β, ϵ_d, [0.0], [0.0]; N=N, μ=μ)
            @test Z ≈ 1 + exp(-β * ϵ_d) atol=1.0e-10
            τs = [(m) * δτ for m in 0:N-1]
            Gan = [-exp(-τ * ϵ_d) / (1 + exp(-β * ϵ_d)) for τ in τs]
            @test G ≈ Gan atol=1.0e-10
        end
    end

    # Each case is (label, spectrum, reference Matsubara cutoff, fine N): the
    # discrete spectra sum poles directly and are cheap even at n = 4000, while
    # the semicircular continuous spectrum needs one quadgk per Matsubara
    # frequency, so a smaller n = 1000 is used (converged to ~1e-5 there).  The
    # Grassmann-tensor ADT cost grows sharply with N (2N legs), so the fine N is
    # chosen as small as the 0.05 accuracy threshold allows: N = 5 for the
    # discrete spectra, N = 6 for the semicircular one.
    cases = (("1 mode", DiscreteSpectrum([1.0], [0.5] .^ 2), 4000, 5),
             ("2 modes", DiscreteSpectrum([1.0, 2.0], [0.5, 0.3] .^ 2), 4000, 5),
             ("semicircular", bounded(ϵ -> 0.4 * semicircular(1.0)(ϵ), -1.0, 1.0), 1000, 6))
    for (label, spec, nref, Nfine) in cases
        @testset "$label" begin
            Ncoarse = 3   # coarse N for the convergence check

            Z, G, δτ = quapi_rlm_gt(β, ϵ_d, spec; N=Nfine, μ=μ)
            τs = [(m) * δτ for m in 0:Nfine-1]
            # ImpurityModelBase's toulouse_Gτ uses the raw-expectation convention
            # G = +⟨d(τ)d†(0)⟩, the negative of the standard convention returned
            # by quapi_rlm_gt, so compare with the sign flipped.
            Gan = [-real(toulouse_Gτ(spec, τ; β=β, ϵ_d=ϵ_d, μ=μ, n=nref)) for τ in τs]
            @test maximum(abs.(G .- Gan)) < 0.05

            # convergence: finer discretization (smaller δτ) ⇒ smaller error
            _, Gc, _ = quapi_rlm_gt(β, ϵ_d, spec; N=Ncoarse, μ=μ)
            τsc = [(m) * β / Ncoarse for m in 0:Ncoarse-1]
            Ganc = [-real(toulouse_Gτ(spec, τ; β=β, ϵ_d=ϵ_d, μ=μ, n=nref)) for τ in τsc]
            @test maximum(abs.(G .- Gan)) < maximum(abs.(Gc .- Ganc))
        end
    end
end
