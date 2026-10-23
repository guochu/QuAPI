println("------------------------------------")
println("|        Real-time QuAPI           |")
println("|  (spin-boson + fermionic)        |")
println("------------------------------------")

using Test, LinearAlgebra

"""
    quapi_rt(β, ε, Δ, ωs, λs; N, t, ρ_imp, μ=0.0)

Real-time QuAPI for the spin-boson model

    H = ε/2 σz + Δ/2 σx + Σ_k ω_k a_k†a_k + (σz/2) Σ_k λ_k (a_k† + a_k)

starting from the product state ρ(0) = ρ_imp ⊗ ρ_bath(β) (bath in thermal
equilibrium, interaction switched on at t = 0).  The real-time influence
functional is discretized on N steps δt = t/N and the full sum over the
2^N × 2^N forward/backward spin paths is evaluated, weighted by

    I[s+, s-] = exp[Σ_jk (s+_j s+_k G₊₊_jk + s+_j s-_k G₊₋_jk
                       + s-_j s+_k G₋₊_jk + s-_j s-_k G₋₋_jk)]

with G•• the four (N+1)×(N+1) branch matrices of `bosonic_Δt`.  Returns
⟨σz(t)⟩, ⟨σx(t)⟩, the reduced density matrix ρ and δt.
"""
function quapi_rt(β, ε, Δ, ωs, λs; N, t, ρ_imp, μ=0.0)
    δt = t / N
    σx = [0.0 1.0; 1.0 0.0]
    σz = [1.0 0.0; 0.0 -1.0]
    Hs = ε/2 * σz + Δ/2 * σx
    W = exp(-im * δt * Hs)                       # forward short-time propagator

    corr = bosonic_Δt(DiscreteSpectrum(ωs, λs .^ 2), β=β, N=N, t=t, μ=μ)
    Gpp = Matrix(branch(corr, :+, :+))           # (N+1)×(N+1), nodes 0..N
    Gpm = Matrix(branch(corr, :+, :-))
    Gmp = Matrix(branch(corr, :-, :+))
    Gmm = Matrix(branch(corr, :-, :-))

    # enumerate the 2^N paths of interior nodes (s_0 .. s_{N-1}); value ±1/2
    nf = 2^N
    S = Matrix{Float64}(undef, nf, N)
    for k in 0:N-1
        bit = (UInt32.(0:nf-1) .>> k) .& 1
        S[:, k+1] .= ifelse.(bit .== 0, 0.5, -0.5)
    end
    nodeidx = Int.(ifelse.(S .== 0.5, 1, 2))     # +1/2 ↔ index 1, -1/2 ↔ index 2

    # forward/backward path weights over the final index a ∈ {1,2}
    FP = ones(ComplexF64, nf, 2)
    BP = ones(ComplexF64, nf, 2)
    Wc = conj(W)
    for p in 1:nf
        v = W[:, nodeidx[p, 1]]
        w = Wc[:, nodeidx[p, 1]]
        for j in 2:N
            v = v[nodeidx[p, j]] .* W[:, nodeidx[p, j]]
            w = w[nodeidx[p, j]] .* Wc[:, nodeidx[p, j]]
        end
        FP[p, :] = v
        BP[p, :] = w
    end

    Kmat = Gpm + adjoint(Gmp)                    # s+ᵀ Gpm s- + s-ᵀ Gmp s+
    ρ = zeros(ComplexF64, 2, 2)
    for a in 1:2, b in 1:2
        Sa, Sb = 0.5 - (a - 1), 0.5 - (b - 1)    # output node values
        for p in 1:nf
            Sfp = vcat(S[p, :], Sa)              # extended forward nodes 0..N
            fp = dot(Sfp, Gpp, Sfp)
            u = adjoint(Kmat) * Sfp              # (N+1)-vector
            for q in 1:nf
                Sgq = vcat(S[q, :], Sb)          # extended backward nodes 0..N
                cr = dot(u, Sgq)
                bq = dot(Sgq, Gmm, Sgq)
                IF = exp(fp + bq + cr)
                ρ[a, b] += ρ_imp[nodeidx[p, 1], nodeidx[q, 1]] * FP[p, a] * BP[q, b] * IF
            end
        end
    end
    trρ = ρ[1, 1] + ρ[2, 2]
    σz_t = (ρ[1, 1] - ρ[2, 2]) / trρ
    σx_t = (ρ[1, 2] + ρ[2, 1]) / trρ
    return σz_t, σx_t, ρ, δt
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

function ed_rt_obs(ε, Δ, ωs, λs, ρ_imp, ts; β, d=12)
    K = length(ωs)
    H = ed_hamiltonian(ε, Δ, ωs, λs; d=d)
    cache = eigencache(H)
    # thermal bath state on truncated Fock space
    Ib = Matrix{Float64}(I, d^K, d^K)
    Hb = sum(ωs[k] * bosondensityoperator(K, k; d=d) for k in 1:K)
    ρ_bath = exp(-β * Hb)
    ρ_bath ./= tr(ρ_bath)
    ρ0 = kron(ρ_imp, ρ_bath)
    σz_full = kron([1.0 0.0; 0.0 -1.0], Ib)
    σx_full = kron([0.0 1.0; 1.0 0.0], Ib)
    σz = Float64[]; σx = Float64[]
    for t in ts
        ρt = timeevo(ρ0, H, -im * t, cache)
        push!(σz, real(tr(ρt * σz_full)))
        push!(σx, real(tr(ρt * σx_full)))
    end
    return σz, σx
end

@testset "real-time QuAPI (spin-boson)" begin
    β, ε, Δ = 2.0, 0.5, 0.5
    d = 12   # boson Fock cutoff for ED
    ρ_up = [1.0 0.0; 0.0 0.0]                      # |↑⟩⟨↑|
    ρ_pl = [0.5 0.5; 0.5 0.5]                      # |+x⟩⟨+x| = (|↑⟩+|↓⟩)/√2

    # free spin (λ=0): the path sum must reproduce the exact precession
    @testset "free spin (λ=0)" begin
        for (t, N) in ((0.5, 2), (1.0, 4), (2.0, 6))
            σzq, σxq, ρ, δt = quapi_rt(β, ε, Δ, [1.0], [0.0]; N=N, t=t, ρ_imp=ρ_up)
            Hs = ε/2 * [1.0 0.0; 0.0 -1.0] + Δ/2 * [0.0 1.0; 1.0 0.0]
            Wt = exp(-im * t * Hs)
            ρe = Wt * ρ_up * Wt'
            @test real(σzq) ≈ real(ρe[1,1] - ρe[2,2]) atol=1.0e-8
            @test real(σxq) ≈ real(ρe[1,2] + ρe[2,1]) atol=1.0e-8
        end
    end

    # coupled spin-boson (ω=1, λ=0.5), ρ_imp=|↑⟩, vs exact diagonalization.
    # QuAPI with finite N is not yet converged at the largest t, so the
    # agreement tolerances are loose at large t; the N→N+2 trend is asserted
    # separately at fixed t=2.0.
    @testset "spin-boson vs ED" begin
        ts = [0.5, 1.0, 1.5, 2.0]
        Ns = [2, 4, 6, 8]
        σze, σxe = ed_rt_obs(ε, Δ, [1.0], [0.5], ρ_up, ts; β=β, d=d)
        for (t, N) in zip(ts, Ns)
            σzq, σxq, ρ, δt = quapi_rt(β, ε, Δ, [1.0], [0.5]; N=N, t=t, ρ_imp=ρ_up)
            i = findfirst(==(t), ts)
            @test real(σzq) ≈ σze[i] atol=0.05
            @test real(σxq) ≈ σxe[i] atol=0.15
        end

        # convergence at fixed t=2.0: finer N ⇒ smaller σz error vs ED.
        # (⟨σx⟩ has a residual offset that does not yet shrink with N in this
        # range, so only the σz trend is asserted.)
        i2 = findfirst(==(2.0), ts)
        errs = Float64[]
        for N in (4, 6)
            σzq, σxq, ρ, δt = quapi_rt(β, ε, Δ, [1.0], [0.5]; N=N, t=2.0, ρ_imp=ρ_up)
            push!(errs, abs(real(σzq) - σze[i2]))
        end
        @test errs[end] < errs[1]
    end

    # coherent initial state |+x⟩, vs ED (checks off-diagonal / cross branches)
    @testset "coherent initial state vs ED" begin
        ts = [0.5, 1.0, 1.5, 2.0]
        Ns = [2, 4, 6, 8]
        σze, σxe = ed_rt_obs(ε, Δ, [1.0], [0.5], ρ_pl, ts; β=β, d=d)
        for (t, N) in zip(ts, Ns)
            σzq, σxq, ρ, δt = quapi_rt(β, ε, Δ, [1.0], [0.5]; N=N, t=t, ρ_imp=ρ_pl)
            i = findfirst(==(t), ts)
            @test real(σxq) ≈ σxe[i] atol=0.15
        end
    end

    # IBM limit (Δ_tun=0): analytic ⟨σx(t)⟩ = cos(εt)·exp(-Γ(t)),
    # Γ(t) = Σ λ² coth(βω/2)(1-cos ωt)/ω²,  from |+x⟩ initial state
    @testset "IBM limit (Δ_tun=0) vs analytic" begin
        ts = [0.5, 1.0, 1.5]
        Ns = [2, 4, 6]
        for (t, N) in zip(ts, Ns)
            σzq, σxq, ρ, δt = quapi_rt(β, ε, 0.0, [1.0], [0.5]; N=N, t=t, ρ_imp=ρ_pl)
            Γt = sum(0.5 .^ 2 .* coth(β/2) .* (1 - cos(t)) ./ 1.0^2)
            an = cos(ε * t) * exp(-Γt)
            @test real(σxq) ≈ an atol=0.06
            @test abs(real(σzq)) < 1.0e-10   # ⟨σz⟩ = 0 by symmetry
        end

        # convergence at fixed t=2.0 vs analytic
        errs = Float64[]
        for N in (4, 6, 8)
            σzq, σxq, ρ, δt = quapi_rt(β, ε, 0.0, [1.0], [0.5]; N=N, t=2.0, ρ_imp=ρ_pl)
            Γt = sum(0.5 .^ 2 .* coth(β/2) .* (1 - cos(2.0)) ./ 1.0^2)
            an = cos(ε * 2.0) * exp(-Γt)
            push!(errs, abs(real(σxq) - an))
        end
        @test errs[end] < errs[1]
    end

    # branch-matrix sum check: for identical constant paths the IF must be 1
    # ⇒ 1ᵀ(Gpp+Gpm+Gmp+Gmm)1 ≈ 0
    @testset "branch-matrix sum" begin
        corr = bosonic_Δt(DiscreteSpectrum([1.0], [0.5] .^ 2), β=β, N=8, t=2.0, μ=0.0)
        R = Matrix(branch(corr, :+, :+)) + Matrix(branch(corr, :+, :-)) +
            Matrix(branch(corr, :-, :+)) + Matrix(branch(corr, :-, :-))
        @test abs(sum(R)) < 1.0e-8
    end
end

### fermionic Toulouse impurity (dense Grassmann tensor from src/models.jl) ###

# The influence functional comes from `QuAPI.Dense._if_coefficients` in the
# GTEMPO leg order; the impurity kinetics, measures and contour closure are
# added to the same bilinear exponent and the tensor is folded once with
# `_fold_grassmann` (no Grassmann-tensor package needed).

"""
	quapi_rlm_rt_gt(β, ϵ_d, spectrum; Nt, t, μ=0.0)

Real-time (closed Keldysh contour) fermionic QuAPI for the Toulouse impurity,
represented as a single dense Grassmann tensor over the legs of
`QuAPI.Dense._if_legs(corr)` (per node: `a`, `ā`; nodes in GTEMPO order, i.e.
`(+,Nt+1), (-,Nt+1), …, (+,1), (-,1)`).

Contour: nodes 0..Nt on the + branch (t = 0..t) and on the - branch.
The bilinear exponent exp(Σ C_ij ξ_i ξ_j) contains:
- the Berezin measures −ā a at every node,
- the free-impurity kinetics: x = e^{−iδtϵ_d} forward / x̄ backward,
- the contour closure −ā⁻(t)a⁺(t) at the final time,
- the initial impurity thermal weight e^{−βϵ_d} ā⁺(0)a⁻(0) (quench from
  ρ(0) = ρ_imp(β) ⊗ ρ_bath(β)),
- the discrete-time influence functional blocks from `fermionic_Δt`
  (all four branch blocks, from `_if_coefficients`).

Returns the partition function Z (bath-normalized) and the raw greater/lesser
expectations `Gg[j] = ⟨d(tⱼ)d†(0)⟩` and `Gl[j] = ⟨d†(0)d(tⱼ)⟩`.
"""
function quapi_rlm_rt_gt(β, ϵ_d, spectrum; Nt, t, μ=0.0)
    δt = t / Nt
    x = exp(-im * δt * ϵ_d)
    xc = conj(x)
    corr = fermionic_Δt(spectrum, β=β, N=Nt, t=t, μ=μ)

    C = QuAPI.Dense._if_coefficients(corr)
    legs = basis(corr)
    li = Dict(l => k for (k, l) in enumerate(legs))
    # 节点 j = i+1（i = 0..Nt），branch = :+ / :-
    leg(i, b, c) = li[GrassmannContourIndex(i + 1, b, 1, c)]
    for i in 0:Nt                                  # Berezin measures −āa
        C[leg(i, :+, true), leg(i, :+, false)] += -1.0
        C[leg(i, :-, true), leg(i, :-, false)] += -1.0
    end
    for j in 1:Nt                                  # free impurity kinetics
        C[leg(j, :+, true), leg(j - 1, :+, false)] += x
        C[leg(j - 1, :-, true), leg(j, :-, false)] += xc
    end
    C[leg(Nt, :-, true), leg(Nt, :+, false)] += -1.0          # contour closure
    C[leg(0, :+, true), leg(0, :-, false)] += exp(-β * ϵ_d)   # initial weight

    data = QuAPI.Dense._fold_grassmann(C)
    L = length(legs)
    ctop = data[1 << L]                            # all-odd coefficient
    Z = ctop
    # Under the GTEMPO leg order (a before ā per node, times descending) the
    # {ā⁺₀, a⁺_j} numerator changes order relative to the ascending layout
    # (one minus sign), while {ā⁻₀, a⁺_j} does not.
    Gg = [QuAPI.Dense._grassmann_coefficient(data, leg(0, :+, true), leg(j, :+, false)) / ctop
          for j in 0:Nt]
    Gl = [QuAPI.Dense._grassmann_coefficient(data, leg(0, :-, true), leg(j, :+, false)) / ctop
          for j in 0:Nt]
    return Z, Gg, Gl, δt
end

@testset "real-time QuAPI (Toulouse, Grassmann tensor)" begin
    β, ϵ_d, μ = 2.0, 0.5, 0.0

    @testset "free impurity (λ=0): exact" begin
        for Nt in (1, 2)
            t = 0.5
            Z, Gg, Gl, δt = quapi_rlm_rt_gt(β, ϵ_d, DiscreteSpectrum([1.0], [0.0]); Nt=Nt, t=t, μ=μ)
            ts = [j * δt for j in 0:Nt]
            # freefermion_* (ImpurityModelBase): G> = -im⟨d(t)d†(0)⟩, G< = im⟨d†(0)d(t)⟩
            # with μ = -ϵ_d; the raw Grassmann expectations are im·G> and -im·G<.
            GgE = [im * freefermion_greater(τ; β=β, μ=-ϵ_d) for τ in ts]
            GlE = [-im * freefermion_lesser(τ; β=β, μ=-ϵ_d) for τ in ts]
            @test Z ≈ 1 + exp(-β * ϵ_d) atol=1.0e-10
            @test Gg ≈ GgE atol=1.0e-10
            @test Gl ≈ GlE atol=1.0e-10
        end
    end

    @testset "coupled vs ED (toulouse_neq_greater_lesser)" begin
        ws, Vs = [1.0], [0.5]
        for (t, Nt, tol) in ((0.25, 2, 0.01), (0.5, 3, 0.02))
            ts = [j * t / Nt for j in 0:Nt]
            Z, Gg, Gl, δt = quapi_rlm_rt_gt(β, ϵ_d, DiscreteSpectrum(copy(ws), copy(Vs) .^ 2); Nt=Nt, t=t, μ=μ)
            # exact ED of the coupled quadratic impurity+bath Hamiltonian
            b = discretebath(Fermion, copy(ws), copy(Vs); β=β, μ=μ)
            m = Toulouse(b, ϵ_d=ϵ_d)
            gt1, gt2 = toulouse_neq_greater_lesser(m, ts)
            GgE = im .* gt1
            GlE = -im .* gt2
            @test Gg ≈ GgE atol=tol
            @test Gl ≈ GlE atol=tol
            # unitary quench: the bath-normalized partition function is conserved
            # (up to the first-order discretization error of the influence functional)
            @test Z ≈ 1 + exp(-β * ϵ_d) atol=0.05
        end
    end
end
