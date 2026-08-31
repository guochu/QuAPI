println("------------------------------------")
println("|        Mixed-time QuAPI          |")
println("|     (Kadanoff / L-shaped)        |")
println("------------------------------------")

using Test, LinearAlgebra

"""
    quapi_mixed(β, ε, Δ, ωs, λs; Nt, Nτ, t, μ=0.0)

Kadanoff (L-shaped) contour QuAPI for the spin-boson model

    H = ε/2 σz + Δ/2 σx + Σ_k ω_k a_k†a_k + (σz/2) Σ_k λ_k (a_k† + a_k)

starting from the full thermal state e^{-βH}/Z.  The contour has three
branches: real-time forward (+: t=0→t), real-time backward (-: t→0) and
imaginary-time (τ: 0→β).  At the junction the forward node s+_0 (end of the
τ branch, τ=β) and the backward node s-_0 (start of the τ branch, τ=0) close
the loop, and the imaginary-time propagator chain runs from s-_0, through the
interior τ nodes c_2..c_Nτ, back to s+_0.  The influence functional is built
from the nine branch blocks (η/ξ/ζ) of `bosonic_Δm` assembled into one
block matrix G, with weight exp(sᵀ G s) for the spin path vector s.

The spin path is discretized on Nt real steps δt = t/Nt and Nτ imaginary
steps δτ = β/Nτ (requires Nτ ≥ 2), and the full sum over the 2^(2Nt+Nτ-1)
free spin configurations is evaluated.  Returns ⟨σz(t)⟩, ⟨σx(t)⟩, the
two-time correlation ⟨σz(t)σz(0)⟩ (with σz inserted at the junction s+_0),
the unnormalized reduced density matrix ρ = e^{-iHt} e^{-βH} e^{iHt} (traced
over the bath), and the step sizes δt, δτ.
"""
function quapi_mixed(β, ε, Δ, ωs, λs; Nt, Nτ, t, μ=0.0)
    δt = t / Nt
    δτ = β / Nτ
    σx = [0.0 1.0; 1.0 0.0]
    σz = [1.0 0.0; 0.0 -1.0]
    Hs = ε/2 * σz + Δ/2 * σx
    Wp = exp(-im * δt * Hs)     # forward real-time propagator
    Wm = exp(+im * δt * Hs)     # backward real-time propagator
    Wτ = exp(-δτ * Hs)          # imaginary-time propagator

    corr = bosonic_Δm(DiscreteSpectrum(ωs, λs .^ 2), β=β, Nτ=Nτ, t=t, Nt=Nt, μ=μ)
    Gpp = Matrix(branch(corr, b1=:+, b2=:+))   # (Nt+1)×(Nt+1), nodes 0..Nt
    Gpm = Matrix(branch(corr, b1=:+, b2=:-))
    Gmp = Matrix(branch(corr, b1=:-, b2=:+))
    Gmm = Matrix(branch(corr, b1=:-, b2=:-))
    Gpt = Matrix(branch(corr, b1=:+, b2=:τ))   # ×(Nτ), τ nodes s-_0, c_2..c_Nτ
    Gmt = Matrix(branch(corr, b1=:-, b2=:τ))
    Gtp = Matrix(branch(corr, b1=:τ, b2=:+))
    Gtm = Matrix(branch(corr, b1=:τ, b2=:-))
    Gtt = Matrix(branch(corr, b1=:τ, b2=:τ))   # (Nτ)×(Nτ)
    G = [Gpp Gpm Gpt; Gmp Gmm Gmt; Gtp Gtm Gtt]

    # enumerate free spin configurations (values ±1/2):
    #   s+_0..s+_{Nt-1}, s-_0..s-_{Nt-1}, c_2..c_Nτ
    nf = 2^(2Nt + Nτ - 1)
    ρ = zeros(ComplexF64, 2, 2)
    ρ0 = zeros(ComplexF64, 2, 2)    # with σz(s+_0) inserted at the junction
    σzval(v) = ifelse(v == 0.5, 1.0, -1.0)
    @inbounds for a in 1:2, b in 1:2
        sa = ifelse(a == 1, 0.5, -0.5)   # fixed endpoint value s+_Nt
        sb = ifelse(b == 1, 0.5, -0.5)   # fixed endpoint value s-_Nt
        for m in 0:nf-1
            sp = [ifelse(((m >> j) & 1) == 0, 0.5, -0.5) for j in 0:Nt-1]
            sm = [ifelse(((m >> (Nt + j)) & 1) == 0, 0.5, -0.5) for j in 0:Nt-1]
            ct = [ifelse(((m >> (2Nt + j)) & 1) == 0, 0.5, -0.5) for j in 0:Nτ-2]
            idx(v) = v == 0.5 ? 1 : 2
            # forward path ⟨s+_Nt|e^{-iHt}|s+_0⟩
            Pplus = 1.0 + 0im
            for j in 0:Nt-1
                snext = j == Nt - 1 ? sa : sp[j+2]
                scur  = sp[j+1]
                Pplus *= Wp[idx(snext), idx(scur)]
            end
            # backward path ⟨s-_Nt|e^{+iHt}|s-_0⟩
            Pminus = 1.0 + 0im
            for j in 0:Nt-1
                scur  = j == 0 ? sm[1] : sm[j+1]
                snext = j == Nt - 1 ? sb : sm[j+2]
                Pminus *= Wm[idx(snext), idx(scur)]
            end
            # imaginary propagator ⟨s+_0|e^{-βH}|s-_0⟩ via c_1 = s-_0, c_2..c_Nτ
            Pτ = Wτ[idx(ct[1]), idx(sm[1])]
            for j in 2:Nτ-1
                Pτ *= Wτ[idx(ct[j]), idx(ct[j-1])]
            end
            Pτ *= Wτ[idx(sp[1]), idx(ct[end])]
            # influence functional
            Sfp = vcat(sp, sa)              # extended + nodes 0..Nt
            Sgq = vcat(sm, sb)              # extended - nodes 0..Nt
            Sτ = vcat(sm[1], ct)            # τ nodes s-_0, c_2..c_Nτ
            v = vcat(Sfp, Sgq, Sτ)
            IF = exp(dot(v, G, v))
            w = Pplus * Pminus * Pτ * IF
            ρ[a, b] += w
            ρ0[a, b] += σzval(sp[1]) * w
        end
    end
    trρ = ρ[1, 1] + ρ[2, 2]
    σz_t = (ρ[1, 1] - ρ[2, 2]) / trρ
    σx_t = (ρ[1, 2] + ρ[2, 1]) / trρ
    Czz = (ρ0[1, 1] - ρ0[2, 2]) / trρ      # σz(t)=+1 for |↑⟩, -1 for |↓⟩
    return σz_t, σx_t, Czz, ρ, δt, δτ
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

function ed_mixed_obs(ε, Δ, ωs, λs, ts; β, d=12)
    K = length(ωs)
    H = ed_hamiltonian(ε, Δ, ωs, λs; d=d)
    cache = eigencache(H)
    Ib = Matrix{Float64}(I, d^K, d^K)
    ρ0 = exp(-β * H)
    Z = tr(ρ0)
    ρ0n = ρ0 / Z
    σz_full = kron([1.0 0.0; 0.0 -1.0], Ib)
    σx_full = kron([0.0 1.0; 1.0 0.0], Ib)
    σz = Float64[]; σx = Float64[]; Czz = ComplexF64[]
    for t in ts
        ρt = timeevo(ρ0n, H, -im * t, cache)
        push!(σz, real(tr(ρt * σz_full)))
        push!(σx, real(tr(ρt * σx_full)))
        # ⟨σz(t)σz(0)⟩ = tr[ρ0 σz(t) σz], σz(t) = e^{iHt} σz e^{-iHt}
        U = cache.U * Diagonal(exp.(-im .* t .* cache.λs)) * adjoint(cache.U)
        szH = U' * σz_full * U
        push!(Czz, tr(ρ0n * szH * σz_full))
    end
    return σz, σx, Czz
end

@testset "mixed-time QuAPI (Kadanoff contour, spin-boson)" begin
    β, ε, Δ = 2.0, 0.5, 0.5
    d = 12   # boson Fock cutoff for ED

    # free spin (λ=0): the path sum must reproduce the exact thermal
    # expectation values and the two-time correlation exactly (IF = 1)
    @testset "free spin (λ=0)" begin
        Hs = ε/2 * [1.0 0.0; 0.0 -1.0] + Δ/2 * [0.0 1.0; 1.0 0.0]
        σz = [1.0 0.0; 0.0 -1.0]
        σx = [0.0 1.0; 1.0 0.0]
        ρth = exp(-β * Hs)
        ρth ./= tr(ρth)
        for (t, Nt, Nτ) in ((0.5, 2, 3), (1.0, 4, 4), (2.0, 6, 6))
            σzq, σxq, Czzq, ρ, δt, δτ = quapi_mixed(β, ε, Δ, [1.0], [0.0]; Nt=Nt, Nτ=Nτ, t=t)
            Ut = exp(im * t * Hs)
            ρe = Ut' * ρth * Ut
            @test real(σzq) ≈ real(tr(ρe * σz)) atol=1.0e-8
            @test real(σxq) ≈ real(tr(ρe * σx)) atol=1.0e-8
            # ⟨σz(t)σz(0)⟩ = tr[σz e^{-iHt} σz ρ_th e^{iHt}] (cyclicity of the trace)
            CzzE = tr(σz * Ut' * σz * ρth * Ut)
            @test abs(Czzq - CzzE) < 1.0e-8
        end
    end

    # coupled spin-boson (ω=1, λ=0.5) from the thermal state, vs ED.  The
    # discretized influence functional converges as Nt, Nτ grow; the largest
    # error occurs at the largest t, so tolerances are loose there, and the
    # N → N+1 trend is asserted separately at fixed t=2.0.
    @testset "spin-boson vs ED" begin
        ts = [0.5, 1.0, 1.5, 2.0]
        Ns = [(2, 2), (3, 3), (4, 4), (5, 5)]
        σze, σxe, CzzE = ed_mixed_obs(ε, Δ, [1.0], [0.5], ts; β=β, d=d)
        for (t, (Nt, Nτ)) in zip(ts, Ns)
            σzq, σxq, Czzq, ρ, δt, δτ = quapi_mixed(β, ε, Δ, [1.0], [0.5]; Nt=Nt, Nτ=Nτ, t=t)
            i = findfirst(==(t), ts)
            @test real(σzq) ≈ σze[i] atol=0.05
            @test real(σxq) ≈ σxe[i] atol=0.06
            @test abs(Czzq - CzzE[i]) < 0.06
        end
    end

    # convergence at fixed t=2.0: finer discretization ⇒ smaller error vs ED
    @testset "convergence at fixed t=2.0" begin
        σze, σxe, CzzE = ed_mixed_obs(ε, Δ, [1.0], [0.5], [2.0]; β=β, d=d)
        errsz = Float64[]; errsx = Float64[]; errCzz = Float64[]
        for (Nt, Nτ) in ((4, 4), (5, 5), (6, 6))
            σzq, σxq, Czzq, ρ, δt, δτ = quapi_mixed(β, ε, Δ, [1.0], [0.5]; Nt=Nt, Nτ=Nτ, t=2.0)
            push!(errsz, abs(real(σzq) - σze[1]))
            push!(errsx, abs(real(σxq) - σxe[1]))
            push!(errCzz, abs(Czzq - CzzE[1]))
        end
        @test errsz[end] < errsz[1]
        @test errsx[end] < errsx[1]
        @test errCzz[end] < errCzz[1]
    end

    # branch-matrix sum checks.  Unlike the closed real-time contour, on the
    # L-shaped contour the constant path does not give IF = 1: the imaginary
    # branch contributes the static (displaced-oscillator) weight.  The exact
    # identities that hold are:
    #   Σ_{i,j} Gττ[i,j] = β Σ_k λ_k²/ω_k          (τ-loop sum)
    #   G⁺τ = -G⁻τ and Gτ⁺ = -Gτ⁻ (mixed blocks cancel pairwise)
    @testset "branch-matrix sum" begin
        corr = bosonic_Δm(DiscreteSpectrum([1.0], [0.5] .^ 2), β=β, Nτ=6, t=2.0, Nt=6, μ=0.0)
        Gττ = Matrix(branch(corr, b1=:τ, b2=:τ))
        @test sum(Gττ) ≈ β * sum([0.5] .^ 2 ./ [1.0]) atol=1.0e-10
        Gpt = Matrix(branch(corr, b1=:+, b2=:τ))
        Gmt = Matrix(branch(corr, b1=:-, b2=:τ))
        Gtp = Matrix(branch(corr, b1=:τ, b2=:+))
        Gtm = Matrix(branch(corr, b1=:τ, b2=:-))
        @test sum(Gpt + Gmt) ≈ 0 atol=1.0e-12
        @test sum(Gtp + Gtm) ≈ 0 atol=1.0e-12
    end
end
