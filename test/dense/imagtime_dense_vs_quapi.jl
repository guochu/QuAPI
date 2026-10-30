println("------------------------------------")
println("|   dense vs QuAPI: imag time     |")
println("------------------------------------")

# Dense folding vs analytic closed form (free limit): two consumers of the
# same fermionic_Δτ data.
using Test, LinearAlgebra
using QuAPI
using QuAPI.Dense
using ImpurityModelBase

function quapi_gt_free(β, ϵ_d, N)
    δτ = β / N
    x = exp(-δτ * ϵ_d)
    corr = fermionic_Δτ(DiracDelta(ω=1.0, α=0.0); β=β, N=N, μ=0.0)
    h = GrassmannHybridizationMatrix(corr)
    legs = h.basis
    li = Dict(l => k for (k, l) in enumerate(legs))
    C = copy(h.Δ)
    leg(i, c) = li[GrassmannContourIndex(i, :τ, 1, c)]
    for i in 1:N
        C[leg(i, true), leg(i, false)] += -1.0
    end
    for i in 2:N
        C[leg(i, true), leg(i-1, false)] += x
    end
    C[leg(1, true), leg(N, false)] += -x
    gif = GrassmannIF(GrassmannHybridizationMatrix(legs, C))
    ctop = Zvalue(gif)
    # gif.basis elements are NOT identical objects to legs elements — build the
    # position dict on the SAME array used for complementindex lookups
    idxof(i, c) = GrassmannContourIndex(i, :τ, 1, c)
    pos = Dict(idxof(i, c) => k for (k, i) in enumerate(1:N) for c in (false, true))
    msign(ia, iab) = ia < iab ? 1.0 : -1.0
    G = [msign(pos[idxof(1, true)], pos[idxof(m, false)]) *
         gif[complementindex(gif.basis, legs[leg(1, true)],
                             legs[leg(m, false)])...] / ctop
         for m in 1:N]        # τ_m = (m−1)δτ
    return G, δτ   # msign 已吸收提取符号，G 即标准约定 −⟨d(τ)d†⟩/Z 归一值
end

@testset "dense vs QuAPI (fermionic imag time)" begin
    β, ϵ_d = 2.0, 0.5
    for N in (4, 6)
        G, δτ = quapi_gt_free(β, ϵ_d, N)
        ref = [-exp(-((m - 1) * δτ) * ϵ_d) / (1 + exp(-β * ϵ_d)) for m in 1:N]
        @test G ≈ ref atol=1.0e-12
    end
end
