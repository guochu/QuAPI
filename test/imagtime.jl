println("------------------------------------")
println("|      Imaginary-time Δτ           |")
println("------------------------------------")



@testset "fermionic Δτ" begin
	β = 5
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			corr1 = fermionic_Δτ(f, β=β, N=10, μ=μ)
			corr2 = fermionic_Δτ2(f, β=β, N=10, μ=μ)
			@test corr1 ≈ corr2 atol=1.0e-8
		end
	end
end


@testset "bosonic Δτ" begin
	β = 5
	for f in (Leggett(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0)
			corr1 = bosonic_Δτ(f, β=β, N=10, μ=μ)
			corr2 = bosonic_Δτ2(f, β=β, N=10, μ=μ)
			@test corr1 ≈ corr2 atol=1.0e-8
		end
	end
end