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

			corr3 = transpose(corr1)
			for i in 1:size(corr3, 1), j in 1:size(corr3, 2)
				@test index(corr1, i, j) == index(corr3, j, i)
			end
		end
	end
	f1, f2 = semicircular(), DiracDelta(ω=0.5)
	for μ in (-0.3, 0, 0.3)
		corr1 = fermionic_Δτ(f1, β=β, N=5, μ=μ)
		corr2 = fermionic_Δτ(f2, β=β, N=5, μ=μ)

		corr3 = corr1 + corr2
		for i in 1:size(corr3, 1), j in 1:size(corr3, 2)
			@test index(corr3, i, j) ≈ index(corr1, i, j) + index(corr2, i, j)  
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

			corr3 = transpose(corr1)
			for i in 1:size(corr3, 1), j in 1:size(corr3, 2)
				@test index(corr1, i, j) == index(corr3, j, i)
			end
		end
	end
	f1, f2 = Leggett(), DiracDelta(ω=0.5)

	for μ in (-0.3, 0)
		corr1 = bosonic_Δτ(f1, β=β, N=10, μ=μ)
		corr2 = bosonic_Δτ(f2, β=β, N=10, μ=μ)

		corr3 = corr1 - corr2
		for i in 1:size(corr3, 1), j in 1:size(corr3, 2)
			@test index(corr3, i, j) ≈ index(corr1, i, j) - index(corr2, i, j) 
		end
	end	
end