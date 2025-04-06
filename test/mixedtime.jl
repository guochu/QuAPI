println("------------------------------------")
println("|        Mixed-time Δτ             |")
println("------------------------------------")



@testset "fermionic Δτ" begin
	β = 5
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			x = fermionic_Δm(f, β=β, t=2, Nt=10, Nτ=5, μ=μ)
			y = fermionic_Δm2(f, β=β, t=2, Nt=10, Nτ=5, μ=μ)
			for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
				@test branch(x, b1=b1, b2=b2) ≈ branch(y, b1=b1, b2=b2) atol=1.0e-8
			end
			@test x ≈ y atol=1.0e-8
		end
	end
end


@testset "bosonic Δτ" begin
	β = 5
	for f in (Leggett(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0)
			x = bosonic_Δm(f, β=β, t=2, Nt=10, Nτ=5, μ=μ)
			y = bosonic_Δm2(f, β=β, t=2, Nt=10, Nτ=5, μ=μ)
			for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
				@test branch(x, b1=b1, b2=b2) ≈ branch(y, b1=b1, b2=b2) atol=1.0e-8
			end
			@test x ≈ y atol=1.0e-8
		end
	end
end