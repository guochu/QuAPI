println("------------------------------------")
println("|        Real-time Δt              |")
println("------------------------------------")


@testset "fermionic Δt" begin
	β = 5
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			x = fermionic_Δt(f, β=β, t=2, N=10, μ=μ)

			z = transpose(x)
			for i in 1:size(z, 1), b1 in (:+, :-)
				for j in 1:size(z, 2), b2 in (:+, :-)
					@test index(x, i, j, b1=b1, b2=b2) == index(z, j, i, b1=b2, b2=b1)
				end
			end
		end
	end

	f1, f2 = semicircular(), DiracDelta(ω=0.5)

	for μ in (-0.3, 0, 0.4)
		x = fermionic_Δt(f1, β=β, t=2, N=10, μ=μ)
		y = fermionic_Δt(f2, β=β, t=2, N=10, μ=μ)

		z = x + y
		for i in 1:size(z, 1), b1 in (:+, :-)
			for j in 1:size(z, 2), b2 in (:+, :-)
				@test index(z, i, j, b1=b1, b2=b2) ≈ index(x, i, j, b1=b1, b2=b2) + index(y, i, j, b1=b1, b2=b2) 
			end
		end
	end	
end


@testset "bosonic Δt" begin
	β = 5
	for f in (Leggett(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0)
			x = bosonic_Δt(f, β=β, t=2, N=5, μ=μ)

			z = transpose(x)
			for i in 1:size(z, 1), b1 in (:+, :-)
				for j in 1:size(z, 2), b2 in (:+, :-)
					@test index(x, i, j, b1=b1, b2=b2) == index(z, j, i, b1=b2, b2=b1)
				end
			end
		end
	end

	f1, f2 = Leggett(), DiracDelta(ω=0.5)

	for μ in (-0.3, 0)
		x = bosonic_Δt(f1, β=β, t=2, N=5, μ=μ)
		y = bosonic_Δt(f2, β=β, t=2, N=5, μ=μ)

		z = x - y
		for i in 1:size(z, 1), b1 in (:+, :-)
			for j in 1:size(z, 2), b2 in (:+, :-)
				@test index(z, i, j, b1=b1, b2=b2) ≈ index(x, i, j, b1=b1, b2=b2) - index(y, i, j, b1=b1, b2=b2)
			end
		end
	end	
end