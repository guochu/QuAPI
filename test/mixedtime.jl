println("------------------------------------")
println("|        Mixed-time Δτ             |")
println("------------------------------------")



@testset "fermionic Δm" begin
	β = 2
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			x = fermionic_Δm(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)
			y = fermionic_Δm2(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)
			for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
				@test branch(x, b1=b1, b2=b2) ≈ branch(y, b1=b1, b2=b2) atol=1.0e-8
			end
			@test x ≈ y atol=1.0e-8

			z = transpose(x)
			for b1 in (:+, :-, :τ)
				k1 = ifelse(b1 == :τ, isize(z), rsize(z))
				for b2 in (:+, :-, :τ)
					k2 = ifelse(b2 == :τ, isize(z), rsize(z))
					for i in 1:k1, j in 1:k2
						@test index(x, i, j, b1=b1, b2=b2) == index(z, j, i, b1=b2, b2=b1)
					end
				end
			end
		end
	end

	f1, f2 = semicircular(t=1), semicircular(t=2)
	μ = 0.23
	x = fermionic_Δm(f1, β=β, t=2, Nt=4, Nτ=5, μ=μ)
	y = fermionic_Δm(f2, β=β, t=2, Nt=4, Nτ=5, μ=μ)

	z = x - y

	for b1 in (:+, :-, :τ)
		k1 = ifelse(b1 == :τ, isize(z), rsize(z))
		for b2 in (:+, :-, :τ)
			k2 = ifelse(b2 == :τ, isize(z), rsize(z))
			for i in 1:k1, j in 1:k2
				@test index(z, i, j, b1=b1, b2=b2) ≈ index(x, i, j, b1=b1, b2=b2) - index(y, i, j, b1=b1, b2=b2)
			end
		end
	end
end


@testset "bosonic Δm" begin
	β = 2
	for f in (Leggett(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0)
			x = bosonic_Δm(f, β=β, t=1, Nt=2, Nτ=5, μ=μ)
			y = bosonic_Δm2(f, β=β, t=1, Nt=2, Nτ=5, μ=μ)
			for b1 in (:+, :-, :τ), b2 in (:+, :-, :τ)
				@test branch(x, b1=b1, b2=b2) ≈ branch(y, b1=b1, b2=b2) atol=1.0e-8
			end
			@test x ≈ y atol=1.0e-8

			z = transpose(x)
			for b1 in (:+, :-, :τ)
				k1 = ifelse(b1 == :τ, isize(z), rsize(z))
				for b2 in (:+, :-, :τ)
					k2 = ifelse(b2 == :τ, isize(z), rsize(z))
					for i in 1:k1, j in 1:k2
						@test index(x, i, j, b1=b1, b2=b2) == index(z, j, i, b1=b2, b2=b1)
					end
				end
			end
		end
	end

	f1, f2 = Leggett(), Leggett(d=2)
	μ = -0.2
	x = bosonic_Δm(f1, β=β, t=1, Nt=2, Nτ=5, μ=μ)
	y = bosonic_Δm(f2, β=β, t=1, Nt=2, Nτ=5, μ=μ)
	z = x + y
	for b1 in (:+, :-, :τ)
		k1 = ifelse(b1 == :τ, isize(z), rsize(z))
		for b2 in (:+, :-, :τ)
			k2 = ifelse(b2 == :τ, isize(z), rsize(z))
			for i in 1:k1, j in 1:k2
				@test index(z, i, j, b1=b1, b2=b2) ≈ index(x, i, j, b1=b1, b2=b2) + index(y, i, j, b1=b1, b2=b2) 
			end
		end
	end
end