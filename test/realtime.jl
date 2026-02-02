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

@testset "infinite fermionic Δt" begin
	tol = 1.0e-2
	β = 5
	t = 10.
	δt = 0.1
	Nt = round(Int, t/δt)
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			bath = fermionicbath(f, β=β, μ=μ)
			corr1 = Δt(bath, N=Nt, t=t)
			x = infinite_Δt(bath, δt=δt, t=100.)
			corr2 = RealCorrelationFunction(branch(x, :+, :+), branch(x, :+, :-), branch(x, :-, :+), branch(x, :-, :-))
			for b1 in (:+, :-), b2 in (:+, :-)
				m1 = branch(corr1, :+, :-)[1:Nt, 1:Nt]
				m2 = branch(corr2, :+, :-)[1:Nt, 1:Nt]
				@test norm(m1-m2)/norm(m1) < tol
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


@testset "infinite bosonic Δt" begin
	tol = 1.0e-2
	β = 1
	t = 10.
	δt = 0.1
	Nt = round(Int, t/δt)
	for f in (Leggett(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0)
			bath = bosonicbath(f, β=β, μ=μ)
			corr1 = Δt(bath, N=Nt, t=t)
			x = infinite_Δt(bath, δt=δt, t=100.)
			corr2 = RealCorrelationFunction(branch(x, :+, :+), branch(x, :+, :-), branch(x, :-, :+), branch(x, :-, :-))
			for b1 in (:+, :-), b2 in (:+, :-)
				m1 = branch(corr1, :+, :-)[1:Nt, 1:Nt]
				m2 = branch(corr2, :+, :-)[1:Nt, 1:Nt]
				@test norm(m1-m2)/norm(m1) < tol
			end
		end
	end
end