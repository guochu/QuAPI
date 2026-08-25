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

@testset "infinite fermionic Δτ" begin
	tol = 2.0e-2
	β = 100
	δτ = 0.1
	N = round(Int, β/δτ)
	Nh = div(N, 2)
	spec = semicircular()
	for μ in (-0.3, 0, 0.5)
		bath1 = fermionicbath(spec, β=β, μ=μ)
		bath2 = FermionicVacuum(spec, μ=μ)
		corr1 = Δτ(bath1, N=N)
		corr2 = infinite_Δτ(bath2, δτ=δτ, β=β)
		corr2 = ImagCorrelationFunction(CorrelationMatrix{Float64}(corr2.ηⱼₖ, corr2.ηₖⱼ))
		
		m1 = corr1.data[1:Nh, 1:Nh]
		m2 = corr2.data[1:Nh, 1:Nh]

		@test norm(m1 - m2) / norm(m1) < tol
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

@testset "infinite bosonic Δτ" begin
	tol = 2.0e-2
	β = 100
	δτ = 0.1
	N = round(Int, β/δτ)
	Nh = div(N, 2)
	spec = Leggett()

	bath1 = bosonicbath(spec, β=β)
	bath2 = BosonicVacuum(spec)
	corr1 = Δτ(bath1, N=N)
	corr2 = infinite_Δτ(bath2, δτ=δτ, β=β)
	corr2 = ImagCorrelationFunction(CorrelationMatrix{Float64}(corr2.ηⱼₖ, corr2.ηₖⱼ))
	
	m1 = corr1.data[1:Nh, 1:Nh]
	m2 = corr2.data[1:Nh, 1:Nh]

	@test norm(m1 - m2) / norm(m1) < tol

end

@testset "CorrelationMatrix" begin
	A = CorrelationMatrix([1.0+0im, 0.5, 0.2], [0.5+0im, 0.3, 0.1])
	B = CorrelationMatrix([1.0+0im, 0.5, 0.2], [0.5+0im, 0.3, 0.1])
	C = CorrelationMatrix([1.0+0im, 0.5, 0.2], [0.5+0im, 0.3, 0.1+1.0e-9im])
	@test A == B
	@test A != C
	@test A ≈ B
	@test A ≈ C atol=1.0e-8
	@test !isapprox(A, C, atol=1.0e-12)
	@test iszero(A - A)
	@test !iszero(A)
	@test transpose(A)[1, 2] == A[2, 1]
	@test (A + B)[1, 3] ≈ 2*A[1, 3]
	@test (-A)[2, 1] ≈ -A[2, 1]
end

@testset "Δiw_to_Δτ" begin
	β = 5.0
	N = 20
	ω0 = 1.0
	nmax = 2000

	# analytic benchmark: for a single fermionic level at ω0+μ the
	# Matsubara-domain hybridization Δ(iωₙ) = 1/(iωₙ-(ω0+μ)) reproduces
	# fermionic_Δτ(DiracDelta(ω=ω0); μ) via the doubly-integrated Matsubara
	# sum.  The two storage conventions split the diagonal differently, so the
	# materialized Δ[i,j] matrices are compared.
	for μ in (0.0, 0.3)
		Δiw = [1/(im*ω - (ω0+μ)) for ω in ifrequencies(β, nmax)]
		m1 = Matrix(Δiw_to_Δτ(Δiw, β=β, N=N).data)
		m2 = Matrix(fermionic_Δτ(DiracDelta(ω=ω0), β=β, N=N, μ=μ).data)
		@test m1 ≈ m2 atol=1.0e-8
	end

	# an odd number of Matsubara frequencies is rejected
	@test_throws ArgumentError Δiw_to_Δτ([1.0, 2.0, 3.0], β=β, N=N)

	# zero hybridization gives the zero correlation matrix
	@test iszero(Δiw_to_Δτ(zeros(2), β=β, N=N))

	# homogeneity: scaling the hybridization scales the output linearly
	Δiw = [1/(im*ω - ω0) for ω in ifrequencies(β, nmax)]
	c1 = Δiw_to_Δτ(3 .* Δiw, β=β, N=N)
	c2 = Δiw_to_Δτ(Δiw, β=β, N=N)
	for i in 1:N, j in 1:N
		@test index(c1, i, j) ≈ 3*index(c2, i, j)
	end

	# a real, even hybridization gives a symmetric Toeplitz matrix
	Δiw = [2ω0/(ω0^2 + ω^2) for ω in ifrequencies(β, nmax)]
	c3 = Δiw_to_Δτ(Δiw, β=β, N=N)
	@test transpose(c3) == c3
end