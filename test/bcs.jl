println("------------------------------------")
println("|           BCS Δτ Δt Δm           |")
println("------------------------------------")


@testset "Imaginary time" begin
	β = 5
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			corr1 = fermionic_Δτ(f, β=β, N=10, μ=μ)
			corr2 = bcs_Δτ(f, β=β, N=10, μ=μ)
			@test corr1 ≈ corr2[true, false] atol=1.0e-8
			@test corr1 ≈ -transpose(corr2[false, true]) atol=1.0e-8
			@test iszero(corr2[true, true])
			@test iszero(corr2[false, false])
		end
	end
end



@testset "Real time" begin
	β = 5
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			x = fermionic_Δt(f, β=β, t=2, N=10, μ=μ)
			y = bcs_Δt(f, β=β, t=2, N=10, μ=μ)

			@test x ≈ y[true, false] atol=1.0e-8
			@test x ≈ -transpose(y[false, true]) atol=1.0e-8
			@test iszero(y[true, true])
			@test iszero(y[false, false])

		end
	end
end


@testset "Mixed time" begin
	β = 2
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			x = fermionic_Δm(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)
			y = bcs_Δm(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)

			@test x ≈ y[true, false] atol=1.0e-8
			@test x ≈ -transpose(y[false, true]) atol=1.0e-8
			@test iszero(y[true, true])
			@test iszero(y[false, false])

		end	
	end
end
