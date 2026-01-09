println("------------------------------------")
println("|           BCS Δτ Δt Δm           |")
println("------------------------------------")


@testset "Imaginary time" begin
	atol = 1.0e-8
	β = 5
	μ = 0
	for f in (Leggett(), DiracDelta(ω=0.5))
		corr1 = fermionic_Δτ(f, β=β, N=10, μ=μ)
		corr2 = bcs_Δτ(f, β=β, N=10, μ=μ, Δ=complex(0))
		@test corr2[true, true] === corr2[:↑, :↓] === corr2.ud === corr2.cc
		@test corr2[true, false] === corr2[:↑, :↑] === corr2.uu === corr2.cn
		@test corr2[false, true] === corr2[:↓, :↓] === corr2.dd === corr2.nc
		@test corr2[false, false] === corr2[:↓, :↑] === corr2.du === corr2.nn
		@test corr1 ≈ corr2[true, false] atol=atol
		@test corr1 ≈ -transpose(corr2[false, true]) atol=atol
		@test iszero(corr2[true, true])
		@test iszero(corr2[false, false])
	end
end



@testset "Real time" begin
	atol = 1.0e-8
	β = 5
	μ = 0
	for f in (semicircular(), DiracDelta(ω=0.5))
		x = fermionic_Δt(f, β=β, t=2, N=10, μ=μ)
		y = bcs_Δt(f, β=β, t=2, N=10, μ=μ)

		@test x ≈ y[true, false] atol=atol
		@test x ≈ -transpose(y[false, true]) atol=atol
		@test iszero(y[true, true])
		@test iszero(y[false, false])
	end
end


@testset "Mixed time" begin
	atol = 1.0e-8
	β = 2
	μ = 0
	for f in (semicircular(), DiracDelta(ω=0.5))
		x = fermionic_Δm(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)
		y = bcs_Δm(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)

		@test x ≈ y[true, false] atol=atol
		@test x ≈ -transpose(y[false, true]) atol=atol
		@test iszero(y[true, true])
		@test iszero(y[false, false])
	end
end