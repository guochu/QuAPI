println("------------------------------------")
println("|        Mixed-time Δτ             |")
println("------------------------------------")



@testset "fermionic Δm" begin
	β = 2
	for f in (semicircular(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0, 0.3)
			x = fermionic_Δm(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)
			y = fermionic_Δm_deprecated(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)
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
						@test branch(x; b1=b1, b2=b2)[i, j] == branch(z; b1=b2, b2=b1)[j, i]
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
				@test branch(z; b1=b1, b2=b2)[i, j] ≈ branch(x; b1=b1, b2=b2)[i, j] - branch(y; b1=b1, b2=b2)[i, j]
			end
		end
	end
end


@testset "bosonic Δm" begin
	β = 2
	for f in (Leggett(), DiracDelta(ω=0.5))
		for μ in (-0.3, 0)
			x = bosonic_Δm(f, β=β, t=1, Nt=2, Nτ=5, μ=μ)
			y = bosonic_Δm_deprecated(f, β=β, t=1, Nt=2, Nτ=5, μ=μ)
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
						@test branch(x; b1=b1, b2=b2)[i, j] == branch(z; b1=b2, b2=b1)[j, i]
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
				@test branch(z; b1=b1, b2=b2)[i, j] ≈ branch(x; b1=b1, b2=b2)[i, j] + branch(y; b1=b1, b2=b2)[i, j] 
			end
		end
	end
end

@testset "Δm bath convenience" begin
	β = 2
	fb = fermionicbath(semicircular(), β=β, μ=0)
	@test Δm(fb, Nt=4, Nτ=5, t=2) ≈ fermionic_Δm(semicircular(), β=β, Nt=4, Nτ=5, t=2, μ=0) atol=1.0e-8

	bb = bosonicbath(Leggett(), β=β, μ=0)
	@test Δm(bb, Nt=2, Nτ=5, t=1) ≈ bosonic_Δm(Leggett(), β=β, Nt=2, Nτ=5, t=1, μ=0) atol=1.0e-8
end
println("------------------------------------")
println("|   Δm ↔ Δτ/Δt branch alignment    |")
println("------------------------------------")

# Cross-contour consistency of the frozen kernels: the imaginary branch of the
# mixed (Kadanoff-Baym) contour correlation must reproduce the pure
# imaginary-time kernel, and the four real-time branch combinations must
# reproduce the pure real-time kernel blocks (up to the documented Keldysh
# relative signs, which differ between statistics).  All equalities hold to
# quadrature determinism — the same cell weights are evaluated on both sides.
# μ = 0 only: at μ ≠ 0 the τ sector uses different (but individually
# consistent) μ conventions in fermionic_Δτ (spectrumshift) and the Δm τ
# branch (direct μ in the kernel), so the two are not expected to coincide.
@testset "fermionic Δm ↔ Δτ/Δt alignment" begin
	β = 2
	for f in (semicircular(), DiracDelta(ω=0.5))
		begin
			μ = 0
			dm = fermionic_Δm(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)
			dτ = fermionic_Δτ(f, β=β, N=5, μ=μ)
			dt = fermionic_Δt(f, β=β, t=2, N=4, μ=μ)
			# τ branch ≡ imaginary-time kernel
			ξ = branch(dm, b1=:τ, b2=:τ)
			@test ξ.ηⱼₖ ≈ dτ.data.ηⱼₖ atol=1.0e-12
			@test ξ.ηₖⱼ ≈ dτ.data.ηₖⱼ atol=1.0e-12
			# real branches ≡ real-time kernel blocks — unified convention on
			# both sides, all four branches match directly (the historical
			# relative minus of the fermionic cross branches is now inside
			# the Δt blocks themselves)
			@test branch(dm, b1=:+, b2=:+).ηⱼₖ ≈ branch(dt, :+, :+).ηⱼₖ atol=1.0e-12
			@test branch(dm, b1=:+, b2=:+).ηₖⱼ ≈ branch(dt, :+, :+).ηₖⱼ atol=1.0e-12
			@test branch(dm, b1=:+, b2=:-).ηⱼₖ ≈ branch(dt, :+, :-).ηⱼₖ atol=1.0e-12
			@test branch(dm, b1=:+, b2=:-).ηₖⱼ ≈ branch(dt, :+, :-).ηₖⱼ atol=1.0e-12
			@test branch(dm, b1=:-, b2=:+).ηⱼₖ ≈ branch(dt, :-, :+).ηⱼₖ atol=1.0e-12
			@test branch(dm, b1=:-, b2=:+).ηₖⱼ ≈ branch(dt, :-, :+).ηₖⱼ atol=1.0e-12
			@test branch(dm, b1=:-, b2=:-).ηⱼₖ ≈ branch(dt, :-, :-).ηⱼₖ atol=1.0e-12
			@test branch(dm, b1=:-, b2=:-).ηₖⱼ ≈ branch(dt, :-, :-).ηₖⱼ atol=1.0e-12
		end
	end
end

@testset "bosonic Δm ↔ Δτ/Δt alignment" begin
	β = 2
	for f in (Leggett(), DiracDelta(ω=0.5))
		begin
			μ = 0
			dm = bosonic_Δm(f, β=β, t=2, Nt=4, Nτ=5, μ=μ)
			dτ = bosonic_Δτ(f, β=β, N=5, μ=μ)
			dt = bosonic_Δt(f, β=β, t=2, N=4, μ=μ)
			ξ = branch(dm, b1=:τ, b2=:τ)
			@test ξ.ηⱼₖ ≈ dτ.data.ηⱼₖ atol=1.0e-12
			@test ξ.ηₖⱼ ≈ dτ.data.ηₖⱼ atol=1.0e-12
			# unified convention on both sides: all four branches match
			@test branch(dm, b1=:+, b2=:+).ηⱼₖ ≈ branch(dt, :+, :+).ηⱼₖ atol=1.0e-12
			@test branch(dm, b1=:+, b2=:+).ηₖⱼ ≈ branch(dt, :+, :+).ηₖⱼ atol=1.0e-12
			@test branch(dm, b1=:+, b2=:-).ηⱼₖ ≈ branch(dt, :+, :-).ηⱼₖ atol=1.0e-12
			@test branch(dm, b1=:+, b2=:-).ηₖⱼ ≈ branch(dt, :+, :-).ηₖⱼ atol=1.0e-12
			@test branch(dm, b1=:-, b2=:+).ηⱼₖ ≈ branch(dt, :-, :+).ηⱼₖ atol=1.0e-12
			@test branch(dm, b1=:-, b2=:+).ηₖⱼ ≈ branch(dt, :-, :+).ηₖⱼ atol=1.0e-12
			@test branch(dm, b1=:-, b2=:-).ηⱼₖ ≈ branch(dt, :-, :-).ηⱼₖ atol=1.0e-12
			@test branch(dm, b1=:-, b2=:-).ηₖⱼ ≈ branch(dt, :-, :-).ηₖⱼ atol=1.0e-12
		end
	end
end
