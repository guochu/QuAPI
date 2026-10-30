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

@testset "BCS bath convenience" begin
	atol = 1.0e-8
	β = 5
	μ = 0
	Δ = 0.5
	for f in (semicircular(), DiracDelta(ω=0.5))
		bath = bcsbath(f, β=β, μ=μ, Δ=Δ)

		x = Δτ(bath, N=10)
		y = bcs_Δτ(f, β=β, N=10, μ=μ, Δ=Δ)
		@test x.uu ≈ y.uu atol=atol
		@test x.ud ≈ y.ud atol=atol
		@test x.du ≈ y.du atol=atol
		@test x.dd ≈ y.dd atol=atol

		x = Δt(bath, N=10, t=2)
		y = bcs_Δt(f, β=β, N=10, t=2, μ=μ, Δ=Δ)
		@test x.uu ≈ y.uu atol=atol
		@test x.ud ≈ y.ud atol=atol
		@test x.du ≈ y.du atol=atol
		@test x.dd ≈ y.dd atol=atol

		x = Δm(bath, Nt=4, Nτ=5, t=2)
		y = bcs_Δm(f, β=β, Nt=4, Nτ=5, t=2, μ=μ, Δ=Δ)
		@test x.uu ≈ y.uu atol=atol
		@test x.ud ≈ y.ud atol=atol
		@test x.du ≈ y.du atol=atol
		@test x.dd ≈ y.dd atol=atol
	end
end
println("------------------------------------")
println("|   BCS Δm ↔ Δτ/Δt alignment       |"
)
println("------------------------------------")

# Cross-contour consistency for the superconducting bath: per flavor block, the
# imaginary branch of bcs_Δm must reproduce bcs_Δτ and the real branch
# combinations must reproduce bcs_Δt (fermionic Keldysh relative signs).  The
# complex-Δ case guards the un-conjugated backward-sector storage: with a
# conjugated backward route the anomalous (ud/du) components of the (+,−),
# (−,+) and (−,−) combinations would pick up conj(weight) and fail here.
@testset "bcs Δm ↔ Δτ/Δt alignment" begin
	atol = 1.0e-12
	β, μ = 2, 0
	for Δ in (0.5, 0.2 + 0.45im)
		for f in (semicircular(), DiracDelta(ω=0.5))
			dm = bcs_Δm(f, β=β, Nτ=5, t=2, Nt=4, μ=μ, Δ=Δ)
			dτ = bcs_Δτ(f, β=β, N=5, μ=μ, Δ=Δ)
			dt = bcs_Δt(f, β=β, N=4, t=2, μ=μ, Δ=Δ)
			for nm in (:uu, :ud, :du, :dd)
				mblk, τblk, tblk = getproperty(dm, nm), getproperty(dτ, nm), getproperty(dt, nm)
				# τ branch ≡ imaginary-time kernel
				ξ = branch(mblk, b1=:τ, b2=:τ)
				@test ξ.ηⱼₖ ≈ τblk.data.ηⱼₖ atol=atol
				@test ξ.ηₖⱼ ≈ τblk.data.ηₖⱼ atol=atol
				# real branches ≡ real-time kernel blocks (fermionic Keldysh signs)
				@test branch(mblk, b1=:+, b2=:+).ηⱼₖ ≈ branch(tblk, :+, :+).ηⱼₖ atol=atol
				@test branch(mblk, b1=:+, b2=:+).ηₖⱼ ≈ branch(tblk, :+, :+).ηₖⱼ atol=atol
				@test branch(mblk, b1=:+, b2=:-).ηⱼₖ ≈ branch(tblk, :+, :-).ηⱼₖ atol=atol
				@test branch(mblk, b1=:+, b2=:-).ηₖⱼ ≈ branch(tblk, :+, :-).ηₖⱼ atol=atol
				@test branch(mblk, b1=:-, b2=:+).ηⱼₖ ≈ branch(tblk, :-, :+).ηⱼₖ atol=atol
				@test branch(mblk, b1=:-, b2=:+).ηₖⱼ ≈ branch(tblk, :-, :+).ηₖⱼ atol=atol
				@test branch(mblk, b1=:-, b2=:-).ηⱼₖ ≈ branch(tblk, :-, :-).ηⱼₖ atol=atol
				@test branch(mblk, b1=:-, b2=:-).ηₖⱼ ≈ branch(tblk, :-, :-).ηₖⱼ atol=atol
			end
		end
	end
end
