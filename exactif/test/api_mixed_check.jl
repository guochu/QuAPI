# 混合时间 API 的锚定检验
include("exactif_api.jl")
using LinearAlgebra, Test
function main()
    Γ0, Wb, β = 0.15, 3.0, 2.0
    f = bounded(ω -> Γ0, -Wb, Wb)
    Nτ, Nt, t = 6, 4, 0.4
    Δm = exact_fermionic_Δm(f; β=β, Nτ=Nτ, t=t, Nt=Nt, M=48)
    Δt_, _ = exact_fermionic_Δt(f; β=β, N=Nt, t=t, M=48)
    Δτ = exact_fermionic_Δτ(f; β=β, N=Nτ, M=48)
    @testset "mixed API 锚定" begin
        # 1) Keldysh 交叉块 = 实时版（一般 Matrix 存储，直接矩阵对比）
        @test maximum(abs.(Δm.G₊₋ - Δt_.G₊₋)) < 1e-13   # 同构造（求积顺序差 ~1e-16）
        @test maximum(abs.(Δm.G₋₊ - Δt_.G₋₊)) < 1e-13
        # 2) ττ 块 = 虚时版
        @test Δm.Gττ.ηⱼₖ == Δτ.data.ηⱼₖ
        @test Δm.Gττ.ηₖⱼ == Δτ.data.ηₖⱼ
        # 3) μ=0 费米反对称：Gτ₊ = −G₊τᵀ
        @test maximum(abs.(Δm.Gτ₊ .+ transpose(Δm.G₊τ))) < 1e-10
        # 4) branch 规则：(−,τ) = −G₊τ；(τ,−) = −Gτ₊
        @test branch(Δm; b1=:-, b2=:τ) == .-Δm.G₊τ
        @test branch(Δm; b1=:τ, b2=:-) == .-Δm.Gτ₊
        # 5) 尺寸
        @test size(Δm.G₊τ) == (Nt + 1, Nτ)
        @test size(Δm.Gτ₊) == (Nτ, Nt + 1)
        # 6) V→0：τ-交叉 → ½（杂质自线，系统腿——调用方测度抵消）
        Δ0 = exact_fermionic_Δm(bounded(ω -> 0.0, -Wb, Wb); β=β, Nτ=3, t=0.2, Nt=2, M=8)
        @test maximum(abs.(Δ0.G₊τ .- 0.5)) < 1e-12
        @test maximum(abs.(Δ0.Gτ₊ .+ 0.5)) < 1e-12
    end
    println("mixed API 锚定全部通过 ✓")
end
main()
