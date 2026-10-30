# API 契约 + 物理检验（模式离散版，M = 求积节点数）
include(joinpath(@__DIR__, "exactif_api.jl"))
using LinearAlgebra, Test
function main()
    Γ0, Wb, β = 0.15, 3.0, 2.0
    f = bounded(ω -> Γ0, -Wb, Wb)
    @testset "exactif_api QuAPI 契约" begin
        # 虚时：反周期闭合逐位 + 接触值
        Δτ = exact_fermionic_Δτ(f; β=β, N=8)
        ηⱼₖ = Δτ.data.ηⱼₖ
        @test length(ηⱼₖ) == 8
        @test all(imag.(ηⱼₖ) .== 0)
        @test maximum(abs.(Δτ.data.ηₖⱼ[2:end] .+ reverse(ηⱼₖ[2:end]))) == 0
        # G(0⁻) = −⟨d†d⟩ = −Σv²f(λ)：对称谱 ⟹ −½（精确），GL 求积偏差由
        # M-收敛控制
        @test isapprox(real(Δτ.data.ηₖⱼ[1]), -0.5; atol=5e-2)   # GL 端点求积
        Δτ96 = exact_fermionic_Δτ(f; β=β, N=8, M=96)
        @test isapprox(real(Δτ.data.ηₖⱼ[1]), real(Δτ96.data.ηₖⱼ[1]); atol=1e-2)
        # 实时：V→0 ⟹ Δ=0（纯浴语义）
        _, b0v = exact_fermionic_Δt(bounded(ω -> 0.0, -Wb, Wb); β=β, N=2, t=0.4)
        @test maximum(abs.(b0v.Gpp)) < 1e-12
        # M-收敛（求积节点 → 连续谱）
        _, bA = exact_fermionic_Δt(f; β=β, N=2, t=0.4, M=48)
        _, bB = exact_fermionic_Δt(f; β=β, N=2, t=0.4, M=96)
        @test maximum(abs.(bA.Gpp - bB.Gpp)) < 1e-3
        # μ 平移等价（谱整体平移；GL 节点随移 ⟹ 求积容差）
        Δμ, _ = exact_fermionic_Δt(f; β=β, N=2, t=0.4, μ=0.3, M=32)
        Δh, _ = exact_fermionic_Δt(bounded(ω -> Γ0, -Wb+0.3, Wb+0.3); β=β, N=2, t=0.4, μ=0.0, M=32)
        @test maximum(abs.(Δμ.G₊₋ - Δh.G₊₋)) < 5e-2
    end
    println("api_consistency 全部通过 ✓")
end
main()
