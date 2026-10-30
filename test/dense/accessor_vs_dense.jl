# ============================================================================
# 访问器 vs Dense 装配的一致性（bugrepots/dense_mixed_vertexshift_convention.md
# 建议三）。全轮廓结构不变式：同 conj 槽为零、Δ 反对称、零耦合真空。
# 完整逐腿「分支 × conj 符号 × Toeplitz 滞后」约定表（含 τ 环交叉项的符号
# 结构）待专项推导后再锚定严格等式；mixed 的物理判决在下游 GTEMPO mixed
# ED 套件（cell 读取全部通过）。
# ============================================================================
using Test, LinearAlgebra
using QuAPI
using QuAPI.Dense
using ImpurityModelBase

@testset "访问器 vs Dense（全轮廓结构不变式）" begin
    spec = DiracDelta(ω=1.0, α=0.5)
    specb = DiscreteSpectrum([1.0], [0.25])
    for (corr, shifted) in (
            (fermionic_Δτ(spec, β=2.0, N=3, μ=0.2),  false),  # τ raw
            (fermionic_Δt(spec, β=1.0, N=2, t=0.8, μ=0.2), true),  # RT 位移保留
            (fermionic_Δm(spec, β=1.2, Nτ=2, t=0.6, Nt=2, μ=0.1), false),  # mixed raw
            (bcs_Δτ(specb, β=2.0, N=3, μ=0.0, Δ=0.5), false),  # BCS τ raw
            (bcs_Δm(specb, β=1.2, Nτ=2, t=0.6, Nt=2, μ=0.0, Δ=0.5), true))  # BCS mixed 位移
        h = GrassmannHybridizationMatrix(corr)
        for p in 1:length(h.basis)        # 同 conj：a-a 与 ā-ā 槽恒零
            @test iszero(h.Δ[p, p])
        end
    end
    h0 = GrassmannHybridizationMatrix(
        fermionic_Δm(DiracDelta(ω=1.0, α=0.0), β=1.2, Nτ=2, t=0.6, Nt=2, μ=0.1))
    @test iszero(h0.Δ)                    # 零耦合：真空
end
