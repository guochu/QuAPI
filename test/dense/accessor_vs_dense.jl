# ============================================================================
# 访问器 vs Dense 装配的一致性（bugrepots/dense_mixed_vertexshift_convention.md
# 建议三）。全轮廓结构不变式：同 conj 槽为零、零耦合真空；逐元素对齐：
# Dense 的 Δ 放置（(conj,annihil) 规范位型 + RT 跨支 Keldysh 负号）必须与
# 访问器读数逐位相等——访问器是 ED 判决后的唯一约定基准（2026-10-30 范围
# 收缩、2026-11-04 BCS RT 对齐后全轮廓一致：仅实时 Keldysh 位移）。
# ============================================================================
using Test, LinearAlgebra
using QuAPI
using QuAPI.Dense
using ImpurityModelBase

@testset "访问器 vs Dense（全轮廓结构不变式）" begin
    spec = DiracDelta(ω=1.0, α=0.5)
    specb = DiscreteSpectrum([1.0], [0.25])
    for corr in (fermionic_Δτ(spec, β=2.0, N=3, μ=0.2),
                 fermionic_Δt(spec, β=1.0, N=2, t=0.8, μ=0.2),
                 fermionic_Δm(spec, β=1.2, Nτ=2, t=0.6, Nt=2, μ=0.1),
                 bcs_Δτ(specb, β=2.0, N=3, μ=0.0, Δ=0.5),
                 bcs_Δt(specb, β=1.0, N=2, t=0.8, μ=0.0, Δ=0.5),
                 bcs_Δm(specb, β=1.2, Nτ=2, t=0.6, Nt=2, μ=0.0, Δ=0.5))
        h = GrassmannHybridizationMatrix(corr)
        for p in 1:length(h.basis)        # 同 conj：a-a 与 ā-ā 槽恒零
            @test iszero(h.Δ[p, p])
        end
    end
    h0 = GrassmannHybridizationMatrix(
        fermionic_Δm(DiracDelta(ω=1.0, α=0.0), β=1.2, Nτ=2, t=0.6, Nt=2, μ=0.1))
    @test iszero(h0.Δ)                    # 零耦合：真空
end

@testset "访问器 vs Dense（逐元素对齐）" begin
    spec = DiracDelta(ω=1.0, α=0.5)
    specb = DiscreteSpectrum([1.0], [0.25])
    # 普通费米 bath：Dense 在 (conj,annihil) 规范位型放置块，与访问器 a-ā 读数
    # 同位型同符号，其余位型恒零。跨支 Keldysh 负号只属于纯实时轮廓
    # （(branches == (:+,:-))；混合轮廓按 ED 校准的历史约定直接消费原块，
    # 不带任何附加符号——与 partialif 的读数条件一致）。
    for corr in (fermionic_Δτ(spec, β=2.0, N=3, μ=0.2),            # τ 环 raw
                 fermionic_Δt(spec, β=1.0, N=2, t=0.8, μ=0.2),     # RT 位移
                 fermionic_Δm(spec, β=1.2, Nτ=2, t=0.6, Nt=2, μ=0.1))  # mixed raw
        h = GrassmannHybridizationMatrix(corr)
        rt = QuAPI.branches(corr) == (:+, :-)
        for p in 1:length(h.basis), q in 1:length(h.basis)
            lp, lq = h.basis[p], h.basis[q]
            expected = lp.conj && !lq.conj ?
                corr[GrassmannContourIndex(lp.j, lp.branch, lp.band, false),
                     GrassmannContourIndex(lq.j, lq.branch, lq.band, true)] *
                (rt && lp.branch !== lq.branch ? -1 : 1) :
                zero(eltype(h.Δ))
            @test h.Δ[p, q] ≈ expected atol=1e-13
        end
    end
    # BCS bath：访问器为有向指数系数（(ā↑,a↑)/(ā↓,a↓)/(ā↑,ā↓)/(a↓,a↑) 四类
    # 规范位型，其余零），Dense 的 _delta_terms 同一分解，逐位相等。
    for corrb in (bcs_Δτ(specb, β=2.0, N=3, μ=0.0, Δ=0.5),             # τ 环 raw
                  bcs_Δt(specb, β=1.0, N=2, t=0.8, μ=0.0, Δ=0.5),      # BCS RT 位移
                  bcs_Δm(specb, β=1.2, Nτ=2, t=0.6, Nt=2, μ=0.0, Δ=0.5))  # BCS mixed raw
        h = GrassmannHybridizationMatrix(corrb)
        for p in 1:length(h.basis), q in 1:length(h.basis)
            @test h.Δ[p, q] ≈ corrb[h.basis[p], h.basis[q]] atol=1e-13
        end
    end
end
