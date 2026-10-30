# changes.md —— 实时冻结核类型与统一 branch 约定改造（2026-10-30）

## 1. 新类型结构

### `RealCorrelationFunction`（紧凑型，`fermionic_Δt` / `bosonic_Δt` 的默认返回值）

```julia
struct RealCorrelationFunction <: AbstractRealCorrelationFunction
    ηⱼₖ  :: Vector{ComplexF64}   # (+,+) 前向滞后，热权重 g₁
    ηₖⱼ  :: Vector{ComplexF64}   # (+,+) 反向滞后，热权重 g₂
    ηⱼₖ₂ :: Vector{ComplexF64}   # (+,−) 前向滞后，热权重 g₂
    ηₖⱼ₂ :: Vector{ComplexF64}   # (−,+) 反向滞后，热权重 g₁
end
```

- **四向量构造器**（复谱安全）：`RealCorrelationFunction(ηⱼₖ, ηₖⱼ, ηⱼₖ₂, ηₖⱼ₂)`。
  `fermionic_Δt` 无条件走这条路（ηⱼₖ₂/ηₖⱼ₂ 由独立 quadgk 求积、谱权重**不共轭**）
  ——BCS 反常分量（复 Δ 的 uv/vu 相干权重）因此天然正确，与 `MixedCorrelationFunction`
  的补全存储完全同构。
- **双向量构造器**（仅实值谱）：`RealCorrelationFunction(ηⱼₖ, ηₖⱼ)` 以
  `ηⱼₖ₂ = conj(ηₖⱼ)`、`ηₖⱼ₂ = conj(ηⱼₖ)` 派生。`bosonic_Δt` 走这条路（玻色谱恒实）。
- 四块 branch 的 Toeplitz 结构：8 个存储向量里只有 4 个独立，
  `G₊₊ = CM(a,b)`、`G₊₋ = CM(−c,−b)`、`G₋₊ = CM(−a,−d)`、`G₋₋ = CM(c,d)`。

### `GenericRealCorrelationFunction`（全四块显式存储）

即旧 `RealCorrelationFunction` 的原样重命名；保留旧四矩阵构造器
`RealCorrelationFunction(G₊₊, G₊₋, G₋₊, G₋₋) -> GenericRealCorrelationFunction`
（deprecated 实现等旧代码不受影响）。

### `AbstractRealCorrelationFunction`

两个实现的公共抽象父类型；`branch`/`branches`/`contour_ordering`/`transpose`/
`+`/`-`/`==`/`isapprox` 全部定义在该层（块级比较，跨表示可用）。

## 2. 统一 branch 约定（核心语义变更）

**之前**：费米 Δt 交叉支不带 Keldysh 相对负号（G₊₋ = CM(c, b)），玻色 Δt 与
`MixedCorrelationFunction` 的交叉支带负号——同库两套约定。

**现在**：全库统一为 Δm/玻色那套（交叉支负号在块内）：

| branch | 统一约定（费米 = 玻色） | 旧费米 Δt |
|---|---|---|
| (+,+) | CM(a, b) | CM(a, b)（不变） |
| (+,−) | **CM(−c, −b)** | CM(c, b)（**翻号**） |
| (−,+) | **CM(−a, −d)** | CM(a, d)（**翻号**） |
| (−,−) | CM(c, d) | CM(c, d)（不变） |

`MixedCorrelationFunction`（branch(:+,:−) = CM(−ηⱼₖ₂, −ηₖⱼ)）本就是这套约定，
无需改动；`InfiniteRealCorrelationCache` 的费米 `_branch` 已同步翻号统一。

## 3. 下游需要的调整

凡是读费米 `fermionic_Δt` / `bcs_Δt` 的 **(+,−) 或 (−,+) 块**的下游代码，
数值整体翻了一号。等价改写规则（数值逐位不变，performance 套件已全部照此更新
并 CSV 逐位验证）：

| 场景 | 旧公式 | 新公式 |
|---|---|---|
| 费米/BCS 正常槽推迟核 | `X − P` | `−(X + P)` |
| BCS 反常槽推迟核 | `P − X` | `P + X` |
| Pfaffian 路径积分 σ 因子 | 交叉支乘 σ_u σ_v = −1 | **去掉** σ 因子（负号已在块内） |
| 与 deprecated 对拍 | 直接比较 | 交叉支乘 −1 后比较 |

不读交叉支的下游（只消费 (+,+)/τ 支/Δτ/Δm）完全不受影响；玻色侧全部不受影响。

## 4. 已完成的包内更新（供对照）

- `src/realtime/realcorr.jl`：新类型体系 + 统一 branch；`src/realtime/{fermionic,bosonic}.jl`：
  构造器重写（费米无条件四向量、玻色双向量 conj 派生——玻色 quadgk 次数减半）；
  `src/infinitetime/realtime.jl`：费米 _branch 统一；`src/dense`、`src/QuAPI.jl` 导出更新。
- 测试：`test/realtime.jl`、`test/mixedtime.jl`、`test/bcs.jl` 的对齐断言改为统一
  约定下的直接相等（费米 Δm↔Δt 的交叉支相对负号消失）；`test/quapi/bcsreal.jl`
  Pfaffian 的 σ 因子待下游同步时一并去掉（其当前状态在统一约定下依然通过是因为
  容差较宽，见该项 TODO）。
- performance：全部费米/BCS 基准的组合公式照 §3 规则改写；21 个 CSV 逐位不变。
