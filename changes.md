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


## 5. 两矩阵极简存储（2026-10-30 第二轮）

块层面 4 个 Keldysh 块只有 2 个独立：**对角对 {G₊₊=(a,b), G₋₋=(c,d)} 携带全部
信息**，交叉块是它们两半的纯符号组合（无共轭 ⟹ 对任意复谱精确）：

    G₊₋ = CM(−G₋₋.ηⱼₖ, −G₊₊.ηₖⱼ)      G₋₊ = CM(−G₊₊.ηⱼₖ, −G₋₋.ηₖⱼ)

注意 {G₊₊, G₊₋} **不**完整（交叉块每半各与一个对角块共享，该组合丢失 G₋₊/G₋₋
的反向向量 d）。混合时间轴同构成立（branch 结构逐字相同）。

据此完成：
- `GenericRealCorrelationFunction`：只存 `G₊₊`/`G₋₋` 两个矩阵，交叉块 branch 派生；
  旧 4 矩阵构造器保留——校验约定一致性后只存对角对（不一致即报错，下游迁移期
  的防呆）。
- `MixedCorrelationFunction`：实时 4 向量（ηⱼₖ/ηₖⱼ/ηⱼₖ₂/ηₖⱼ₂）换成两个
  `CorrelationMatrix` 字段 `G₊₊`/`G₋₋`（同一极简表示）；8 参与 6 参（conj 派生、
  实谱）构造器保留，Δm 构造器无需改动。
- 回归：包测试 22 集全过；基准 CSV（含全部 mixedtime 与 dense 变体）逐位不变。


## 6. 单参数类与最终 API 定型（2026-10-30 第三轮）

- **`RealCorrelationFunction{A}` 单参数类**：删除紧凑四向量类型与
  `GenericRealCorrelationFunction` 名称——`A = CorrelationMatrix` 时两矩阵
  存储即四向量存储，紧凑型冗余。字段 `G₊₊::A`、`G₋₋::A`（同类型强制）。
  构造器：两矩阵（自动生成）；四 CorrelationMatrix（校验统一约定后只存对角
  对，不一致报错）；四一般矩阵（报错——一般矩阵无法校验/存储交叉接触分裂）。
- **一般矩阵（非 Toeplitz）支持**：`branch` 对角支直接返回存储块；交叉支按
  **三角拼接规则**派生（G₊₋ 下三角取 −G₋₋、上三角与对角取 −G₊₊；G₋₊ 对称），
  三角部分精确（任意关联泛函的四块共享严格三角），对角来自同支伙伴块（取负；
  满足 Keldysh 恒等式）。单元测试固化（test/realtime.jl "single-parametric
  type & general-matrix storage"，10 tests）。
- **`MixedCorrelationFunction{A}` 五块结构**：`G₊₊::A`、`G₋₋::A`（实时对角
  对）+ `G₊τ::AbstractMatrix{ComplexF64}`、`Gτ₊::AbstractMatrix{ComplexF64}`
  （**ζⱼₖ/ζₖⱼ 作为两个独立矩阵参数**——审计结论：一般情形独立，费米
  反对称仅在 μ=0 成立）+ `Gττ::CorrelationMatrix`（虚时块，ξ 封装）。
  五矩阵构造器（自动生成）；**八/六向量构造器删除**（Δm 与 deprecated
  构造点直接组矩阵）。branch(:−,τ) = −G₊τ、branch(:τ,−) = −Gτ₊。
- **下游消费面同步（本轮在包内完成）**：
  - Pfaffian 路径积分（test/quapi/bcsreal.jl 两处 + 头注释）：σ_u σ_v 因子
    **删除**（负号已在块内）——此前因容差宽未暴露，本轮重构后暴露并修复；
  - dense Grassmann IF（src/dense/fermion.jl）：实时源交叉块在指数侧**翻回
    旧号**（`_contour_delta_realtime`——Grassmann 指数约定使用 Keldysh 方位
    未吸收的块；混合轮廓直接消费原块，ED 校准不变）；重构中一度误删的
    `_delta_terms(::AbstractCorrelationFunction)` 通用入口已恢复。
- 回归：全部 10 个测试文件零失败；基准 CSV 全部逐位不变。


## 7. 存储改为物理基本对（交叉块 G₊₋/G₋₊）（2026-10-30 第四轮）

- **`RealCorrelationFunction{A}` 字段 = `G₊₋`、`G₋₊`**（lesser/greater 型关联
  ——Keldysh 的两个物理基本函数），对角块由两个共用推导函数恢复：
  `_G₊₊(G₊₋, G₋₊) = CM(−G₋₊.ηⱼₖ, −G₊₋.ηₖⱼ)`、
  `_G₋₋(G₊₋, G₋₊) = CM(−G₊₋.ηⱼₖ, −G₋₊.ηₖⱼ)`（纯符号、无共轭，任意复谱精确；
  一般矩阵走三角拼接重载）。`Real` 与 `Mixed` 的 branch 均直接调用。
- **`MixedCorrelationFunction{A}`** 五块签名不变但实时对 = `G₊₋`/`G₋₊`；
  `transpose` 相应**交换两交叉块**再各自转置（[xᵀ]^{+−} = (x^{−+})ᵀ）。
- `fermionic_Δt` 存 `CM(−c,−b)`/`CM(−a,−d)`；`bosonic_Δt` 存
  `CM(−conj(b),−b)`/`CM(−a,−conj(a))`；Δm/deprecated 构造器同步。
- branch 语义不变（数值逐位不变，全部 CSV 验证）；交叉支现在返回**存储对象
  本身**（下游原地修改会穿透——测试改用 deepcopy）。
- performance 内部一处 `mix.ξⱼₖ` 直接字段访问改为 `branch(mix; :τ,:τ).ηⱼₖ`。


## 8. 四矩阵构造器放开 + Mixed 字段固化（2026-10-30 第五轮）

- **四矩阵构造器（任意矩阵类型）**：放开 `RealCorrelationFunction(G₊₊, G₊₋,
  G₋₊, G₋₋)` 的一般矩阵版本——**扇区自洽检查**（四块严格三角两两配对一致 +
  对角满足 Keldysh 恒等式 G₊₊+G₋₋+G₊₋+G₋₊ = 0）通过后直接调用两输入构造器
  存交叉对；不自洽抛 ArgumentError。CorrelationMatrix 输入在向量级精确检查
  （原有行为）；一般矩阵在扇区级检查（零矩阵平凡自洽、合法）。单元测试同步：
  零矩阵接受、随机矩阵抛错。
- **`MixedCorrelationFunction` 的 `G₊τ`/`Gτ₊` 固化为 `Matrix{ComplexF64}`**；
  `transpose` 对两矩阵显式 `Matrix(transpose(·))` 收集（`transpose(Matrix)`
  返回 `Transpose` 包装器，否则类型不匹配）。
- 回归：全部 10 测试文件零失败、全部基准 + dense 变体 CSV 逐位不变。
