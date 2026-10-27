# QuAPI

准绝热传播子路径积分（**Qu**asi-**A**diabatic **P**ropagator path integral）——
把 Feynman-Vernon 影响泛函的杂化函数离散为虚时（Matsubara）、实时（Keldysh）
和混合时（Kadanoff L 形回路）网格上的 Toeplitz 矩阵，供自旋-玻色、
Toulouse（共振能级）与 BCS 模型等做路径积分模拟。

**符号约定**：影响泛函指数中的负号被吸收进杂化矩阵的定义，
离散路径权重为 ``I = \exp[+\sum_{ij}\bar\psi_i\,\Delta_{ij}\,\psi_j]``
（费米）或 ``I=\exp[+\sum_{ij}s_i\,\Delta_{ij}s_j]``（玻色自旋路径），
而非常见的 ``\exp[-\sum\bar\psi\,\Delta\,\psi]``。

## 功能

- **虚时** `fermionic_Δτ` / `bosonic_Δτ`：β=∞ 的无穷温窗版本
  `infinite_Δτ`；Matsubara 频域入口 `Δiw_to_Δτ`；
- **实时**（闭合 Keldysh 回路，四块 G⁺⁺/G⁺⁻/G⁻⁺/G⁻⁻）`fermionic_Δt` /
  `bosonic_Δt` / `infinite_Δt`；
- **混合时**（L 形 Kadanoff-Baym 回路，九块）`fermionic_Δm` /
  `bosonic_Δm`；
- **BCS 超导 bath**（Bogoliubov 相干因子，正常 + 反常块）`bcs_Δτ` /
  `bcs_Δt` / `bcs_Δm`；
- 所有离散杂化矩阵以 O(N) 的 Toeplitz `CorrelationMatrix` 存储，支持
  `+ - transpose` 及逐元素索引。

## 安装与环境

本包依赖 [ImpurityModelBase](../ImpurityModelBase)（谱函数、bath 类型与
严格对角化解），并与它并列开发。克隆后保持布局一致即可直接使用：

```
Missile/
├── QuAPI/
├── ImpurityModelBase/
```

> `Manifest.toml` 不入库：每个环境（根目录 / `test/` / `docs/`）首次使用时
> 运行 `julia --project=. -e 'using Pkg; Pkg.resolve()'` 生成各自的
> Manifest。本地包已用相对路径记录在 `Project.toml` 的 `[deps]`/dev 源中。

加载：

```julia
using QuAPI
```

## 快速上手

谱函数与 bath 类型来自 ImpurityModelBase：

```julia
using QuAPI

β, μ = 5.0, 0.3
f = semicircular(t = 1.0)              # 费米谱（Bethe）
j = Leggett(α = 0.1, ωc = 1.0)         # 玻色谱（Ohmic）

fbath = fermionicbath(f; β = β, μ = μ)
bbath = bosonicbath(j; β = β, μ = 0.0)
bcs   = bcsbath(f; β = β, Δ = 0.5, μ = μ)
```

### 虚时（Matsubara）

```julia
corr = Δτ(fbath; N = 10)               # δτ = β/N，bath 便捷包装
corr.data[1, 1]                         # 等时对角元
corr.data[3, 1]                         # 平移不变矩阵元
```

### 实时（Keldysh）

```julia
corr = Δt(bbath; N = 20, t = 2.0)      # δt = t/N，节点 0..N
Gpp = branch(corr; b1 = :+, b2 = :+)   # 时序块 (N+1)×(N+1)
Gpm = branch(corr; b1 = :+, b2 = :-)   # lesser 型
```

### 混合时（Kadanoff L 形回路）

```julia
corr = Δm(bbath; Nτ = 10, Nt = 4, t = 2.0)
Gττ = branch(corr; b1 = :τ, b2 = :τ)   # 虚时块 Nτ×Nτ
Gpt = branch(corr; b1 = :+, b2 = :τ)   # 混合块 (Nt+1)×Nτ
```

### BCS

```julia
corr = bcs_Δτ(f; β = β, N = 10, Δ = 0.5, μ = μ)
corr.uu    # 正常块 ā↑Δ↑↑a↑（= cn）
corr.ud    # 反常块 ā↑Δ↑↓ā↓（= cc）
```

`Δ = 0` 时 BCS 结果退化为普通费米结果。

## 数学定义

Δτ、Δt、Δm 及 BCS 杂化函数的完整数学定义（Toeplitz 结构、时间核的封闭
形式、每条生成向量的谱积分表达式、Keldysh/九块装配）见
[`docs/src/hybridization.md`](docs/src/hybridization.md)，或用
[Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) 构建 HTML：

```bash
julia --project=docs docs/make.jl
```

本地实时预览（改 `docs/src/*.md` 后自动重建、浏览器自动刷新；首次需先
`julia --project=docs -e 'using Pkg; Pkg.add("LiveServer")'`）：

```bash
julia --project=docs docs/liveserver.jl
```

## 运行测试

```bash
# 直接运行（绕过 Pkg.test 的环境解析）
julia --project=test test/runtests.jl

# 或标准入口
julia --project=. -e 'using Pkg; Pkg.test()'
```

测试覆盖杂化矩阵与旧实现（`src/deprecated/`）的交叉验证、BCS 退化极限，
以及参考路径积分与严格对角化 / 解析解的对比。

## License

见 [LICENSE](LICENSE)。
