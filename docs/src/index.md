# QuAPI

QuAPI（Quasi-Adiabatic Propagator Path Integral）用于玻色子和费米子杂化函数的
离散化，以及在此基础上的准绝热传播子路径积分。文档结构：

- **[API](api.md)**（`src/`，导出）：把连续谱的杂化函数（influence
  functional 的核）离散为虚时、实时、混合时（Kadanoff L 形回路）三种时间
  网格上的 `CorrelationMatrix`，支持玻色、费米以及 BCS 超导 bath。
- **[Hybridization Functions](hybridization.md)**：Δτ、Δt、Δm 及 BCS
  杂化函数的完整数学定义——Toeplitz 结构、每个生成元素的谱积分表达式与
  时间核的封闭形式。

## 全局符号约定

影响泛函（IF）指数中的负号被吸收进杂化矩阵的定义（见
`src/QuAPI.jl` 模块注释），即离散路径权重为

```
I = exp[ +Σᵢⱼ ψ̄ᵢ Δᵢⱼ ψⱼ ]        (费米)
I = exp[ +Σᵢⱼ sᵢ  Δᵢⱼ sⱼ ]        (玻色自旋路径)
```

而非常见的 `exp[−Σ ψ̄ Δ ψ]`。使用外部文献公式时注意相差一个负号。

## 安装

QuAPI 依赖 [ImpurityModelBase](https://github.com/guochu/ImpurityModelBase)
（谱函数、bath 类型、ED 参考解）。本地开发布局要求几个包并列放置：

```
Missile/
├── QuAPI/
├── ImpurityModelBase/
├── GrassmannTensors/     # 测试间接依赖（可选）
└── HadamardProduct/      # 测试间接依赖（可选）
```

`Project.toml` / `Manifest.toml` 中的本地依赖用**相对路径**记录
（如 `path = "../ImpurityModelBase"`），因此在任何机器上克隆后布局一致即可
直接使用，无需修改路径。

加载包：

```julia
using QuAPI
```

## 目录结构

```
src/
├── QuAPI.jl              # 模块入口与全部 export
├── util.jl               # Fermi/Bose 因子（私有）
├── corrmat.jl            # CorrelationMatrix（Toeplitz）
├── abstractcorr.jl       # AbstractCorrelationFunction / index / branch
├── imagtime/             # 虚时：fermionic/bosonic Δτ、Δiw_to_Δτ
├── realtime/             # 实时：fermionic/bosonic Δt（Keldysh 四块）
├── mixedtime/            # 混合时：fermionic/bosonic Δm（九块）
├── infinitetime/         # 无穷温窗：infinite_Δτ / infinite_Δt
├── bcs/                  # BCS：bcs_Δτ / bcs_Δt / bcs_Δm
│   └── infinitetime/     # 占位（未来 BCS 无穷温窗）
└── deprecated/           # 旧实现（交叉验证用，未导出）
test/
├── runtests.jl           # 入口
├── imagtime/realtime/mixedtime.jl, bcs.jl    # src 层测试
└── quapi/                # QuAPI 算法实现 + ED 对比
docs/
└── src/                  # 本文档源
```

## 运行测试

```bash
# 方式一：直接运行（推荐，绕过 Pkg.test 的沙箱环境解析）
julia --project=test test/runtests.jl

# 方式二：标准入口
julia --project=. -e 'using Pkg; Pkg.test()'
```

当前 25 个 testset 全部通过，覆盖：虚/实/混合时费米与玻色 Δ、无穷温窗、
`CorrelationMatrix` 代数、`Δiw_to_Δτ`、BCS 退化极限、以及全部 QuAPI 算法
的 ED / 解析解对比。
