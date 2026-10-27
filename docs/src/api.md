# API

本页列出 `QuAPI` 的全部导出 API。函数签名与说明直接取自源码 docstring。

## 核心类型与接口

```@docs
CorrelationMatrix
AbstractCorrelationFunction
branch
branches
```

## 虚时（Matsubara）

```@docs
ImagCorrelationFunction
Δτ
fermionic_Δτ
bosonic_Δτ
Δiw_to_Δτ
```

## 实时（Keldysh）

```@docs
RealCorrelationFunction
Δt
fermionic_Δt
bosonic_Δt
```

## 混合时（Kadanoff L 形回路）

```@docs
AbstractMixedCorrelationFunction
MixedCorrelationFunction
Δm
fermionic_Δm
bosonic_Δm
```

## 无穷温窗

```@docs
infinite_Δτ
infinite_Δt
QuAPI.InfiniteImagCorrelationCache
QuAPI.InfiniteRealCorrelationCache
QuAPI.current_size
QuAPI.current_temperature
QuAPI.current_time
QuAPI.compute!
QuAPI.compute_next!
```

## BCS 超导 bath

```@docs
BCSCorrelationFunction
bcs_Δτ
bcs_Δt
bcs_Δm
```

## 稠密 IF（QuAPI.Dense）

```@docs
QuAPI.Dense.MAXTENSORDIM
QuAPI.Dense.GrassmannADT
QuAPI.Dense.GrassmannHybridizationMatrix
QuAPI.Dense.GrassmannIF
QuAPI.Dense.BosonicADT
QuAPI.Dense.BosonicPT
QuAPI.Dense.BosonicSplitDelta
QuAPI.Dense.BosonicHybridizationMatrix
QuAPI.Dense.BosonicIF
QuAPI.Dense.Zvalue
QuAPI.Dense.margin
QuAPI.Dense.complementindex
QuAPI.Dense.contour_ordering
QuAPI.Dense.basis
```

## 废弃实现（未导出）

旧版本内核仅用于交叉验证（新实现与之一致性由测试保证），命名带
`_deprecated` 后缀：`fermionic_Δτ_deprecated`、`bosonic_Δτ_deprecated`、
`bosonic_Δt_deprecated`、`fermionic_Δm_deprecated`、`bosonic_Δm_deprecated`
及 bath 便捷方法 `Δτ_deprecated` / `Δt_deprecated` / `Δm_deprecated`。
使用时需显式引入：

```julia
using QuAPI: fermionic_Δτ_deprecated
```

```@docs
QuAPI.fermionic_Δτ_deprecated
QuAPI.bosonic_Δτ_deprecated
QuAPI.bosonic_Δt_deprecated
QuAPI.fermionic_Δm_deprecated
QuAPI.bosonic_Δm_deprecated
```
