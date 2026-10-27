# Dense IF Assembly

本页讲三种稠密影响泛函——费米 [`QuAPI.Dense.GrassmannADT`](@ref)、玻色对角
[`QuAPI.Dense.BosonicADT`](@ref) 与玻色非对角 [`QuAPI.Dense.BosonicPT`](@ref)——如何从
[Hybridization Functions](@ref "The Hybridization Functions") 章的离散杂化函数
**组装**出来：腿/节点的排列约定、「杂化矩阵元 ⟶ 张量元」的映射公式，以及
组装引入的误差标度。记号：谱强度记作 ``\Gamma`` 或 ``\alpha``（即 ``V^{2}``，
``V`` 为系统–库耦合项的系数；谱宽记 ``D``）；误差标度以耦合系数 ``V`` 计——
``O(V^{2})`` 即「线性于谱强度」阶。

## 1. 三种稠密表示

装配都是两步：先**杂化矩阵**（把 `CorrelationMatrix` 的 Toeplitz 生成向量铺
到 contour 腿/节点上，`GrassmannHybridizationMatrix` /
`BosonicHybridizationMatrix`），再**折叠成稠密张量**（`GrassmannIF` /
`BosonicIF`）：

| 类型 | 耦合 | `basis` | `data` | 构造 |
|---|---|---|---|---|
| [`QuAPI.Dense.GrassmannADT`](@ref) | 费米（含 BCS） | `L` 条 contour 腿 | `2^L` 个 Grassmann 单项式系数 | `GrassmannIF` ∘ `GrassmannHybridizationMatrix` |
| [`QuAPI.Dense.BosonicADT`](@ref) | 玻色对角（`hyb` 对角化） | `n` 个节点 | `d^n` 个路径权重 | `BosonicIF(corr, λ)`，`λ` = 本征值 |
| [`QuAPI.Dense.BosonicPT`](@ref) | 玻色非对角（任意 `hyb`） | `n` 个节点 | `d^{2n}` 个矩阵元 | `BosonicIF(corr, hyb)`（内部 `split=true`） |

```julia
using QuAPI, QuAPI.Dense

h  = GrassmannHybridizationMatrix(corr)   # 腿基 basis + 配对矩阵 Δ
I  = GrassmannIF(h)                       # GrassmannADT
Z  = Zvalue(I)                            # Berezin 积分（全腿单项式系数）

adt = BosonicIF(bc, [0.5, -0.5])          # BosonicADT（对角耦合）
pt  = BosonicIF(bc, X)                    # BosonicPT（非对角耦合）
```

### 1.1 装配流水线与腿基约定（示例 ``N_\tau=2,\ N_t=2``）

两种 HybridizationMatrix 与 correlation function 携带**同一信息**：后者存分支
生成向量（``\eta_{jk}/\eta_{kj},\ \xi_{jk}/\xi_{kj},\ \zeta``），前者按
[`QuAPI.Dense.contour_ordering`](@ref) 链序把它们铺成**单个**矩阵 ``\Delta``，并把
`basis` 显式列出来。稠密 IF 的折叠/组装因此**只消费 `basis` + `Δ`**（玻色
PT 另加耦合算符 `hyb`——它是系统侧输入，不属于 bath 数据）：

```text
费米:  CorrelationFunction ─铺腿 _if_coefficients→ C ─顶点移位 _vertex_shifted_IF→ D
       (basis, D) ＝ GrassmannHybridizationMatrix
       GrassmannIF(h) ＝ _fold_grassmann(h.Δ)            # 有限次单项式乘法，只读 basis + Δ
玻色:  CorrelationFunction ─按 contour_ordering 组装→ (basis, Δ) ＝ BosonicHybridizationMatrix
       split=false: Δ = M1+M2             → BosonicIF(h, λ)    # ADT，对角嵌入
       split=true:  Δ = SplitDelta(M1,M2) → BosonicIF(h, hyb)  # PT，显式 exp(H)
```

即本页第 2 节各公式的 ``C`` 指**移位前**矩阵（`_if_coefficients` 的输出，
顶点在 cell 边界的滞后约定），``D`` 指移位后矩阵——**``D`` 就是
`GrassmannHybridizationMatrix.Δ`**。以混合轮廓 ``N_\tau=2,\ N_t=2`` 为例
（8 节点、单 flavor 费米 16 腿；下图约定已用 `DiscreteSpectrum` 实际构造
逐位核对）：

```@raw html
<div style="max-width:760px;margin:12px 0">
<svg viewBox="0 0 720 330" width="100%" height="auto" role="img" aria-label="mixed contour basis convention" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <marker id="mdoc1" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="8" markerHeight="8" markerUnits="userSpaceOnUse" orient="auto">
      <path d="M1 1 L7 4 L1 7 Z" fill="#59666E"/>
    </marker>
  </defs>
  <g font-family="system-ui, 'Segoe UI', 'PingFang SC', 'Helvetica Neue', Arial, sans-serif">
  <text x="20" y="26" font-size="16" font-weight="500" fill="#1F2328">混合轮廓 Nτ=2, Nt=2：contour_ordering 链序与两种 HybridizationMatrix 的 basis</text>
  <path d="M50 118 H180" stroke="#59666E" stroke-width="1.5" fill="none" stroke-linecap="round" marker-end="url(#mdoc1)"/>
  <circle cx="85" cy="118" r="3.5" fill="#59666E"/><circle cx="150" cy="118" r="3.5" fill="#59666E"/>
  <text x="85" y="106" font-size="11" fill="#59666E" text-anchor="middle">τ₁</text>
  <text x="150" y="106" font-size="11" fill="#59666E" text-anchor="middle">τ₂</text>
  <text x="115" y="140" font-size="11" fill="#59666E" text-anchor="middle">虚时 τ（2 cells）</text>
  <path d="M185 118 V76" stroke="#59666E" stroke-width="1.5" fill="none" stroke-linecap="round"/>
  <path d="M185 76 H340" stroke="#59666E" stroke-width="1.5" fill="none" stroke-linecap="round" marker-end="url(#mdoc1)"/>
  <circle cx="215" cy="76" r="3.5" fill="#59666E"/><circle cx="275" cy="76" r="3.5" fill="#59666E"/><circle cx="335" cy="76" r="3.5" fill="#59666E"/>
  <text x="215" y="64" font-size="11" fill="#59666E" text-anchor="middle">+₁</text>
  <text x="275" y="64" font-size="11" fill="#59666E" text-anchor="middle">+₂</text>
  <text x="335" y="64" font-size="11" fill="#59666E" text-anchor="middle">+₃</text>
  <text x="285" y="46" font-size="11" fill="#59666E" text-anchor="middle">+ 支（t₁ → t₃，Nt=2 cells）</text>
  <path d="M345 76 V158" stroke="#59666E" stroke-width="1.5" fill="none" stroke-linecap="round"/>
  <path d="M345 158 H190" stroke="#59666E" stroke-width="1.5" fill="none" stroke-linecap="round" marker-end="url(#mdoc1)"/>
  <circle cx="215" cy="158" r="3.5" fill="#59666E"/><circle cx="275" cy="158" r="3.5" fill="#59666E"/><circle cx="335" cy="158" r="3.5" fill="#59666E"/>
  <text x="215" y="176" font-size="11" fill="#59666E" text-anchor="middle">−₁</text>
  <text x="275" y="176" font-size="11" fill="#59666E" text-anchor="middle">−₂</text>
  <text x="335" y="176" font-size="11" fill="#59666E" text-anchor="middle">−₃</text>
  <text x="285" y="194" font-size="11" fill="#59666E" text-anchor="middle">− 支（t₃ → t₁）</text>
  <path d="M185 158 V118" stroke="#59666E" stroke-width="1" stroke-dasharray="4 4" fill="none"/>
  <g>
    <rect x="20" y="210" width="78" height="70" rx="8" fill="#F6F8FA" stroke="#D0D7DE"/>
    <text x="59" y="227" font-size="11" fill="#59666E" text-anchor="middle">p = 1</text>
    <text x="59" y="245" font-size="13" font-weight="600" fill="#1F2328" text-anchor="middle">(τ,1)</text>
    <text x="59" y="261" font-size="11" fill="#59666E" text-anchor="middle">G腿 1:a 2:ā</text>
    <text x="59" y="276" font-size="11" fill="#59666E" text-anchor="middle">B节点 1</text>
    <rect x="106" y="210" width="78" height="70" rx="8" fill="#F6F8FA" stroke="#D0D7DE"/>
    <text x="145" y="227" font-size="11" fill="#59666E" text-anchor="middle">p = 2</text>
    <text x="145" y="245" font-size="13" font-weight="600" fill="#1F2328" text-anchor="middle">(τ,2)</text>
    <text x="145" y="261" font-size="11" fill="#59666E" text-anchor="middle">G腿 3:a 4:ā</text>
    <text x="145" y="276" font-size="11" fill="#59666E" text-anchor="middle">B节点 2</text>
    <rect x="192" y="210" width="78" height="70" rx="8" fill="#F6F8FA" stroke="#D0D7DE"/>
    <text x="231" y="227" font-size="11" fill="#59666E" text-anchor="middle">p = 3</text>
    <text x="231" y="245" font-size="13" font-weight="600" fill="#1F2328" text-anchor="middle">(+,1)</text>
    <text x="231" y="261" font-size="11" fill="#59666E" text-anchor="middle">G腿 5:a 6:ā</text>
    <text x="231" y="276" font-size="11" fill="#59666E" text-anchor="middle">B节点 3</text>
    <rect x="278" y="210" width="78" height="70" rx="8" fill="#F6F8FA" stroke="#D0D7DE"/>
    <text x="317" y="227" font-size="11" fill="#59666E" text-anchor="middle">p = 4</text>
    <text x="317" y="245" font-size="13" font-weight="600" fill="#1F2328" text-anchor="middle">(−,1)</text>
    <text x="317" y="261" font-size="11" fill="#59666E" text-anchor="middle">G腿 7:a 8:ā</text>
    <text x="317" y="276" font-size="11" fill="#59666E" text-anchor="middle">B节点 4</text>
    <rect x="364" y="210" width="78" height="70" rx="8" fill="#F6F8FA" stroke="#D0D7DE"/>
    <text x="403" y="227" font-size="11" fill="#59666E" text-anchor="middle">p = 5</text>
    <text x="403" y="245" font-size="13" font-weight="600" fill="#1F2328" text-anchor="middle">(+,2)</text>
    <text x="403" y="261" font-size="11" fill="#59666E" text-anchor="middle">G腿 9:a 10:ā</text>
    <text x="403" y="276" font-size="11" fill="#59666E" text-anchor="middle">B节点 5</text>
    <rect x="450" y="210" width="78" height="70" rx="8" fill="#F6F8FA" stroke="#D0D7DE"/>
    <text x="489" y="227" font-size="11" fill="#59666E" text-anchor="middle">p = 6</text>
    <text x="489" y="245" font-size="13" font-weight="600" fill="#1F2328" text-anchor="middle">(−,2)</text>
    <text x="489" y="261" font-size="11" fill="#59666E" text-anchor="middle">G腿 11:a 12:ā</text>
    <text x="489" y="276" font-size="11" fill="#59666E" text-anchor="middle">B节点 6</text>
    <rect x="536" y="210" width="78" height="70" rx="8" fill="#F6F8FA" stroke="#D0D7DE"/>
    <text x="575" y="227" font-size="11" fill="#59666E" text-anchor="middle">p = 7</text>
    <text x="575" y="245" font-size="13" font-weight="600" fill="#1F2328" text-anchor="middle">(+,3)</text>
    <text x="575" y="261" font-size="11" fill="#59666E" text-anchor="middle">G腿 13:a 14:ā</text>
    <text x="575" y="276" font-size="11" fill="#59666E" text-anchor="middle">B节点 7</text>
    <rect x="622" y="210" width="78" height="70" rx="8" fill="#F6F8FA" stroke="#D0D7DE"/>
    <text x="661" y="227" font-size="11" fill="#59666E" text-anchor="middle">p = 8</text>
    <text x="661" y="245" font-size="13" font-weight="600" fill="#1F2328" text-anchor="middle">(−,3)</text>
    <text x="661" y="261" font-size="11" fill="#59666E" text-anchor="middle">G腿 15:a 16:ā</text>
    <text x="661" y="276" font-size="11" fill="#59666E" text-anchor="middle">B节点 8</text>
  </g>
  <text x="20" y="302" font-size="12" fill="#59666E">Grassmann basis：每节点 2 条腿，band 内 (a, ā) 相邻 —— GrassmannContourIndex(j, branch, band=1, conj)；a 占奇位、ā 占偶位</text>
  <text x="20" y="320" font-size="12" fill="#59666E">Bosonic basis：每节点 1 个路径指标 —— BosonicContourIndex(j, branch)。行/列序都由 contour_ordering 给出（τ 升序 → 实时按时刻，每时刻 + 先 − 后）</text>
  </g>
</svg>
</div>
```

```@raw html
<div style="max-width:760px;margin:12px 0">
<svg viewBox="0 0 720 400" width="100%" height="auto" role="img" aria-label="Grassmann hybridization matrix D structure" xmlns="http://www.w3.org/2000/svg">
  <g font-family="system-ui, 'Segoe UI', 'PingFang SC', 'Helvetica Neue', Arial, sans-serif">
  <text x="20" y="26" font-size="16" font-weight="500" fill="#1F2328">GrassmannHybridizationMatrix.Δ ＝ 文档的 D（16 腿 × 16 腿，Nτ=2, Nt=2）</text>
  <rect x="84" y="60" width="272" height="272" rx="2" fill="#F6F8FA" stroke="#D0D7DE"/>
  <path d="M101 60V332M118 60V332M135 60V332M152 60V332M169 60V332M186 60V332M203 60V332M220 60V332M237 60V332M254 60V332M271 60V332M288 60V332M305 60V332M322 60V332M339 60V332" stroke="#E1E6EB" stroke-width="0.5" fill="none"/>
  <path d="M84 77H356M84 94H356M84 111H356M84 128H356M84 145H356M84 162H356M84 179H356M84 196H356M84 213H356M84 230H356M84 247H356M84 264H356M84 281H356M84 298H356M84 315H356" stroke="#E1E6EB" stroke-width="0.5" fill="none"/>
  <g fill="#D3D4DA" opacity="0.55">
    <rect x="84" y="77" width="17" height="17"/><rect x="118" y="77" width="17" height="17"/><rect x="152" y="77" width="17" height="17"/><rect x="186" y="77" width="17" height="17"/><rect x="220" y="77" width="17" height="17"/><rect x="254" y="77" width="17" height="17"/><rect x="288" y="77" width="17" height="17"/><rect x="322" y="77" width="17" height="17"/>
    <rect x="84" y="315" width="17" height="17"/><rect x="118" y="315" width="17" height="17"/><rect x="152" y="315" width="17" height="17"/><rect x="186" y="315" width="17" height="17"/><rect x="220" y="315" width="17" height="17"/><rect x="254" y="315" width="17" height="17"/><rect x="288" y="315" width="17" height="17"/><rect x="322" y="315" width="17" height="17"/>
    <rect x="186" y="111" width="17" height="17"/><rect x="288" y="111" width="17" height="17"/>
    <rect x="186" y="145" width="17" height="17"/><rect x="288" y="145" width="17" height="17"/>
    <rect x="186" y="179" width="17" height="17"/><rect x="288" y="179" width="17" height="17"/>
    <rect x="186" y="213" width="17" height="17"/><rect x="288" y="213" width="17" height="17"/>
    <rect x="186" y="247" width="17" height="17"/><rect x="288" y="247" width="17" height="17"/>
    <rect x="186" y="281" width="17" height="17"/><rect x="288" y="281" width="17" height="17"/>
  </g>
  <g fill="#85FABF">
    <rect x="84" y="111" width="17" height="17"/><rect x="118" y="111" width="17" height="17"/><rect x="254" y="111" width="17" height="17"/><rect x="322" y="111" width="17" height="17"/>
    <rect x="84" y="145" width="17" height="17"/><rect x="118" y="145" width="17" height="17"/><rect x="254" y="145" width="17" height="17"/><rect x="322" y="145" width="17" height="17"/>
    <rect x="84" y="179" width="17" height="17"/><rect x="118" y="179" width="17" height="17"/><rect x="254" y="179" width="17" height="17"/><rect x="322" y="179" width="17" height="17"/>
    <rect x="84" y="247" width="17" height="17"/><rect x="118" y="247" width="17" height="17"/><rect x="254" y="247" width="17" height="17"/><rect x="322" y="247" width="17" height="17"/>
    <rect x="152" y="213" width="17" height="17"/><rect x="220" y="213" width="17" height="17"/>
    <rect x="152" y="281" width="17" height="17"/><rect x="220" y="281" width="17" height="17"/>
  </g>
  <g fill="#F87454">
    <rect x="152" y="111" width="17" height="17"/><rect x="220" y="111" width="17" height="17"/>
    <rect x="152" y="145" width="17" height="17"/><rect x="220" y="145" width="17" height="17"/>
    <rect x="152" y="179" width="17" height="17"/><rect x="220" y="179" width="17" height="17"/>
    <rect x="152" y="247" width="17" height="17"/><rect x="220" y="247" width="17" height="17"/>
    <rect x="84" y="213" width="17" height="17"/><rect x="118" y="213" width="17" height="17"/><rect x="254" y="213" width="17" height="17"/><rect x="322" y="213" width="17" height="17"/>
    <rect x="84" y="281" width="17" height="17"/><rect x="118" y="281" width="17" height="17"/><rect x="254" y="281" width="17" height="17"/><rect x="322" y="281" width="17" height="17"/>
  </g>
  <g font-size="11" fill="#59666E" text-anchor="middle">
    <text x="70" y="72">a</text><text x="70" y="89">ā</text><text x="70" y="106">a</text><text x="70" y="123">ā</text><text x="70" y="140">a</text><text x="70" y="157">ā</text><text x="70" y="174">a</text><text x="70" y="191">ā</text><text x="70" y="208">a</text><text x="70" y="225">ā</text><text x="70" y="242">a</text><text x="70" y="259">ā</text><text x="70" y="276">a</text><text x="70" y="293">ā</text><text x="70" y="310">a</text><text x="70" y="327">ā</text>
    <text x="50" y="81" text-anchor="end">τ₁</text><text x="50" y="115" text-anchor="end">τ₂</text><text x="50" y="149" text-anchor="end">+₁</text><text x="50" y="183" text-anchor="end">−₁</text><text x="50" y="217" text-anchor="end">+₂</text><text x="50" y="251" text-anchor="end">−₂</text><text x="50" y="285" text-anchor="end">+₃</text><text x="50" y="319" text-anchor="end">−₃</text>
    <text x="92.5" y="348">a</text><text x="109.5" y="348">ā</text><text x="126.5" y="348">a</text><text x="143.5" y="348">ā</text><text x="160.5" y="348">a</text><text x="177.5" y="348">ā</text><text x="194.5" y="348">a</text><text x="211.5" y="348">ā</text><text x="228.5" y="348">a</text><text x="245.5" y="348">ā</text><text x="262.5" y="348">a</text><text x="279.5" y="348">ā</text><text x="296.5" y="348">a</text><text x="313.5" y="348">ā</text><text x="330.5" y="348">a</text><text x="347.5" y="348">ā</text>
    <text x="101" y="364">τ₁</text><text x="135" y="364">τ₂</text><text x="169" y="364">+₁</text><text x="203" y="364">−₁</text><text x="237" y="364">+₂</text><text x="271" y="364">−₂</text><text x="305" y="364">+₃</text><text x="339" y="364">−₃</text>
  </g>
  <text x="16" y="196" font-size="11" fill="#59666E" text-anchor="middle" transform="rotate(-90 16 196)">行：指标 k 后紧跟 ā_k（k 奇行位、ā_k 偶行位）</text>
  <text x="220" y="380" font-size="11" fill="#59666E" text-anchor="middle">列：指标 k 后紧跟 ā_k——行列均按节点 k 分组（同 Bosonic 基），每节点占 2 行/2 列</text>
  <text x="220" y="394" font-size="11" fill="#59666E" text-anchor="middle">非零元只在（行 ā_k × 列 a_j）</text>
  <text x="396" y="60" font-size="12" font-weight="600" fill="#59666E">装配流水线</text>
  <text x="396" y="80" font-size="12" fill="#59666E">corr ─生成向量铺腿→ C（顶点在 cell 边界）</text>
  <text x="396" y="98" font-size="12" fill="#59666E">C ─顶点移位→ D ＝ h.Δ（本图的矩阵）</text>
  <text x="396" y="116" font-size="12" fill="#59666E">GrassmannIF(h)：折叠只读 basis + Δ</text>
  <text x="396" y="146" font-size="13" font-weight="500" fill="#1F2328" font-family="ui-monospace, Menlo, Consolas, monospace">D[2p, 2q−1] = s · C[ā_ρ(p), a_κ(q)]</text>
  <text x="396" y="162" font-size="11" fill="#59666E">ā_p 在腿 2p，a_q 在腿 2q−1；普通 bath 只有 (ā,a) 型项</text>
  <text x="396" y="190" font-size="12" font-weight="600" fill="#59666E">源映射（· = 无源 → 全零）</text>
  <text x="396" y="208" font-size="11" fill="#59666E">节点   τ₁  τ₂  +₁  −₁  +₂  −₂  +₃  −₃</text>
  <text x="396" y="224" font-size="11" fill="#59666E">行源 ρ  ·   τ₁  τ₂  −₁  +₁  −₂  +₂  ·</text>
  <text x="396" y="240" font-size="11" fill="#59666E">列源 κ  τ₁  τ₂  +₁  ·   +₂  −₁  ·   −₂</text>
  <text x="396" y="268" font-size="12" font-weight="600" fill="#59666E">符号 s（代码 _mx_maps / _vs_mixed_acc）</text>
  <text x="396" y="286" font-size="12" fill="#59666E">s = −1 ⟺ 行源、列源恰一个在 + 支</text>
  <text x="396" y="302" font-size="12" fill="#59666E">行源在 + 支：p = +₂,+₃（腿 10,14）</text>
  <text x="396" y="318" font-size="12" fill="#59666E">列源在 + 支：q = +₁,+₂（腿 5,9）</text>
  <rect x="396" y="332" width="12" height="12" rx="3" fill="#85FABF"/><text x="414" y="342" font-size="12" fill="#59666E">s = +1（24 元）</text>
  <rect x="396" y="352" width="12" height="12" rx="3" fill="#F87454"/><text x="414" y="362" font-size="12" fill="#59666E">s = −1（12 元）</text>
  <rect x="500" y="332" width="12" height="12" rx="3" fill="#D3D4DA" opacity="0.55"/><text x="518" y="342" font-size="12" fill="#59666E">无源零元：行 τ₁/−₃，</text>
  <text x="518" y="362" font-size="12" fill="#59666E">列 −₁/+₃</text>
  </g>
</svg>
</div>
```

普通 bath 的费米 ``D`` 的行/列同玻色基一样**按节点 ``k`` 分组**，每节点占
2 行/2 列——指标 ``k``（``a_k``，奇行位/奇列位）后紧跟 ``\bar k``（``\bar a_k``，
偶位）；非零元只在 ``(行\ \bar a_k,\ 列\ a_j)``：
``D[2k,\,2j-1] \;=\; s\cdot C[\bar a_{\rho(k)},\, a_{\kappa(j)}]``。无源
（``\rho`` 或 ``\kappa`` 为空）的行/列整行/整列为零；符号 ``s=-1`` 当且仅当
行源、列源**恰一个**在 ``+`` 支。BCS 的反常块 ``(ā↓,a↑)/(ā↑,ā↓)`` 逐块做
同构移位。装配规则的统一数学描述见 1.2 节。

费米侧同样的生成元来源分解（取 ``N_t=1`` 减少指标数；``N_\tau=2``，6 节点
12 腿）。注意顶点移位符号 ``s`` 已折入各表值：

```@raw html
<div style="max-width:760px;margin:12px 0">
<svg viewBox="0 0 860 410" width="100%" height="auto" role="img" aria-label="Grassmann hybridization matrix generator provenance (fermionic, Nt=1)" xmlns="http://www.w3.org/2000/svg">
  <g font-family="system-ui, 'Segoe UI', 'PingFang SC', 'Helvetica Neue', Arial, sans-serif">
  <text x="20" y="26" font-size="16" font-weight="500" fill="#1F2328">GrassmannHybridizationMatrix.Δ 的生成元来源（费米，Nτ=2, Nt=1）</text>
  <rect x="96" y="64" width="264" height="264" rx="2" fill="#F6F8FA" stroke="#D0D7DE"/>
  <path d="M118 64V328M140 64V328M162 64V328M184 64V328M206 64V328M228 64V328M250 64V328M272 64V328M294 64V328M316 64V328M338 64V328" stroke="#E1E6EB" stroke-width="0.5" fill="none"/>
  <path d="M96 86H360M96 108H360M96 130H360M96 152H360M96 174H360M96 196H360M96 218H360M96 240H360M96 262H360M96 284H360M96 306H360" stroke="#E1E6EB" stroke-width="0.5" fill="none"/>
  <g fill="#D3D4DA" opacity="0.55">
    <rect x="96" y="64" width="22" height="22"/><rect x="140" y="64" width="22" height="22"/><rect x="184" y="64" width="22" height="22"/><rect x="228" y="64" width="22" height="22"/><rect x="272" y="64" width="22" height="22"/><rect x="316" y="64" width="22" height="22"/>
    <rect x="96" y="306" width="22" height="22"/><rect x="140" y="306" width="22" height="22"/><rect x="184" y="306" width="22" height="22"/><rect x="228" y="306" width="22" height="22"/><rect x="272" y="306" width="22" height="22"/><rect x="316" y="306" width="22" height="22"/>
    <rect x="228" y="130" width="22" height="22"/><rect x="228" y="152" width="22" height="22"/><rect x="228" y="174" width="22" height="22"/><rect x="228" y="196" width="22" height="22"/>
    <rect x="272" y="130" width="22" height="22"/><rect x="272" y="152" width="22" height="22"/><rect x="272" y="174" width="22" height="22"/><rect x="272" y="196" width="22" height="22"/>
  </g>
  <g fill="#80C1FF">
    <rect x="96" y="130" width="22" height="22"/><rect x="140" y="130" width="22" height="22"/>
    <rect x="96" y="174" width="22" height="22"/><rect x="140" y="174" width="22" height="22"/>
  </g>
  <g fill="#02743B">
    <rect x="184" y="130" width="22" height="22"/><rect x="316" y="130" width="22" height="22"/>
    <rect x="184" y="174" width="22" height="22"/><rect x="316" y="174" width="22" height="22"/>
    <rect x="96" y="218" width="22" height="22"/><rect x="140" y="218" width="22" height="22"/>
    <rect x="96" y="262" width="22" height="22"/><rect x="140" y="262" width="22" height="22"/>
  </g>
  <g fill="#85FABF">
    <rect x="184" y="218" width="22" height="22"/><rect x="316" y="218" width="22" height="22"/>
    <rect x="184" y="262" width="22" height="22"/><rect x="316" y="262" width="22" height="22"/>
  </g>
  <g font-size="11" fill="#59666E" text-anchor="middle">
    <text x="84" y="80">a</text><text x="84" y="102">ā</text><text x="84" y="124">a</text><text x="84" y="146">ā</text><text x="84" y="168">a</text><text x="84" y="190">ā</text><text x="84" y="212">a</text><text x="84" y="234">ā</text><text x="84" y="256">a</text><text x="84" y="278">ā</text><text x="84" y="300">a</text><text x="84" y="322">ā</text>
    <text x="66" y="80" text-anchor="end">τ₁</text><text x="66" y="124" text-anchor="end">τ₂</text><text x="66" y="168" text-anchor="end">+₁</text><text x="66" y="212" text-anchor="end">−₁</text><text x="66" y="256" text-anchor="end">+₂</text><text x="66" y="300" text-anchor="end">−₂</text>
    <text x="107" y="342">a</text><text x="129" y="342">ā</text><text x="151" y="342">a</text><text x="173" y="342">ā</text><text x="195" y="342">a</text><text x="217" y="342">ā</text><text x="239" y="342">a</text><text x="261" y="342">ā</text><text x="283" y="342">a</text><text x="305" y="342">ā</text><text x="327" y="342">a</text><text x="349" y="342">ā</text>
    <text x="107" y="358">τ₁</text><text x="151" y="358">τ₂</text><text x="195" y="358">+₁</text><text x="239" y="358">−₁</text><text x="283" y="358">+₂</text><text x="327" y="358">−₂</text>
  </g>
  <text x="20" y="382" font-size="12" fill="#59666E">行/列按节点 k 分组（k 后紧跟 ā_k）。零行 τ₁、−₂（行源 ρ=∅）；零列 +₂、−₁（列源 κ=∅）。颜色 = 来源块：蓝 G^ττ、深绿 ζ 跨块、绿 G^σσ</text>
  <text x="470" y="70" font-size="12" font-weight="600" fill="#59666E">非零 4×4 子块（移位符号 s 已折入）</text>
  <g font-size="11" fill="#59666E" text-anchor="middle">
    <text x="566" y="104">a(τ₁)</text><text x="662" y="104">a(τ₂)</text><text x="758" y="104">a(+₁)</text><text x="818" y="104">a(−₂)</text>
  </g>
  <g font-size="11" fill="#59666E" text-anchor="end">
    <text x="462" y="138">ā(τ₂)</text><text x="462" y="178">ā(+₁)</text><text x="462" y="218">ā(−₁)</text><text x="462" y="258">ā(+₂)</text>
  </g>
  <g fill="#F6F8FA" stroke="#D0D7DE">
    <rect x="470" y="112" width="192" height="40"/><rect x="666" y="112" width="192" height="40"/>
    <rect x="470" y="152" width="192" height="40"/><rect x="666" y="152" width="192" height="40"/>
    <rect x="470" y="192" width="192" height="40"/><rect x="666" y="192" width="192" height="40"/>
    <rect x="470" y="232" width="192" height="40"/><rect x="666" y="232" width="192" height="40"/>
  </g>
  <g font-size="11" fill="#1F2328" text-anchor="middle">
    <text x="566" y="136">−ξ₊₁−ξ₋₁</text><text x="762" y="136">−ξ₋₂</text>
    <text x="566" y="176">−ξ₊₂</text><text x="762" y="176">−ξ₊₁−ξ₋₁</text>
  </g>
  <g font-size="11" fill="#02743B" text-anchor="middle">
    <text x="518" y="136">iζ₊[1,1]</text><text x="714" y="136">iζ₊[1,1]</text>
    <text x="518" y="176">iζ₊[2,1]</text><text x="714" y="176">iζ₊[2,1]</text>
    <text x="566" y="216">iζ₋[1,1]</text><text x="762" y="216">iζ₋[1,2]</text>
    <text x="566" y="256">iζ₋[1,1]</text><text x="762" y="256">iζ₋[1,2]</text>
  </g>
  <g font-size="11" fill="#3D7A34" text-anchor="middle">
    <text x="518" y="216">η₊₁+η₊₁*</text><text x="714" y="216">η₋₁*+η₊₁*</text>
    <text x="518" y="256">η₊₁+η₋₁</text><text x="714" y="256">η₋₁*+η₋₁</text>
  </g>
  <text x="470" y="300" font-size="12" fill="#59666E">表值 = s·C[ā_ρ(k), a_κ(j)]（顶点移位符号已折入）；</text>
  <text x="470" y="320" font-size="12" fill="#59666E">等时和 = M₁（T₋/S₋ 半胞）+ M₂（T₊/S₊ 半胞），</text>
  <text x="470" y="340" font-size="12" fill="#59666E">split=true 时分给 X_νY_ν / Y_νX_ν（§4）</text>
  </g>
</svg>
</div>
```

```@raw html
<div style="max-width:760px;margin:12px 0">
<svg viewBox="0 0 860 600" width="100%" height="auto" role="img" aria-label="Bosonic hybridization matrix block provenance" xmlns="http://www.w3.org/2000/svg">
  <g font-family="system-ui, 'Segoe UI', 'PingFang SC', 'Helvetica Neue', Arial, sans-serif">
  <text x="20" y="26" font-size="16" font-weight="500" fill="#1F2328">BosonicHybridizationMatrix.Δ（8 节点 × 8 节点，Nτ=2, Nt=2）—— 逐元生成元来源</text>
  <text x="20" y="46" font-size="11" fill="#59666E">记号：η₊ₙ ≡ η_jk[n]（前向/行靠后核）、η₋ₙ ≡ η_kj[n]（反向核），ξ 同理；ζ₊ₖᵢ ≡ ζ_jk[k,i]（行 τ 列 t）、ζ₋ᵢₖ ≡ ζ_kj[i,k]；* 为复共轭</text>
  <rect x="104" y="64" width="496" height="496" rx="2" fill="#F6F8FA" stroke="#D0D7DE"/>
  <path d="M166 64V560M228 64V560M290 64V560M352 64V560M414 64V560M476 64V560M538 64V560" stroke="#E1E6EB" stroke-width="0.5" fill="none"/>
  <path d="M104 126H600M104 188H600M104 250H600M104 312H600M104 374H600M104 436H600M104 498H600" stroke="#E1E6EB" stroke-width="0.5" fill="none"/>
  <g fill="#80C1FF">
    <rect x="104" y="64" width="62" height="62"/><rect x="166" y="64" width="62" height="62"/>
    <rect x="104" y="126" width="62" height="62"/><rect x="166" y="126" width="62" height="62"/>
  </g>
  <g fill="#02743B">
    <rect x="228" y="64" width="62" height="62"/><rect x="290" y="64" width="62" height="62"/><rect x="352" y="64" width="62" height="62"/><rect x="414" y="64" width="62" height="62"/><rect x="476" y="64" width="62" height="62"/><rect x="538" y="64" width="62" height="62"/>
    <rect x="228" y="126" width="62" height="62"/><rect x="290" y="126" width="62" height="62"/><rect x="352" y="126" width="62" height="62"/><rect x="414" y="126" width="62" height="62"/><rect x="476" y="126" width="62" height="62"/><rect x="538" y="126" width="62" height="62"/>
    <rect x="104" y="188" width="62" height="62"/><rect x="166" y="188" width="62" height="62"/><rect x="104" y="250" width="62" height="62"/><rect x="166" y="250" width="62" height="62"/><rect x="104" y="312" width="62" height="62"/><rect x="166" y="312" width="62" height="62"/><rect x="104" y="374" width="62" height="62"/><rect x="166" y="374" width="62" height="62"/><rect x="104" y="436" width="62" height="62"/><rect x="166" y="436" width="62" height="62"/><rect x="104" y="498" width="62" height="62"/><rect x="166" y="498" width="62" height="62"/>
  </g>
  <g fill="#85FABF">
    <rect x="228" y="188" width="62" height="62"/><rect x="352" y="188" width="62" height="62"/><rect x="476" y="188" width="62" height="62"/>
    <rect x="290" y="250" width="62" height="62"/><rect x="414" y="250" width="62" height="62"/><rect x="538" y="250" width="62" height="62"/>
    <rect x="228" y="312" width="62" height="62"/><rect x="352" y="312" width="62" height="62"/><rect x="476" y="312" width="62" height="62"/>
    <rect x="290" y="374" width="62" height="62"/><rect x="414" y="374" width="62" height="62"/><rect x="538" y="374" width="62" height="62"/>
    <rect x="228" y="436" width="62" height="62"/><rect x="352" y="436" width="62" height="62"/><rect x="476" y="436" width="62" height="62"/>
    <rect x="290" y="498" width="62" height="62"/><rect x="414" y="498" width="62" height="62"/><rect x="538" y="498" width="62" height="62"/>
  </g>
  <g fill="#0FDC78">
    <rect x="290" y="188" width="62" height="62"/><rect x="414" y="188" width="62" height="62"/><rect x="538" y="188" width="62" height="62"/>
    <rect x="228" y="250" width="62" height="62"/><rect x="352" y="250" width="62" height="62"/><rect x="476" y="250" width="62" height="62"/>
    <rect x="290" y="312" width="62" height="62"/><rect x="414" y="312" width="62" height="62"/><rect x="538" y="312" width="62" height="62"/>
    <rect x="228" y="374" width="62" height="62"/><rect x="352" y="374" width="62" height="62"/><rect x="476" y="374" width="62" height="62"/>
    <rect x="290" y="436" width="62" height="62"/><rect x="414" y="436" width="62" height="62"/><rect x="538" y="436" width="62" height="62"/>
    <rect x="228" y="498" width="62" height="62"/><rect x="352" y="498" width="62" height="62"/><rect x="476" y="498" width="62" height="62"/>
  </g>
  <g font-size="11" fill="#1F2328" text-anchor="middle">
    <text x="135" y="99">−ξ₊₁−ξ₋₁</text><text x="197" y="99">−ξ₋₂</text><text x="383" y="99">η₊₁+η₋₁</text><text x="445" y="99">−η₋₁*−η₋₁</text><text x="507" y="99">η₋₂</text><text x="569" y="99">−η₋₃</text>
    <text x="135" y="161">−ξ₊₂</text><text x="197" y="161">−ξ₊₁−ξ₋₁</text><text x="383" y="161">η₊₂</text><text x="445" y="161">−η₋₂*</text><text x="507" y="161">η₋₂</text><text x="569" y="161">−η₋₃</text>
    <text x="259" y="223">η₊₁+η₋₁</text><text x="321" y="223">−η₋₁*−η₋₁</text><text x="383" y="223">η₋₂</text><text x="445" y="223">−η₋₂</text><text x="507" y="223">η₋₃</text><text x="569" y="223">−η₋₃</text>
    <text x="259" y="285">−η₊₁−η₊₁*</text><text x="321" y="285">η₋₁*+η₊₁*</text><text x="383" y="285">−η₊₂*</text><text x="445" y="285">η₊₂*</text><text x="507" y="285">−η₊₃*</text><text x="569" y="285">η₊₃*</text>
    <text x="259" y="347">η₊₂</text><text x="321" y="347">−η₋₂*</text><text x="383" y="347">η₊₁+η₋₁</text><text x="445" y="347">−η₋₁*−η₋₁</text><text x="507" y="347">η₋₂</text><text x="569" y="347">−η₋₂</text>
    <text x="259" y="409">−η₊₂</text><text x="321" y="409">η₋₂*</text><text x="383" y="409">−η₊₁−η₊₁*</text><text x="445" y="409">η₋₁*+η₊₁*</text><text x="507" y="409">−η₊₂*</text><text x="569" y="409">η₊₂*</text>
    <text x="259" y="471">η₊₃</text><text x="321" y="471">−η₋₃*</text><text x="383" y="471">η₊₂</text><text x="445" y="471">−η₋₂*</text><text x="507" y="471">η₊₁+η₋₁</text><text x="569" y="471">−η₋₁*−η₋₁</text>
    <text x="259" y="533">−η₊₃</text><text x="321" y="533">η₋₃*</text><text x="383" y="533">−η₊₂</text><text x="445" y="533">η₋₂*</text><text x="507" y="533">−η₊₁−η₊₁*</text><text x="569" y="533">η₋₁*+η₊₁*</text>
  </g>
  <g font-size="11" fill="#FFFFFF" text-anchor="middle">
    <text x="259" y="99">−iζ₊₁₁</text><text x="321" y="99">iζ₊₁₁</text><text x="383" y="99">−iζ₊₁₂</text><text x="445" y="99">iζ₊₁₂</text><text x="507" y="99">−iζ₊₁₃</text><text x="569" y="99">iζ₊₁₃</text>
    <text x="259" y="161">−iζ₊₂₁</text><text x="321" y="161">iζ₊₂₁</text><text x="383" y="161">−iζ₊₂₂</text><text x="445" y="161">iζ₊₂₂</text><text x="507" y="161">−iζ₊₂₃</text><text x="569" y="161">iζ₊₂₃</text>
    <text x="135" y="223">−iζ₋₁₁</text><text x="197" y="223">−iζ₋₁₂</text>
    <text x="135" y="285">iζ₋₁₁</text><text x="197" y="285">iζ₋₁₂</text>
    <text x="135" y="347">−iζ₋₂₁</text><text x="197" y="347">−iζ₋₂₂</text>
    <text x="135" y="409">iζ₋₂₁</text><text x="197" y="409">iζ₋₂₂</text>
    <text x="135" y="471">−iζ₋₃₁</text><text x="197" y="471">−iζ₋₃₂</text>
    <text x="135" y="533">iζ₋₃₁</text><text x="197" y="533">iζ₋₃₂</text>
  </g>
  <g font-size="12" fill="#59666E" text-anchor="middle">
    <text x="135" y="56">τ₁</text><text x="197" y="56">τ₂</text><text x="259" y="56">+₁</text><text x="321" y="56">−₁</text><text x="383" y="56">+₂</text><text x="445" y="56">−₂</text><text x="507" y="56">+₃</text><text x="569" y="56">−₃</text>
  </g>
  <g font-size="12" fill="#59666E" text-anchor="end">
    <text x="96" y="99">τ₁</text><text x="96" y="161">τ₂</text><text x="96" y="223">+₁</text><text x="96" y="285">−₁</text><text x="96" y="347">+₂</text><text x="96" y="409">−₂</text><text x="96" y="471">+₃</text><text x="96" y="533">−₃</text>
  </g>
  <g font-size="12" fill="#59666E">
    <rect x="104" y="576" width="12" height="12" rx="3" fill="#80C1FF"/><text x="122" y="586">G^ττ ← ξ</text>
    <rect x="220" y="576" width="12" height="12" rx="3" fill="#85FABF"/><text x="238" y="586">G^σσ ← η</text>
    <rect x="336" y="576" width="12" height="12" rx="3" fill="#0FDC78"/><text x="354" y="586">G^±∓ ← η* 组合</text>
    <rect x="500" y="576" width="12" height="12" rx="3" fill="#02743B"/><text x="518" y="586">G^τσ / G^στ ← ∓iζ</text>
  </g>
  <text x="622" y="70" font-size="12" font-weight="600" fill="#59666E">对角元的 split 规则</text>
  <text x="622" y="94" font-size="12" fill="#59666E">η 型对角 = η₊₁ + η₋₁</text>
  <text x="622" y="114" font-size="12" fill="#59666E">η₊₁ = η_jk[1]（T₋/S₋ 半胞）→ M₁</text>
  <text x="622" y="134" font-size="12" fill="#59666E">η₋₁ = η_kj[1]（T₊/S₊ 半胞）→ M₂</text>
  <text x="622" y="164" font-size="12" fill="#59666E">split=false：两项相加 = Δ 对角（ADT）</text>
  <text x="622" y="184" font-size="12" fill="#59666E">split=true：分给 X_νY_ν 与 Y_νX_ν</text>
  <text x="622" y="204" font-size="12" fill="#59666E">两种算符积（§4 等时公式）</text>
  <text x="622" y="234" font-size="12" fill="#59666E">ζ / 普通块对角整体入 M₁；</text>
  <text x="622" y="254" font-size="12" fill="#59666E">跨支、τ×± 用求和核 (M₁+M₂)，</text>
  <text x="622" y="274" font-size="12" fill="#59666E">不受 split 影响</text>
  <text x="622" y="304" font-size="12" font-weight="600" fill="#59666E">读法示例</text>
  <text x="622" y="324" font-size="12" fill="#59666E">Δ[+₁,−₁] = −η₋₁*−η₋₁：同一时刻</text>
  <text x="622" y="344" font-size="12" fill="#59666E">跨支配对，读 (+,−) 块等时对角</text>
  <text x="622" y="374" font-size="12" font-weight="600" fill="#59666E">两个构造入口（只读 basis+Δ）</text>
  <text x="622" y="394" font-size="12" fill="#59666E">split=false → BosonicIF(h, λ)</text>
  <text x="622" y="414" font-size="12" fill="#59666E">split=true → BosonicIF(h, hyb)，</text>
  <text x="622" y="434" font-size="12" fill="#59666E">H 按分支规则组装，I = exp(H)</text>
  </g>
</svg>
</div>
```

玻色 ``Δ`` 的行列序同样是链序——九个分支块被 ``\pm`` 交错**打散**进矩阵，
不是按分支拼成 ``3\times3`` 连续块。混合装配用到的行/列源映射与符号在代码里
由 `Dense._mx_maps`（装配）与 `_vs_mixed_acc`（腿指标访问器）给出，二者逐位
一致。

### 1.2 装配的统一公式：不等时项与等时项

费米与玻色 Dense IF 的装配在「**核槽位 × 有序二算符积**」这一层是同一套数学。
记 contour 节点 ``\nu=1,\dots,n``（chain 序），每个节点有一对**顶点半腿**——
出腿 ``F_\nu`` 与入腿 ``G_\nu``：

- 费米：``F_\nu=\bar a_\nu``（行腿）、``G_\nu=a_\nu``（列腿），非零槽只有
  ``(\bar a_\nu, a_\mu)``（BCS 另有反常槽，见第 2 节）；
- 玻色 PT：``(F_\nu,G_\nu)=(X_\nu, Y_\nu)``（``+``、``\tau`` 支）或
  ``(Y_\nu, X_\nu)``（``-`` 支），``X_\nu, Y_\nu`` 为第 4 节的分支行腿矩阵
  （下标 ``\nu=(\sigma_\nu, j_\nu)`` 即节点的分支–时间标签）；
  玻色 ADT：``X_\nu=Y_\nu=s_\nu``（对角嵌入）；
- 核拆分 ``M_1+M_2=\Delta``：按**块内矩阵位置**拆——行时刻靠后（``i>j``）的元
  入 ``M_1``、行靠前（``i<j``）入 ``M_2``；各块的上/下三角具体读哪条生成向量
  见第 5 节（如 ``-{-}`` 块上三角读 ``\eta_{kj}^*``）。等时对角：``\eta`` 型块
  两半各入 ``M_1``/``M_2``（``\eta_{jk}[1]``/``\eta_{kj}[1]``）；``\zeta`` 等
  普通矩阵块对角全入 ``M_1``（求和核规则不受影响）。

#### 不等时（``\nu\neq\mu``）

同支节点对（``\nu`` 时刻靠后、``\mu`` 靠前；``+``/``-``/``\tau`` 三支同式）：

```math
H \;\mathrel{+}=\; M_1[\nu,\mu]\,F_\nu G_\mu \;+\; M_2[\mu,\nu]\,G_\nu F_\mu .
```

- 费米：``F_\nu G_\mu=\bar a_\nu a_\mu`` ⟹ 两项恰为 Toeplitz 的下、上三角两个
  槽位；``i<j`` 槽位的 Grassmann 反对易符号在 `_fold_grassmann` 折叠时计入；
- 玻色 ADT：``F=G`` ⟹ 两项合并为 ``(M_1{+}M_2)[\nu,\mu]\,s_\nu s_\mu``；
- 玻色 PT：``+``/``\tau`` 支 ``(F,G)=(X_1,X_2)``、``-`` 支交换，即第 4 节的
  「同支不同节点」规则。

跨支与 ``\tau\times\pm`` 用**求和核**（核留原槽位，配对算符标签按支互换）：

```math
(+\nu,-\mu):\quad
H \;\mathrel{+}=\; (M_1{+}M_2)[\nu,\mu]\,X_1^\nu X_1^\mu
            + (M_1{+}M_2)[\mu,\nu]\,X_2^\nu X_2^\mu ,
```

```math
(\tau\nu,\pm\mu):\quad
H \;\mathrel{+}=\; (M_1{+}M_2)[\nu,\mu]\,X_1^\nu G_\mu
            + (M_1{+}M_2)[\mu,\nu]\,X_2^\nu F_\mu .
```

费米没有标签交换问题：跨支槽位仍是 ``(ā_\nu,a_\mu)``，只是核值取对应跨支块
（``G^{\sigma\sigma'}``、``\zeta``）；混合轮廓再经顶点移位把这些槽位的核统一
改写为 ``D[\bar a_\nu,a_\mu] = s\,C[\bar a_{\rho(\nu)},\,a_{\kappa(\mu)}]``
（``\rho/\kappa`` 与 ``s`` 见第 2 节）。

#### 等时（``\nu=\mu``）

同一 cell 上的双重积分拆成两种算符次序的三角核 ``T_\mp``（虚时）/``S_\mp``
（实时），分别存 ``M_1[\nu,\nu]=\eta_{jk}[1]``、``M_2[\nu,\nu]=\eta_{kj}[1]``。
两种统计对「同一顶点上的两种次序」处理不同：

- **玻色 PT**：两算符挤在同一节点的 ``d\times d`` 空间，矩阵积次序显式保留：

```math
H \;\mathrel{+}=\;
\begin{cases}
M_1[\nu,\nu]\,X_2^\nu X_1^\nu + M_2[\nu,\nu]\,X_1^\nu X_2^\nu,
  & \sigma_\nu=+ ,\\[4pt]
M_1[\nu,\nu]\,X_1^\nu X_2^\nu + M_2[\nu,\nu]\,X_2^\nu X_1^\nu,
  & \sigma_\nu=-\ \text{或}\ \tau .
\end{cases}
```

  ``op_1=op_2``（厄米）时两行相同 ⟹ 合并为 ``\Delta[\nu,\nu]``（对角耦合下
  无次序信息可丢）。

- **玻色 ADT**：``X_1X_2=X_2X_1=s_\nu^2`` ⟹
  ``H \mathrel{+}= \Delta[\nu,\nu]\,s_\nu^2``，与不等时公式同形（即第 3 节
  路径权重的对角项）。

- **费米**：``\bar a_\nu\neq a_\nu`` 是不同 Grassmann 变量 ⟹ 两种次序核
  **相加**进同一槽位，无次序拆分。顶点移位前
  ``D[\bar a_\nu,a_\nu]=\eta_{jk}[1]+\eta_{kj}[1]``；移位后等时内容搬到
  **相邻节点的「动能位」**：``\tau`` 环
  ``D[\bar a_k,a_{k-1}]=c_1=\eta_{jk}[1]+\eta_{kj}[1]``（``k=1`` 环绕取负），
  对角位放 lag ``-1`` 核 ``c_0=\eta_{kj}[2]``；实时
  ``D[\bar a^+_j,a^+_{j-1}]=d_+``、``D[\bar a^-_j,a^-_{j+1}]=d_-``
  （``d_\pm=G^{\sigma\sigma}_{00}`` 即该支等时和）；混合 ``+`` 支首节点行源为
  ``\tau_{N_\tau}``，``\beta\to0`` 边界的等时和成为 τ→t 的「接头」项。

## 2. 费米：GrassmannADT

### 腿基与配对矩阵

腿的排列（`basis`）：节点按 contour 序（虚时/实时时间升序，实时每时刻
`+` 先 `−` 后，混合 τ 先 ± 后），每节点 band 内 `(a, ā)` 相邻。IF 是腿上的
Grassmann 高斯指数

```math
I[\xi_1,\dots,\xi_L] \;=\; \exp\Big(\sum_{p,q=1}^{L} C_{pq}\,\xi_p\,\xi_q\Big),
```

配对矩阵 ``C`` 的非零结构由 `_delta_terms` 决定（本式的配对矩阵即移位后的
``D`` ＝ `GrassmannHybridizationMatrix.Δ`，见 1.1 节的记号约定）：普通 bath 只有
``(\bar a, a)`` 块 ``C[\bar a_p, a_q] = \Delta_{pq}``（``\Delta`` 的行列序即
`contour_ordering`）；BCS 另有第二 flavor 的正常块与反常块
``(\bar a_\uparrow,\bar a_\downarrow)``、``(a_\downarrow,a_\uparrow)``。
由于 ``\xi^2=0``，指数在总阶数 ``\le L`` 后自动终止——**折叠（`_fold_grassmann`）
是有限次的单项式乘法，无任何近似**。

### 张量元、求值与边缘化

`data` 按腿 1 最快展开单项式系数，单腿占据态用位掩码索引：

```math
I \;=\; \sum_{m=0}^{2^L-1} \texttt{data}[m{+}1]\;\xi^m,
\qquad \xi^m = \prod_{k:\,m_k=1}\xi_k,
```

单项式系数带**反对易符号**（按 `basis` 升序重排计入；腿重复返回零）。
配对展开为

```math
[\xi_{i_1}\cdots\xi_{i_r}]\,I \;=\; \sum_{\text{配对方案 } M}
\operatorname{sgn}(M)\,\prod_{(\bar p,\,q)\in M} C[\bar p, q],
\qquad i_1<\cdots<i_r,
```

全腿单项式系数即 Berezin 积分（配分函数）：

```math
Z \;=\; \int\!\mathrm{d}\xi_L\cdots\mathrm{d}\xi_1\; I
\;=\; [\xi_1\cdots\xi_L]\,I \;=\; \texttt{data}[2^L],
```

即 `Zvalue`。[`QuAPI.Dense.margin`](@ref) 对未列出的腿占据态求和（Grassmann
边缘化）；两个 IF 的乘法 `A * B` 把共享腿按 Grassmann 融合、外腿并入
`basis`。

### 顶点移位后的 ``C``：τ 环

冻结核的生成向量按「顶点在 cell 边界」的滞后约定存储；稠密组装把顶点
移到节点——生成值重组为 **skew-circulant**（单 flavor，``N`` 节点）：

```math
c_0 = \eta_{kj}[2], \qquad
c_1 = \eta_{jk}[1] + \eta_{kj}[1] \;(\text{等时和}), \qquad
c_k = \eta_{jk}[k] \;(k\ge 2),
```

```math
C[\bar a_i, a_j] \;=\;
\begin{cases}
c_{\,i-j}, & i > j \\[2pt]
c_0, & i = j \\[2pt]
-\,c_{\,N+i-j}, & i < j
\end{cases}
```

（下三角为正向传播核；上三角取反周期边界。旧装配把 ``c_1`` 放在
``(i, i{-}1)`` 位——残差因此停在 ``O(V^{2})``，见第 5 节。）

BCS 把同一移位施加到四个 ``(row, col)`` 块上（正常
``(\bar a_\uparrow,a_\uparrow)``、``(\bar a_\downarrow,a_\downarrow)`` 都读
`uu`；反常 ``(\bar a_\uparrow,\bar a_\downarrow)`` 读 `ud`、
``(a_\downarrow,a_\uparrow)`` 读 `du`）。一般块移位公式（``r_i``/``s_j`` 为
块内腿位，``n`` 为节点数）：

```math
\gamma_1 = C[r_1,s_1], \quad
\gamma_k = C[r_k,s_1]\;(k\ge2), \quad
\gamma_0 = C[r_1,s_2],
```

```math
D[r_i,s_j] \;=\;
\begin{cases}
\gamma_0, & i=j \\
\gamma_{\,i-j}, & i>j \\
-\gamma_{\,n+i-j}, & i<j
\end{cases}
```

### 实时 Keldysh 块

节点 ``j_0 = 0,\dots,N``，腿位
``\mathrm{pos}(j_0,\sigma,\bar a) = 4j_0 + [\sigma=+]\cdot 0 + [\sigma=-]\cdot 2 + [\bar a]\cdot 1``。
记 ``d_{\pm} = C[\mathrm{pos}(0,\pm,\bar a),\,\mathrm{pos}(0,\pm,a)]``（Toeplitz
对角常值）。四个 ``(ā, a)`` 块的移位：

```math
\begin{aligned}
D[\bar a^{+}_{j_1}, a^{-}_{j_2}] &= C[\bar a^{+}_{j_1}, a^{-}_{j_2}],
&& \text{窗口：}\ \bar a^{+}_0\ \text{行与}\ a^{-}_0\ \text{列置零},\\[2pt]
D[\bar a^{-}_{j_1}, a^{+}_{j_2}] &= C[\bar a^{-}_{j_1}, a^{+}_{j_2}],
&& \text{窗口：}\ a^{+}_N\ \text{列与}\ \bar a^{-}_N\ \text{行置零},\\[2pt]
D[\bar a^{+}_{j_1}, a^{+}_{j_2}] &=
\begin{cases} d_{+}, & j_2 = j_1-1 \\ C[\bar a^{+}_{j_1},\, a^{+}_{j_2+1}], & \text{否则} \end{cases}
&& \text{（列移位 } a_j\!\to\! a_{j+1}\text{，动能位取本支对角）},\\[2pt]
D[\bar a^{-}_{j_1}, a^{-}_{j_2}] &=
\begin{cases} d_{-}, & j_2 = j_1+1 \\ C[\bar a^{-}_{j_1},\, a^{-}_{j_2-1}], & \text{否则} \end{cases}
&& \text{（镜像移位 } a_j\!\to\! a_{j-1}\text{）}.
\end{aligned}
```

（交叉块 ``(ā⁺,a⁻)``, ``(ā⁻,a⁺)`` 数值不变，仅把顶点窗口收紧到 cell
``1..N``；``(ā^+,a^+)``/``(ā^-,a^-)`` 的列移位把 a 算符顶点从 cell 边界搬到
节点，动能位 ``e^{-\delta t\,\epsilon_d}`` 由对角元给出。）

### 混合轮廓

``\tau`` 段与 ``\pm`` 段的顶点源映射（行源 ``\rho``、列源 ``\kappa``；无源 ⟹ 该
元为零；``N_\tau`` 为虚时节点数、``N_t`` 为实时 cell 数）。符号
``s(\nu,\mu) = -1`` 当且仅当行源、列源**恰一个**在 ``+`` 支（``+_1`` 行例外：
其源为 ``\tau_{N_\tau}``，行侧不算 + 支）：

```math
\rho:\quad
\tau_j \to \begin{cases}\tau_{j-1}, & j\ge 2\\ \varnothing, & j=1\end{cases}
\qquad
+_j \to \begin{cases}\tau_{N_\tau}, & j=1\\ +_{j-1}, & j\ge2\end{cases}
\qquad
-_j \to \begin{cases}-_j, & j\le N_t\\ \varnothing, & j=N_t{+}1\end{cases}
```

```math
\kappa:\quad
\tau_k \to \tau_k
\qquad
+_k \to \begin{cases}+_k, & k\le N_t\\ \varnothing, & k=N_t{+}1\end{cases}
\qquad
-_k \to \begin{cases}\varnothing, & k=1\\ -_{k-1}, & k\ge2\end{cases}
```

```math
D[\nu, \mu] \;=\; s(\nu,\mu)\; C[\rho(\nu),\, \kappa(\mu)].
```

（数值标定：混合移位对 `exact_fermionic_Δm` 严格参照做 Richardson 外推
至 ``O(V^{2})`` 逐位一致。）

## 3. 玻色对角：BosonicADT

对角耦合 ``S`` 可对角化（本征值 ``\lambda_1,\dots,\lambda_d``）时，IF 是路径
权重的**对角嵌入**：路径 ``\nu`` 处算符取本征值 ``s_\nu\in\{\lambda_p\}``，

```math
I[s_1,\dots,s_n] \;=\;
\exp\Big(\sum_{\nu,\mu=1}^{n} s_\nu\,\Delta_{\nu\mu}\,s_\mu\Big),
```

`data` 按节点 1 最快存全部 ``d^n`` 个路径权重；配分函数是全部路径之和

```math
Z \;=\; \sum_{s_1,\dots,s_n} I[s_1,\dots,s_n] \;=\; \texttt{sum(data)}
```

（`Zvalue`）。链上相邻 ADT 的乘法收缩共享节点。**这一形式对任意
`split=false` 核是精确的等式**（第 5 节的负对照：对严格 cell 参照残差
``\sim 10^{-15}`` 平线）。

## 4. 玻色非对角：BosonicPT

``hyb`` 不可对角化（或本征基外）时，IF 是作用在系统路径上的**算符**。
代码取 ``\mathrm{op}_1 = S``、``\mathrm{op}_2 = S^{\dagger}``（`op2 = adjoint(hyb)`）。
节点 ``\nu``（下标即其分支–时间标签 ``(\sigma_\nu, j_\nu)``）的行腿矩阵随分支取
转置 / 伴随 / 原矩阵（张量指标为
``(\mathrm{out},\mathrm{in})``；+ 支元素 ``\mathrm{op}[a,m]``、− 支为
``\mathrm{op}^{\dagger}`` 的元素 ``\overline{\mathrm{op}[b,m]}``、τ 支同 + 支
但腿序相反）：

```math
X_{\nu} = \begin{cases}\mathrm{op}_1^{\mathsf T} & (\sigma_\nu=+)\\ \mathrm{op}_1^{\dagger} & (\sigma_\nu=-)\\ \mathrm{op}_1 & (\sigma_\nu=\tau)\end{cases},
\qquad
Y_{\nu} = \begin{cases}\mathrm{op}_2^{\mathsf T} & (\sigma_\nu=+)\\ \mathrm{op}_2^{\dagger} & (\sigma_\nu=-)\\ \mathrm{op}_2 & (\sigma_\nu=\tau)\end{cases}.
```

``M_1 + M_2 = \Delta`` 为等时算符次序拆分（``\eta_{jk}``/``\eta_{kj}`` 侧，
第 9 节）。``H`` 按节点对 ``(\nu,\mu)`` 的分支关系组装（``j_\nu`` 为 contour
步数，``\nu`` 靠后指 ``j_\nu>j_\mu``）：

```math
\boxed{\;
\begin{aligned}
\text{同支不同节点}\;\; H &\mathrel{+}= M_1[\nu,\mu]\,X_{\nu}Y_{\mu} + M_2[\mu,\nu]\,Y_{\nu}X_{\mu},
&& \sigma_\nu=\sigma_\mu=+;\\
&\mathrel{+}= M_1[\nu,\mu]\,Y_{\nu}X_{\mu} + M_2[\mu,\nu]\,X_{\nu}Y_{\mu},
&& \sigma_\nu=\sigma_\mu=-\ (\tau\ \text{同}\ +\ \text{支});\\[4pt]
\text{同节点等时}\;\; H &\mathrel{+}= M_1[\nu\nu]\,Y_{\nu}X_{\nu} + M_2[\nu\nu]\,X_{\nu}Y_{\nu},
&& \sigma_\nu=+;\\
&\mathrel{+}= M_1[\nu\nu]\,X_{\nu}Y_{\nu} + M_2[\nu\nu]\,Y_{\nu}X_{\nu},
&& \sigma_\nu=-\ \text{或}\ \tau;\\[4pt]
\text{跨支}\ (+\nu,-\mu)\;\; H &\mathrel{+}= (M_1{+}M_2)[\nu,\mu]\,X_{\nu}X_{\mu}
+ (M_1{+}M_2)[\mu,\nu]\,Y_{\nu}Y_{\mu},
&& \text{从 } + \text{ 侧处理一次};\\[4pt]
\tau\nu\times\pm\mu\;\; H &\mathrel{+}= (M_1{+}M_2)[\nu,\mu]\,X_{\nu}\,Z_{\mu}
+ (M_1{+}M_2)[\mu,\nu]\,Y_{\nu}\,Z'_{\mu},
&& \text{从 } \tau \text{ 侧处理一次},
\end{aligned}
\;}
```

其中 ``\tau\times\pm`` 混合配对的 ``\pm`` 侧算符标签与跨支相反：

```math
Z_{\mu} = \begin{cases}Y_{\mu}, & \mu\text{ 为 } + \text{ 节点}\\ X_{\mu}, & \mu\text{ 为 } - \text{ 节点}\end{cases},
\qquad
Z'_{\mu} \text{ 取另一个}.
```

（物理来源：``b@\tau \leftrightarrow b^{\dagger}@\pm`` 的收缩把跨支的
「求和核配同算符」翻成「核留原槽位、± 侧换配对算符」。混合轮廓中
``\pm`` 节点之间的对仍用纯实时规则。）

稠密张量 = 显式矩阵指数（**不含任何 Prony/指数展开近似**）：

```math
I \;=\; \exp(H)\Big|_{\text{vec}} ,
\qquad \text{指标 } 2p{-}1 = \text{节点 } p \text{ 的 ket 路径},\ 2p = \text{bra 路径},
```

对 ``\mathrm{op}_1=\mathrm{op}_2``（厄米 ``S``）或对角 ``S``，以上退化为
``H = \sum M_1\,\mathrm{op}_1\mathrm{op}_2 + M_2\,\mathrm{op}_2\mathrm{op}_1``，
再退化为第 3 节的 ADT 对角嵌入。

## 5. 误差标度

顶点移位后，三种表示的核内容 = 严格 B-cell IF（同网格、同腿基的
exact 参照：费米为 `exact_fermionic_Δτ/Δt/Δm`，玻色为「精确 cell 演化
算符 + 单一 bath 迹」）的 **``O(V^{2})`` 系数逐位一致**（``V`` 为耦合
系数；以谱强度计即线性阶 ``O(\Gamma)``）。据此分三种情形：

1. **玻色对角 / ``\mathrm{op}_1=\mathrm{op}_2``（厄米）/ ``[\mathrm{op}_1,\mathrm{op}_2]=0``**：
   严格 IF 恰为 Gaussian ⟹ 组装**机器精确**（任意 ``\delta t``；负对照实测
   ``\sim10^{-15}`` 平线）。
2. **玻色 ``\mathrm{op}_1\neq\mathrm{op}_2``（如 ``S=\sigma^{+}``）**：严格 IF 含
   ``O(V^{4})`` 非高斯累积量（Magnus 对易子项），组装与严格 IF 的
   生成元一致到 ``O(V^{2})``，同网格残差按 ``V^{4}`` 标度。
3. **费米（冻结核）**：核闭式只保证 ``O(V^{2})`` 系数；同网格残差 =
   胞内两配对过程的 ``V^{4}`` 截断。整体 ``Z/\rho`` 另带一个
   ``O(V^{2})`` 归一化常数（「跳去又跳回」泄漏 = exact 参照的
   `prefactor`），迹归一化观测量中自动消去。

fixed-``N`` 标度实测（`performance/assembly_scaling.jl`，残差为逐 cell
局部量 ⟹ 同网格差值 = ``N\times`` 单胞截断，方案误差完全抵消）：

| 检验 | 残差比 |
|---|---|
| ``V`` 减半（固定网格；τ / RT / PT） | ``\times\,15\text{–}17``（``V^{4}``） |
| ``\delta t`` 加倍（fixed ``N``；费米 RT，``\beta`` 固定无漂移） | ``\times\,15\text{–}16``（首倍增） |
| ``\delta\tau`` 加倍（fixed ``N``；小 ``\omega`` 抑热漂移） | ``\times\,\sim14`` |
| 玻色 ``\sigma_x``（op1 = op2）负对照 | 机器精确平线 |

即**首阶同网格残差 ``\approx N\cdot V^{4}\delta t^{4}``/cell**（以谱强度
计 ``\propto\Gamma^{2}\delta t^{4}``）——塌缩参数为 ``V^{2}\delta t^{2}``
（单次胞内配对积分
``\int_{\text{cell}}\!\!\int_{\text{cell}} \Delta \propto V^{2}\delta t^{2}``）；
大 ``\delta t`` 端（``V^{2}\delta t\gtrsim0.1``）出现高阶修正使比值回落。
注意 fixed-``N`` 扫描测的是装配/核截断的 ``\delta t`` 标度；对连续极限的
全局收敛是方案误差 ``O(V^{2}\delta t^{3})/\text{cell}``（exact 核），
需固定 ``T`` 变 ``N`` 另测，不进同网格残差。

验证脚本：`performance/fermionicG4.jl`（费米 ``V`` 平台）、
`performance/assembly_scaling.jl`（fixed-``N`` ``\delta t/\delta\tau`` 扫描）、
`performance/nondiagcoupling.jl`（玻色 PT 严格 IF 对照）。
