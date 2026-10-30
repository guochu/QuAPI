# =============================================================================
# benchmark/bosons/make_results.jl
#
# 生成误差 ~ 离散化步长 的对数-对数图与数据表，写入 benchmark/bosons/results/。
#   · imagtime_error_vs_dtau.{csv,svg} : 标准 Trotter（1/2/4 阶）与冻结顶点规则
#   · realtime_error_vs_dt.{csv,svg}   : 标准 Trotter（1/2/4 阶）
#
# 运行：julia --startup-file=no --project=. benchmark/bosons/make_results.jl
# =============================================================================
include(joinpath(@__DIR__, "..", "..", "performance", "common", "common.jl"))

const RESDIR = joinpath(@__DIR__, "results")

# ------------------------------ 极简 log-log SVG ------------------------------ #
struct Series
    label::String
    xs::Vector{Float64}
    ys::Vector{Float64}
    color::String
    marker::Symbol
end

const MARKERPATH = Dict(
    :circle   => "M-3.4,0a3.4,3.4 0 1,0 6.8,0a3.4,3.4 0 1,0 -6.8,0",
    :square   => "M-3.2,-3.2h6.4v6.4h-6.4z",
    :triangle => "M0,-3.8L3.6,3.0H-3.6z",
    :diamond  => "M0,-3.8L3.8,0L0,3.8L-3.8,0z",
    :cross    => "M-3,-3L3,3M-3,3L3,-3",
)

const SUPDIG = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                    '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹', '-' => '⁻')
"把整数指数写成 10⁰ / 10⁻² 形式。"
function deglabel(k::Int)
    s = k < 0 ? "⁻" : ""
    for c in string(abs(k))
        s *= SUPDIG[c]
    end
    return "10" * s
end

function svg_loglog(path; title, xlabel, ylabel, series::Vector{Series}, note="")
    W, H = 820, 500
    ml, mr, mt, mb = 86, 252, 58, 64
    pw, ph = W - ml - mr, H - mt - mb
    allx = vcat([s.xs for s in series]...)
    ally = vcat([s.ys for s in series]...)
    lxmin = Int(floor(minimum(log10.(allx)))); lxmax = Int(ceil(maximum(log10.(allx))))
    lymin = Int(floor(minimum(log10.(ally)))); lymax = Int(ceil(maximum(log10.(ally))))
    (lxmax - lxmin) < 1 && (lxmax = lxmin + 1)
    (lymax - lymin) < 1 && (lymax = lymin + 1)
    X(x) = ml + (log10(x) - lxmin) / (lxmax - lxmin) * pw
    Y(y) = mt + (lymax - log10(y)) / (lymax - lymin) * ph

    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" viewBox="0 0 $W $H">""")
    println(io, """<rect width="$W" height="$H" fill="#ffffff"/>""")
    println(io, """<g font-family="system-ui, 'Segoe UI', 'PingFang SC', 'Helvetica Neue', Arial, sans-serif">""")
    println(io, """<text x="$(W/2)" y="32" font-size="17" font-weight="600" fill="#1F2328" text-anchor="middle">$(title)</text>""")

    # 网格与刻度（十进位）
    for k in lxmin:lxmax
        px = X(10.0^k)
        println(io, """<line x1="$px" y1="$mt" x2="$px" y2="$(mt+ph)" stroke="#EAEEF2" stroke-width="1"/>""")
        println(io, """<text x="$px" y="$(mt+ph+22)" font-size="12" fill="#59666E" text-anchor="middle">$(deglabel(k))</text>""")
    end
    for k in lymin:lymax
        py = Y(10.0^k)
        println(io, """<line x1="$ml" y1="$py" x2="$(ml+pw)" y2="$py" stroke="#EAEEF2" stroke-width="1"/>""")
        println(io, """<text x="$(ml-10)" y="$(py+4)" font-size="12" fill="#59666E" text-anchor="end">$(deglabel(k))</text>""")
    end
    println(io, """<rect x="$ml" y="$mt" width="$pw" height="$ph" fill="none" stroke="#59666E" stroke-width="1.2"/>""")
    println(io, """<text x="$(ml+pw/2)" y="$(H-22)" font-size="14" fill="#1F2328" text-anchor="middle">$(xlabel)</text>""")
    println(io, """<text transform="translate(26,$(mt+ph/2)) rotate(-90)" font-size="14" fill="#1F2328" text-anchor="middle">$(ylabel)</text>""")
    note != "" && println(io, """<text x="$ml" y="$(H-6)" font-size="11" fill="#8A949C">$(note)</text>""")

    # 曲线 + 名义斜率虚线
    for s in series
        pts = join(["$(X(x)),$(Y(y))" for (x, y) in zip(s.xs, s.ys)], " ")
        println(io, """<polyline points="$pts" fill="none" stroke="$(s.color)" stroke-width="2"/>""")
        for (x, y) in zip(s.xs, s.ys)
            println(io, """<path transform="translate($(X(x)),$(Y(y)))" d="$(MARKERPATH[s.marker])" fill="$(s.color)" stroke="$(s.color)" stroke-width="1.4"/>""")
        end
        # 过最细网格点、斜率 = 实测局部阶 的参考虚线
        n = length(s.xs)
        if n >= 2 && s.ys[end] > 0 && s.ys[end-1] > 0
            p = log(s.ys[end] / s.ys[end-1]) / log(s.xs[end] / s.xs[end-1])
            x0, y0 = s.xs[end], s.ys[end]
            x1, y1 = s.xs[1], y0 * (s.xs[1] / x0)^p
            println(io, """<line x1="$(X(x0))" y1="$(Y(y0))" x2="$(X(x1))" y2="$(Y(y1))" stroke="$(s.color)" stroke-width="1" stroke-dasharray="5 5" opacity="0.75"/>""")
        end
    end

    # 图例（右侧）
    for (i, s) in enumerate(series)
        yb = mt + 16 + 30 * (i - 1)
        n = length(s.xs)
        p = n >= 2 && s.ys[end] > 0 && s.ys[end-1] > 0 ?
            log(s.ys[end] / s.ys[end-1]) / log(s.xs[end] / s.xs[end-1]) : NaN
        println(io, """<line x1="$(ml+pw+22)" y1="$(yb-4)" x2="$(ml+pw+52)" y2="$(yb-4)" stroke="$(s.color)" stroke-width="2.4"/>""")
        println(io, """<path transform="translate($(ml+pw+37),$(yb-4))" d="$(MARKERPATH[s.marker])" fill="$(s.color)" stroke="$(s.color)" stroke-width="1.2"/>""")
        lab = isnan(p) ? s.label : @sprintf("%s   (p≈%.2f)", s.label, p)
        println(io, """<text x="$(ml+pw+62)" y="$(yb)" font-size="12" fill="#1F2328">$(lab)</text>""")
    end
    println(io, """<text x="$(ml+pw+22)" y="$(mt+30*length(series)+24)" font-size="11" fill="#8A949C">虚线：过最细网格点、斜率=实测阶</text>""")

    println(io, "</g></svg>")
    write(path, String(take!(io)))
end

# ------------------------------ 扫描 ------------------------------ #
function imag_sweep(mdl, Ns, fracs)
    Gex = exact_Gτ(mdl)
    rows = Vector{NamedTuple}()
    for N in Ns
        ms = round.(Int, fracs .* N)
        ex = [Gex(mod(m * mdl.β / N, mdl.β)) for m in ms]
        sc = maximum(abs.(ex))
        errs = zeros(3)
        for (k, order) in enumerate((1, 2, 4))
            vals, _, _ = chain_Gτ(mdl, order, N; ms=ms)
            errs[k] = maximum(abs.(vals .- ex)) / sc
        end
        exx = [Gex(mod(m * mdl.β / N, mdl.β)) for m in 1:N-1]
        efr = zeros(2)
        for (k, p) in enumerate((1, 4))
            vals, _ = frozen_dyson_Gτ(mdl, N; cellrule=p)
            efr[k] = maximum(abs.(vals .- exx)) / maximum(abs.(exx))
        end
        push!(rows, (; N, dtau=mdl.β / N, trotter1=errs[1], trotter2=errs[2],
                     trotter4=errs[3], frozen=efr[1], frozen4=efr[2]))
    end
    return rows
end

function real_sweep(mdl, Ns, T, fracs, γ)
    gl = exact_greater_lesser(mdl)
    rows = Vector{NamedTuple}()
    for N in Ns
        ms = round.(Int, fracs .* N)
        ts = ms .* (T / N)
        gt_ex = [gl(t)[1] for t in ts]
        lt_ex = [gl(t)[2] for t in ts]
        sc = maximum(abs.(gt_ex))
        errs = zeros(3)
        for (k, order) in enumerate((1, 2, 4))
            gt, lt, _ = chain_greater_lesser(mdl, order, N, T; γ=γ, ms=ms)
            errs[k] = maximum(vcat(abs.(gt .- gt_ex), abs.(lt .- lt_ex))) / sc
        end
        push!(rows, (; N, dt=T / N, trotter1=errs[1], trotter2=errs[2], trotter4=errs[3]))
    end
    return rows
end

function main()
    isdir(RESDIR) || mkpath(RESDIR)
    mdl = boson_toulouse(; Lmodes=48)
    Ns = [8, 16, 32, 64, 128, 256]
    fracs = [0.25, 0.5, 0.75]

    # ---------------- 虚时 ----------------
    # 等时项算符恒等式自检（分开保留 ⟺ 系数相加，只差整体常数 e^{ηkj1}）
    cs = bosonic_Δτ(mdl.bath.spectrum, β=mdl.β, N=16, μ=mdl.μ).data
    d_id, r_id = equal_time_identity_check(cs.ηⱼₖ[1], cs.ηₖⱼ[1])
    println(@sprintf("等时项自检：max|Osplit − e^{ηkj1}Osum| = %.2e，归一化期望偏差 = %.2e",
                     d_id, abs(r_id)))
    rows = imag_sweep(mdl, Ns, fracs)
    open(joinpath(RESDIR, "imagtime_error_vs_dtau.csv"), "w") do f
        println(f, "N,dtau,trotter1,trotter2,trotter4,frozen,frozen_p4")
        for r in rows
            @printf(f, "%d,%.8e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                    r.N, r.dtau, r.trotter1, r.trotter2, r.trotter4, r.frozen, r.frozen4)
        end
    end
    svg_loglog(joinpath(RESDIR, "imagtime_error_vs_dtau.svg");
        title="虚时：G(τ) 相对误差 vs δτ（β=4，ϵ_d=1，48 模玻色浴）",
        xlabel="δτ", ylabel="max |G(τ)−G_exact| / max|G_exact|",
        series=[
            Series("Trotter 1 阶 (Lie, 严格 cell)", [r.dtau for r in rows], [r.trotter1 for r in rows], "#1f77b4", :circle),
            Series("Trotter 2 阶 (Strang, 严格 cell)", [r.dtau for r in rows], [r.trotter2 for r in rows], "#2ca02c", :square),
            Series("Trotter 4 阶 (Suzuki, 严格 cell)", [r.dtau for r in rows], [r.trotter4 for r in rows], "#d62728", :triangle),
            Series("冻结顶点 (节点重构)", [r.dtau for r in rows], [r.frozen for r in rows], "#9467bd", :diamond),
            Series("冻结顶点 (三次重构)", [r.dtau for r in rows], [r.frozen4 for r in rows], "#ff7f0e", :cross),
        ],
        note="Δτ 越小越准（左低右高）；虚线为实测局部斜率。冻结顶点两条曲线（节点/三次重构）都停在 1 阶；真 2/4 阶来自严格 cell。")

    # ---------------- 实时 ----------------
    γ = exact_gamma(mdl)
    T = 2.0
    rows2 = real_sweep(mdl, Ns, T, fracs, γ)
    open(joinpath(RESDIR, "realtime_error_vs_dt.csv"), "w") do f
        println(f, "N,dt,trotter1,trotter2,trotter4")
        for r in rows2
            @printf(f, "%d,%.8e,%.6e,%.6e,%.6e\n", r.N, r.dt, r.trotter1, r.trotter2, r.trotter4)
        end
    end
    svg_loglog(joinpath(RESDIR, "realtime_error_vs_dt.svg");
        title="实时：G>、G< 相对误差 vs δt（T=2，β=4，ϵ_d=1，48 模玻色浴）",
        xlabel="δt", ylabel="max |G>,G< − G_exact| / max|G>_exact|",
        series=[
            Series("Trotter 1 阶 (Lie, 严格 cell)", [r.dt for r in rows2], [r.trotter1 for r in rows2], "#1f77b4", :circle),
            Series("Trotter 2 阶 (Strang, 严格 cell)", [r.dt for r in rows2], [r.trotter2 for r in rows2], "#2ca02c", :square),
            Series("Trotter 4 阶 (Suzuki, 严格 cell)", [r.dt for r in rows2], [r.trotter4 for r in rows2], "#d62728", :triangle),
        ],
        note="平衡初态取严格 γ（等价于 Δm 虚时段充分收敛）；实时传播用单粒子链。")

    # ---------------- 实时（实轴、直积初态）冻结核 ----------------
    nsys = 0.3
    Aex = exact_return_amplitude(mdl)
    Ns3 = [16, 32, 64, 128]
    fracs3 = [0.25, 0.5, 0.75]
    rows3 = Vector{NamedTuple}()
    for N in Ns3
        ms = round.(Int, fracs3 .* N)
        sc = maximum(abs.(Aex(m * T / N)) for m in ms)
        fA, _ = frozen_retarded_A(mdl, N, T)
        ef = maximum(abs.(fA[ms .+ 1] .- [Aex(m * T / N) for m in ms])) / sc
        eo = zeros(3)
        for (k, order) in enumerate((1, 2, 4))
            v, _ = chain_return_amplitude(mdl, order, N, T)
            eo[k] = maximum(abs.(v[ms .+ 1] .- [Aex(m * T / N) for m in ms])) / sc
        end
        push!(rows3, (; dt=T / N, frozen=ef, o1=eo[1], o2=eo[2], o4=eo[3]))
    end
    open(joinpath(RESDIR, "realtime_frozen_error_vs_dt.csv"), "w") do f
        println(f, "dt,frozen_quapi_Dt_retarded,trotter1,trotter2,trotter4")
        for r in rows3
            @printf(f, "%.8e,%.6e,%.6e,%.6e,%.6e\n", r.dt, r.frozen, r.o1, r.o2, r.o4)
        end
    end
    svg_loglog(joinpath(RESDIR, "realtime_frozen_error_vs_dt.svg");
        title="实时（实轴、直积初态 n_sys=0.3）：A(t) 相对误差 vs δt",
        xlabel="δt", ylabel="max |A−A_exact| / max|A_exact|",
        series=[
            Series("冻结核 (bosonic_Δt, −Δ⁺⁺−Δ⁺₋)", [r.dt for r in rows3], [r.frozen for r in rows3], "#9467bd", :diamond),
            Series("Trotter 1 阶 (严格 cell)", [r.dt for r in rows3], [r.o1 for r in rows3], "#1f77b4", :circle),
            Series("Trotter 2 阶 (严格 cell)", [r.dt for r in rows3], [r.o2 for r in rows3], "#2ca02c", :square),
            Series("Trotter 4 阶 (严格 cell)", [r.dt for r in rows3], [r.o4 for r in rows3], "#d62728", :triangle),
        ],
        note="初态为直积态，参考解 toulouse_neq_greater_lesser；冻结核为真 1 阶。")

    # ---------------- 实时（实轴、直积初态）冻结核 ----------------
    nsys = 0.3
    Aex = exact_return_amplitude(mdl)
    Ns3 = [16, 32, 64, 128]
    fracs3 = [0.25, 0.5, 0.75]
    rows3 = Vector{NamedTuple}()
    for N in Ns3
        ms = round.(Int, fracs3 .* N)
        sc = maximum(abs.(Aex(m * T / N)) for m in ms)
        fA, _ = frozen_retarded_A(mdl, N, T)
        ef = maximum(abs.(fA[ms .+ 1] .- [Aex(m * T / N) for m in ms])) / sc
        eo = zeros(3)
        for (k, order) in enumerate((1, 2, 4))
            v, _ = chain_return_amplitude(mdl, order, N, T)
            eo[k] = maximum(abs.(v[ms .+ 1] .- [Aex(m * T / N) for m in ms])) / sc
        end
        push!(rows3, (; dt=T / N, frozen=ef, o1=eo[1], o2=eo[2], o4=eo[3]))
    end
    open(joinpath(RESDIR, "realtime_frozen_error_vs_dt.csv"), "w") do f
        println(f, "dt,frozen_quapi_Dt_retarded,trotter1,trotter2,trotter4")
        for r in rows3
            @printf(f, "%.8e,%.6e,%.6e,%.6e,%.6e\n", r.dt, r.frozen, r.o1, r.o2, r.o4)
        end
    end
    svg_loglog(joinpath(RESDIR, "realtime_frozen_error_vs_dt.svg");
        title="实时（实轴、直积初态 n_sys=0.3）：A(t) 相对误差 vs δt",
        xlabel="δt", ylabel="max |A−A_exact| / max|A_exact|",
        series=[
            Series("冻结核 (bosonic_Δt, −Δ⁺⁺−Δ⁺₋)", [r.dt for r in rows3], [r.frozen for r in rows3], "#9467bd", :diamond),
            Series("Trotter 1 阶 (严格 cell)", [r.dt for r in rows3], [r.o1 for r in rows3], "#1f77b4", :circle),
            Series("Trotter 2 阶 (严格 cell)", [r.dt for r in rows3], [r.o2 for r in rows3], "#2ca02c", :square),
            Series("Trotter 4 阶 (严格 cell)", [r.dt for r in rows3], [r.o4 for r in rows3], "#d62728", :triangle),
        ],
        note="初态为直积态，参考解 toulouse_neq_greater_lesser；冻结核为真 1 阶。")

    # ---------------- 控制台摘要 ----------------
    println("写入：")
    for f in ("imagtime_error_vs_dtau.csv", "imagtime_error_vs_dtau.svg",
              "realtime_error_vs_dt.csv", "realtime_error_vs_dt.svg",
              "realtime_frozen_error_vs_dt.csv", "realtime_frozen_error_vs_dt.svg")
        println("  ", joinpath(RESDIR, f))
    end
    println("\n虚时数据（δτ, 1阶, 2阶, 4阶, 冻结节点, 冻结三次）：")
    for r in rows
        @printf("  %.4e  %.3e  %.3e  %.3e  %.3e  %.3e\n",
                r.dtau, r.trotter1, r.trotter2, r.trotter4, r.frozen, r.frozen4)
    end
    println("\n实时数据（δt, 1阶, 2阶, 4阶）：")
    for r in rows2
        @printf("  %.4e  %.3e  %.3e  %.3e\n", r.dt, r.trotter1, r.trotter2, r.trotter4)
    end
    println("\n实时（直积初态）冻结核数据（δt, 冻结, 1阶, 2阶, 4阶）：")
    for r in rows3
        @printf("  %.4e  %.3e  %.3e  %.3e  %.3e\n", r.dt, r.frozen, r.o1, r.o2, r.o4)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
