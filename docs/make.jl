using Documenter
using QuAPI

makedocs(;
    modules = [QuAPI],
    sitename = "QuAPI",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = nothing,
        edit_link = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => "guide.md",
        "API" => "api.md",
        "Hybridization Functions" => "hybridization.md",
        "IF OSEE" => "osee.md",
    ],
    # Dense / contourindex 的 docstring 尚未全部收录进 api.md，降级为警告
    warnonly = [:docs_block, :missing_docs],
)

# ---------------------------------------------------------------------------
# 本地自托管 KaTeX：把 HTML 里的 cdnjs 地址改写为 build 内的本地资源。
# cdnjs 首次加载实测可超过 RequireJS 的 7 s 等待上限，公式会停留在未渲染的
# 原始 LaTeX 状态；本地资源在 file:// 下也能正常加载。KaTeX 资源
# （0.16.8，MIT）位于 docs/src/assets/katex/，makedocs 会拷到 build/assets/。
# ---------------------------------------------------------------------------
const KATEX_CDN = "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.8"
const BUILD = joinpath(@__DIR__, "build")
for (root, _, files) in walkdir(BUILD), f in files
    endswith(f, ".html") || continue
    path = joinpath(root, f)
    # 页面相对 build 根的前缀：根目录页为 ""，目录式页（CI prettyurls）为 "../"
    prefix = let r = relpath(BUILD, root)
        r == "." ? "" : r * "/"
    end
    s = read(path, String)
    s = replace(s,
        "$(KATEX_CDN)/katex.min.css" => "$(prefix)assets/katex/katex.min.css",
    )
    write(path, s)
end

# KaTeX 的 JS 由 documenter.js 在运行时经 RequireJS 加载。页面以
# `data-main="assets/documenter.js"` 引导，RequireJS 的 baseUrl 即 assets/
# 目录，因此路径写作 "katex/..."（RequireJS 自动补 .js）。
docjs = joinpath(BUILD, "assets", "documenter.js")
s = read(docjs, String)
s = replace(s,
    "$(KATEX_CDN)/katex.min" => "katex/katex.min",
    "$(KATEX_CDN)/contrib/auto-render.min" => "katex/auto-render.min",
)
write(docjs, s)

deploydocs(;
    repo = "github.com/guochu/QuAPI.git",
    devbranch = "master",
    push_preview = true,
)
