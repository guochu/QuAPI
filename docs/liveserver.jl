# ============================================================================
# liveserver.jl —— 文档本地实时预览
#
#   julia --project=docs docs/liveserver.jl [port]
#
# 首次自动构建 docs/build，然后监听 docs/src/*.md：保存后自动重建，
# LiveServer 检测到 build 目录变化后浏览器自动刷新。用 IDE 内置浏览器或
# 任意浏览器打开终端提示的地址（默认 http://localhost:8000/）即可预览。
#
# 首次使用前需安装 LiveServer（只需一次）：
#   julia --project=docs -e 'using Pkg; Pkg.add("LiveServer")'
#
# make.jl 末尾的 deploydocs 在非 CI 环境会被 Documenter 自动跳过
# （"Documenter could not auto-detect the building environment" 警告属预期）。
# ============================================================================

using LiveServer

const DOCS = @__DIR__
const SRC = joinpath(DOCS, "src")
const BUILD = joinpath(DOCS, "build")

function newest_src_mtime()
    t = 0.0
    for f in readdir(SRC; join = true)
        endswith(f, ".md") || continue
        t = max(t, stat(f).mtime)
    end
    return t
end

# 首次构建
include(joinpath(DOCS, "make.jl"))

# 监听源文件：变化 → 重建（构建失败保留上一次 build，修复后保存即可）
@async begin
    t = newest_src_mtime()
    while true
        sleep(1)
        t2 = newest_src_mtime()
        t2 > t || continue
        t = t2
        @info "检测到 docs/src 变化，重新构建…"
        try
            include(joinpath(DOCS, "make.jl"))
        catch err
            @warn "构建失败，浏览器中仍是上一次的版本" err
        end
    end
end

port = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 8000
@info "开始预览：浏览器打开 http://localhost:$(port)/ （Ctrl-C 退出）"
serve(dir = BUILD; port = port)
