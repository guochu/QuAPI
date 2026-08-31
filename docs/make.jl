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
    ],
    warnonly = false,
)

deploydocs(;
    repo = "github.com/guochu/QuAPI.git",
    devbranch = "master",
    push_preview = true,
)
