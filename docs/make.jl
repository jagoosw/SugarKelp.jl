using Documenter
include("../src/SugarKelp.jl")
Documenter.makedocs(
    root = "docs",
    source = "src",
    build = "build",
    clean = true,
    doctest = true,
    modules = Module[SugarKelp],
    repo = "https://github.com/jagoosw/SugarKelp.jl",
    highlightsig = true,
    sitename = "SugarKelp.jl Documentation",
    expandfirst = [],
    pages = [
        "Index" => "index.md",
        "Function Documentation" => "functions.md"
    ],
    format = Documenter.HTML(prettyurls = false)
)

deploydocs(
    repo = "github.com/jagoosw/SugarKelp.jl.git",
    devbranch = "main"
)