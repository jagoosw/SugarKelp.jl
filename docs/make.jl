using Documenter
include("../src/Kelp.jl")
Documenter.makedocs(
    root = "docs",
    source = "src",
    build = "build",
    clean = true,
    doctest = true,
    modules = Module[Kelp],
    repo = "https://github.com/jagoosw/Kelp.jl",
    highlightsig = true,
    sitename = "Kelp.jl Documentation",
    expandfirst = [],
    pages = [
        "Index" => "index.md",
        "Function Documentation" => "functions.md"
    ],
    format = Documenter.HTML(prettyurls = false)
)

deploydocs(
    repo = "github.com/jagoosw/Kelp.jl.git",
    devbranch = "main"
)