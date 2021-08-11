using Documenter
include("../src/Kelp.jl")
Documenter.makedocs(
    root = "./",
    source = "src",
    build = "build",
    clean = true,
    doctest = true,
    modules = Module[Kelp],
    repo = "https://github.com/jagoosw/Kelp",
    highlightsig = true,
    sitename = "Kelp Documentation",
    expandfirst = [],
    pages = [
        "Index" => "index.md",
    ],
    format = Documenter.HTML(prettyurls = false)
)
