using Condor
using Documenter

DocMeta.setdocmeta!(Condor, :DocTestSetup, :(using Condor); recursive=true)

makedocs(;
    modules=[Condor],
    authors="Yusuke Takase",
    repo="https://github.com/yusuke-takase/Condor.jl/blob/{commit}{path}#{line}",
    sitename="Condor.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://yusuke-takase.github.io/Condor.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yusuke-takase/Condor.jl",
)
