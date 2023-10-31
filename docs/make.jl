using RubberBearings
using Documenter

DocMeta.setdocmeta!(RubberBearings, :DocTestSetup, :(using RubberBearings); recursive=true)

makedocs(;
    modules=[RubberBearings],
    authors="Eduardo J. Montalto",
    repo="https://github.com/eduardomontalto/RubberBearings.jl/blob/{commit}{path}#{line}",
    sitename="RubberBearings.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://eduardomontalto.github.io/RubberBearings.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/eduardomontalto/RubberBearings.jl",
    devbranch="master",
)
