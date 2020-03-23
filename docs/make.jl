using Documenter, SMILES

makedocs(;
    modules=[SMILES],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/caseykneale/SMILES.jl/blob/{commit}{path}#L{line}",
    sitename="SMILES.jl",
    authors="Casey Kneale",
    assets=String[],
)

deploydocs(;
    repo="github.com/caseykneale/SMILES.jl",
)
