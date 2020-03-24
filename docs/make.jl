using Documenter, OpenSMILES

makedocs(;
    modules=[SMILES],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/caseykneale/OpenSMILES.jl/blob/{commit}{path}#L{line}",
    sitename="OpenSMILES.jl",
    authors="Casey Kneale",
    assets=String[],
)

deploydocs(;
    repo="github.com/caseykneale/OpenSMILES.jl",
)
