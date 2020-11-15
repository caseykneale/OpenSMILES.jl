module OpenSMILES
    using SimpleWeightedGraphs, LightGraphs
    import Base: ==
    include("Elements.jl")
    export bondoperators, specialoperators, operators, bonds, valence,
            aliphatics, aromatics, bracket, bracket_aromatic, isoperator,
            isbondoperator, isspecialoperator, ispm, FindNumerics, FindPMs,
            Element, GraphElement, abbreviation, implicitH, explicitH, H,
            charge, isotope, ==

    include("ReadFns.jl")
    export ReadNextElement, ReadNextNumeric, ReadNextCharge

    include("Brackets.jl")
    export BracketState, BracketParseException, ParseBracket

    include("SMILESParser.jl")
    export SMILESParseException, ParseSMILES

    include("Utils.jl")
    export WeightedToSimple, countitems

    include("deprecations.jl")
end
