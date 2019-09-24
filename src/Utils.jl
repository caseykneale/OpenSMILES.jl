WeightedToSimple(graph) = LightGraphs.Graph( SimpleWeightedGraphs.adjacency_matrix( graph ) )

function countitems(x)
    un = unique(x)
    counts = [ sum(x .== i) for i in un ]
    return un .=> counts
end

function EmpiricalFormula(data::Array{Element, 1} )
    Hydrogens = sum( H.( data ) )
    countem = Dict(countitems( abbreviation.( data ) ))
    if "H" in keys( countem )
        countem["H"] += Hydrogens
    else
        countem["H"] = Hydrogens
    end
    sortem = sort( collect( keys( countem ) ) )
    return join([ ele * string( countem[ele] ) for ele in sortem])
end
