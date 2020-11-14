WeightedToSimple(graph) = LightGraphs.Graph( SimpleWeightedGraphs.adjacency_matrix( graph ) )

function countitems(x)
    un = unique(x)
    counts = [ sum(x .== i) for i in un ]
    return un .=> counts
end

function EmpiricalFormula(data::Array{Element, 1} )
    countem = Dict(countitems( abbreviation.( data ) ))
    nhydrogens = sum( H.( data ) )
    if nhydrogens > 0
        countem["H"] = get(countem, "H", 0) + nhydrogens
    end
    sortem = sort( collect( keys( countem ) ) )
    return join([ (countem[ele] == 1) ? ele : ele * string( countem[ele] ) for ele in sortem])
end
