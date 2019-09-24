WeightedToSimple(graph) = LightGraphs.Graph( SimpleWeightedGraphs.adjacency_matrix( graph ) )

function countitems(x)
    un = unique(x)
    counts = [ sum(x .== i) for i in un ]
    return un .=> counts
end
