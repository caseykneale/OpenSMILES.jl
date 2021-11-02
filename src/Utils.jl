WeightedToSimple(graph) = Graphs.Graph( SimpleWeightedGraphs.adjacency_matrix( graph ) )

function countitems(x)
    un = unique(x)
    counts = [ sum(x .== i) for i in un ]
    return un .=> counts
end

function EmpiricalFormula(data::Array{Element, 1} )
    countem = Dict(countitems( abbreviation.( data ) ))
    nhydrogens = sum(H, data)
    if nhydrogens > 0
        countem["H"] = get(countem, "H", 0) + nhydrogens
    end
    sortem = sort( collect( keys( countem ) ) )
    return join([ (countem[ele] == 1) ? ele : ele * string( countem[ele] ) for ele in sortem])
end

"""
    gH, atomnodesH = instantiate_hydrogens(g, atomnodes)

Given a molecule, where `g` is the connectivity graph and `atomnodes` contains additional
information about the vertices, return a new representation `gH, atomnodesH` where all hydrogens
in `g`, `atomnodes` have been instantiated as full nodes in `gH, atomnodesH`.
"""
function instantiate_hydrogens(g::AbstractGraph, atomnodes::AbstractVector{Element})
    gH, atomnodesH = copy(g), [copy(atom) for atom in atomnodes]
    nhydrogens = sum(H, atomnodes)
    nhydrogens == 0 && return gH, atomnodesH
    lastH = nv(g)
    add_vertices!(gH, nhydrogens)
    for (i, atom) in enumerate(atomnodesH)
        for _ = 1:H(atom)
            add_edge!(gH, i, lastH += 1)
        end
        atom.explicithydrogens = atom.implicithydrogens = 0
    end
    append!(atomnodesH, [Element("H") for _ in 1:nhydrogens])
    return gH, atomnodesH
end
