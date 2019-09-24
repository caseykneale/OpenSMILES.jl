WeightedToSimple(graph) = LightGraphs.Graph( SimpleWeightedGraphs.adjacency_matrix( graph ) )

#immutable struct to hold elemental information - less memory then "Element"
struct GraphElement
    symbol::String
    isotope::Union{Int16, Nothing}
    aromatic::Bool
    hydrogens::Int8
    charge::Int8
end

GraphElement(E::Element) = GraphElement(E.symbol, E.isotope, E.aromatic,
                                        E.hydrogens, E.charge)
