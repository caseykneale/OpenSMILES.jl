push!(LOAD_PATH, "/home/caseykneale/Desktop/SMILES.jl/")
import SMILES
using Pkg
Pkg.add("Compose");
Pkg.add("GraphPlot");
import Compose, GraphPlot, LightGraphs


Graph, Data = SMILES.ParseSMILES("C1CC12CC2")
GraphPlot.gplot( LightGraphs.Graph( LightGraphs.adjacency_matrix( Graph ) ) )



SMILES.ParseBracket("CH3-")
SMILES.ParseBracket("CH2--") == SMILES.ParseBracket("CH2-2")

function Base.:(==)(a::SMILES.Element, b::SMILES.Element)
    return all( [ a.symbol == b.symbol, a.isotope == b.isotope,
                a.aromatic == b.aromatic, a.ringID == b.ringID,
                a.explicithydrogens == b.explicithydrogens,
                a.charge == b.charge])
end
SMILES.ParseBracket("CH2--") == SMILES.ParseBracket("CH2-2")
