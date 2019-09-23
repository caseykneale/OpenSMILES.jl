push!(LOAD_PATH, "/home/caseykneale/Desktop/SMILES.jl/")
import SMILES
using Pkg
Pkg.add("Compose");
Pkg.add("GraphPlot");
import Compose, GraphPlot, LightGraphs


Graph, Data = SMILES.ParseSMILES("C1CC12CC2")
GraphPlot.gplot( LightGraphs.Graph( LightGraphs.adjacency_matrix( Graph ) ) )
