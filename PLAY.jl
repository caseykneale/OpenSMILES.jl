push!(LOAD_PATH, "/home/caseykneale/Desktop/SMILES.jl/");
import SMILES
import Compose, GraphPlot, LightGraphs
#Tryptophan
#Graph, Data = SMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")

#Bowtie
#Graph, Data = SMILES.ParseSMILES("C1CC12CC2")

#Anthracene
Graph, Data = SMILES.ParseSMILES("C1=CC=C2C=C3C=CC=CC3=CC2=C1")
GraphPlot.gplot( LightGraphs.Graph( LightGraphs.adjacency_matrix( Graph ) ) )

Data
