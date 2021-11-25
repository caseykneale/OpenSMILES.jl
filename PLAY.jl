# push!(LOAD_PATH, "/home/caseykneale/Desktop/OpenSMILES.jl/");
# import OpenSMILES
using Pkg
Pkg.activate(".")

import GraphPlot, LightGraphs#, Compose, PeriodicTable
#Tryptophan
Graph, Data = OpenSMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
#Bowtie
#Graph, Data = SMILES.ParseSMILES("C1CC12CC2")

GraphPlot.gplot( OpenSMILES.WeightedToSimple( Graph ), nodelabel = OpenSMILES.abbreviation.( Data ) )

SMILES.EmpiricalFormula( Data )

#Data
