push!(LOAD_PATH, "/home/caseykneale/Desktop/SMILES.jl/");
import SMILES
using Pkg
import Compose, GraphPlot, LightGraphs, PeriodicTable
#Tryptophan
Graph, Data = SMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
#Bowtie
#Graph, Data = SMILES.ParseSMILES("C1CC12CC2")

GraphPlot.gplot( SMILES.WeightedToSimple( Graph ), nodelabel = SMILES.abbreviation.( Data ) )

SMILES.EmpiricalFormula( Data )

#Data
