push!(LOAD_PATH, "/home/caseykneale/Desktop/OpenSMILES.jl/");
import SMILES
using Pkg
import Compose, GraphPlot, Graphs, PeriodicTable
#Tryptophan
Graph, Data = SMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
#Bowtie
#Graph, Data = SMILES.ParseSMILES("C1CC12CC2")

GraphPlot.gplot( SMILES.WeightedToSimple( Graph ), nodelabel = SMILES.abbreviation.( Data ) )

SMILES.EmpiricalFormula( Data )

#Data
git remote set-url origin git@github.com:caseykneale/OpenSMILES.jl.git
