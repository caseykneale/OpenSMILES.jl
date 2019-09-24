push!(LOAD_PATH, "/home/caseykneale/Desktop/SMILES.jl/");
import SMILES
using Pkg
import Compose, GraphPlot, LightGraphs, PeriodicTable
#Tryptophan
#Graph, Data = SMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
#Bowtie
#Graph, Data = SMILES.ParseSMILES("C1CC12CC2")
#Anthracene
#Graph, Data = SMILES.ParseSMILES("C1=CC=C2C=C3C=CC=CC3=CC2=C1")
#Lysergic Acid Diethylamide
#Graph, Data = SMILES.ParseSMILES("CCN(CC)C(=O)C1CN(C2CC3=CNC4=CC=CC(=C34)C2=C1)C")
#Firefly Luciferin
#Graph, Data = SMILES.ParseSMILES("O=C(O)[CH]1N=C(SC1)c2sc3cc(O)ccc3n2")

Graph, Data = SMILES.ParseSMILES("CCCC(C)C")
GraphPlot.gplot( SMILES.WeightedToSimple( Graph ), nodelabel = SMILES.abbreviation.( Data ) )

Data

SMILES.H.(Data)

countitems(SMILES.H.(Data))
