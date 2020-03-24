using OpenSMILES
using GraphPlot
using Plots
using Cairo, Compose
# Tryptophan
Graph, Data = OpenSMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
g = GraphPlot.gplot( OpenSMILES.WeightedToSimple( Graph ), nodelabel = OpenSMILES.abbreviation.( Data ) )
OpenSMILES.EmpiricalFormula( Data )
draw( PNG("/home/caseykneale/.julia/dev/SMILES/output/Tryptophan.png"), g )

# Bowtie ( not a real moleculer :P )
Graph, Data = OpenSMILES.ParseSMILES("C1CC12CC2")
g = GraphPlot.gplot( OpenSMILES.WeightedToSimple( Graph ), nodelabel = OpenSMILES.abbreviation.( Data ) )
OpenSMILES.EmpiricalFormula( Data )

draw( PNG("/home/caseykneale/.julia/dev/SMILES/output/Bowtie.png"), g )
