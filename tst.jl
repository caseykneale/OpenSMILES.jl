using OpenSMILES
using GraphPlot
using Plots
using Cairo, Compose
# Tryptophan
Graph, Data = OpenSMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
g = GraphPlot.gplot( Graph, nodelabel = OpenSMILES.abbreviation.( Data ) )
OpenSMILES.EmpiricalFormula( Data )
draw( PNG("output/Tryptophan.png"), g )

# Bowtie ( not a real molecule :P )
Graph, Data = OpenSMILES.ParseSMILES("C1CC12CC2")
g = GraphPlot.gplot( Graph , nodelabel = OpenSMILES.abbreviation.( Data ) )
OpenSMILES.EmpiricalFormula( Data )

draw( PNG("output/Bowtie.png"), g )
