# OpenSMILES.jl
| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://caseykneale.github.io/OpenSMILES.jl/stable) | [![Build Status](https://travis-ci.com/caseykneale/OpenSMILES.jl.svg?branch=master)](https://travis-ci.com/caseykneale/OpenSMILES.jl)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://caseykneale.github.io/OpenSMILES.jl/dev)
[![Codecov](https://codecov.io/gh/caseykneale/OpenSMILES.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/caseykneale/OpenSMILES.jl)



This is a SMILES parser in Julia that aims to follow the [OpenSMILES format](http://opensmiles.org/opensmiles.html) (to the best of my ability). Theres probably bugs, this isn't inventive its just a parser that turns SMILES into a weighted Graphs.jl graph. Contributions welcome!

# Notice
Although this package does mostly what it intends too (still missing chiral support, etc), an excellent package[MolecularGraph](https://github.com/mojaie/MolecularGraph.jl) has been released. It is highly reccommended users contribute to that package instead of this one, unless they specifically want Graphs.jl integration.

# Examples

## Tryptophan
```Julia
using OpenSMILES, GraphPlot

# Tryptophan
Graph, Data = OpenSMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
GraphPlot.gplot( OpenSMILES.WeightedToSimple( Graph ), nodelabel = OpenSMILES.abbreviation.( Data ) )
OpenSMILES.EmpiricalFormula( Data ) # C11H12N2O2
```
![Tryptophan](https://raw.githubusercontent.com/caseykneale/OpenSMILES.jl/master/output/Tryptophan.png)

```Julia
# Bowtie ( not a real molecule :P )
Graph, Data = OpenSMILES.ParseSMILES("C1CC12CC2")
GraphPlot.gplot( OpenSMILES.WeightedToSimple( Graph ), nodelabel = OpenSMILES.abbreviation.( Data ) )
OpenSMILES.EmpiricalFormula( Data ) #C5H8
```
![Bowtie](https://raw.githubusercontent.com/caseykneale/OpenSMILES.jl/master/output/Bowtie.png)

Cool? Enjoy!
