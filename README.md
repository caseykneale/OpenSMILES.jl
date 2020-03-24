# OpenSMILES.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://caseykneale.github.io/OpenSMILES.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://caseykneale.github.io/OpenSMILES.jl/dev)
[![Build Status](https://travis-ci.com/caseykneale/OpenSMILES.jl.svg?branch=master)](https://travis-ci.com/caseykneale/OpenSMILES.jl)
[![Codecov](https://codecov.io/gh/caseykneale/OpenSMILES.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/caseykneale/OpenSMILES.jl)
[![Coveralls](https://coveralls.io/repos/github/caseykneale/OpenSMILES.jl/badge.svg?branch=master)](https://coveralls.io/github/caseykneale/OpenSMILES.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/caseykneale/OpenSMILES.jl.svg)](https://cirrus-ci.com/github/caseykneale/OpenSMILES.jl)


This is a SMILES parser in Julia following the OpenSMILES format (to the best of my ability). Theres probably bugs, this isn't inventive its just a parser that turns SMILES into a weighted LightGraphs graph. Contributions welcome!

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
