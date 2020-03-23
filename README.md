# SMILES

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://caseykneale.github.io/SMILES.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://caseykneale.github.io/SMILES.jl/dev)
[![Build Status](https://travis-ci.com/caseykneale/SMILES.jl.svg?branch=master)](https://travis-ci.com/caseykneale/SMILES.jl)
[![Codecov](https://codecov.io/gh/caseykneale/SMILES.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/caseykneale/SMILES.jl)
[![Coveralls](https://coveralls.io/repos/github/caseykneale/SMILES.jl/badge.svg?branch=master)](https://coveralls.io/github/caseykneale/SMILES.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/caseykneale/SMILES.jl.svg)](https://cirrus-ci.com/github/caseykneale/SMILES.jl)


This is a SMILES parser in Julia following the OpenSMILES format (to the best of my ability). Theres probably bugs, this isn't inventive its just a parser that turns SMILES into a weighted LightGraphs graph. Contributions welcome!

# Examples

## Tryptophan
```Julia
using SMILES, GraphPlot

# Tryptophan
Graph, Data = SMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
GraphPlot.gplot( SMILES.WeightedToSimple( Graph ), nodelabel = SMILES.abbreviation.( Data ) )
SMILES.EmpiricalFormula( Data ) # C11H12N2O2
```
![Tryptophan](https://raw.githubusercontent.com/caseykneale/SMILES.jl/master/output/Tryptophan.png)

```Julia
# Bowtie ( not a real molecule :P )
Graph, Data = SMILES.ParseSMILES("C1CC12CC2")
GraphPlot.gplot( SMILES.WeightedToSimple( Graph ), nodelabel = SMILES.abbreviation.( Data ) )
SMILES.EmpiricalFormula( Data ) #C5H8
```
![Bowtie](https://raw.githubusercontent.com/caseykneale/SMILES.jl/master/output/Bowtie.png)

Enjoy!
