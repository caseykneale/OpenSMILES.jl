push!(LOAD_PATH, "/home/caseykneale/Desktop/SMILES.jl/")
import SMILES
#using Pkg
#Pkg.add("Compose");
#Pkg.add("GraphPlot");
import Compose, GraphPlot, LightGraphs


Graph, Data = SMILES.ParseSMILES("C1CC12CC2")
GraphPlot.gplot( LightGraphs.Graph( LightGraphs.adjacency_matrix( Graph ) ) )

sum( Graph.weights[2,:].nzval )

Graph.weights[2,:]
sum(Graph.weights[1,:].value)
#Now we need to flesh out the number of hydrogens in each bond
for ( i, atom ) in enumerate( Data )
    #Find all edges with this atom.
    implicitH = 0
    if atom.explicithydrogens == 0
        #Lookup valence
        if atom.symbol in keys( SMILES.valence )
            Valence = SMILES.valence[ atom.symbol ]
            BondedElectrons = sum(Graph.weights[:,i].nzval)
            Aromaticity = (isa(atom.aromatic, Nothing) ? 0 : atom.aromatic)
            implicitH = Valence - BondedElectrons - Aromaticity
            if implicitH < 0
                #Wiggle to find proper valence
                @warn("Illegal number of implicit hydrogens in $(atom.symbol) (Atom #$i). Defaulting to 0 hydrogens.")
                implicitH = 0
            end
        else
            #Valence not known to package :/
            @warn("The valence for $(atom.symbol) isn't known, implicit hydrogens not determined.")
        end
    else #User specified hydrogens explicitly - use them.
        implicitH = atom.explicithydrogens
    end

    Data[i].implicithydrogens = implicitH
end

Data
