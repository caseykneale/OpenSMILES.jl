struct SMILESParseException <: Exception
    msg::String
    smiles::String
    charloc::Int
end
Base.showerror(io::IO, e::SMILESParseException) = print(io, "Failed to parse SMILES on character ", e.charloc, "!")

function ParseSMILES(s::AbstractString, calculate_implicit_hydrogens = true )
    MoleculeGraph = SimpleWeightedGraph();
    MolecularData = Element[]

    len = length(s)

    chainstack = Int16[]

    nextedgestart = 0
    weight = 1
    lastcursor = 0
    RingClosures = Dict{Int,Vector{Int}}()

    idx = 1
    while idx <= len
        moiety = nothing
        cursor = s[idx]

        if isdigit(cursor) || cursor == '%' #Handle rings
            if isdigit(cursor)
                # Single-digit ring id
                id = Int16(cursor - '0')
                idx = nextind(s, idx)
            else
                # Two-digit ring id
                id = parse(Int16, s[idx+1:idx+2])
                idx += 3
            end
            push!( MolecularData[end].ringID, id )
            if id in keys( RingClosures )
                push!( RingClosures[ id ], length(MolecularData) )
            else
                RingClosures[ id ] = [ length(MolecularData) ]
            end
        elseif isletter( cursor )
            symbollist = isuppercase( cursor ) ? aliphatics : aromatics
            aromatic = islowercase( cursor )
            symbol, idx = tryparseelement(s, idx, symbollist)
            if isa(symbol, Nothing)
                throw(SMILESParseException("Invalid SMILES. Unrecognized chemical symbol.", s, idx))
            end
            moiety = Element(symbol)
        elseif isspecialoperator(cursor)    #Handle operators
            if cursor == '['
                idxclose = findnext(']', s, idx)
                idxclose === nothing && throw(SMILESParseException("Invalid SMILES. Bracket is either empty, or does not end.", s, idx))
                moiety = ParseBracket(SubString(s, idx+1, idxclose-1))
                idx = nextind(s, idxclose)
            end
            if cursor == '('
                if lastcursor == 0
                    push!( chainstack, length( MolecularData )  )
                else
                    push!( chainstack, lastcursor  )
                end
                idx = nextind(s, idx)
            end
            if cursor == ')'
                isempty(chainstack) && throw(SMILESParseException("Invalid SMILES. Parenthesis/chain does not have a beginning.", s, idx))
                idx = nextind(s, idx)
                lastcursor = nextedgestart = pop!( chainstack )
            end
        elseif isbondoperator(cursor) #Handle bonds
            weight = bonds[ cursor ]
            idx = nextind(s, idx)
        end
        #New atom/moiety was parsed
        if isa( moiety, Element )
            lastcursor = 0
            add_vertex!( MoleculeGraph )
            push!( MolecularData, moiety )
            l = length( MolecularData )
            if (l > 1)
                if nextedgestart == 0
                    add_edge!(MoleculeGraph, l - 1, l, weight )
                else
                    add_edge!(MoleculeGraph, nextedgestart, l, weight )
                    nextedgestart = 0
                end
                if weight > 1
                    weight = 1
                end
            end
        end
    end
    #Close rings!
    for (k, v) in RingClosures
        for id in 1 : ( length( v ) - 1 )
            if ( v[id+1] - v[id] ) != -1
                add_edge!(MoleculeGraph, v[id], v[id+1])
            end
        end
    end
    # Add implicit hydrogens to graph structure
    if calculate_implicit_hydrogens
        for ( i, atom ) in enumerate( MolecularData )
            #Find all edges with this atom.
            implicitH = 0
            #Lookup valence
            if atom.symbol in keys( valence )
                Valence = valence[ atom.symbol ]
                BondedElectrons = sum( MoleculeGraph.weights[:,i].nzval )
                Aromaticity = (isa(atom.aromatic, Nothing) ? 0 : atom.aromatic)
                implicitH = Valence .- BondedElectrons .- Aromaticity
                #Look I'm not keeping track of non-bonding e- pairs
                #All I'm doing here is seeing if this is likely invalid structure...
                if length(implicitH) > 1
                    if all( implicitH .< 0 )
                        #Wiggle to find proper valence
                        @warn("Illegal number of implicit hydrogens in $(atom.symbol) (Atom #$i). Defaulting to 0 hydrogens.")
                        implicitH = 0
                    else
                        implicitH = implicitH[ implicitH .> 0 ]
                        implicitH = reduce(min, implicitH .> 0)
                    end
                end
            else
                #Valence not known to package :/
                @warn("The valence for $(atom.symbol) isn't known, implicit hydrogens not determined.")
            end
            # implicitH >= atom.explicithydrogens || error("more hydrogens were supplied than supported by the valence in $Sorig (at $S)")
            MolecularData[i].implicithydrogens = max(0, implicitH - atom.explicithydrogens)
        end
    end

    return MoleculeGraph, MolecularData
end
