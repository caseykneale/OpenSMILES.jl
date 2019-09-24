struct SMILESParseException <: Exception
    charloc::Int64
end
Base.showerror(io::IO, e::SMILESParseException) = print(io, "Failed to parse SMILES on character ", e.charloc, "!")

function ParseSMILES( S::String, calculate_implicit_hydrogens = true )
    MoleculeGraph = SimpleWeightedGraph();
    MolecularData = Element[]

    Sorig = deepcopy(S)
    origlen = length(S)
    curlen = deepcopy(origlen)
    lastlen = deepcopy(origlen)

    chainstack = Int16[]

    nextedgestart = 0
    weight = 1
    lastcursor = 0
    RingClosures = Dict()

    while curlen > 0
        moiety = nothing
        cursor = S[1]

        if isdigit( cursor )  #Handle rings
            id = parse(Int16, cursor)
            push!( MolecularData[end].ringID, id )
            if id in keys( RingClosures )
                push!( RingClosures[ id ], length(MolecularData) )
            else
                RingClosures[ id ] = [ length(MolecularData) ]
            end
            S = S[ 2 : end]
        elseif cursor == '%'#2 decimal ring
            push!(MolecularData[end].ringID, parse(Int16, S[ (BracketClose+1) : (BracketClose+2) ] ) )
            S = S[ 4 : end]
        elseif isletter( cursor )
            symbollist = isuppercase( cursor ) ? aliphatics : aromatics
            aromatic = islowercase( cursor )
            symbol, S = ReadNextElement( S, symbollist )
            if isa(symbol, Nothing)
                @warn("Invalid SMILES. Unrecognized chemical symbol.")
            else
                moiety = Element(symbol)
            end
        elseif isspecialoperator(cursor)    #Handle operators
            if cursor == '['
                BracketClose = findfirst( [ s == ']' for s in S ] ) - 1
                if (BracketClose > 0) && !isa(BracketClose, Nothing)
                    moiety = ParseBracket( S[ 2 : BracketClose ] )
                    S = S[ (BracketClose+2) : end]
                else
                    @warn("Invalid SMILES. Bracket is either empty, or does not end.")
                end
            end
            if cursor == '('
                if lastcursor == 0
                    push!( chainstack, length( MolecularData )  )
                else
                    push!( chainstack, lastcursor  )
                end
                S = S[ 2 : end ]
            end
            if cursor == ')'
                if length(chainstack) > 0
                    S = S[ 2 : end ]
                    nextedgestart = pop!( chainstack )
                    lastcursor = nextedgestart
                else
                    @warn("Invalid SMILES. Parenthesis/chain does not have a beginning.")
                end
            end
        elseif isbondoperator(cursor) #Handle bonds
            weight = bonds[ string(cursor) ]
            S = S[ 2 : end ]
        end
        #New atom/moiety was parsed
        if isa( moiety, Element )
            lastcursor = 0
            add_vertex!( MoleculeGraph )
            push!( MolecularData, moiety )
            len = length( MolecularData )
            if (len > 1)
                if nextedgestart == 0
                    add_edge!(MoleculeGraph, len - 1, len, weight )
                else
                    add_edge!(MoleculeGraph, nextedgestart, len, weight )
                    nextedgestart = 0
                end
                if weight > 1
                    weight = 1
                end
            end
        end
        #Ensure we don't get stuck in an infinite loop
        curlen = length( S )
        if lastlen == curlen
            throw(SMILESParseException( (origlen - curlen + 1) ))
        end
        lastlen = curlen
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
            if atom.explicithydrogens == 0
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
            else #User specified hydrogens explicitly - use them.
                implicitH = atom.explicithydrogens
            end

            MolecularData[i].implicithydrogens = implicitH
        end
    end

    return MoleculeGraph, MolecularData
end
