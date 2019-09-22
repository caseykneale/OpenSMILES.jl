bondoperators = ["@", "-", "\\", "/", "=", "#", "\$"];
specialoperators = [ ".", "]", "[", "(", ")"];
operators = [ specialoperators;  bondoperators ];

bonds = Dict( "-" => 1, "\\" => 1, "/" => 1, "=" => 2, "#" => 3, "\$" => 4 );
valence = Dict( "B" => 3, "C" => 4, "N" => (3, 5), "O" => 2, "P" => (3, 5),
                "S" => (2, 4, 6), "F" => 1, "Cl" => 1, "Br" => 1, "I" => 1 );

aliphatics  = ["B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"];
aromatics   = ["b", "c", "n", "o", "s", "p"];

#bracket_atom ::= "[" isotope? symbol chiral? hcount? charge? class? "]"
bracket = [ "H", "He" ,"Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
            "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
            "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
            "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
            "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Hf",
            "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
            "Po", "At", "Rn", "Fr", "Ra", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
            "Ds", "Rg", "Cn", "Fl", "Lv", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
            "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th",
            "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
            "No", "Lr", "Mg"];

bracket_aromatic = ["b", "c", "n", "o", "p", "s", "se", "as"];

isoperator(x) = any( x .== [operators; Symbol.(operators); first.(operators) ] );
isbondoperator(x) = any( x .== [bondoperators; Symbol.(bondoperators); first.(bondoperators) ] );
isspecialoperator(x) = any( x .== [specialoperators; Symbol.(specialoperators); first.(specialoperators) ] );
ispm(x) = any( x .== ["+", "-", :+, :-, '+', '-'] )
FindNumerics( S::Vector{Char} ) = findall( isnumeric.( S ) );
FindPMs( S::Vector{Char} ) = findall( ispm.( S ) )

#Returns next element symbol and returns SMILES string one element fewer.
function ReadNextElement( S::String , List::Array{ String, 1 } )
    NextElement = nothing
    if length( S ) > 1
        #See if the second char is lowercase
        if islowercase( S[ 2 ] )
            if S[ 1:2 ] in List
                NextElement = S[1:2]
                S = (length(S) > 2) ? S[3:end] : ""
            end
        else
            #Get list items 1 Char in length
            AvailableList = List[ length.( List ) .== 1 ]
            if string(S[ 1 ]) in AvailableList
                NextElement = string(S[ 1 ])
                S = S[2:end]
            end
        end
    elseif length(S) == 1
        #Get list items 1 Char in length
        AvailableList = List[ length.( List ) .== 1 ]
        if S in AvailableList
            NextElement = S
            S = ""
        end
    end
    return NextElement, S
end

ReadNextElement( "ScC", bracket )
ReadNextElement( "ScAs", bracket )
ReadNextElement( "Sc", bracket )
ReadNextElement( "CO", bracket )
ReadNextElement( "C", bracket )

#Fn for parsing the inside of a bracket...
function ReadNextNumeric(S::String)
    len = length(S)
    Sc = Vector{Char}(S)
    isotope = nothing
    if isdigit( Sc[ 1 ] )
        offset = findfirst( isdigit.( Sc ) .== false )
        if isa(offset, Nothing)
            isotope = parse( Int16, S[ 1:end ] )
            S = ""
        else
            isotope = parse( Int16, S[ 1 : ( offset - 1 ) ] )
            S = S[ offset : end ]
        end
    end
    return isotope, S
end

ReadNextNumeric("123UrC")
ReadNextNumeric("UrC")
ReadNextNumeric("12CH4")


function ReadNextCharge(S::String)
    len = length(S)
    charge = nothing
    if len == 1
        charge = Int( "+" == S ) - Int( "-" == S )
    else
        Sc = Vector{Char}( S )
        if all( [ isdigit( s ) for s in S[ 2:end ] ] )
            charge = parse( Int16, S[ 1:end ] )
        elseif all( [ ispm( s ) for s in S ] )
            charge = sum( '+' .== Sc ) - sum( '-' .== Sc )
        else
            @warn("Invalid SMILES. Charge representation should be final attribute in a bracket. Or unacceptable mixed charge notation.")
        end
    end
    return charge, ""
end

@enum BracketState begin
    ISOTOPE = Int8(0)
    SYMBOL = Int8(1)
    MULTIPLICITY = Int8(2)
    CHARGE = Int8(3)
    FINISH = Int8(4)
end

struct BracketParseException <: Exception
    bracketcontents::String
    charloc::Int64
end
Base.showerror(io::IO, e::BracketParseException) = print(io, "Failed to parse bracket(\"[$(e.bracketcontents)]\") on character ", e.charloc, "!")

mutable struct Element
    symbol::String
    isotope::Union{Int16, Nothing}
    aromatic::Bool
    ringID::Array{Int16,1}
    explicithydrogens::Int8
    charge::Int8
end

function Element(symbol::String)
    if islowercase(symbol[1])
        return Element(uppercasefirst(symbol), nothing, true, [], 0, 0 )
    else
        return Element(symbol, nothing, false, [], 0, 0 )
    end
end

function ParseBracket(S)
    Sorig = deepcopy(S)
    origlen = length(S)
    curlen = deepcopy(origlen)
    lastlen = deepcopy(origlen)

    state = ISOTOPE
    symbol, isotope, aromatic, ringID, hydrogens, charge = nothing, nothing, false, [], 0, 0

    while curlen > 0
        cursor = S[1]
        #Handle isotope, explicit hydrogen multiplicity
        if isdigit( cursor )
            if state == ISOTOPE
                isotope, S = ReadNextNumeric( S )
                state = SYMBOL
            elseif state == MULTIPLICITY
                hydrogens, S = ReadNextNumeric( S )
                state = CHARGE
            else
                @warn("Invalid SMILES. Charge must be denoted with +/-.")
            end
        #Handle Element symbol parsing
        elseif isletter( cursor )
            if ( state == SYMBOL ) || (state == ISOTOPE)
                if isa( symbol, Nothing )
                    symbollist = isuppercase( cursor ) ? bracket : bracket_aromatic
                    aromatic = islowercase( cursor )
                    symbol, S = ReadNextElement( S, symbollist )
                    if isa(symbol, Nothing)
                        @warn("Invalid SMILES. Unrecognized chemical symbol.")
                    end
                elseif S[1] == 'H'
                    state = MULTIPLICITY
                    hydrogens = 1
                    S = (length(S) > 1) ? S[2:end] : ""
                else
                    @warn("Invalid SMILES. Bracket should contain one element (with or without Hydrogen).")
                end
            else
                @warn("Invalid SMILES. Bracket should contain one element.")
            end
        #Handle charge
        elseif ispm(cursor)
            if ( state == CHARGE ) || ( state == MULTIPLICITY ) || ( state == SYMBOL )
                charge, S = ReadNextCharge( S )
                state = SYMBOL
            else
                @warn("Invalid SMILES. Charge should follow element symbol.")
            end
        #Handle operators
        elseif isspecialoperator(cursor)
            @warn("Invalid SMILES. Special operators(" * join(specialoperators, ",") * ") should not be contained in a bracket.")
        elseif isbondoperator(cursor)
            @warn("Invalid SMILES. Bond operators(" * join(bondoperators, ",") * ") should not be contained in a bracket.")
        end
        #Ensure we don't get stuck in an infinite loop
        curlen = length( S )
        if lastlen == curlen
            throw(BracketParseException( Sorig, (origlen - curlen + 1) ))
        end
        lastlen = curlen
    end
    return Element( symbol, isotope, aromatic, ringID, hydrogens, charge )
end

ParseBracket("22NaH")

struct SMILESParseException <: Exception
    charloc::Int64
end
Base.showerror(io::IO, e::BracketParseException) = print(io, "Failed to parse bracket on character ", e.charloc, "!")

using SimpleWeightedGraphs
using LightGraphs
using GraphPlot
using Compose
#Parse SMILES - whew here goes...
MoleculeGraph = SimpleGraph()#SimpleWeightedGraph();
MolecularData = Element[]

S = "C1CCCCC1C1CCCCC1"
#"C[CH4+](OC(OCC)CC)CC"#"C12(C(CCC)CCCC1)CCCCC2"#"C[CH4+](OC(OCC)CC)CC"
Sorig = deepcopy(S)
origlen = length(S)
curlen = deepcopy(origlen)
lastlen = deepcopy(origlen)

chainstack = Int16[]

nextedgestart = 0
weight = 1
symbol, isotope, aromatic, ringID, hydrogens, charge = nothing, nothing, false, nothing, 0, 0

while curlen > 0
    moiety = nothing
    cursor = S[1]

    if isdigit( cursor )  #Handle rings
        push!(MolecularData[end].ringID, parse(Int16, cursor))
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
            push!( chainstack, length( MolecularData ) + 1 )
            S = S[ 2 : end ]
        end
        if cursor == ')'
            S = S[ 2 : end ]
            nextedgestart = pop!( chainstack )
        end
    elseif isbondoperator(cursor) #Handle bonds
        weight = bonds[ string(cursor) ]
        S = S[ 2 : end ]
    end
    #New atom/moiety was parsed
    if isa( moiety, Element )
        add_vertex!( MoleculeGraph )
        push!( MolecularData, moiety )
        len = length( MolecularData )
        if (len > 1)
            if nextedgestart == 0
                add_edge!(MoleculeGraph, len - 1, len)#, weight )
            else
                add_edge!(MoleculeGraph, nextedgestart-1, len)#, weight )
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

RingClosures = Dict()
for ( i, molecule) in enumerate(MolecularData)
    if length(molecule.ringID) > 0
        for id in molecule.ringID
            if id in keys( RingClosures )
                push!( RingClosures[ id ], i)
            else
                RingClosures[ id ] = [i]
            end
        end
    end
end

RingClosures
for (k, v) in RingClosures
    for id in 2:length(v)
        add_edge!(MoleculeGraph, v[1], v[id])
    end
end

length(MolecularData)

#Close rings



Sorig


gplot( MoleculeGraph )
println("Wuh")
