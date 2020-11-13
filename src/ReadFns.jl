#Returns next element symbol and returns SMILES string one element fewer.
function ReadNextElement( S::String , List::Array{ String, 1 } )
    NextElement = nothing
    if length( S ) > 1
        #See if the second char is lowercase
        if islowercase( S[ 2 ] ) && isuppercase( S[ 1 ] )
            if S[ 1:2 ] in List   # might fail if second letter is aromatic
                return S[1:2], S[3:end]
            end
        end
        #Get list items 1 Char in length
        AvailableList = List[ length.( List ) .== 1 ]
        if string(S[ 1 ]) in AvailableList
            NextElement = string(S[ 1 ])
            S = S[2:end]
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

function tryparseelement(s::AbstractString, idx::Int, elements)
    len = length(s)
    idx > len && return nothing, idx
    c = s[idx]
    isletter(c) || return nothing, idx
    idxnext = nextind(s, idx)
    cnext = idxnext <= len ? s[idxnext] : '1'
    if islowercase(cnext)
        # Try a 2-letter element
        element = c*cnext
        element ∈ elements && return element, nextind(s, idxnext)
    end
    # Try a 1-letter element
    element = string(c)
    element ∈ elements && return element, idxnext
    return nothing, idx
end

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

function tryparseint16(s::AbstractString, idx::Int)
    len = length(s)
    lastidx = prevind(s, idx)
    while lastidx < len
        i = nextind(s, lastidx)
        isdigit(s[i]) || break
        lastidx = i
    end
    return lastidx >= idx ? (parse(Int16, SubString(s, idx, lastidx)), nextind(s, lastidx)) : (nothing, idx)
end

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

function tryparsecharge(s::AbstractString, idx::Int)
    len = length(s)
    idx > len && return nothing, idx
    c = s[idx]
    ispm(c) || return nothing, idx
    sgn = c == '+' ? 1 : -1
    idxnext = nextind(s, idx)
    count = 1
    while idxnext <= len
        cnext = s[idxnext]
        if cnext == c
            count += 1
            @warn("++ and -- are deprecated in the OpenSMILES specification, please use + or - followed by multiplicity.", maxlog=1)
            idxnext = nextind(s, idxnext)
        elseif ispm(cnext)
            error("+- and -+ are not allowed for charge in SMILES")
        elseif isdigit(cnext)
            count == 1 || error("cannot use ++ or -- in charge specification followed by multiplicity")
            break
        else
            return sgn*count, idxnext
        end
    end
    idxnext > len && return sgn*count, idxnext
    mult, idxnext = tryparseint16(s, idxnext)
    mult === nothing && return sgn, idxnext
    @assert count == 1
    return sgn*mult, idxnext
end

function tryparsechirality(s::AbstractString, idx::Int)
    m = match(r"(@TH[12]|@AL[12]|@SP[1-3]|@TB\d\d?|@OH\d\d?|@@?)", s, idx)
    (m === nothing || m.offset > idx) && return nothing, idx
    chirality = m.captures[1]::AbstractString
    return chirality, m.offset + lastindex(chirality)
end
