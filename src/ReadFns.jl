"""
    symbol, idxnext = tryparseelement(s, idx, elements)

Extract the element symbol starting at position `idx` in `s`. It must be one of the items in `elements`.

`symbol` is the element symbol (a string), or `nothing` if no valid element can be found.
`idxnext` is the first position in `s` after the parsed element.
"""
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

"""
    i, idxnext = tryparseint16(s, idx)

Parse an `Int16` starting at position `idx` in `s`.

`i` is the `Int16` value, or `nothing` if `s` does not contain an integer starting at `idx`.
`idxnext` is the first position in `s` after the parsed integer.
"""
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

"""
    z, idxnext = tryparseint16(s, idx)

Parse charge state (as an `Int`) starting at position `idx` in `s`.

`z` is the charge, or `nothing` if `s` does not contain a charge starting at `idx`.
`idxnext` is the first position in `s` after the parsed charge.
"""
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
