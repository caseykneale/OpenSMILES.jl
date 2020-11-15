if Base.VERSION >= v"1.6.0-DEV.843"
    depwarn(msg, funcsym) = Base.depwarn(msg, funcsym; force=true)
else
    depwarn(msg, funcsym) = Base.depwarn(msg, funcsym)
end

function ReadNextElement( S::String , List::Array{ String, 1 } )
    depwarn("`ReadNextElement(s, elements)` is deprecated, please use `tryparseelement(s, idx, elements)` instead.", :ReadNextElement)
    symbol, idxnext = tryparseelement(S, 1, List)
    return symbol, S[idxnext:end]
end

function ReadNextNumeric( S::String)
    depwarn("`ReadNextNumeric(s)` is deprecated, please use `tryparseint16(s, idx)` instead.", :ReadNextNumeric)
    i, idxnext = tryparseint16(S, 1)
    return i, S[idxnext:end]
end

function ReadNextCharge( S::String)
    depwarn("`ReadNextCharge(s)` is deprecated, please use `tryparsecharge(s, idx)` instead.", :ReadNextCharge)
    c, idxnext = tryparsecharge(S, 1)
    return c, S[idxnext:end]
end
