@enum BracketState begin
    ISOTOPE = Int8(0)
    SYMBOL = Int8(1)
    MULTIPLICITY = Int8(2)
    CHARGE = Int8(3)
    FINISH = Int8(4)
end

struct BracketParseException <: Exception
    msg::String
    bracketcontents::String
    charloc::Int
end
Base.showerror(io::IO, e::BracketParseException) = print(io, "Failed to parse bracket(\"[$(e.bracketcontents)]\") on character ", e.charloc, "!\n", e.msg)

function ParseBracket(s::AbstractString)
    len = length(s)

    state = ISOTOPE
    symbol, isotope, aromatic, ringID, hydrogens, charge = nothing, nothing, false, [], 0, 0

    idx = 1
    while idx <= len
        cursor = s[idx]
        #Handle isotope, explicit hydrogen multiplicity
        if isdigit(cursor)
            if state == ISOTOPE
                isotope, idx = tryparseint16(s, idx)
                state = SYMBOL
            elseif state == MULTIPLICITY
                hydrogens, idx = tryparseint16(s, idx)
                state = CHARGE
            else
               throw(BracketParseException("Invalid SMILES. Charge must be denoted with +/-.", s, idx))
            end
        #Handle Element symbol parsing
        elseif isletter(cursor)
            if state == SYMBOL || state == ISOTOPE
                if isa(symbol, Nothing)
                    symbollist = isuppercase(cursor) ? bracket : bracket_aromatic
                    aromatic = islowercase(cursor)
                    symbol, idx = tryparseelement(s, idx, symbollist)
                    if isa(symbol, Nothing)
                        throw(BracketParseException("Invalid SMILES. Unrecognized chemical symbol.", s, idx))
                    end
                elseif cursor == 'H'
                    state = MULTIPLICITY
                    hydrogens = 1
                    idx = nextind(s, idx)
                else
                    throw(BracketParseException("Invalid SMILES. Bracket should contain one element (with or without Hydrogen).", s, idx))
                end
            else
                throw(BracketParseException("Invalid SMILES. Bracket should contain one element.", s, idx))
            end
        #Handle charge
        elseif ispm(cursor)
            if state == CHARGE || state == MULTIPLICITY || state == SYMBOL || state == ISOTOPE
                charge, idx = tryparsecharge(s, idx)
                state = SYMBOL
            else
                throw(BracketParseException("Invalid SMILES. Charge should follow element symbol.", s, idx))
            end
        #Handle operators
        elseif isspecialoperator(cursor)
            throw(BracketParseException("Invalid SMILES. Special operators (" * join(specialoperators, ",") * ") should not be contained in a bracket.", s, idx))
        elseif isbondoperator(cursor)
            if cursor == '@'
                @warn("Chirality is not yet supported, ignoring @ operator", maxlog=1)
                _, idx = tryparsechirality(s, idx)
            else
                throw(BracketParseException("Invalid SMILES. Bond operators(" * join(setdiff(bondoperators, '@'), ",") * ") should not be contained in a bracket.", s, idx))
            end
        end
    end
    return Element( symbol, isotope, aromatic, ringID, hydrogens, 0, charge )
end
