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
            if ( state == CHARGE ) || ( state == MULTIPLICITY ) || ( state == SYMBOL ) || ( state == ISOTOPE )
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
