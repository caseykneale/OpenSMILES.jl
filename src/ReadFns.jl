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
