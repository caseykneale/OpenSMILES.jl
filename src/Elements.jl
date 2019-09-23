bondoperators = ["@", "-", "\\", "/", "=", "#", "\$"];
specialoperators = [ ".", "]", "[", "(", ")"];
operators = [ specialoperators;  bondoperators ];

bonds = Dict( "-" => 1, "\\" => 1, "/" => 1, "=" => 2, "#" => 3, "\$" => 4 );
valence = Dict( "B" => 3, "C" => 4, "N" => (3, 5), "O" => 2, "P" => (3, 5),
                "S" => (2, 4, 6), "F" => 1, "Cl" => 1, "Br" => 1, "I" => 1 );

aliphatics  = ["B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"];
aromatics   = ["b", "c", "n", "o", "s", "p"];

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

mutable struct Element
    symbol::String
    isotope::Union{Int16, Nothing}
    aromatic::Bool
    ringID::Array{Int16,1}
    explicithydrogens::Int8
    implicithydrogens::Int8
    charge::Int8
end

function Element(symbol::String)
    if islowercase(symbol[1])
        return Element(uppercasefirst(symbol), nothing, true, [], 0, 0, 0 )
    else
        return Element(symbol, nothing, false, [], 0, 0, 0 )
    end
end


function Base.:(==)(a::Element, b::Element)
    return all( [ a.symbol == b.symbol, a.isotope == b.isotope,
                a.aromatic == b.aromatic, a.ringID == b.ringID,
                a.explicithydrogens == b.explicithydrogens,
                a.charge == b.charge])
end
