const bondoperators = ['@', '-', '\\', '/', '=', '#', '\$']
const specialoperators = [ '.', ']', '[', '(', ')']
const operators = [ specialoperators;  bondoperators ]

const bonds = Dict( '-' => 1, '\\' => 1, '/' => 1, '=' => 2, '#' => 3, '\$' => 4 )
const valence = Dict( "B" => 3, "C" => 4, "N" => 3, "O" => 2, "P" => 3,
                      "S" => 2, "F" => 1, "Cl" => 1, "Br" => 1, "I" => 1,
                      "H" => 1, "Po" => 4, "As" => 5, "Co" => 6)

const aliphatics  = ["B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
const aromatics   = ["b", "c", "n", "o", "s", "p"]

const bracket = [ "H", "He" ,"Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
                  "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
                  "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
                  "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                  "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Hf",
                  "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
                  "Po", "At", "Rn", "Fr", "Ra", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
                  "Ds", "Rg", "Cn", "Fl", "Lv", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
                  "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th",
                  "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
                  "No", "Lr", "Mg"]

const bracket_aromatic = ["b", "c", "n", "o", "p", "s", "se", "as"]

isoperator(x::AbstractChar) = x ∈ operators
isbondoperator(x::AbstractChar) = x ∈ bondoperators
isspecialoperator(x::AbstractChar) = x ∈ specialoperators
ispm(x::AbstractChar) = x ∈ ('+', '-')
FindNumerics( S::Vector{Char} ) = findall( isnumeric.( S ) )
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

Element(e::Element) = Element(e.symbol, e.isotope, e.aromatic, copy(e.ringID), e.explicithydrogens, e.implicithydrogens, e.charge)

function Base.:(==)(a::Element, b::Element)
    return all( [ a.symbol == b.symbol, a.isotope == b.isotope,
                a.aromatic == b.aromatic, a.ringID == b.ringID,
                a.explicithydrogens == b.explicithydrogens,
                a.charge == b.charge])
end

Base.copy(e::Element) = Element(e)

#immutable struct to hold elemental information - less memory then "Element"
struct GraphElement
    symbol::String
    isotope::Union{Int16, Nothing}
    aromatic::Bool
    hydrogens::Int8
    charge::Int8
end

GraphElement(E::Element) = GraphElement(E.symbol, E.isotope, E.aromatic,
                                        E.hydrogens, E.charge)

abbreviation(E::Union{GraphElement,Element}) = E.symbol
implicitH(E::Element) = E.implicithydrogens
explicitH(E::Element) = E.explicithydrogens
H(E::Element) = E.implicithydrogens + E.explicithydrogens
Hydrogens(E::Element) = E.implicithydrogens + E.explicithydrogens
charge(E::Union{GraphElement,Element}) = E.charge
isotope(E::Union{GraphElement,Element}) = E.isotope
