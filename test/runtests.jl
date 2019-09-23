using SMILES
using Test

@testset "ReadNextElement" begin
    @test SMILES.ReadNextElement( "ScC", SMILES.bracket ) == ("Sc", "C")
    @test SMILES.ReadNextElement( "ScAs", SMILES.bracket ) == ("Sc", "As")
    @test SMILES.ReadNextElement( "Sc", bracket ) == ("Sc", "")
    @test SMILES.ReadNextElement( "CO", bracket ) == ("C", "O")
    @test SMILES.ReadNextElement( "C", bracket ) == ("C", "")
end

@testset "ReadNextNumeric" begin
    @test SMILES.ReadNextNumeric("123UrC") == (123, "UrC")
    @test SMILES.ReadNextNumeric("UrC") == (nothing, "UrC")
    @test SMILES.ReadNextNumeric("12CH4") == (12, "CH4")
end

@testset "ReadNextCharge" begin
    @test SMILES.ReadNextCharge("+") == (1, "")
    @test SMILES.ReadNextCharge("-") == (-1, "")
    @test SMILES.ReadNextCharge("+++") == (3, "")
    @test SMILES.ReadNextCharge("---") == (-3, "")
    @test SMILES.ReadNextCharge("+2") == (2, "")
    @test SMILES.ReadNextCharge("-2") == (-2, "")
end

@testset "Brackets" begin
    @test SMILES.ParseBracket("22NaH") == SMILES.Element("Na", 22, false, Int16[], 1, 0)
    @test SMILES.ParseBracket("O-") == SMILES.Element("O", nothing, false, Int16[], 0, -1)
    @test SMILES.ParseBracket("CH4") == SMILES.Element("C", nothing, false, Int16[], 4, 0)
    @test SMILES.ParseBracket("CH3-") == SMILES.Element("C", nothing, false, Int16[], 3, -1)
    @test SMILES.ParseBracket("CH2--") == SMILES.ParseBracket("CH2-2")
    @test SMILES.ParseBracket("CH1---") == SMILES.ParseBracket("CH1-3")
    @test SMILES.ParseBracket("Fe++") == SMILES.ParseBracket("Fe+2")
    @test SMILES.ParseBracket("Fe+++") == SMILES.ParseBracket("Fe+3")
end
