using SMILES
using Test

@testset "Read Next Element" begin
    @test SMILES.ReadNextElement( "ScC", SMILES.bracket ) == ("Sc", "C")
    @test SMILES.ReadNextElement( "ScAs", SMILES.bracket ) == ("Sc", "As")
    @test SMILES.ReadNextElement( "Sc", SMILES.bracket ) == ("Sc", "")
    @test SMILES.ReadNextElement( "CO", SMILES.bracket ) == ("C", "O")
    @test SMILES.ReadNextElement( "C", SMILES.bracket ) == ("C", "")
    @test SMILES.ReadNextElement( "s", SMILES.aromatics ) == ("s", "")
end

@testset "Read Next Numeric" begin
    @test SMILES.ReadNextNumeric("123UrC") == (123, "UrC")
    @test SMILES.ReadNextNumeric("UrC") == (nothing, "UrC")
    @test SMILES.ReadNextNumeric("12CH4") == (12, "CH4")
end

@testset "Read Next Charge" begin
    @test SMILES.ReadNextCharge("+") == (1, "")
    @test SMILES.ReadNextCharge("-") == (-1, "")
    @test SMILES.ReadNextCharge("+++") == (3, "")
    @test SMILES.ReadNextCharge("---") == (-3, "")
    @test SMILES.ReadNextCharge("+2") == (2, "")
    @test SMILES.ReadNextCharge("-2") == (-2, "")
end

@testset "Brackets Known" begin
    @test SMILES.ParseBracket("22NaH") == SMILES.Element("Na", 22, false, Int16[], 1, 0, 0)
    @test SMILES.ParseBracket("O-") == SMILES.Element("O", nothing, false, Int16[], 0, 0, -1)
    @test SMILES.ParseBracket("CH4") == SMILES.Element("C", nothing, false, Int16[], 4, 0, 0)
    @test SMILES.ParseBracket("CH3-") == SMILES.Element("C", nothing, false, Int16[], 3, 0, -1)
end

@testset "Brackets Equivalences" begin
    @test SMILES.ParseBracket("CH2--") == SMILES.ParseBracket("CH2-2")
    @test SMILES.ParseBracket("CH1---") == SMILES.ParseBracket("CH1-3")
    @test SMILES.ParseBracket("Fe++") == SMILES.ParseBracket("Fe+2")
    @test SMILES.ParseBracket("Fe+++") == SMILES.ParseBracket("Fe+3")
end

@testset "ParseSMILES Checking Implicit Hydrogens on Simple Molecules" begin
    _, Data = SMILES.ParseSMILES("CCCC(C)C")
    checkhydrogens = Dict( SMILES.countitems( SMILES.H.( Data ) ) )
    @test checkhydrogens[1] == 1
    @test checkhydrogens[2] == 2
    @test checkhydrogens[3] == 3

    #Check a Ring
    _, Data = SMILES.ParseSMILES("C1CCC(C)C1")
    checkhydrogens = Dict( SMILES.countitems( SMILES.H.( Data ) ) )
    @test checkhydrogens[1] == 1
    @test checkhydrogens[2] == 4
    @test checkhydrogens[3] == 1
end

@testset "ParseSMILES Empirical Formulas of Complicated Molecules" begin
    #Anthracene
    _, Data = SMILES.ParseSMILES("C1=CC=C2C=C3C=CC=CC3=CC2=C1")
    @test SMILES.EmpiricalFormula( Data ) == "C14H10"
    #Tryptophan
    _, Data = SMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
    @test SMILES.EmpiricalFormula( Data ) == "C11H12N2O2"
    #Lysergic Acid Diethylamide
    _, Data = SMILES.ParseSMILES("CCN(CC)C(=O)C1CN(C2CC3=CNC4=CC=CC(=C34)C2=C1)C")
    @test SMILES.EmpiricalFormula( Data ) == "C20H25N3O"
    #Firefly Luciferin
    _, Data = SMILES.ParseSMILES("O=C(O)[CH]1N=C(SC1)c2sc3cc(O)ccc3n2")
    @test SMILES.EmpiricalFormula( Data ) == "C11H8N2O3S2"
end
