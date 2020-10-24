using OpenSMILES
using LightGraphs
using Test

@testset "Read Next Element" begin
    @test OpenSMILES.ReadNextElement( "ScC", OpenSMILES.bracket ) == ("Sc", "C")
    @test OpenSMILES.ReadNextElement( "ScAs", OpenSMILES.bracket ) == ("Sc", "As")
    @test OpenSMILES.ReadNextElement( "Sc", OpenSMILES.bracket ) == ("Sc", "")
    @test OpenSMILES.ReadNextElement( "CO", OpenSMILES.bracket ) == ("C", "O")
    @test OpenSMILES.ReadNextElement( "C", OpenSMILES.bracket ) == ("C", "")
    @test OpenSMILES.ReadNextElement( "s", OpenSMILES.aromatics ) == ("s", "")
    @test OpenSMILES.ReadNextElement( "Oc", [OpenSMILES.aliphatics; OpenSMILES.aromatics] ) == ("O", "c")
end

@testset "Read Next Numeric" begin
    @test OpenSMILES.ReadNextNumeric("123UrC") == (123, "UrC")
    @test OpenSMILES.ReadNextNumeric("UrC") == (nothing, "UrC")
    @test OpenSMILES.ReadNextNumeric("12CH4") == (12, "CH4")
end

@testset "Read Next Charge" begin
    @test OpenSMILES.ReadNextCharge("+") == (1, "")
    @test OpenSMILES.ReadNextCharge("-") == (-1, "")
    @test OpenSMILES.ReadNextCharge("+++") == (3, "")
    @test OpenSMILES.ReadNextCharge("---") == (-3, "")
    @test OpenSMILES.ReadNextCharge("+2") == (2, "")
    @test OpenSMILES.ReadNextCharge("-2") == (-2, "")
end

@testset "Brackets Known" begin
    @test OpenSMILES.ParseBracket("22NaH") == OpenSMILES.Element("Na", 22, false, Int16[], 1, 0, 0)
    @test OpenSMILES.ParseBracket("O-") == OpenSMILES.Element("O", nothing, false, Int16[], 0, 0, -1)
    @test OpenSMILES.ParseBracket("CH4") == OpenSMILES.Element("C", nothing, false, Int16[], 4, 0, 0)
    @test OpenSMILES.ParseBracket("CH3-") == OpenSMILES.Element("C", nothing, false, Int16[], 3, 0, -1)
end

@testset "Brackets Equivalences" begin
    @test OpenSMILES.ParseBracket("CH2--") == OpenSMILES.ParseBracket("CH2-2")
    @test OpenSMILES.ParseBracket("CH1---") == OpenSMILES.ParseBracket("CH1-3")
    @test OpenSMILES.ParseBracket("Fe++") == OpenSMILES.ParseBracket("Fe+2")
    @test OpenSMILES.ParseBracket("Fe+++") == OpenSMILES.ParseBracket("Fe+3")
end

@testset "ParseOpenSMILES Checking Implicit Hydrogens on Simple Molecules" begin
    _, Data = OpenSMILES.ParseSMILES("CCCC(C)C")
    checkhydrogens = Dict( OpenSMILES.countitems( OpenSMILES.H.( Data ) ) )
    @test checkhydrogens[1] == 1
    @test checkhydrogens[2] == 2
    @test checkhydrogens[3] == 3

    #Check a Ring
    _, Data = OpenSMILES.ParseSMILES("C1CCC(C)C1")
    checkhydrogens = Dict( OpenSMILES.countitems( OpenSMILES.H.( Data ) ) )
    @test checkhydrogens[1] == 1
    @test checkhydrogens[2] == 4
    @test checkhydrogens[3] == 1
end

@testset "Chirality" begin
    # Chirality is not fully supported, but we shouldn't error
    G1, Data = OpenSMILES.ParseSMILES("N[C@](Br)(O)C")
    @test OpenSMILES.EmpiricalFormula( Data ) == "BrC2H6NO"
    G2, Data = OpenSMILES.ParseSMILES("N[C@@](Br)(C)O")
    @test OpenSMILES.EmpiricalFormula( Data ) == "BrC2H6NO"
    @test G1 == G2
end

@testset "ParseOpenSMILES Empirical Formulas of Complicated Molecules" begin
    #Phenol
    g, Data = OpenSMILES.ParseSMILES("Oc1ccccc1")
    @test OpenSMILES.EmpiricalFormula( Data ) == "C6H6O"
    for i = 1:6
        @test has_edge(g, i+1, mod1(i-1, 6)+1)   # +1 is because O is the first atom
        @test has_edge(g, i+1, mod1(i+1, 6)+1)   #              "
    end
    #Anthracene
    _, Data = OpenSMILES.ParseSMILES("C1=CC=C2C=C3C=CC=CC3=CC2=C1")
    @test OpenSMILES.EmpiricalFormula( Data ) == "C14H10"
    #Tryptophan
    _, Data = OpenSMILES.ParseSMILES("C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N")
    @test OpenSMILES.EmpiricalFormula( Data ) == "C11H12N2O2"
    #Lysergic Acid Diethylamide
    _, Data = OpenSMILES.ParseSMILES("CCN(CC)C(=O)C1CN(C2CC3=CNC4=CC=CC(=C34)C2=C1)C")
    @test OpenSMILES.EmpiricalFormula( Data ) == "C20H25N3O"
    #Firefly Luciferin
    _, Data = OpenSMILES.ParseSMILES("O=C(O)[CH]1N=C(SC1)c2sc3cc(O)ccc3n2")
    @test OpenSMILES.EmpiricalFormula( Data ) == "C11H8N2O3S2"
end
