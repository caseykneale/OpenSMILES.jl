using OpenSMILES
using Graphs
using Test

@testset "Read Next Element" begin
    @test OpenSMILES.tryparseelement( "ScC", 1, OpenSMILES.bracket ) == ("Sc", 3)
    @test OpenSMILES.tryparseelement( "ScAs", 1, OpenSMILES.bracket ) == ("Sc", 3)
    @test OpenSMILES.tryparseelement( "ScAs", 3, OpenSMILES.bracket ) == ("As", 5)
    @test OpenSMILES.tryparseelement( "Sc", 1, OpenSMILES.bracket ) == ("Sc", 3)
    @test OpenSMILES.tryparseelement( "CO", 1, OpenSMILES.bracket ) == ("C", 2)
    @test OpenSMILES.tryparseelement( "C", 1, OpenSMILES.bracket ) == ("C", 2)
    @test OpenSMILES.tryparseelement( "s", 1, OpenSMILES.aromatics ) == ("s", 2)
    @test OpenSMILES.tryparseelement( "Oc", 1, [OpenSMILES.aliphatics; OpenSMILES.aromatics] ) == ("O", 2)
end

@testset "Read Next Numeric" begin
    @test OpenSMILES.tryparseint16("123UrC", 1) == (123, 4)
    @test OpenSMILES.tryparseint16("UrC", 1) == (nothing, 1)
    @test OpenSMILES.tryparseint16("12CH4", 1) == (12, 3)
    @test OpenSMILES.tryparseint16("12CH4", 5) == (4, 6)
end

@testset "Read Next Charge" begin
    @test OpenSMILES.tryparsecharge("+", 1) == (1, 2)
    @test OpenSMILES.tryparsecharge("-", 1) == (-1, 2)
    @test OpenSMILES.tryparsecharge("+++", 1) == (3, 4)
    @test OpenSMILES.tryparsecharge("---", 1) == (-3, 4)
    @test OpenSMILES.tryparsecharge("+2", 1) == (2, 3)
    @test OpenSMILES.tryparsecharge("-2", 1) == (-2, 3)
    @test OpenSMILES.tryparsecharge("Na+Cl-", 3) == ( 1, 4)
    @test OpenSMILES.tryparsecharge("Na+Cl-", 6) == (-1, 7)
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

    # Atoms with multiple valence states (here, S)
    _, Data = OpenSMILES.ParseSMILES("OS(=O)(=O)O")
    @test OpenSMILES.EmpiricalFormula( Data ) == "H2O4S"
end

@testset "Chirality" begin
    # Chirality is not fully supported, but we shouldn't error
    # Tetrahedral
    for smiles in ("N[C@](Br)(O)C", "N[C@TH1](Br)(O)C", "N[C@@](Br)(O)C", "N[C@TH2](Br)(O)C")
        _, Data = OpenSMILES.ParseSMILES(smiles)
        @test OpenSMILES.EmpiricalFormula( Data ) == "BrC2H6NO"
    end
    # Allenal
    for smiles in ("NC(Br)=[C@]=C(O)C", "NC(Br)=[C@AL1]=C(O)C", "NC(Br)=[C@@]=C(O)C", "NC(Br)=[C@AL2]=C(O)C")
        _, Data = OpenSMILES.ParseSMILES(smiles)
        @test OpenSMILES.EmpiricalFormula( Data ) == "BrC4H6NO"
    end
    # Square-planar
    for smiles in ("F[Po@SP1](Cl)(Br)I", "F[Po@SP2](Cl)(Br)I", "F[Po@SP3](Cl)(Br)I")
        _, Data = OpenSMILES.ParseSMILES(smiles)
        @test OpenSMILES.EmpiricalFormula( Data ) == "BrClFIPo"
    end
    # Trigonal-bipyramidal
    for smiles in ("S[As@TB1](F)(Cl)(Br)N", "S[As@TB2](Br)(Cl)(F)N", "S[As@TB5](F)(N)(Cl)Br", "F[As@TB10](S)(Cl)(N)Br",
                   "F[As@TB15](Cl)(S)(Br)N", "Br[As@TB20](Cl)(S)(F)N")
        _, Data = OpenSMILES.ParseSMILES(smiles)
        @test OpenSMILES.EmpiricalFormula( Data ) == "AsBrClFH3NS"
    end
    # Octahedral
    for smiles in ("C[Co@](F)(Cl)(Br)(I)S", "F[Co@@](S)(I)(C)(Cl)Br", "S[Co@OH5](F)(I)(Cl)(C)Br",
                   "Br[Co@OH9](C)(S)(Cl)(F)I", "Br[Co@OH12](Cl)(I)(F)(S)C", "Cl[Co@OH15](C)(Br)(F)(I)S",
                   "Cl[Co@OH19](C)(I)(F)(S)Br", "I[Co@OH27](Cl)(Br)(F)(S)C")
        _, Data = OpenSMILES.ParseSMILES(smiles)
        @test OpenSMILES.EmpiricalFormula( Data ) == "BrCClCoFH4IS"
    end

    # Hydrogens as bracket atoms are equivalent to Hcount hydrogens
    strh = "N[C@H](O)C"                 # hcount
    strb = "N[C@]([H])(O)C"             # bracket
    _, Data = OpenSMILES.ParseSMILES(strh)
    @test OpenSMILES.EmpiricalFormula( Data ) == "C2H7NO"
    _, Data = OpenSMILES.ParseSMILES(strb)
    @test OpenSMILES.EmpiricalFormula( Data ) == "C2H7NO"
end

@testset "Ring identifiers" begin
    # spiro[5.5]undecane, https://pubchem.ncbi.nlm.nih.gov/#query=spiro%5B5.5%5Dundecane
    g1, elms1 = ParseSMILES("C12(CCCCC1)CCCCC2")     # single-digit ring identifiers
    g2, elms2 = ParseSMILES("C%012(CCCCC%01)CCCCC2") # two-digit ring identifiers
    @test g1 == g2 && elms1 == elms2
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

@testset "instantiate_hydrogens" begin
    g, atomnodes = OpenSMILES.ParseSMILES("O")
    gH, atomnodesH = instantiate_hydrogens(g, atomnodes)
    # Ensure the input is not modified
    @test nv(g) == 1
    @test length(atomnodes) == 1 && atomnodes[1].symbol == "O"
    @test OpenSMILES.EmpiricalFormula(atomnodes) == "H2O"
    # Check the instantiated version
    @test nv(gH) == 3
    @test length(atomnodesH) == 3
    @test atomnodesH[1].symbol == "O"
    @test atomnodesH[2].symbol == "H"
    @test atomnodesH[3].symbol == "H"
    @test has_edge(gH, 1, 2) && has_edge(gH, 1, 3) && !has_edge(gH, 2, 3)
    @test OpenSMILES.EmpiricalFormula(atomnodesH) == "H2O"

    g, atomnodes = OpenSMILES.ParseSMILES("[H]O[H]")
    gH, atomnodesH = instantiate_hydrogens(g, atomnodes)
    @test g == gH && atomnodes == atomnodesH
    @test g !== gH && atomnodes !== atomnodesH    # it's a copy, not the same object
end
