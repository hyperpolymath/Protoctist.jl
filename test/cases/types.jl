@testset "TaxPath — ladders, ranks, truncation" begin
    p = TaxPath("Eukaryota;TSAR;Alveolata;;Ciliophora")
    @test p.ladder === :pr2
    @test at_rank(p, :domain)      == "Eukaryota"
    @test at_rank(p, :division)    == "Alveolata"
    @test at_rank(p, :subdivision) === nothing   # blank is not a claim
    @test at_rank(p, :class)       == "Ciliophora"
    @test at_rank(p, :genus)       === nothing   # not reached
    @test at_rank(p, :phylum)      === nothing   # not on the PR2 ladder

    @test truncate_to(p, :division) == TaxPath("Eukaryota;TSAR;Alveolata")

    s = TaxPath("Bacteria;Proteobacteria"; ladder = :silva)
    @test at_rank(s, :phylum) == "Proteobacteria"
    @test s != TaxPath("Bacteria;Proteobacteria")  # different ladder, not equal

    @test_throws ArgumentError TaxPath("a;b"; ladder = :greengenes)
    @test_throws ArgumentError TaxPath(fill("x", 20))          # longer than ladder
    @test_throws ArgumentError truncate_to(p, :phylum)         # rank off-ladder
end

@testset "EpiStatus and residuals" begin
    @test FACTIVE !== BELIEF
    r = ResidualInfo(3, false, true)
    @test r.named_candidates == 3
    @test r.novel_ok == false
    @test r.contaminated == true
end
