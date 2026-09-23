@testset "taxonomy tree — build, blanks, rollup" begin
    paths = ["Eukaryota;TSAR;Alveolata;;Ciliophora",
             "Eukaryota;TSAR;Rhizaria;;Cercozoa",
             "Eukaryota;Excavata;Metamonada;;Fornicata;Diplomonadida;Giardia"]
    t = build_taxonomy_tree(paths)

    labs = [n.label for n in t.nodes]
    # A blank subdivision must NOT truncate the lineage: PR2 leaves it empty
    # for lineages that do reach genus, and these are the target taxa.
    @test "Cercozoa" in labs
    @test "Giardia"  in labs
    @test "Ciliophora" in labs

    # shared prefixes merge — one Eukaryota, one TSAR
    @test count(==("Eukaryota"), labs) == 1
    @test count(==("TSAR"), labs) == 1

    rows = [(path = TaxPath(p), count = 3, status = FACTIVE,
             residual = ResidualInfo(1, false, false)) for p in paths]
    a = rollup_counts(t, rows)
    @test a.cum[1] == 9                 # root carries every read
    @test a.f[1]   == 9                 # all factive
    @test a.contam[1] == 0

    # mixing ladders in one tree is refused, not coerced
    @test_throws ArgumentError build_taxonomy_tree(
        [TaxPath("Bacteria;Proteobacteria"; ladder = :silva)]; ladder = :pr2)
end

@testset "annotated Extended Newick" begin
    t = build_taxonomy_tree(["Eukaryota;TSAR;Alveolata"])
    a = rollup_counts(t, [(path = TaxPath("Eukaryota;TSAR;Alveolata"), count = 5)])
    s = annotate_newick(t, a)
    @test endswith(s, ";")
    @test occursin("&&NHX", s)
    @test occursin("cum=5", s)
    @test occursin("Alveolata", s)
    @test count(==('('), s) == count(==(')'), s)

    # a reserved character in a label must be quoted, not emitted raw
    t2 = build_taxonomy_tree([TaxPath(["Eukaryota", "group, incertae sedis"])])
    s2 = annotate_newick(t2, nothing)
    @test occursin("'group, incertae sedis'", s2)
end
