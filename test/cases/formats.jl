@testset "ECSV round-trip" begin
    mktempdir() do dir
        p = joinpath(dir, "t.ecsv")
        write(
            p,
            """
   # %ECSV 1.0
   # ---
   # datatype: [{name: lineage}, {name: count}]
   lineage,count
   "Eukaryota;TSAR;Alveolata",12
   "Eukaryota;Excavata;Metamonada,odd",3
   """,
        )
        t = Protoctist.Formats.read_ecsv(p)
        @test t.names == ["lineage", "count"]
        @test length(t) == 2
        # a comma inside a quoted field must survive
        @test t["lineage"][2] == "Eukaryota;Excavata;Metamonada,odd"
        @test length(t.meta) == 3          # header preserved

        q = joinpath(dir, "out.ecsv")
        Protoctist.Formats.write_ecsv(t, q)
        t2 = Protoctist.Formats.read_ecsv(q)
        @test t2.names == t.names
        @test t2.columns == t.columns
        @test t2.meta == t.meta
    end
end

@testset "jplace round-trip" begin
    mktempdir() do dir
        p = joinpath(dir, "t.jplace")
        write(
            p,
            """
   {
     "tree": "((A:0.1{0},B:0.2{1}):0.3{2});",
     "version": 3,
     "fields": ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"],
     "placements": [
       {"p": [[1, -100.5, 0.9, 0.01, 0.02]], "n": ["query1"]},
       {"p": [[2, -110.25, 0.8, 0.03, 0.04]], "n": ["query2"]}
     ]
   }
   """,
        )
        b = Protoctist.Formats.read_jplace(p)
        @test b.version == 3
        @test occursin("{0}", b.tree)
        @test length(b.placements) == 2
        @test b.placements[1].name == "query1"
        @test b.placements[1].edge_num == 1
        @test b.placements[1].like_weight_ratio ≈ 0.9
        @test b.placements[2].name == "query2"

        q = joinpath(dir, "out.jplace")
        Protoctist.Formats.write_jplace(b, q)
        b2 = Protoctist.Formats.read_jplace(q)
        @test b2.tree == b.tree
        @test length(b2.placements) == length(b.placements)
        @test b2.placements[1].name == "query1"
        @test b2.placements[2].edge_num == 2
    end
end

@testset "jplace keeps every candidate placement" begin
    mktempdir() do dir
        p = joinpath(dir, "multi.jplace")
        # query1 lists its WEAKER edge first: selection must go by
        # like_weight_ratio, never by position. query2 names two identical
        # sequences via "nm"; query3 has a null field.
        write(
            p,
            """
   {
     "tree": "((A:0.1{0},B:0.2{1}):0.3{2});",
     "version": 3,
     "fields": ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"],
     "placements": [
       {"p": [[0, -120.0, 0.35, 0.01, 0.02], [1, -119.5, 0.65, 0.02, 0.03]], "n": ["query1"]},
       {"p": [[2, -110.0, 1.0, 0.05, 0.06]], "nm": [["query2a", 3], ["query2b", 1]]},
       {"p": [[1, null, 0.5, 0.0, 0.1]], "n": ["query3"]}
     ]
   }
   """,
        )
        b = Protoctist.Formats.read_jplace(p)
        @test length(b.placements) == 5
        q1 = filter(x -> x.name == "query1", b.placements)
        @test [x.edge_num for x in q1] == [0, 1]
        @test [x.like_weight_ratio for x in q1] ≈ [0.35, 0.65]
        @test sum(x.like_weight_ratio for x in q1) ≈ 1.0
        @test Set(x.name for x in b.placements) ==
              Set(["query1", "query2a", "query2b", "query3"])
        q3 = only(filter(x -> x.name == "query3", b.placements))
        @test isnan(q3.likelihood)
        @test q3.like_weight_ratio ≈ 0.5           # columns after the null stay aligned
        @test q3.pendant_length ≈ 0.1

        best = Protoctist.Formats.best_placements(b)
        @test [x.name for x in best] == ["query1", "query2a", "query2b", "query3"]
        @test best[1].edge_num == 1                # highest LWR, not first listed
        @test best[1].like_weight_ratio ≈ 0.65

        q = joinpath(dir, "out.jplace")
        Protoctist.Formats.write_jplace(b, q)
        out = read(q, String)
        @test !occursin("NaN", out)                # JSON has no NaN
        @test count("\"n\":", out) == 4            # one object per query
        b2 = Protoctist.Formats.read_jplace(q)
        @test length(b2.placements) == length(b.placements)
        @test [(x.name, x.edge_num) for x in b2.placements] == [(x.name, x.edge_num) for x in b.placements]
        @test isnan(only(filter(x -> x.name == "query3", b2.placements)).likelihood)
    end
end

@testset "iTOL bundle + run pipeline" begin
    mktempdir() do dir
        p = joinpath(dir, "run.ecsv")
        write(
            p,
            """
   # %ECSV 1.0
   lineage,count
   Eukaryota;Excavata;Metamonada;;Fornicata;Diplomonadida;Giardia,40
   Eukaryota;TSAR;Rhizaria;;Cercozoa,10
   """,
        )
        run = load_protist_run(p)
        @test length(run.table) == 2
        @test run.ladder === :pr2

        cum = clade_cumulus(run)
        @test first(cum).first == "root"
        @test first(cum).second == 50
        @test ("Giardia" => 40) in cum

        out = joinpath(dir, "itol")
        files = export_iqtree_annotated(run, out)
        @test length(files) == 4
        @test all(isfile, files)
        @test occursin("Giardia", read(joinpath(out, "tree.newick"), String))
        @test occursin("DATASET_SIMPLEBAR", read(joinpath(out, "itol_counts.txt"), String))
    end
end

@testset "counts are never silently zeroed; TSV loads" begin
    mktempdir() do dir
        # Integral floats (as pandas/R write them) are whole counts.
        p = joinpath(dir, "float.csv")
        write(
            p,
            "lineage,count\nEukaryota;TSAR;Rhizaria,12.0\nEukaryota;TSAR;Alveolata,3\n",
        )
        run = load_protist_run(p)
        @test first(clade_cumulus(run)) == ("root" => 15)

        # A count that is not a whole number is an error naming the row,
        # never a zero that drops the reads.
        for bad in ("NA", "", "12.5")
            q = joinpath(dir, "bad.csv")
            write(
                q,
                "lineage,count\nEukaryota;TSAR;Rhizaria,5\nEukaryota;TSAR;Alveolata,$bad\n",
            )
            err = try
                load_protist_run(q)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("row 2", sprint(showerror, err))
        end

        # Tab-separated, as DADA2 / QIIME 2 export.
        t = joinpath(dir, "run.tsv")
        write(t, "Taxon\treads\nEukaryota;TSAR;Rhizaria;;Cercozoa\t7\n")
        run = load_protist_run(t; sep = '\t')
        @test ("Cercozoa" => 7) in clade_cumulus(run)
    end
end

@testset "a field containing the separator round-trips" begin
    mktempdir() do dir
        # QIIME 2 lineages contain spaces; Astropy ECSV is space-delimited.
        t = Protoctist.Formats.ECSVTable(
            ["taxonomy", "reads"],
            [["d__Eukaryota; p__Cercozoa", "a\tb"], ["4", "2"]],
            ["# %ECSV 1.0"],
        )
        for sep in (' ', '\t', ',')
            p = joinpath(dir, "sep.ecsv")
            Protoctist.Formats.write_ecsv(t, p; sep = sep)
            @test Protoctist.Formats.read_ecsv(p; sep = sep).columns == t.columns
        end
    end
end

@testset "jplace: pretty-printed placements are read" begin
    mktempdir() do dir
        p = joinpath(dir, "pretty.jplace")
        # Whitespace and newlines between the brackets of "p", as a JSON
        # pretty-printer writes them. The old scanner returned 0 placements.
        write(
            p,
            """
   {
     "tree": "((A:0.1{0},B:0.2{1}):0.3{2});",
     "version": 3,
     "fields": ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"],
     "placements": [
       {
         "p": [
           [ 0, -120.0, 0.35, 0.01, 0.02 ],
           [ 1, -119.5, 0.65, 0.02, 0.03 ]
         ],
         "n": [ "query1" ]
       },
       { "p": [ [2, -110.0, 1.0, 0.05, 0.06] ], "nm": [ [ "query2a", 3 ], [ "query2b", 1 ] ] }
     ]
   }
   """,
        )
        b = Protoctist.Formats.read_jplace(p)
        @test length(b.placements) == 4
        q1 = filter(x -> x.name == "query1", b.placements)
        @test [x.edge_num for x in q1] == [0, 1]
        @test [x.like_weight_ratio for x in q1] ≈ [0.35, 0.65]
        @test Set(x.name for x in b.placements) == Set(["query1", "query2a", "query2b"])
    end
end

@testset "jplace: a null ratio never wins best_placements" begin
    mktempdir() do dir
        p = joinpath(dir, "nan.jplace")
        write(
            p,
            """
   {
     "tree": "((A:0.1{0},B:0.2{1}):0.3{2});",
     "version": 3,
     "fields": ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"],
     "placements": [
       {"p": [[0, -120.0, null, 0.01, 0.02], [1, -119.5, 0.8, 0.02, 0.03]], "n": ["first_null"]},
       {"p": [[0, -120.0, null, 0.01, 0.02], [1, -119.5, null, 0.02, 0.03]], "n": ["all_null"]}
     ]
   }
   """,
        )
        best = Protoctist.Formats.best_placements(Protoctist.Formats.read_jplace(p))
        @test only(filter(x -> x.name == "first_null", best)).edge_num == 1
        # With no number to compare, the first listed is kept.
        @test only(filter(x -> x.name == "all_null", best)).edge_num == 0
    end
end

@testset "jplace: JSON reader, not a pattern scan" begin
    mktempdir() do dir
        hdr = """"fields": ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"]"""
        # A brace and an escaped quote inside a name, the word "placements" and
        # a fake "p" inside the tree, keys in another order, \\u escapes. The
        # old scanner split objects at any brace and searched the raw text.
        p = joinpath(dir, "tricky.jplace")
        write(
            p,
            """
   {"version": 3, $hdr,
    "tree": "((A:0.1{0},B:0.2{1}):0.3{2}); \\"placements\\": [{\\"p\\": [[9]]}]",
    "placements": [
      {"n": ["q{1}", "q\\"2", "caf\\u00e9"], "p": [[1, -1.5e2, 1.0, 0.0, 0.25]]}
    ]}
   """,
        )
        b = Protoctist.Formats.read_jplace(p)
        @test [x.name for x in b.placements] == ["q{1}", "q\"2", "café"]
        @test all(x -> x.edge_num == 1 && x.likelihood == -150.0, b.placements)
        @test occursin("\"placements\"", b.tree)

        # Malformed JSON and a placement without "p" fail loudly with the path,
        # instead of reading as zero placements.
        bad = joinpath(dir, "bad.jplace")
        write(
            bad,
            """{"tree": "(A);", $hdr, "placements": [{"p": [[0, 1.0,]], "n": ["x"]}]}""",
        )
        err = try
            Protoctist.Formats.read_jplace(bad)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin(bad, err.msg) && occursin("byte", err.msg)
        nop = joinpath(dir, "nop.jplace")
        write(nop, """{"tree": "(A);", $hdr, "placements": [{"n": ["x"]}]}""")
        @test_throws ArgumentError Protoctist.Formats.read_jplace(nop)
        dup = joinpath(dir, "dup.jplace")
        write(dup, """{"tree": "(A);", "tree": "(B);", $hdr, "placements": []}""")
        @test_throws ArgumentError Protoctist.Formats.read_jplace(dup)
        @test_throws ArgumentError Protoctist.Formats.read_jplace(
            (write(joinpath(dir, "arr.jplace"), "[1]"); joinpath(dir, "arr.jplace")),
        )

        # A control character in a name is escaped on write and read back.
        q = joinpath(dir, "ctl.jplace")
        b2 = Protoctist.Formats.PlacementBatch(
            "(A);",
            [Protoctist.Formats.Placement("a\tb", 0, -1.0, 1.0, 0.0, 0.0)],
            String[],
            3,
        )
        Protoctist.Formats.write_jplace(b2, q)
        @test only(Protoctist.Formats.read_jplace(q).placements).name == "a\tb"
    end
end

@testset "jplace: the shape is checked, nothing is invented" begin
    mktempdir() do dir
        hdr = """"tree": "(A);", "version": 3, "fields": ["edge_num", "likelihood", "like_weight_ratio"]"""
        bad(name, body) = (p = joinpath(dir, name); write(p, body); p)
        # Each of these used to give a batch with made-up or missing values.
        for body in (
            """{"tree": "(A);", "placements": []}""",                       # no fields
            """{"tree": "(A);", "fields": ["likelihood"], "placements": []}""", # no edge_num
            """{$hdr}""",                                                      # no placements
            """{$hdr, "placements": [{"p": [[1]], "n": ["q"]}]}""",            # short row
            """{$hdr, "placements": [{"p": [[1, -1.0, 0.5, 9]], "n": ["q"]}]}""", # long row
            """{$hdr, "placements": [{"p": [[1.5, -1.0, 0.5]], "n": ["q"]}]}""", # edge 1.5
            """{$hdr, "placements": [{"p": [[null, -1.0, 0.5]], "n": ["q"]}]}""", # edge null
            """{$hdr, "placements": [{"p": [[1, -1.0, 0.5]], "nm": [["q"]]}]}""", # nm no count
            """{$hdr, "placements": [{"p": [[1, -1.0, 0.5]], "nm": [["q", false]]}]}""",
            """{"tree": "(A);", "version": 9223372036854775807, "fields": ["edge_num"], "placements": []}""",
        )
            @test_throws ArgumentError Protoctist.Formats.read_jplace(bad("x.jplace", body))
        end
        # Explicit empty placements and edge 0 are still fine.
        ok = bad(
            "ok.jplace",
            """{$hdr, "placements": [{"p": [[0, -1.0, null]], "nm": [["q", 2]]}]}""",
        )
        b = Protoctist.Formats.read_jplace(ok)
        @test only(b.placements).edge_num == 0
        @test isnan(only(b.placements).like_weight_ratio)
        @test isempty(
            Protoctist.Formats.read_jplace(bad("e.jplace", """{$hdr, "placements": []}""")).placements,
        )

        # Deep nesting is a located error, not a stack overflow.
        deep = bad("deep.jplace", """{$hdr, "x": $("["^10_000)""")
        err = try
            Protoctist.Formats.read_jplace(deep)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("nested deeper", err.msg)

        # An infinity is refused before the file is touched.
        q = bad("keep.jplace", "keep")
        inf = Protoctist.Formats.PlacementBatch(
            "(A);",
            [Protoctist.Formats.Placement("q", 0, -Inf, 1.0, 0.0, 0.0)],
            String[],
            3,
        )
        @test_throws ArgumentError Protoctist.Formats.write_jplace(inf, q)
        @test read(q, String) == "keep"
    end
end
