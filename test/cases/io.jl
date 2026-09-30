@testset "ECSV round-trip" begin
    mktempdir() do dir
        p = joinpath(dir, "t.ecsv")
        write(p, """
        # %ECSV 1.0
        # ---
        # datatype: [{name: lineage}, {name: count}]
        lineage,count
        "Eukaryota;TSAR;Alveolata",12
        "Eukaryota;Excavata;Metamonada,odd",3
        """)
        t = Protoctist.IO.read_ecsv(p)
        @test t.names == ["lineage", "count"]
        @test length(t) == 2
        # a comma inside a quoted field must survive
        @test t["lineage"][2] == "Eukaryota;Excavata;Metamonada,odd"
        @test length(t.meta) == 3          # header preserved

        q = joinpath(dir, "out.ecsv")
        Protoctist.IO.write_ecsv(t, q)
        t2 = Protoctist.IO.read_ecsv(q)
        @test t2.names   == t.names
        @test t2.columns == t.columns
        @test t2.meta    == t.meta
    end
end

@testset "jplace round-trip" begin
    mktempdir() do dir
        p = joinpath(dir, "t.jplace")
        write(p, """
        {
          "tree": "((A:0.1{0},B:0.2{1}):0.3{2});",
          "version": 3,
          "fields": ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"],
          "placements": [
            {"p": [[1, -100.5, 0.9, 0.01, 0.02]], "n": ["query1"]},
            {"p": [[2, -110.25, 0.8, 0.03, 0.04]], "n": ["query2"]}
          ]
        }
        """)
        b = Protoctist.IO.read_jplace(p)
        @test b.version == 3
        @test occursin("{0}", b.tree)
        @test length(b.placements) == 2
        @test b.placements[1].name == "query1"
        @test b.placements[1].edge_num == 1
        @test b.placements[1].like_weight_ratio ≈ 0.9
        @test b.placements[2].name == "query2"

        q = joinpath(dir, "out.jplace")
        Protoctist.IO.write_jplace(b, q)
        b2 = Protoctist.IO.read_jplace(q)
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
        write(p, """
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
        """)
        b = Protoctist.IO.read_jplace(p)
        @test length(b.placements) == 5
        q1 = filter(x -> x.name == "query1", b.placements)
        @test [x.edge_num for x in q1] == [0, 1]
        @test [x.like_weight_ratio for x in q1] ≈ [0.35, 0.65]
        @test sum(x.like_weight_ratio for x in q1) ≈ 1.0
        @test Set(x.name for x in b.placements) == Set(["query1", "query2a", "query2b", "query3"])
        q3 = only(filter(x -> x.name == "query3", b.placements))
        @test isnan(q3.likelihood)
        @test q3.like_weight_ratio ≈ 0.5           # columns after the null stay aligned
        @test q3.pendant_length ≈ 0.1

        best = Protoctist.IO.best_placements(b)
        @test [x.name for x in best] == ["query1", "query2a", "query2b", "query3"]
        @test best[1].edge_num == 1                # highest LWR, not first listed
        @test best[1].like_weight_ratio ≈ 0.65

        q = joinpath(dir, "out.jplace")
        Protoctist.IO.write_jplace(b, q)
        out = read(q, String)
        @test !occursin("NaN", out)                # JSON has no NaN
        @test count("\"n\":", out) == 4            # one object per query
        b2 = Protoctist.IO.read_jplace(q)
        @test length(b2.placements) == length(b.placements)
        @test [(x.name, x.edge_num) for x in b2.placements] ==
              [(x.name, x.edge_num) for x in b.placements]
        @test isnan(only(filter(x -> x.name == "query3", b2.placements)).likelihood)
    end
end

@testset "iTOL bundle + run pipeline" begin
    mktempdir() do dir
        p = joinpath(dir, "run.ecsv")
        write(p, """
        # %ECSV 1.0
        lineage,count
        Eukaryota;Excavata;Metamonada;;Fornicata;Diplomonadida;Giardia,40
        Eukaryota;TSAR;Rhizaria;;Cercozoa,10
        """)
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
        write(p, "lineage,count\nEukaryota;TSAR;Rhizaria,12.0\nEukaryota;TSAR;Alveolata,3\n")
        run = load_protist_run(p)
        @test first(clade_cumulus(run)) == ("root" => 15)

        # A count that is not a whole number is an error naming the row,
        # never a zero that drops the reads.
        for bad in ("NA", "", "12.5")
            q = joinpath(dir, "bad.csv")
            write(q, "lineage,count\nEukaryota;TSAR;Rhizaria,5\nEukaryota;TSAR;Alveolata,$bad\n")
            err = try load_protist_run(q); nothing catch e; e end
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
        t = Protoctist.IO.ECSVTable(["taxonomy", "reads"],
            [["d__Eukaryota; p__Cercozoa", "a\tb"], ["4", "2"]], ["# %ECSV 1.0"])
        for sep in (' ', '\t', ',')
            p = joinpath(dir, "sep.ecsv")
            Protoctist.IO.write_ecsv(t, p; sep = sep)
            @test Protoctist.IO.read_ecsv(p; sep = sep).columns == t.columns
        end
    end
end
