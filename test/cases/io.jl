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
