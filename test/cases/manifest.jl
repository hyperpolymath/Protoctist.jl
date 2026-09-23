@testset "manifest + merge check" begin
    cfg    = (reference = "PR2", ladder = "pr2", version = "5.0.0")
    policy = (factive = 0.95, belief = 0.70)
    man_a  = build_manifest(cfg, policy, Dict("input.fasta" => "abc123"))
    man_b  = build_manifest(cfg, policy, Dict("input.fasta" => "def456"))

    @test man_a["schema"] == "protoctist/manifest/v1"
    # differing INPUTS must still merge — that is what merging is for
    v = validate_merge(man_a, man_b)
    @test v.ok
    @test isempty(v.reasons)

    # a differing POLICY must not merge, and must say why
    man_c = build_manifest(cfg, (factive = 0.99, belief = 0.70), Dict{String,String}())
    v2 = validate_merge(man_a, man_c)
    @test !v2.ok
    @test any(r -> occursin("policy.factive", r), v2.reasons)

    # a differing reference must not merge either
    man_d = build_manifest((reference = "SILVA", ladder = "silva", version = "138"),
                           policy, Dict{String,String}())
    v3 = validate_merge(man_a, man_d)
    @test !v3.ok
    @test length(v3.reasons) >= 2      # every reason, not just the first
end
