# Protoctist — test entry point (julia-library archetype)
#
# Two suites, in this order:
#   1. Aqua package-shape gate — runs FIRST: a package that does not load
#      cleanly, whose deps are not compat-closed, or whose exports are
#      ambiguous is rejected before any behaviour test can pass.
#   2. Behaviour — every test/cases/*.jl is included, in sorted order. There is
#      one file per source module (types, tree, io, manifest); io.jl covers
#      the lineage tables, jplace and the exporters, manifest.jl the
#      manifest and merge check, examples.jl runs examples/*/run.jl. smoke.jl
#      is only the load check.

using Test
using Protoctist

@testset "Protoctist — Aqua package shape" begin
    using Aqua
    Aqua.test_all(Protoctist)
end

@testset "Protoctist — behaviour" begin
    # NOTE: no `continue`/`break` inside a `for` inside @testset — Test's
    # loop-form detection re-parses the loop body at top level, where
    # `continue` is a ParseError. Filter with a guarded call instead.
    # (Measured trap, 2026-09-19 — see the Julia testing guide.)
    cases = joinpath(@__DIR__, "cases")
    if isdir(cases)
        for path in sort(readdir(cases, join = true))
            endswith(path, ".jl") && include(path)
        end
    end
end
