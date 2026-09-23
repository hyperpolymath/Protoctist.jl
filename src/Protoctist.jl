# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>

"""
    Protoctist

The single import point for protistologists.

Protoctist unifies the types, conversions and exports used by PR2/SILVA-driven
protist workflows. It is deliberately small: it **composes** the surrounding
libraries rather than duplicating their logic.

    using Protoctist

    run  = load_protist_run("merged.ecsv")
    cum  = clade_cumulus(run)
    export_iqtree_annotated(run, "out/")

# Design provenance

The API here was recovered from the design thread *"Clades and Other Prompts"*
(2026-09-16). That thread specifies Protoctist as the fourth library of the
BIOCEV epistemic + cladistic stack, sitting on:

- `EpistemicTypes.jl` — receipts, `epi_status`, residual summaries, merge check
- `Cladistics.jl`     — LCA, cumulative rollups, annotated Newick
- `EchoTypes.jl`      — spec tests (test-only; never on a hot path)

Those three are not yet hard dependencies. The evidence vocabulary and the
taxonomy trie are defined here so the package loads and is useful standalone;
when a sibling is present, the matching extension in `ext/` defers to it. That
keeps the honest property that `using Protoctist` works today, without
pretending the siblings are already published.

# Submodules

- [`Protoctist.Types`](@ref)    — `TaxPath`, `EvidenceReceipt`, `EpiStatus`, …
- [`Protoctist.Tree`](@ref)     — taxonomy trees, evidence rollup, Newick
- [`Protoctist.IO`](@ref)       — ECSV, jplace, iTOL bundles
- [`Protoctist.Manifest`](@ref) — run manifests and the merge check
"""
module Protoctist

# --- submodules (order matters: IO uses Tree) ----------------------------
include("types.jl")
include("tree.jl")
include("io.jl")
include("manifest.jl")

using .Types
using .Tree
using .Manifest

# --- public API ----------------------------------------------------------

export Types, Tree, IO, Manifest
export TaxPath, EvidenceReceipt, EpiStatus, ResidualInfo,
       Standpoint, Warrant, ProjectionY,
       PR2_RANKS, SILVA_RANKS, ranks_for, at_rank, truncate_to,
       FACTIVE, BELIEF, COLLAPSED, SANS_FIBRE
export TaxonomyTree, NodeAnnot, build_taxonomy_tree, rollup_counts, annotate_newick
export build_manifest, validate_merge, MergeVerdict
export ProtistRun, load_protist_run, clade_cumulus, export_iqtree_annotated

"""
    ProtistRun(table, tree, annot, ladder)

One loaded protist run: the source table, the taxonomy tree built from its
lineages, and the evidence rolled up over that tree.
"""
struct ProtistRun
    table::IO.ECSVTable
    tree::Tree.TaxonomyTree
    annot::Tree.NodeAnnot
    ladder::Symbol
end

Base.show(io::Base.IO, r::ProtistRun) = print(io,
    "ProtistRun(", length(r.table), " rows, ",
    Tree.nodecount(r.tree), " taxa, ladder=:", r.ladder, ")")

"""
    load_protist_run(path; ladder = :pr2, lineage = nothing, count = nothing) -> ProtistRun

Load an ECSV table, build its taxonomy tree and roll the evidence up — the
three steps every protist analysis starts with, in one call.

`lineage` and `count` name the columns to use. When omitted they are detected
from the usual spellings (`lineage`/`taxonomy`/`taxpath`, `count`/`abundance`/
`reads`). Detection failing for the lineage column is an error, not a guess.
"""
function load_protist_run(path::AbstractString; ladder::Symbol = :pr2,
                          lineage = nothing, count = nothing)
    t = IO.read_ecsv(path)
    lincol = lineage === nothing ? _detect(t, ("lineage", "taxonomy", "taxpath", "taxon")) : String(lineage)
    lincol === nothing && throw(ArgumentError(
        "$path: no lineage column found among $(t.names); pass `lineage=` explicitly"))
    cntcol = count === nothing ? _detect(t, ("count", "abundance", "reads", "n")) : String(count)

    lins = t[lincol]
    paths = [TaxPath(s; ladder = ladder) for s in lins]
    tree = build_taxonomy_tree(paths; ladder = ladder)

    counts = cntcol === nothing ? fill(1, length(paths)) :
             [something(tryparse(Int, s), 0) for s in t[cntcol]]
    rows = [(path = paths[i], count = counts[i]) for i in eachindex(paths)]

    return ProtistRun(t, tree, rollup_counts(tree, rows), ladder)
end

function _detect(t::IO.ECSVTable, candidates)
    for c in candidates
        i = findfirst(n -> lowercase(n) == c, t.names)
        i === nothing || return t.names[i]
    end
    return nothing
end

"""
    clade_cumulus(run::ProtistRun) -> Vector{Pair{String,Int}}

The cumulative count for every clade in `run`, deepest-first — the "clade
cumulus" the design names. Root is included; unnamed nodes are not.
"""
function clade_cumulus(run::ProtistRun)
    out = Pair{String,Int}[]
    for i in eachindex(run.tree.nodes)
        lab = run.tree.nodes[i].label
        isempty(lab) && continue
        push!(out, lab => run.annot.cum[i])
    end
    sort!(out; by = last, rev = true)
    return out
end

"""
    export_iqtree_annotated(run::ProtistRun, outdir) -> Vector{String}

Write the annotated tree and the iTOL bundle for `run` into `outdir`, the
form IQ-TREE's output is carried into iTOL with. Returns the paths written.
"""
export_iqtree_annotated(run::ProtistRun, outdir::AbstractString) =
    IO.export_itol_bundle(run.tree, run.annot, outdir)

end # module Protoctist
