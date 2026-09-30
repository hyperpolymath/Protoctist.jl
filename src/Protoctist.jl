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

# Place in the stack

Protoctist is the integration layer over three sibling libraries:

- `EpistemicTypes.jl` — receipts, `epi_status`, residual summaries, merge check
- `Cladistics.jl`     — LCA, cumulative rollups, annotated Newick
- `EchoTypes.jl`      — spec tests (test-only; never on a hot path)

Version 0.1 has no dependencies at all. The evidence vocabulary and the
taxonomy trie are defined here, so the package loads and works in a bare
environment. Binding to the siblings through package extensions is planned;
no extension ships in this release.

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
    load_protist_run(path; ladder = :pr2, lineage = nothing, count = nothing, sep = ',') -> ProtistRun

Load an ECSV/CSV/TSV table, build its taxonomy tree and roll the evidence
up — the three steps every protist analysis starts with, in one call. Pass
`sep = '\\t'` for the tab-separated tables DADA2 and QIIME 2 export.

`lineage` and `count` name the columns to use. When omitted they are detected
from the usual spellings (`lineage`/`taxonomy`/`taxpath`/`taxon`,
`count`/`abundance`/`reads`/`n`). Detection failing for the lineage column is
an error, not a guess. A count that is not a whole number — blank, `NA`,
`12.5` — is also an error naming the row: reading it as zero would drop
reads without saying so. Integral floats such as `12.0` are accepted.
"""
function load_protist_run(path::AbstractString; ladder::Symbol = :pr2,
                          lineage = nothing, count = nothing, sep::Char = ',')
    t = IO.read_ecsv(path; sep = sep)
    lincol = lineage === nothing ? _detect(t, ("lineage", "taxonomy", "taxpath", "taxon")) : String(lineage)
    lincol === nothing && throw(ArgumentError(
        "$path: no lineage column found among $(t.names); pass `lineage=` explicitly"))
    cntcol = count === nothing ? _detect(t, ("count", "abundance", "reads", "n")) : String(count)

    lins = t[lincol]
    paths = [TaxPath(s; ladder = ladder) for s in lins]
    tree = build_taxonomy_tree(paths; ladder = ladder)

    counts = cntcol === nothing ? fill(1, length(paths)) :
             [_parsecount(s, i, cntcol, path) for (i, s) in enumerate(t[cntcol])]
    rows = [(path = paths[i], count = counts[i]) for i in eachindex(paths)]

    return ProtistRun(t, tree, rollup_counts(tree, rows), ladder)
end

function _parsecount(s::AbstractString, row::Int, col::AbstractString, path::AbstractString)
    v = tryparse(Int, strip(s))
    v === nothing || return v
    f = tryparse(Float64, strip(s))
    (f !== nothing && isfinite(f) && isinteger(f)) && return Int(f)
    throw(ArgumentError("$path: row $row of column \"$col\" is $(repr(s)), " *
                        "not a whole-number count"))
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

The cumulative count for every clade in `run`, largest first — the "clade
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

Write the annotated taxonomy tree and its iTOL dataset files for `run` into
`outdir`. Returns the paths written.

The tree written is the *taxonomy* tree of `run` — ranks, no branch lengths —
not a phylogeny. The dataset files are keyed by taxon label, so they can
also be dropped onto an IQ-TREE phylogeny in iTOL wherever its tip or node
labels use the same taxon names. The name follows the recovered design.
"""
export_iqtree_annotated(run::ProtistRun, outdir::AbstractString) =
    IO.export_itol_bundle(run.tree, run.annot, outdir)

end # module Protoctist
