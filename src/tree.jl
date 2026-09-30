# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# Protoctist.Tree — taxonomy trees, evidence rollup, annotated Newick.
#
# Recovered API (design thread, 2026-09-16):
#   build_taxonomy_tree(paths)::Cladistics.Tree
#   rollup_counts(tree, df, epi_policy)::NodeAnnot
#   annotate_newick(tree, annot)::String
#
# Cladistics.jl owns the general tree algebra. Protoctist does not duplicate
# it: this module builds the *taxonomy* tree, which is a trie over ranked
# lineages. Converting `TaxonomyTree` to Cladistics' `Tree` through a
# package extension is planned and is not present in this release.

module Tree

using ..Types

export TaxonomyTree, TaxNode, NodeAnnot,
       build_taxonomy_tree, rollup_counts, annotate_newick, nodecount, leaves

"""
    TaxNode

One node of the taxonomy trie. `rank` is the ladder rank this node sits at
(`:root` for the root), `label` the taxon name, `children` indices into the
owning tree's node vector.
"""
struct TaxNode
    label::String
    rank::Symbol
    parent::Int
    children::Vector{Int}
end

"""
    TaxonomyTree

A trie over ranked lineages. Node 1 is always the root. Because every path
shares the root, the tree is exactly the set of lineages that were observed,
with shared prefixes merged — which is what makes a cumulative rollup a
single post-order pass.
"""
struct TaxonomyTree
    nodes::Vector{TaxNode}
    ladder::Symbol
end

nodecount(t::TaxonomyTree) = length(t.nodes)
leaves(t::TaxonomyTree) = [i for i in eachindex(t.nodes) if isempty(t.nodes[i].children)]

"""
    build_taxonomy_tree(paths; ladder = :pr2)

Build the taxonomy trie from `paths` (a collection of `TaxPath`, or of
strings parsed against `ladder`). A blank label is an *unfilled* rank, not the end of the lineage: PR2
routinely leaves `subdivision` empty for lineages that do reach genus, so a
blank is skipped and the next filled rank attaches at its own true rank.
Node identity is the pair (label, rank), so skipping cannot collide two
different taxa onto one node.

All paths must share one ladder — a PR2 lineage and a SILVA lineage do not
belong in the same tree, and mixing them is rejected rather than coerced.
"""
function build_taxonomy_tree(paths; ladder::Symbol = :pr2)
    ranks = Types.ranks_for(ladder)
    nodes = [TaxNode("root", :root, 0, Int[])]

    tps = TaxPath[]
    for p in paths
        tp = p isa TaxPath ? p : TaxPath(p; ladder = ladder)
        tp.ladder === ladder || throw(ArgumentError(
            "path ladder :$(tp.ladder) does not match tree ladder :$(ladder)"))
        push!(tps, tp)
    end

    for tp in tps
        cur = 1
        for (i, lab) in enumerate(tp.labels)
            isempty(lab) && continue       # unfilled rank; the lineage continues
            j = findfirst(c -> nodes[c].label == lab && nodes[c].rank == ranks[i],
                          nodes[cur].children)
            if j === nothing
                push!(nodes, TaxNode(lab, ranks[i], cur, Int[]))
                push!(nodes[cur].children, length(nodes))
                cur = length(nodes)
            else
                cur = nodes[cur].children[j]
            end
        end
    end
    return TaxonomyTree(nodes, ladder)
end

"""
    NodeAnnot

Per-node evidence, indexed by node id. `counts` is the node's own tally;
`cum` the cumulative tally over its clade (the node plus all descendants) —
the "clade cumulus" the design names. `f`, `b`, `res`, `novel` and `contam`
are the factive / belief / residual / novel / contaminated tallies, each
cumulative over the clade for the same reason.
"""
struct NodeAnnot
    counts::Vector{Int}
    cum::Vector{Int}
    f::Vector{Int}
    b::Vector{Int}
    res::Vector{Int}
    novel::Vector{Int}
    contam::Vector{Int}
end

NodeAnnot(n::Integer) = NodeAnnot(zeros(Int, n), zeros(Int, n), zeros(Int, n),
                                  zeros(Int, n), zeros(Int, n), zeros(Int, n), zeros(Int, n))

"""
    rollup_counts(tree, rows) -> NodeAnnot

Roll per-row evidence up the taxonomy tree.

`rows` is any iterable of rows exposing `path` (a `TaxPath`), and optionally
`count`, `status` (an `EpiStatus`) and `residual` (a `ResidualInfo`). A row
whose lineage is not in `tree` is attributed to the deepest ancestor the
tree does contain, never dropped silently — losing a read is a worse error
than attributing it coarsely, and the coarse attribution is visible.

The cumulative fields are computed in one post-order pass, so this is linear
in the tree, not quadratic.
"""
function rollup_counts(tree::TaxonomyTree, rows)
    n = nodecount(tree)
    a = NodeAnnot(n)
    ranks = Types.ranks_for(tree.ladder)

    for row in rows
        tp     = _rowpath(row)
        c      = _rowcount(row)
        st     = _rowstatus(row)
        resid  = _rowresidual(row)

        cur = 1
        for (i, lab) in enumerate(tp.labels)
            isempty(lab) && continue
            j = findfirst(k -> tree.nodes[k].label == lab && tree.nodes[k].rank == ranks[i],
                          tree.nodes[cur].children)
            j === nothing && break         # deepest ancestor present wins
            cur = tree.nodes[cur].children[j]
        end

        a.counts[cur] += c
        if st === Types.FACTIVE
            a.f[cur] += c
        elseif st === Types.BELIEF
            a.b[cur] += c
        end
        if resid !== nothing
            resid.named_candidates > 0 && (a.res[cur]    += c)
            resid.novel_ok             && (a.novel[cur]  += c)
            resid.contaminated         && (a.contam[cur] += c)
        end
    end

    _accumulate!(tree, a)
    return a
end

# Post-order accumulation: a node's cumulative tally is its own plus its
# children's. Iterative, so a deep lineage cannot blow the stack.
function _accumulate!(tree::TaxonomyTree, a::NodeAnnot)
    order = Int[]
    stack = [1]
    while !isempty(stack)
        v = pop!(stack)
        push!(order, v)
        append!(stack, tree.nodes[v].children)
    end
    for v in Iterators.reverse(order)
        a.cum[v] += a.counts[v]
        p = tree.nodes[v].parent
        if p != 0
            a.cum[p]    += a.cum[v]
            a.f[p]      += a.f[v]
            a.b[p]      += a.b[v]
            a.res[p]    += a.res[v]
            a.novel[p]  += a.novel[v]
            a.contam[p] += a.contam[v]
        end
    end
    return a
end

_rowpath(r) = r isa TaxPath ? r :
              hasproperty(r, :path) ? getproperty(r, :path) :
              throw(ArgumentError("row has no `path` field"))
_rowcount(r)  = hasproperty(r, :count)    ? Int(getproperty(r, :count)) : 1
_rowstatus(r) = hasproperty(r, :status)   ? getproperty(r, :status)     : nothing
_rowresidual(r) = hasproperty(r, :residual) ? getproperty(r, :residual) : nothing

"""
    annotate_newick(tree, annot; tags = true) -> String

Render `tree` as Extended Newick with NHX-style annotation blocks
(`[&&NHX:key=value:...]`), the form iTOL and the BEAST-lineage tools read.

Labels are quoted when they contain a character Newick reserves, so a taxon
name with a space or a comma round-trips instead of corrupting the string.
"""
function annotate_newick(tree::TaxonomyTree, annot::Union{NodeAnnot,Nothing} = nothing;
                         tags::Bool = true)
    io = IOBuffer()
    _emit(io, tree, 1, annot, tags)
    print(io, ";")
    return String(take!(io))
end

function _emit(io, tree, v, annot, tags)
    kids = tree.nodes[v].children
    if !isempty(kids)
        print(io, "(")
        for (i, c) in enumerate(kids)
            i > 1 && print(io, ",")
            _emit(io, tree, c, annot, tags)
        end
        print(io, ")")
    end
    print(io, _newick_label(tree.nodes[v].label))
    if tags && annot !== nothing
        print(io, "[&&NHX",
              ":rank=",   tree.nodes[v].rank,
              ":count=",  annot.counts[v],
              ":cum=",    annot.cum[v],
              ":f=",      annot.f[v],
              ":b=",      annot.b[v],
              ":res=",    annot.res[v],
              ":novel=",  annot.novel[v],
              ":contam=", annot.contam[v],
              "]")
    end
    return io
end

const _NEWICK_RESERVED = (' ', '\t', '(', ')', '[', ']', ',', ':', ';', '\'')

function _newick_label(s::AbstractString)
    isempty(s) && return ""
    any(c -> c in _NEWICK_RESERVED, s) || return s
    return "'" * replace(s, "'" => "''") * "'"
end

end # module Tree
