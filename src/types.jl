# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# Protoctist.Types — the vocabulary protist workflows share.
#
# Recovered from the design thread "Clades and Other Prompts" (2026-09-16).
# The design names EpistemicTypes.jl as the owner of the evidence vocabulary.
# Protoctist 0.1 defines these shapes itself and takes no dependency on it;
# binding them to EpistemicTypes through a package extension is planned and
# is not present in this release.

module Types

export TaxPath, EvidenceReceipt, EpiStatus, ResidualInfo,
       Standpoint, Warrant, ProjectionY,
       PR2_RANKS, SILVA_RANKS, ranks_for, rank_index, at_rank, truncate_to,
       FACTIVE, BELIEF, COLLAPSED, SANS_FIBRE

"""
    PR2_RANKS

The PR2 (Protist Ribosomal Reference) v5 rank ladder, root-first. PR2 is the
reference the protist 18S pipelines in this stack classify against.
"""
const PR2_RANKS = (:domain, :supergroup, :division, :subdivision,
                   :class, :order, :family, :genus, :species)

"""
    SILVA_RANKS

The SILVA rank ladder, root-first. Shorter than PR2 — it carries no
supergroup/subdivision — so a path is only comparable across the two
ladders at the ranks they share.
"""
const SILVA_RANKS = (:domain, :phylum, :class, :order, :family, :genus, :species)

"""
    ranks_for(ladder::Symbol)

Return the rank tuple for `:pr2` or `:silva`. Throws on an unknown ladder
rather than silently picking one — a lineage read against the wrong ladder
is a wrong lineage, not a warning.
"""
function ranks_for(ladder::Symbol)
    ladder === :pr2   && return PR2_RANKS
    ladder === :silva && return SILVA_RANKS
    throw(ArgumentError("unknown rank ladder $(repr(ladder)); expected :pr2 or :silva"))
end

"""
    EpiStatus

The epistemic standing of a per-row claim, as decided by the warrant gates.

- `FACTIVE`     — warranted at this rank; the receipt verifies.
- `BELIEF`      — supported but below the factive threshold.
- `COLLAPSED`   — the claim collapsed to a coarser rank than requested.
- `SANS_FIBRE`  — no origin witness; the claim carries no receipt at all.
"""
@enum EpiStatus FACTIVE BELIEF COLLAPSED SANS_FIBRE

"""
    Standpoint(source, ladder, version)

Where a claim is being made *from*: which reference, which rank ladder, and
which version of it. Two claims are only comparable from the same standpoint.
"""
struct Standpoint
    source::String
    ladder::Symbol
    version::String
end
Standpoint(source::AbstractString; ladder::Symbol = :pr2, version::AbstractString = "") =
    Standpoint(String(source), ladder, String(version))

"""
    Warrant(score, threshold, gates)

What entitles a claim: the classifier score, the threshold it was judged
against, and the named sample gates that passed. `gates` is ordered so the
canonical form of a warrant is stable, which is what makes a receipt
signable.
"""
struct Warrant
    score::Float64
    threshold::Float64
    gates::Vector{Symbol}
end
Warrant(score::Real, threshold::Real) = Warrant(Float64(score), Float64(threshold), Symbol[])

"""
    ProjectionY(rank, label)

The projection of a claim onto the rank ladder: the rank actually claimed
and the label claimed there. `rank` is a member of the standpoint's ladder.
"""
struct ProjectionY
    rank::Symbol
    label::String
end

"""
    EvidenceReceipt(standpoint, warrant, projection, signature)

A per-row "claim with a receipt" — the fibre binding a visible taxon at a
rank to the origin witness that entitles it. The design aliases this to
`EpistemicTypes.Receipt`; in this release it is a separate type.

`signature` is opaque here: Protoctist carries the receipts it is given and
never mints one, and it does not check signatures itself. Signing and
verification stay with the issuing library.
"""
struct EvidenceReceipt
    standpoint::Standpoint
    warrant::Warrant
    projection::ProjectionY
    signature::String
end
EvidenceReceipt(s::Standpoint, w::Warrant, y::ProjectionY) = EvidenceReceipt(s, w, y, "")

"""
    ResidualInfo(named_candidates, novel_ok, contaminated, note)

What is left over after a claim is made: the named species still in
contention, whether the novelty gates were satisfied, and whether a
contamination gate fired. "Novel" is only ever asserted after the
contamination gates pass — that ordering is the whole point of the field.
"""
struct ResidualInfo
    named_candidates::Int
    novel_ok::Bool
    contaminated::Bool
    note::String
end
ResidualInfo(n::Integer, novel_ok::Bool, contaminated::Bool) =
    ResidualInfo(Int(n), novel_ok, contaminated, "")

"""
    TaxPath(labels; ladder = :pr2)

A ranked lineage. `labels` is root-first and may be shorter than the ladder
(an unresolved lineage), never longer. Missing intermediate ranks are the
empty string, which is distinct from a rank the path does not reach at all.
"""
struct TaxPath
    labels::Vector{String}
    ladder::Symbol

    function TaxPath(labels::AbstractVector{<:AbstractString}; ladder::Symbol = :pr2)
        r = ranks_for(ladder)
        # Trailing blanks are ranks the path does not reach, not unfilled
        # ones, so they are dropped: SILVA and QIIME 2 write a fully resolved
        # lineage with a closing ';', and "a;b;" must equal "a;b".
        labs = String.(collect(labels))
        while !isempty(labs) && isempty(last(labs))
            pop!(labs)
        end
        length(labs) > length(r) && throw(ArgumentError(
            "lineage has $(length(labs)) labels but the $(ladder) ladder has $(length(r)) ranks"))
        new(labs, ladder)
    end
end

TaxPath(s::AbstractString; ladder::Symbol = :pr2, sep = ';') =
    TaxPath(strip.(split(s, sep)); ladder = ladder)

Base.length(p::TaxPath) = length(p.labels)
Base.isempty(p::TaxPath) = isempty(p.labels)
Base.:(==)(a::TaxPath, b::TaxPath) = a.ladder == b.ladder && a.labels == b.labels
Base.hash(p::TaxPath, h::UInt) = hash(p.labels, hash(p.ladder, h))
Base.show(io::IO, p::TaxPath) = print(io, "TaxPath(", join(p.labels, ";"), "; ladder=:", p.ladder, ")")

"""
    rank_index(p::TaxPath, rank::Symbol)

Position of `rank` in this path's ladder, or `nothing` if the ladder has no
such rank. Used to compare a PR2 path against a SILVA one only where the
two ladders actually agree.
"""
function rank_index(p::TaxPath, rank::Symbol)
    r = ranks_for(p.ladder)
    i = findfirst(==(rank), r)
    return i
end

"""
    at_rank(p::TaxPath, rank::Symbol)

The label `p` claims at `rank`, or `nothing` if the path does not reach it.
An empty label at a reached rank returns `nothing` too: an unfilled rank is
not a claim.
"""
function at_rank(p::TaxPath, rank::Symbol)
    i = rank_index(p, rank)
    i === nothing && return nothing
    i > length(p.labels) && return nothing
    lab = p.labels[i]
    return isempty(lab) ? nothing : lab
end

"""
    truncate_to(p::TaxPath, rank::Symbol)

`p` cut back to `rank` — the operation a collapse performs when a claim
cannot be warranted at the rank it was made.
"""
function truncate_to(p::TaxPath, rank::Symbol)
    i = rank_index(p, rank)
    i === nothing && throw(ArgumentError("rank $(repr(rank)) is not on the $(p.ladder) ladder"))
    return TaxPath(p.labels[1:min(i, length(p.labels))]; ladder = p.ladder)
end

end # module Types
