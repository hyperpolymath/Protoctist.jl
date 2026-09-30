# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# Protoctist.IO — the file formats protist pipelines actually exchange.
#
# Recovered API (design thread, 2026-09-16):
#   read_ecsv(path)::DataFrame        write_ecsv(df, path)
#   read_jplace(path)::PlacementBatch write_jplace(batch, path)
#   export_itol_bundle(tree, node_annot, outdir)
#
# The design returns a DataFrame. DataFrames.jl is deliberately NOT a hard
# dependency: the core must load in a bare environment. `read_ecsv` returns a
# column table (a NamedTuple of vectors) which is the Tables.jl shape every
# DataFrame constructor accepts — `DataFrame(read_ecsv(p))` is the design's
# signature, one call away, with no dependency imposed on users who do not
# want it.

module IO

using ..Types
using ..Tree: annotate_newick

export ECSVTable, PlacementBatch, Placement,
       read_ecsv, write_ecsv, read_jplace, write_jplace, best_placements,
       export_itol_bundle

"""
    ECSVTable(names, columns, meta)

An ECSV (Enhanced CSV) table: a CSV body under a YAML-ish comment header
carrying column metadata. `meta` keeps the header lines verbatim so a
round-trip does not quietly discard provenance the pipeline wrote.
"""
struct ECSVTable
    names::Vector{String}
    columns::Vector{Vector{String}}
    meta::Vector{String}
end

Base.length(t::ECSVTable) = isempty(t.columns) ? 0 : length(t.columns[1])
# NOTE: `something(i, throw(...))` would throw unconditionally — Julia
# evaluates call arguments eagerly, so the throw runs before `something`.
function Base.getindex(t::ECSVTable, col::AbstractString)
    i = findfirst(==(String(col)), t.names)
    i === nothing && throw(KeyError(col))
    return t.columns[i]
end
Base.haskey(t::ECSVTable, col::AbstractString) = findfirst(==(String(col)), t.names) !== nothing
Base.keys(t::ECSVTable) = t.names

"""
    columntable(t::ECSVTable)

`t` as a `NamedTuple` of columns — the Tables.jl column-table shape, so
`DataFrame(columntable(t))` gives the DataFrame the design specifies.
"""
columntable(t::ECSVTable) =
    NamedTuple{Tuple(Symbol.(t.names))}(Tuple(t.columns))

"""
    read_ecsv(path; sep = ',') -> ECSVTable

Read an ECSV file. Lines before the body that begin with `#` are header
metadata and are preserved. The first non-comment line is the column header.

Quoted fields are parsed (RFC4180 rules: `""` is a literal quote inside a
quoted field), because a taxon string legitimately contains commas.
"""
function read_ecsv(path::AbstractString; sep::Char = ',')
    meta = String[]
    names = String[]
    rows = Vector{Vector{String}}()
    open(path, "r") do fh
        for line in eachline(fh)
            if startswith(line, "#")
                push!(meta, line)
            elseif isempty(strip(line))
                continue
            elseif isempty(names)
                names = _splitcsv(line, sep)
            else
                push!(rows, _splitcsv(line, sep))
            end
        end
    end
    ncol = length(names)
    cols = [String[] for _ in 1:ncol]
    for (ln, r) in enumerate(rows)
        length(r) == ncol || throw(ArgumentError(
            "row $ln has $(length(r)) fields, header declares $ncol"))
        for j in 1:ncol
            push!(cols[j], r[j])
        end
    end
    return ECSVTable(names, cols, meta)
end

function _splitcsv(line::AbstractString, sep::Char)
    out = String[]
    buf = IOBuffer()
    inq = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if inq
            if c == '"'
                nxt = nextind(line, i)
                if nxt <= lastindex(line) && line[nxt] == '"'
                    print(buf, '"'); i = nxt
                else
                    inq = false
                end
            else
                print(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == sep
            push!(out, String(take!(buf)))
        else
            print(buf, c)
        end
        i = nextind(line, i)
    end
    push!(out, String(take!(buf)))
    return out
end

_csvfield(s::AbstractString) =
    any(c -> c in (',', '"', '\n', '\r'), s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : String(s)

"""
    write_ecsv(t::ECSVTable, path; sep = ',')

Write `t` back out, header metadata first. Round-trips `read_ecsv`.
"""
function write_ecsv(t::ECSVTable, path::AbstractString; sep::Char = ',')
    open(path, "w") do fh
        for m in t.meta
            println(fh, m)
        end
        println(fh, join(_csvfield.(t.names), sep))
        for i in 1:length(t)
            println(fh, join((_csvfield(t.columns[j][i]) for j in eachindex(t.names)), sep))
        end
    end
    return path
end

"""
    Placement(name, edge_num, likelihood, like_weight_ratio, distal_length, pendant_length)

One placement of one query sequence onto one edge of a reference tree.
"""
struct Placement
    name::String
    edge_num::Int
    likelihood::Float64
    like_weight_ratio::Float64
    distal_length::Float64
    pendant_length::Float64
end

"""
    PlacementBatch(tree, placements, fields, version)

A parsed jplace file: the reference tree in Newick with edge numbers, the
placements, and the `fields` order the file declared.
"""
struct PlacementBatch
    tree::String
    placements::Vector{Placement}
    fields::Vector{String}
    version::Int
end

"""
    read_jplace(path) -> PlacementBatch

Read a jplace (phylogenetic placement) file. jplace is JSON; to keep the
core dependency-free this reads the structure with a small focused scanner
rather than pulling in a JSON package. It handles the jplace shape — it is
not a general JSON parser, and says so if the file is not jplace.
"""
function read_jplace(path::AbstractString)
    s = read(path, String)
    tree = _json_string_field(s, "tree")
    tree === nothing && throw(ArgumentError("$path: no \"tree\" field — not a jplace file"))
    fields = _json_string_array(s, "fields")
    ver = something(_json_int_field(s, "version"), 3)
    placements = Placement[]
    for (nm, vals) in _jplace_entries(s, fields)
        push!(placements, Placement(nm,
            Int(get(vals, "edge_num", 0.0)),
            get(vals, "likelihood", 0.0),
            get(vals, "like_weight_ratio", 0.0),
            get(vals, "distal_length", 0.0),
            get(vals, "pendant_length", 0.0)))
    end
    return PlacementBatch(tree, placements, fields, ver)
end

"""
    write_jplace(b::PlacementBatch, path)

Write `b` as a jplace file. Field order follows `b.fields` when it is
populated, otherwise the canonical order.
"""
function write_jplace(b::PlacementBatch, path::AbstractString)
    fields = isempty(b.fields) ?
        ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"] :
        b.fields
    open(path, "w") do fh
        println(fh, "{")
        println(fh, "  \"tree\": ", _json_quote(b.tree), ",")
        println(fh, "  \"version\": ", b.version, ",")
        println(fh, "  \"fields\": [", join(_json_quote.(fields), ", "), "],")
        println(fh, "  \"placements\": [")
        # One object per query, holding all of its candidates, as placers write it.
        order = String[]
        groups = Dict{String,Vector{Placement}}()
        for p in b.placements
            haskey(groups, p.name) || (push!(order, p.name); groups[p.name] = Placement[])
            push!(groups[p.name], p)
        end
        for (i, name) in enumerate(order)
            entries = join(("[" * join((_jplace_value(p, f) for f in fields), ", ") * "]"
                            for p in groups[name]), ", ")
            print(fh, "    {\"p\": [", entries, "], \"n\": [", _json_quote(name), "]}")
            println(fh, i == length(order) ? "" : ",")
        end
        println(fh, "  ]")
        println(fh, "}")
    end
    return path
end

function _jplace_value(p::Placement, f::AbstractString)
    f == "edge_num"          && return string(p.edge_num)
    f == "likelihood"        && return _jnum(p.likelihood)
    f == "like_weight_ratio" && return _jnum(p.like_weight_ratio)
    f == "distal_length"     && return _jnum(p.distal_length)
    f == "pendant_length"    && return _jnum(p.pendant_length)
    return "0"
end

# JSON has no NaN; a missing value read as NaN goes back out as null.
_jnum(x::Float64) = isnan(x) ? "null" : string(x)

_json_quote(s::AbstractString) =
    "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n") * "\""

function _json_string_field(s::AbstractString, key::AbstractString)
    m = match(Regex("\"" * key * "\"\\s*:\\s*\"((?:[^\"\\\\]|\\\\.)*)\""), s)
    m === nothing && return nothing
    return replace(m.captures[1], "\\\"" => "\"", "\\\\" => "\\", "\\n" => "\n")
end

function _json_int_field(s::AbstractString, key::AbstractString)
    m = match(Regex("\"" * key * "\"\\s*:\\s*(-?\\d+)"), s)
    m === nothing && return nothing
    return parse(Int, m.captures[1])
end

function _json_string_array(s::AbstractString, key::AbstractString)
    m = match(Regex("\"" * key * "\"\\s*:\\s*\\[([^\\]]*)\\]"), s)
    m === nothing && return String[]
    return [String(x.captures[1]) for x in eachmatch(r"\"([^\"]*)\"", m.captures[1])]
end

# Scan "placements": [ { "p": [[...],...], "n": [...] or "nm": [[name,mult],...] } ]
#
# One placement object carries EVERY candidate edge for its query (pplacer and
# EPA-ng list several, each with its own like_weight_ratio) and may name
# several identical queries. The result is one entry per (name, candidate):
# nothing is dropped, and no candidate is privileged by its position in the
# file. `null` values (which some placers write for missing fields) become
# NaN in place, so the remaining columns stay aligned with `fields`.
function _jplace_entries(s::AbstractString, fields::Vector{String})
    out = Tuple{String,Dict{String,Float64}}[]
    pm = match(r"\"placements\"\s*:\s*\[", s)
    pm === nothing && return out
    for obj in eachmatch(r"\{[^{}]*\"p\"\s*:\s*\[\[(.*?)\]\][^{}]*\}"s, s[pm.offset:end])
        body = obj.match
        names = _json_string_array(body, "n")
        if isempty(names)
            m2 = match(r"\"nm\"\s*:\s*\[(.*?)\]\s*\]"s, body)
            m2 !== nothing &&
                (names = [String(x.captures[1]) for x in eachmatch(r"\[\s*\"([^\"]*)\"", m2.captures[1])])
        end
        isempty(names) && (names = [""])
        for entry in split(obj.captures[1], r"\]\s*,\s*\[")
            toks = strip.(split(entry, ','))
            vals = Dict{String,Float64}()
            for (i, f) in enumerate(fields)
                i <= length(toks) || break
                vals[f] = toks[i] == "null" ? NaN : parse(Float64, toks[i])
            end
            for name in names
                push!(out, (name, vals))
            end
        end
    end
    return out
end

"""
    best_placements(b::PlacementBatch) -> Vector{Placement}

One placement per query: the candidate with the highest `like_weight_ratio`
(the first listed wins an exact tie). Queries keep the order in which they
first appear. Use this when a downstream step needs a single edge per
query; `b.placements` itself keeps the full distribution.
"""
function best_placements(b::PlacementBatch)
    best = Dict{String,Placement}()
    order = String[]
    for p in b.placements
        cur = get(best, p.name, nothing)
        if cur === nothing
            push!(order, p.name)
            best[p.name] = p
        elseif p.like_weight_ratio > cur.like_weight_ratio
            best[p.name] = p
        end
    end
    return [best[n] for n in order]
end

"""
    export_itol_bundle(tree, annot, outdir) -> Vector{String}

Write an iTOL-ready bundle for `tree`/`annot` into `outdir`:

- `tree.newick`            — the annotated Extended Newick string
- `itol_counts.txt`        — SIMPLEBAR dataset, cumulative clade counts
- `itol_status.txt`        — MULTIBAR dataset, factive / belief / residual
- `itol_contamination.txt` — COLORSTRIP marking contaminated clades

Returns the paths written. iTOL dataset files are plain text with a small
declarative header, so this needs no dependency.
"""
function export_itol_bundle(tree, annot, outdir::AbstractString)
    mkpath(outdir)
    written = String[]

    nwk = joinpath(outdir, "tree.newick")
    write(nwk, annotate_newick(tree, annot))
    push!(written, nwk)

    labels = [n.label for n in tree.nodes]

    p = joinpath(outdir, "itol_counts.txt")
    open(p, "w") do fh
        println(fh, "DATASET_SIMPLEBAR\nSEPARATOR TAB\nDATASET_LABEL\tclade cumulus\nCOLOR\t#1f77b4\nDATA")
        for i in 2:length(labels)
            isempty(labels[i]) && continue
            println(fh, labels[i], "\t", annot.cum[i])
        end
    end
    push!(written, p)

    p = joinpath(outdir, "itol_status.txt")
    open(p, "w") do fh
        println(fh, "DATASET_MULTIBAR\nSEPARATOR TAB\nDATASET_LABEL\tepistemic status")
        println(fh, "FIELD_COLORS\t#2ca02c\t#ff7f0e\t#7f7f7f")
        println(fh, "FIELD_LABELS\tfactive\tbelief\tresidual\nDATA")
        for i in 2:length(labels)
            isempty(labels[i]) && continue
            println(fh, labels[i], "\t", annot.f[i], "\t", annot.b[i], "\t", annot.res[i])
        end
    end
    push!(written, p)

    p = joinpath(outdir, "itol_contamination.txt")
    open(p, "w") do fh
        println(fh, "DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tcontamination\nCOLOR\t#d62728\nDATA")
        for i in 2:length(labels)
            (isempty(labels[i]) || annot.contam[i] == 0) && continue
            println(fh, labels[i], "\t#d62728\tcontaminated")
        end
    end
    push!(written, p)

    return written
end

end # module IO
