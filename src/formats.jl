# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# Protoctist.Formats — the file formats protist pipelines actually exchange.
#
# Recovered API (design thread, 2026-09-16):
#   read_ecsv(path)::DataFrame        write_ecsv(df, path)
#   read_jplace(path)::PlacementBatch write_jplace(batch, path)
#   export_itol_bundle(tree, node_annot, outdir)
#
# The design returns a DataFrame. DataFrames.jl is deliberately NOT a hard
# dependency: the core must load in a bare environment. `read_ecsv` returns an
# `ECSVTable`; `columntable(t)` turns it into a NamedTuple of vectors, which is
# the Tables.jl column-table shape every DataFrame constructor accepts, so
# `DataFrame(Protoctist.Formats.columntable(read_ecsv(p)))` gives the design's
# DataFrame without imposing the dependency on users who do not want it.

module Formats

using ..Types
using ..Tree: annotate_newick, TaxonomyTree, NodeAnnot

export ECSVTable,
    PlacementBatch,
    Placement,
    read_ecsv,
    write_ecsv,
    read_jplace,
    write_jplace,
    best_placements,
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

    # See types.jl's LineageWarrant: explicit inner constructor, no convert().
    ECSVTable(
        names::Vector{String},
        columns::Vector{Vector{String}},
        meta::Vector{String},
    ) = new(names, columns, meta)
end

Base.length(t::ECSVTable) = isempty(t.columns) ? 0 : length(t.columns[1])
# NOTE: `something(i, throw(...))` would throw unconditionally — Julia
# evaluates call arguments eagerly, so the throw runs before `something`.
function Base.getindex(t::ECSVTable, col::AbstractString)
    i = findfirst(==(String(col)), t.names)
    i === nothing && throw(KeyError(col))
    return t.columns[i]
end
Base.haskey(t::ECSVTable, col::AbstractString) =
    findfirst(==(String(col)), t.names) !== nothing
Base.keys(t::ECSVTable) = t.names

"""
    columntable(t::ECSVTable)

`t` as a `NamedTuple` of columns — the Tables.jl column-table shape, so
`DataFrame(columntable(t))` gives the DataFrame the design specifies.
"""
columntable(t::ECSVTable) = NamedTuple{Tuple(Symbol.(t.names))}(Tuple(t.columns))

"""
    read_ecsv(path; sep = ',') -> ECSVTable

Read an ECSV file. Lines before the body that begin with `#` are header
metadata and are preserved. The first non-comment line is the column header.

Quoted fields are parsed (RFC4180 rules: `""` is a literal quote inside a
quoted field), because a taxon string legitimately contains commas.
"""
function read_ecsv(path::AbstractString; sep::Char = ',')
    # Plain loop, no do-block: JET analyzes the bare code cleanly, and no
    # captured-and-reassigned local can widen the ECSVTable constructor's
    # arguments to Any (a Core.Box) the way a closure would.
    meta = String[]
    names = String[]
    rows = Vector{Vector{String}}()
    fh = open(path, "r")
    try
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
    finally
        close(fh)
    end
    ncol = length(names)
    cols = [String[] for _ = 1:ncol]
    for (ln, r) in enumerate(rows)
        length(r) == ncol ||
            throw(ArgumentError("row $ln has $(length(r)) fields, header declares $ncol"))
        for j = 1:ncol
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
                    print(buf, '"')
                    i = nxt
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

# Quote on the separator actually in use: a space-delimited (Astropy) or TSV
# file must quote a field containing a space or tab, not only a comma.
# Quoted via a buffer loop rather than `"\"" * replace(s, ...) * "\""`:
# the three-argument `*` with an AbstractString middle widens under abstract
# interpretation to unions Base's Regex concatenation cannot satisfy.
function _csvfield(s::AbstractString, sep::Char = ',')
    any(c -> c == sep || c in ('"', '\n', '\r'), s) || return String(s)
    buf = IOBuffer()
    print(buf, '"')
    for c in s
        print(buf, c)
        c == '"' && print(buf, c)
    end
    print(buf, '"')
    return String(take!(buf))
end

"""
    write_ecsv(t::ECSVTable, path; sep = ',')

Write `t` back out, header metadata first. Round-trips `read_ecsv`.
"""
function write_ecsv(t::ECSVTable, path::AbstractString; sep::Char = ',')
    open(path, "w") do fh
        for m in t.meta
            println(fh, m)
        end
        println(fh, join(_csvfield.(t.names, sep), sep))
        for i = 1:length(t)
            println(
                fh,
                join((_csvfield(t.columns[j][i], sep) for j in eachindex(t.names)), sep),
            )
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

    # See types.jl's LineageWarrant: explicit inner constructor; the String/Int
    # conversions are explicit and concretely typed, not convert()-ed.
    PlacementBatch(
        tree::AbstractString,
        placements::Vector{Placement},
        fields::Vector{String},
        version::Integer,
    ) = new(String(tree), placements, fields, Int(version))
end

"""
    read_jplace(path) -> PlacementBatch

Read a jplace (phylogenetic placement) file. jplace is JSON; to keep the
core dependency-free, the file is read by a small built-in JSON reader
rather than a JSON package. Invalid JSON, or JSON without the jplace
shape (a `"tree"` string, `"fields"` strings, placement objects with a
`"p"` array), raises an `ArgumentError` that names the file.
"""
function read_jplace(path::AbstractString)
    doc = try
        _json_parse(read(path, String))
    catch e
        e isa ArgumentError ? throw(ArgumentError("$path: $(e.msg)")) : rethrow()
    end
    doc isa Dict{String,Any} ||
        throw(ArgumentError("$path: the top level is not a JSON object, not a jplace file"))
    tree = get(doc, "tree", nothing)
    tree isa String ||
        throw(ArgumentError("$path: no \"tree\" string field, not a jplace file"))
    fields = String[
        f isa String ? f :
        throw(ArgumentError("$path: a \"fields\" entry is not a string")) for
        f in _jplace_array(path, get(doc, "fields", nothing), "fields")
    ]
    # Without edge_num a placement has no edge; do not invent edge 0.
    "edge_num" in fields ||
        throw(ArgumentError("$path: \"fields\" does not declare \"edge_num\""))
    v = get(doc, "version", 3.0)
    (v isa Float64 && isinteger(v) && -2.0^63 <= v < 2.0^63) ||
        throw(ArgumentError("$path: \"version\" is not an integer in range"))
    ver = Int(v)
    placements = Placement[]
    for (nm, vals) in _jplace_entries(path, doc, fields)
        push!(
            placements,
            Placement(
                nm,
                Int(get(vals, "edge_num", 0.0)),
                get(vals, "likelihood", 0.0),
                get(vals, "like_weight_ratio", 0.0),
                get(vals, "distal_length", 0.0),
                get(vals, "pendant_length", 0.0),
            ),
        )
    end
    return PlacementBatch(tree, placements, fields, ver)
end

"""
    write_jplace(b::PlacementBatch, path)

Write `b` as a jplace file. Field order follows `b.fields` when it is
populated, otherwise the canonical order.
"""
function write_jplace(b::PlacementBatch, path::AbstractString)
    fields =
        isempty(b.fields) ?
        ["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"] :
        b.fields
    # JSON has no infinity. Refuse before the destination is opened, so a
    # failed export does not truncate an existing file.
    for p in b.placements
        for x in (p.likelihood, p.like_weight_ratio, p.distal_length, p.pendant_length)
            isinf(x) && throw(
                ArgumentError("cannot write $x for \"$(p.name)\": JSON has no infinity"),
            )
        end
    end
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
            entries = join(
                (
                    "[" * join((_jplace_value(p, f) for f in fields), ", ") * "]" for
                    p in groups[name]
                ),
                ", ",
            )
            print(fh, "    {\"p\": [", entries, "], \"n\": [", _json_quote(name), "]}")
            println(fh, i == length(order) ? "" : ",")
        end
        println(fh, "  ]")
        println(fh, "}")
    end
    return path
end

function _jplace_value(p::Placement, f::AbstractString)
    f == "edge_num" && return string(p.edge_num)
    f == "likelihood" && return _jnum(p.likelihood)
    f == "like_weight_ratio" && return _jnum(p.like_weight_ratio)
    f == "distal_length" && return _jnum(p.distal_length)
    f == "pendant_length" && return _jnum(p.pendant_length)
    return "0"
end

# JSON has no NaN; a missing value read as NaN goes back out as null.
_jnum(x::Float64) = isnan(x) ? "null" : string(x)

# Buffer loop rather than `"\"" * replace(...) * "\""` — see _csvfield.
function _json_quote(s::AbstractString)
    buf = IOBuffer()
    print(buf, '"')
    for c in s
        if c == '\\' || c == '"'
            print(buf, '\\', c)
        elseif c == '\n'
            print(buf, "\\n")
        elseif c < ' '
            # JSON forbids a raw control character in a string.
            print(buf, "\\u", string(UInt32(c); base = 16, pad = 4))
        else
            print(buf, c)
        end
    end
    print(buf, '"')
    return String(take!(buf))
end

# --- JSON reader for jplace -------------------------------------------------
#
# jplace is JSON (Matsen et al. 2012). The core has no dependencies, so this is
# a small reader for the whole of RFC 8259 rather than a JSON package. It used
# to be a set of regular expressions over the raw text; those broke on a brace
# or a quote inside a string, on a key named in the tree, and on any object
# shape they did not expect, and some of those failures read as "no
# placements" instead of an error. A malformed document now fails loudly with
# its byte offset.
#
# Values: object → Dict{String,Any}, array → Vector{Any}, string → String,
# number → Float64, true/false → Bool, null → nothing. A duplicate key is an
# error (I-JSON, RFC 7493).

mutable struct _JReader
    s::String
    i::Int
    depth::Int
end

# RFC 8259 allows a nesting limit; past it, the error names the byte offset
# instead of overflowing the stack.
const _JMAXDEPTH = 512

_jfail(r::_JReader, msg::AbstractString) =
    throw(ArgumentError("invalid JSON at byte $(r.i): $msg"))

function _json_parse(s::String)
    r = _JReader(s, 1, 0)
    v = _jvalue(r)
    _jskip(r)
    r.i <= ncodeunits(s) && _jfail(r, "unexpected data after the top-level value")
    return v
end

function _jskip(r::_JReader)
    n = ncodeunits(r.s)
    while r.i <= n && codeunit(r.s, r.i) in (0x20, 0x09, 0x0a, 0x0d)
        r.i += 1
    end
    return nothing
end

function _jpeek(r::_JReader)
    _jskip(r)
    r.i <= ncodeunits(r.s) || _jfail(r, "unexpected end of input")
    return codeunit(r.s, r.i)
end

function _jexpect(r::_JReader, c::UInt8)
    _jpeek(r) == c || _jfail(r, "expected '$(Char(c))'")
    r.i += 1
    return nothing
end

function _jvalue(r::_JReader)
    c = _jpeek(r)
    if c == UInt8('{') || c == UInt8('[')
        r.depth += 1
        r.depth > _JMAXDEPTH && _jfail(r, "nested deeper than $_JMAXDEPTH levels")
        v = c == UInt8('{') ? _jobject(r) : _jarray(r)
        r.depth -= 1
        return v
    end
    c == UInt8('"') && return _jstring(r)
    (c == UInt8('-') || UInt8('0') <= c <= UInt8('9')) && return _jnumber(r)
    _jliteral(r, "true") && return true
    _jliteral(r, "false") && return false
    _jliteral(r, "null") && return nothing
    _jfail(r, "unexpected character")
end

function _jliteral(r::_JReader, lit::String)
    n = ncodeunits(lit)
    r.i + n - 1 <= ncodeunits(r.s) || return false
    for k = 1:n
        codeunit(r.s, r.i + k - 1) == codeunit(lit, k) || return false
    end
    r.i += n
    return true
end

function _jobject(r::_JReader)
    _jexpect(r, UInt8('{'))
    d = Dict{String,Any}()
    _jpeek(r) == UInt8('}') && (r.i += 1; return d)
    while true
        _jpeek(r) == UInt8('"') || _jfail(r, "expected a string key")
        at = r.i
        k = _jstring(r)
        haskey(d, k) && (r.i = at; _jfail(r, "duplicate key \"$k\""))
        _jexpect(r, UInt8(':'))
        d[k] = _jvalue(r)
        c = _jpeek(r)
        r.i += 1
        c == UInt8('}') && return d
        c == UInt8(',') || (r.i -= 1; _jfail(r, "expected ',' or '}'"))
    end
end

function _jarray(r::_JReader)
    _jexpect(r, UInt8('['))
    a = Any[]
    _jpeek(r) == UInt8(']') && (r.i += 1; return a)
    while true
        push!(a, _jvalue(r))
        c = _jpeek(r)
        r.i += 1
        c == UInt8(']') && return a
        c == UInt8(',') || (r.i -= 1; _jfail(r, "expected ',' or ']'"))
    end
end

function _jhex4(r::_JReader)
    r.i + 3 <= ncodeunits(r.s) || _jfail(r, "truncated \\u escape")
    v = 0
    for k = 0:3
        c = codeunit(r.s, r.i + k)
        d =
            UInt8('0') <= c <= UInt8('9') ? c - UInt8('0') :
            UInt8('a') <= c <= UInt8('f') ? c - UInt8('a') + 10 :
            UInt8('A') <= c <= UInt8('F') ? c - UInt8('A') + 10 :
            _jfail(r, "bad hex digit in \\u escape")
        v = 16v + Int(d)
    end
    r.i += 4
    return v
end

const _JESC = Dict{UInt8,Char}(
    UInt8('"') => '"',
    UInt8('\\') => '\\',
    UInt8('/') => '/',
    UInt8('b') => '\b',
    UInt8('f') => '\f',
    UInt8('n') => '\n',
    UInt8('r') => '\r',
    UInt8('t') => '\t',
)

function _jstring(r::_JReader)
    _jexpect(r, UInt8('"'))
    buf = IOBuffer()
    n = ncodeunits(r.s)
    while true
        r.i <= n || _jfail(r, "unterminated string")
        c = codeunit(r.s, r.i)
        if c == UInt8('"')
            r.i += 1
            s = String(take!(buf))
            isvalid(s) || _jfail(r, "string is not valid UTF-8")
            return s
        elseif c == UInt8('\\')
            r.i += 1
            r.i <= n || _jfail(r, "unterminated string")
            e = codeunit(r.s, r.i)
            r.i += 1
            if e == UInt8('u')
                u = _jhex4(r)
                if 0xD800 <= u <= 0xDBFF
                    _jliteral(r, "\\u") || _jfail(r, "unpaired surrogate in \\u escape")
                    lo = _jhex4(r)
                    0xDC00 <= lo <= 0xDFFF || _jfail(r, "unpaired surrogate in \\u escape")
                    u = 0x10000 + ((u - 0xD800) << 10) + (lo - 0xDC00)
                elseif 0xDC00 <= u <= 0xDFFF
                    _jfail(r, "unpaired surrogate in \\u escape")
                end
                print(buf, Char(u))
            else
                ch = get(_JESC, e, nothing)
                ch === nothing && (r.i -= 1; _jfail(r, "bad escape"))
                print(buf, ch)
            end
        elseif c < 0x20
            _jfail(r, "control character in string")
        else
            write(buf, c)
            r.i += 1
        end
    end
end

function _jnumber(r::_JReader)
    s, n, start = r.s, ncodeunits(r.s), r.i
    isdig(k) = k <= n && UInt8('0') <= codeunit(s, k) <= UInt8('9')
    digits!() = (isdig(r.i) || _jfail(r, "expected a digit"); while isdig(r.i)
        r.i += 1
    end)
    r.i <= n && codeunit(s, r.i) == UInt8('-') && (r.i += 1)
    if r.i <= n && codeunit(s, r.i) == UInt8('0')
        r.i += 1
    else
        digits!()
    end
    r.i <= n && codeunit(s, r.i) == UInt8('.') && (r.i += 1; digits!())
    if r.i <= n && codeunit(s, r.i) in (UInt8('e'), UInt8('E'))
        r.i += 1
        r.i <= n && codeunit(s, r.i) in (UInt8('+'), UInt8('-')) && (r.i += 1)
        digits!()
    end
    return parse(Float64, SubString(s, start, r.i - 1))
end

# --- jplace structure ---------------------------------------------------------

function _jplace_array(path::AbstractString, v, what::AbstractString)
    v isa Vector{Any} ||
        throw(ArgumentError("$path: \"$what\" is missing or not a JSON array"))
    return v
end

function _jplace_names(path::AbstractString, obj::Dict{String,Any})
    if haskey(obj, "n")
        n = obj["n"]
        n isa String && return [n]
        return [
            x isa String ? x :
            throw(ArgumentError("$path: a name in \"n\" is not a string")) for
            x in _jplace_array(path, n, "n")
        ]
    elseif haskey(obj, "nm")
        # "nm": [[name, multiplicity], ...]; the multiplicity is not kept.
        names = String[]
        for e in _jplace_array(path, obj["nm"], "nm")
            (
                e isa Vector{Any} &&
                length(e) == 2 &&
                e[1] isa String &&
                e[2] isa Float64 &&
                e[2] > 0
            ) || throw(ArgumentError("$path: an \"nm\" entry is not [name, multiplicity]"))
            push!(names, e[1])
        end
        return names
    end
    return String[]
end

# One entry per (name, candidate) in every placement object. One object carries
# EVERY candidate edge for its query (pplacer and EPA-ng list several, each with
# its own like_weight_ratio) and may name several identical queries: nothing is
# dropped, and no candidate is privileged by its position in the file. `null`
# (which some placers write for a missing value) becomes NaN in place, so the
# remaining columns stay aligned with `fields`.
function _jplace_entries(
    path::AbstractString,
    doc::Dict{String,Any},
    fields::Vector{String},
)
    out = Tuple{String,Dict{String,Float64}}[]
    for obj in _jplace_array(path, get(doc, "placements", nothing), "placements")
        obj isa Dict{String,Any} ||
            throw(ArgumentError("$path: a placement is not a JSON object"))
        haskey(obj, "p") || throw(ArgumentError("$path: a placement has no \"p\""))
        names = _jplace_names(path, obj)
        isempty(names) && (names = [""])
        for row in _jplace_array(path, obj["p"], "p")
            row = _jplace_array(path, row, "p")
            length(row) == length(fields) || throw(
                ArgumentError(
                    "$path: a \"p\" row has $(length(row)) values for $(length(fields)) fields",
                ),
            )
            vals = Dict{String,Float64}()
            for (f, x) in zip(fields, row)
                f == "edge_num" &&
                    !(x isa Float64 && isinteger(x) && -2.0^63 <= x < 2.0^63) &&
                    throw(ArgumentError("$path: an \"edge_num\" is not an integer"))
                vals[f] =
                    x === nothing ? NaN :
                    x isa Float64 ? x :
                    throw(ArgumentError("$path: \"$f\" in \"p\" is not a number or null"))
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
(the first listed wins an exact tie). A `NaN` ratio, read from `null`, never
wins over a number; if every candidate is `NaN`, the first listed is kept.
Queries keep the order in which they
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
        elseif p.like_weight_ratio > cur.like_weight_ratio ||
               (isnan(cur.like_weight_ratio) && !isnan(p.like_weight_ratio))
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
function export_itol_bundle(tree::TaxonomyTree, annot::NodeAnnot, outdir::AbstractString)
    mkpath(outdir)
    written = String[]

    nwk = joinpath(outdir, "tree.newick")
    write(nwk, annotate_newick(tree, annot))
    push!(written, nwk)

    labels = [n.label for n in tree.nodes]

    p = joinpath(outdir, "itol_counts.txt")
    open(p, "w") do fh
        println(
            fh,
            "DATASET_SIMPLEBAR\nSEPARATOR TAB\nDATASET_LABEL\tclade cumulus\nCOLOR\t#1f77b4\nDATA",
        )
        for i = 2:length(labels)
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
        for i = 2:length(labels)
            isempty(labels[i]) && continue
            println(fh, labels[i], "\t", annot.f[i], "\t", annot.b[i], "\t", annot.res[i])
        end
    end
    push!(written, p)

    p = joinpath(outdir, "itol_contamination.txt")
    open(p, "w") do fh
        println(
            fh,
            "DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tcontamination\nCOLOR\t#d62728\nDATA",
        )
        for i = 2:length(labels)
            (isempty(labels[i]) || annot.contam[i] == 0) && continue
            println(fh, labels[i], "\t#d62728\tcontaminated")
        end
    end
    push!(written, p)

    return written
end

end # module Formats
