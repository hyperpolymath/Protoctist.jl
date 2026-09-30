# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# Protoctist.Manifest — run manifests and the cross-run merge check.
#
# Recovered API (design thread, 2026-09-16):
#   build_manifest(run_cfg, policy, hashes)::Dict
#   validate_merge(manA, manB)  -> uses EpistemicTypes.can_merge
#
# `can_merge` belongs to EpistemicTypes. Protoctist 0.1 takes no dependency on
# it, so `validate_merge` implements the check directly. Any later deferral to
# the sibling must keep the MergeVerdict shape, so callers never see the
# answer change form.

module Manifest

using ..Types

export build_manifest, validate_merge, MergeVerdict

"""
    build_manifest(run_cfg, policy, hashes) -> Dict{String,Any}

Assemble the manifest describing one run: its configuration, the threshold
policy it was judged under, and the content hashes of its inputs.

Hashes are *passed in*, never computed here — the manifest records what the
pipeline hashed, so recomputing them would let this library disagree with
the run it is describing.
"""
function build_manifest(run_cfg, policy, hashes)
    return Dict{String,Any}(
        "schema"   => "protoctist/manifest/v1",
        "run"      => _asdict(run_cfg),
        "policy"   => _asdict(policy),
        "hashes"   => _asdict(hashes),
    )
end

_asdict(x::AbstractDict) = Dict{String,Any}(String(k) => v for (k, v) in x)
_asdict(x::NamedTuple)   = Dict{String,Any}(String(k) => getfield(x, k) for k in keys(x))
function _asdict(x)
    fs = fieldnames(typeof(x))
    isempty(fs) && return Dict{String,Any}("value" => x)
    return Dict{String,Any}(String(f) => getfield(x, f) for f in fs)
end

"""
    MergeVerdict(ok, reasons)

Whether two runs may be merged, and — when they may not — every reason, not
just the first. A caller deciding whether to re-run needs the whole list.
"""
struct MergeVerdict
    ok::Bool
    reasons::Vector{String}
end

Base.show(io::Base.IO, v::MergeVerdict) =
    print(io, v.ok ? "MergeVerdict(ok)" :
              "MergeVerdict(refused: " * join(v.reasons, "; ") * ")")

"""
    validate_merge(manA, manB) -> MergeVerdict

Decide whether two run manifests describe runs that may be pooled.

Two runs merge only when they were judged the same way: the same policy and
the same reference standpoint. Differing *inputs* are expected and fine —
that is what merging is for — so only the judgement surface is compared.
"""
function validate_merge(manA::AbstractDict, manB::AbstractDict)
    reasons = String[]

    sa = get(manA, "schema", nothing); sb = get(manB, "schema", nothing)
    sa == sb || push!(reasons, "schema mismatch: $(repr(sa)) vs $(repr(sb))")

    pa = get(manA, "policy", Dict{String,Any}())
    pb = get(manB, "policy", Dict{String,Any}())
    for k in sort(collect(union(keys(pa), keys(pb))))
        va = get(pa, k, nothing); vb = get(pb, k, nothing)
        va == vb || push!(reasons, "policy.$k differs: $(repr(va)) vs $(repr(vb))")
    end

    ra = get(manA, "run", Dict{String,Any}())
    rb = get(manB, "run", Dict{String,Any}())
    for k in ("reference", "ladder", "version")
        va = get(ra, k, nothing); vb = get(rb, k, nothing)
        (va === nothing && vb === nothing) && continue
        va == vb || push!(reasons, "run.$k differs: $(repr(va)) vs $(repr(vb))")
    end

    return MergeVerdict(isempty(reasons), reasons)
end

end # module Manifest
