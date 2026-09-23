module Epistemic

# Epistemic.jl
# -----------------------------------------------------------------------------
# PURPOSE
#   Per-row "claims with receipts" (a.k.a. fibres):
#   - Bind each visible taxon/rank to its origin witness (standpoint+warrant+projection) with a signature.
#   - Verify/emit receipts (echo:v1?k=&w=&y=&sig=...) that travel in CSV.
#   - Convert raw classifier scores + sample gates into a warranted claim at the right rank:
#       epi_status ∈ { :Factive, :Belief, :Collapsed, :SansFibre }
#   - Summarize residual evidence (named species candidates; "novel" only after contamination gates).
#   - Provide a safe "can_merge" check for cross-run analyses.
#
# PHILOSOPHY
#   This is not a stats/machine-learning module. It is a rigorous, fast,
#   per-row decision engine. Most functions are O(1) per row. Heavy
#   phylogenetics (IQ-TREE, EPA-ng), alignments (MAFFT), and denoising (DADA2)
#   live elsewhere. We annotate/verify their outputs and enforce safe claims.
#
#   Because the hot path here is trivial per row, specialized accelerators
#   (TPU/NPU/VPU) do not provide meaningful benefit. Still, we define a
#   pluggable "backend" interface to batch/parallelize verification where
#   useful, with a robust CPU backend and documented stubs for others.
# -----------------------------------------------------------------------------

using Dates
using JSON3
using SHA
using Base: @kwdef

# -- Utilities ---------------------------------------------------------------

"""
  b64url(x::Vector{UInt8})::String

Base64url without padding, suitable for compact signatures inside CSV.
"""
b64url(x) = replace(base64encode(x), ('+'=>'-','/'=>'_','='=>""))

"""
  canon(k,w,y)::String

Canonical JSON representation (UTF-8, stable key order) that is signed.
We intentionally write an ordered NamedTuple (k,w,y) so keys are stable.
Floats are written as JSON numbers; if you require fixed-precision
rounding, enforce upstream before constructing `Warrant`.
"""
canon(k, w, y) = JSON3.write((k=k, w=w, y=y); allow_inf=false)

# -- Core Types --------------------------------------------------------------

@kwdef struct Standpoint
  tool::String               # "dada2" | "vsearch" | "consensus" | ...
  tool_ver::String           # "1.28.0" | ...
  db_key::String             # "pr2" | "silva" | "unite"
  db_release::String         # "5.0.0" | "138"
  params_hash::String        # "sha256:..." (hash of the taxonomy/step config)
  normalisation::String      # "none" | "rarefy" | "rss" | "css" | ...
end

@kwdef struct Warrant
  # Classifier evidence (one or more fields may be used depending on tool):
  boot_min::Union{Nothing,Float64} = nothing  # DADA2's genus/species bootstrap (0–1)
  id::Union{Nothing,Float64}       = nothing  # vsearch %identity (0–1)
  cov::Union{Nothing,Float64}      = nothing  # vsearch coverage (0–1)
  # Sample/run gates (evidence predicate E):
  depth::Union{Nothing,Int}        = nothing  # per-feature count in sample
  neg_leak::Union{Nothing,Float64} = nothing  # leakage fraction vs negative control
  chimera::Union{Nothing,Bool}     = nothing  # flagged chimera?
end

@kwdef struct ProjectionY
  taxon::String             # "Giardia"
  rank::String              # "Genus" | "Species" | ...
  feature_id::String        # ASV/OTU identifier
  sample_id::String         # if per-sample row
end

"""
Receipt

The portable witness ("fibre") that binds text to provenance and evidence.
Do not put derived status here; keep it raw. Derived fields (epi_status,
warranted_rank, residual info) are computed at load/analysis time.
"""
@kwdef struct Receipt
  v::String            = "echo:v1"  # version tag
  k::Standpoint
  w::Warrant
  y::ProjectionY
  sig::String          # "sha256:..." over canon(k,w,y)
end

# -- Receipt API -------------------------------------------------------------

"""
  make_receipt(k::Standpoint, w::Warrant, y::ProjectionY)::Receipt

Constructs a receipt and computes its signature (sha256 over canon(k,w,y)).
"""
function make_receipt(k::Standpoint, w::Warrant, y::ProjectionY)::Receipt
  payload = canon(k, w, y)
  sig = "sha256:" * b64url(SHA.sha256(payload))
  Receipt(; k, w, y, sig)
end

"""
  verify_receipt(r::Receipt)::Bool

Verifies `r.sig` matches sha256 over canon(r.k, r.w, r.y).
If false, treat row as :SansFibre or flag in loaders.
"""
function verify_receipt(r::Receipt)::Bool
  payload = canon(r.k, r.w, r.y)
  computed = "sha256:" * b64url(SHA.sha256(payload))
  return computed == r.sig
end

"""
  encode_avec_fibre(r::Receipt)::String

URI-like compact value for a CSV column:
echo:v1?k=<json>&w=<json>&y=<json>&sig=<sha256_base64url>
Note: values are percent-encoded when embedded into CSV.
"""
function encode_avec_fibre(r::Receipt)::String
  qk = JSON3.write(r.k)
  qw = JSON3.write(r.w)
  qy = JSON3.write(r.y)
  return "echo:v1?k=$(qk)&w=$(qw)&y=$(qy)&sig=$(r.sig)"
end

"""
  parse_avec_fibre(str)::Receipt

Parses a compact `avec_fibre` string back into a Receipt. Verifies version tag.
Note: We do not verify the signature here; call `verify_receipt` explicitly.
"""
function parse_avec_fibre(str::AbstractString)::Receipt
  startswith(str, "echo:v1?") || error("Bad fibre version or missing 'echo:v1' prefix")
  qs = split(str[9:end], '&')  # after "echo:v1?"
  kv = Dict{String,String}()
  for p in qs
    i = findfirst(==('='), p); i === nothing && continue
    k = p[1:i-1]; v = p[i+1:end]
    kv[k] = v
  end
  haskey(kv, "k") || error("fibre missing k")
  haskey(kv, "w") || error("fibre missing w")
  haskey(kv, "y") || error("fibre missing y")
  haskey(kv, "sig") || error("fibre missing sig")

  k = JSON3.read(kv["k"], Standpoint)
  w = JSON3.read(kv["w"], Warrant)
  y = JSON3.read(kv["y"], ProjectionY)
  return Receipt(; k, w, y, sig=kv["sig"])
end

# -- Policy and Gates --------------------------------------------------------

"""
ThresholdPolicy

Rank thresholds for tools. You can construct this from YAML/JSON at startup.
Example defaults follow the common practice for DADA2/vsearch (tune per DB/locus).
"""
@kwdef struct ThresholdPolicy
  dada2_by_rank::Dict{String,Float64} = Dict("Species"=>0.98, "Genus"=>0.88, "Class"=>0.75, "Division"=>0.70)
  vsearch_by_rank::Dict{String,NamedTuple} = Dict(
    "Species"=>(id=0.99, cov=0.95),
    "Genus"  =>(id=0.95, cov=0.90),
    "Division"=>(id=0.85, cov=0.80)
  )
  depth_min::Int = 100            # minimal per-feature count for any positive claim
  neg_leak_max::Float64 = 0.005   # maximal allowed negative-control leakage (fraction)
end

"""
  passes_gates(w::Warrant, pol::ThresholdPolicy)::Bool

Implements evidence predicate E: sample-level gates (depth/leak/chimera).
Extend here to include replicates, blacklists, placement sanity, etc.
"""
function passes_gates(w::Warrant, pol::ThresholdPolicy)::Bool
  (w.depth !== nothing && w.depth ≥ pol.depth_min) &&
  (w.neg_leak === nothing || w.neg_leak ≤ pol.neg_leak_max) &&
  (w.chimera === nothing || w.chimera == false)
end

"""
  highest_warranted_rank(k::Standpoint, w::Warrant, pol::ThresholdPolicy)::Union{Nothing,String}

Finds the finest (lowest) rank at which the classifier meets the configured threshold.
Note: Does not check gates; combine with `passes_gates` before promoting a claim.
"""
function highest_warranted_rank(k::Standpoint, w::Warrant, pol::ThresholdPolicy)::Union{Nothing,String}
  if k.tool == "dada2" && w.boot_min !== nothing
    for r in ("Species","Genus","Class","Division")
      thr = get(pol.dada2_by_rank, r, 2.0)
      if w.boot_min ≥ thr; return r; end
    end
  elseif k.tool == "vsearch" && w.id !== nothing && w.cov !== nothing
    for r in ("Species","Genus","Division")
      thr = get(pol.vsearch_by_rank, r, nothing)
      if thr !== nothing && w.id ≥ thr.id && w.cov ≥ thr.cov; return r; end
    end
  end
  return nothing
end

# -- Epistemic Status and Residuals -----------------------------------------

"""
  epi_status(r::Receipt, pol::ThresholdPolicy)::Symbol

Returns one of:
  :Factive     - positive claim warranted at (possibly demoted) rank
  :Belief      - measurable evidence but not promoted (below thresholds)
  :Collapsed   - insufficient or failed gates; residue exists but no positive claim
  :SansFibre   - receipt cannot be verified; treat as unverifiable
"""
function epi_status(r::Receipt, pol::ThresholdPolicy)::Symbol
  if !verify_receipt(r); return :SansFibre; end
  if !passes_gates(r.w, pol); return :Collapsed; end
  wr = highest_warranted_rank(r.k, r.w, pol)
  return wr === nothing ? :Belief : :Factive
end

"""
  ZeroKind

true_absence   - adequate power/gates to support absence
undetected     - inadequate power (e.g., low depth) to conclude absence
"""
const ZeroKind = Union{Val{:true_absence},Val{:undetected}}

"""
  disambiguate_zero(w::Warrant, pol::ThresholdPolicy)::ZeroKind

Maps a zero cell (no reads) to either true_absence or undetected based on gates.
"""
function disambiguate_zero(w::Warrant, pol::ThresholdPolicy)::ZeroKind
  if w.depth !== nothing && w.depth ≥ pol.depth_min
    return Val(:true_absence)
  else
    return Val(:undetected)
  end
end

@kwdef struct ResidualInfo
  warranted_rank::Union{Nothing,String} = nothing
  status::Symbol                        = :SansFibre
  residual_named_count::Int             = 0
  residual_has_novel::Bool              = false
  contam_flag::Bool                     = false
  novel_reasons::Vector{String}         = String[]
  contam_reasons::Vector{String}        = String[]
end

"""
  residual_info(r::Receipt; named_candidates_at_species::Int=0, gates_ok_for_novel::Bool=false)::ResidualInfo

Residual summary at the next finer rank:
- residual_named_count: how many named species candidates are plausible (e.g., from vsearch top hits above a consideration floor)
- residual_has_novel: only true if parent clade is warranted (epi_status :Factive at Genus+), no named species clears species threshold,
                      and contamination gates pass (represented here by `gates_ok_for_novel`)
"""
function residual_info(r::Receipt; named_candidates_at_species::Int=0, gates_ok_for_novel::Bool=false, pol::ThresholdPolicy=ThresholdPolicy())
  st = epi_status(r, pol)
  wr = highest_warranted_rank(r.k, r.w, pol)
  has_novel = false
  novel_reasons = String[]
  contam_flag = false
  contam_reasons = String[]

  # "novel" is permitted only if:
  # - parent clade is warranted (Factive at some rank >= Genus),
  # - species threshold fails (wr != "Species"),
  # - and contamination gates/heuristics pass (e.g., neg-control, prevalence, chimera, replicates, blacklist).
  if st == :Factive && wr != "Species"
    if named_candidates_at_species == 0 && gates_ok_for_novel
      has_novel = true
      push!(novel_reasons, "Species threshold failed; parent clade warranted; contamination checks passed.")
    elseif named_candidates_at_species == 0 && !gates_ok_for_novel
      # treat as underdetermined or contaminant based on caller flags
      contam_flag = true
      push!(contam_reasons, "Failed contamination gates; no named species clears species threshold.")
    end
  end

  ResidualInfo(; warranted_rank=wr, status=st, residual_named_count=named_candidates_at_species,
               residual_has_novel=has_novel, contam_flag, novel_reasons, contam_reasons)
end

# -- Merge Checker -----------------------------------------------------------

"""
  can_merge(manA::Dict, manB::Dict)::NamedTuple

Decides if two runs are safe to analyze together.
Returns one of:
  (status="OK",)
  (status="OK_COLLAPSE", collapse_to="Division", reason="db_release_mismatch")
  (status="REFUSE", reason="normalisation_mismatch", details="rarefy vs rss")
Assumes manifests contain:
  man["normalisation"]          :: String
  man["db"]["key"]              :: String
  man["db"]["release"]          :: String
  man["db"]["ranks"]            :: Vector{String} (ordered, finest-last)
"""
function can_merge(manA::Dict, manB::Dict; allow_relative_abundance_if_norm_equal::Bool=true)
  # Normalisation mismatch blocks relative-abundance analysis
  if manA["normalisation"] != manB["normalisation"] && allow_relative_abundance_if_norm_equal
    return (status="REFUSE", reason="normalisation_mismatch", details="$(manA["normalisation"]) vs $(manB["normalisation"])")
  end
  # DB-key mismatch: collapse to lowest common rank
  if manA["db"]["key"] != manB["db"]["key"]
    common = lowest_common_rank(manA["db"]["ranks"], manB["db"]["ranks"])
    return common === nothing ? (status="REFUSE", reason="no_common_rank", details="") :
                                (status="OK_COLLAPSE", collapse_to=common, reason="db_key_mismatch")
  end
  # Same DB key, different release: collapse
  if manA["db"]["release"] != manB["db"]["release"]
    common = lowest_common_rank(manA["db"]["ranks"], manB["db"]["ranks"])
    return (status="OK_COLLAPSE", collapse_to=common, reason="db_release_mismatch")
  end
  # Otherwise OK
  return (status="OK",)
end

"""
  lowest_common_rank(a::Vector{String}, b::Vector{String})::Union{Nothing,String}

Given two ordered rank ladders (coarse->fine or fine->coarse; we treat as a simple intersection),
return the finest rank that exists in both (or nothing).
"""
function lowest_common_rank(a::Vector{String}, b::Vector{String})::Union{Nothing,String}
  s = Set(b)
  for r in a
    if r in s; return r; end
  end
  return nothing
end

# -- Backend Abstraction (lightweight) --------------------------------------
# Most of these functions are O(1) per row and trivially parallelizable on CPU.
# Specialized accelerators (TPU/NPU/VPU) do not fit well: this is not dense
# matmul or deep-learning inference. Still, we define a minimal interface to
# support batch verification on CPU and document where GPU/accelerators *might*
# be considered (e.g., GPU hashing), although not recommended.

abstract type ComputeBackend end

"""
  CPUBackend

Default backend. Uses multi-threading where helpful. This is the only
backend implemented; others are stubs with explanatory errors.
"""
struct CPUBackend <: ComputeBackend end

"""
  bulk_verify!(backend::CPUBackend, receipts::Vector{Receipt})::Vector{Bool}

Batch-verifies signatures. Uses Threads.@threads for large vectors.
"""
function bulk_verify!(::CPUBackend, receipts::Vector{Receipt})::Vector{Bool}
  n = length(receipts)
  out = Vector{Bool}(undef, n)
  # simple threaded loop; in practice, signature verification is already fast
  Threads.@threads for i in 1:n
    out[i] = verify_receipt(receipts[i])
  end
  return out
end

# --- Notional stubs for TPU/NPU/VPU backends -------------------------------
# These are deliberately unimplemented. Exposing them prevents future devs
# from assuming "hardware magic" is appropriate here. If someone proposes a
# GPU hashing kernel or other batched accelerator, add it as an *optional*
# specialized backend with careful determinism notes.

struct TPUBackend   <: ComputeBackend end
struct NPUBackend   <: ComputeBackend end
struct VPUBackend   <: ComputeBackend end

bulk_verify!(::TPUBackend, receipts::Vector{Receipt}) =
  error("TPU/NPU/VPU acceleration is not applicable: per-row hashes and O(1) checks are CPU-fast and not dense-matmul workloads.")

bulk_verify!(::NPUBackend, receipts::Vector{Receipt}) =
  error("TPU/NPU/VPU acceleration is not applicable: per-row hashes and O(1) checks are CPU-fast and not dense-matmul workloads.")

bulk_verify!(::VPUBackend, receipts::Vector{Receipt}) =
  error("TPU/NPU/VPU acceleration is not applicable: per-row hashes and O(1) checks are CPU-fast and not dense-matmul workloads.")

# -- End Module --------------------------------------------------------------

end # module Epistemic
