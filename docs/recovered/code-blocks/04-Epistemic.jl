# src/core/Epistemic.jl
# Faithful executable shadow of the Agda Echo + Epistemic + ResidualEvidence specification
# See: hyperpolymath/echo-types and residual-evidence-types for the formal spec

module Epistemic

export Standpoint, Warrant, EchoFiber, ResidualEvidence, Evidence,
       recover, attach_evidence!, is_factive, is_belief,
       present_in_every_admissible_world, filter_factive_only,
       merge_if_compatible

using DataFrames, UUIDs, Dates, SHA

# ─────────────────────────────────────────────────────────────
# 1. Standpoints (κ) — matches Agda standpoint-indexed modalities
# ─────────────────────────────────────────────────────────────
abstract type Standpoint end

struct DADA2_PR2 <: Standpoint
    release::String
    error_profile_hash::String
end

struct VSEARCH_SILVA <: Standpoint
    release::String
    identity_threshold::Float64
end

# Add more later (Nanopore, consensus, etc.)

# ─────────────────────────────────────────────────────────────
# 2. Warrants and Modalities
# ─────────────────────────────────────────────────────────────
abstract type Warrant end
struct BootstrapWarrant <: Warrant
    value::Float64   # 0–100
    rank::String
end

struct EchoFiber
    admissible_candidates::Vector{String}   # the "cloud" / homotopy fiber
    proof_hash::String                      # witness
end

struct ResidualEvidence
    count::Int
    modality::Symbol      # :factive, :belief, :presence_only, :sans_fibre
    worlds::Vector{Dict{Symbol,Any}}
end

struct Evidence{S<:Standpoint,W<:Warrant}
    standpoint::S
    warrant::W
    fiber::EchoFiber
    residual::ResidualEvidence
    uuid::UUID
    timestamp::DateTime
end

# ─────────────────────────────────────────────────────────────
# 3. Recovery function — direct translation of your Agda Recover
# ─────────────────────────────────────────────────────────────
function recover(taxon::String, bootstrap::Float64, standpoint::Standpoint)
    # This is the place where the mathematician should ensure exact match with Agda
    candidates = bootstrap ≥ 85.0 ? [taxon] : [taxon, taxon * " spp.", taxon * " complex"]

    hash_input = sha256(string(taxon, bootstrap, nameof(typeof(standpoint)), length(candidates)))
    proof_hash = bytes2hex(hash_input)[1:16]

    modality = bootstrap ≥ 90 ? :factive :
               bootstrap ≥ 60 ? :belief : :presence_only

    worlds = [Dict(:taxon => c, :bootstrap => bootstrap) for c in candidates]

    fiber = EchoFiber(candidates, proof_hash)
    residual = ResidualEvidence(length(candidates), modality, worlds)

    Evidence(standpoint, BootstrapWarrant(bootstrap, "unknown"), fiber, residual, uuid4(), now())
end

# ─────────────────────────────────────────────────────────────
# 4. Attachment and analysis helpers
# ─────────────────────────────────────────────────────────────
function attach_evidence!(row::Dict)
    taxon = row["assigned_taxon"]
    boot = get(row, "bootstrap", 50.0)
    standpoint = DADA2_PR2("PR2_v5.0", "default_error_hash")

    ev = recover(taxon, boot, standpoint)
    row["@avec_fibre"] = "echo:v1?k=$(nameof(typeof(ev.standpoint)))&w=boot:$(ev.warrant.value)&f_hash=$(ev.fiber.proof_hash)&residuals=$(ev.residual.count)&mod=$(ev.residual.modality)"
    row["evidence"] = ev
    return row
end

is_factive(ev::Evidence) = ev.residual.modality === :factive
is_belief(ev::Evidence)  = ev.residual.modality === :belief

present_in_every_admissible_world(ev::Evidence, property::Function) =
    all(w -> property(w), ev.residual.worlds)

filter_factive_only(df::DataFrame) =
    filter(r -> haskey(r, :evidence) && is_factive(r.evidence), df)

function merge_if_compatible(ev1::Evidence, ev2::Evidence)
    ev1.standpoint == ev2.standpoint || error("Incompatible standpoints - cannot merge")
    return true
end

end  # module Epistemic
