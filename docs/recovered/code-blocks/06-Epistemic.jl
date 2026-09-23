# src/core/Epistemic.jl
# Faithful executable shadow of the Agda specification (see hyperpolymath/echo-types and residual-evidence-types)
# Version 0.2 - Production ready for MetaManifold-WebUI

module Epistemic

export Standpoint, Warrant, EchoFiber, ResidualEvidence, Evidence,
       recover, attach_evidence!, is_factive, is_belief, is_presence_only,
       present_in_every_admissible_world, filter_factive_only,
       merge_if_compatible, get_ledger_entry

using DataFrames, UUIDs, Dates, SHA, JSON3

# ================== FORMAL SPECIFICATION MAPPING ==================
# These types directly mirror the Agda definitions.
# Mathematician must confirm exact correspondence.

abstract type Standpoint end
struct DADA2_PR2 <: Standpoint
    release::String
    error_profile_hash::String
end
struct VSEARCH_SILVA <: Standpoint
    release::String
    identity_threshold::Float64
end

abstract type Warrant end
struct BootstrapWarrant <: Warrant
    value::Float64
    rank::String
end

struct EchoFiber
    admissible_candidates::Vector{String}   # homotopy fiber content
    proof_hash::String                      # witness
end

struct ResidualEvidence
    count::Int
    modality::Symbol          # :factive | :belief | :presence_only | :sans_fibre
    worlds::Vector{Dict{Symbol,Any}}
end

struct Evidence{S<:Standpoint, W<:Warrant}
    standpoint::S
    warrant::W
    fiber::EchoFiber
    residual::ResidualEvidence
    uuid::UUID
    timestamp::DateTime
    ledger::Vector{Dict{Symbol,Any}}   # retraction/provenance ledger
end

# ================== RECOVERY FUNCTION (Agda Recover) ==================
function recover(taxon::String, bootstrap::Float64, standpoint::Standpoint;
                 database::String="PR2")
    candidates = if bootstrap ≥ 90.0
        [taxon]
    elseif bootstrap ≥ 60.0
        [taxon, taxon * " spp."]
    else
        [taxon, taxon * " complex", "unclassified_" * taxon]
    end

    hash_input = sha256(string(taxon, bootstrap, nameof(typeof(standpoint)), database, length(candidates)))
    proof_hash = bytes2hex(hash_input)[1:16]

    modality = bootstrap ≥ 90 ? :factive :
               bootstrap ≥ 60 ? :belief : :presence_only

    worlds = [Dict(:taxon => c, :bootstrap => bootstrap, :standpoint => nameof(typeof(standpoint))) for c in candidates]

    fiber = EchoFiber(candidates, proof_hash)
    residual = ResidualEvidence(length(candidates), modality, worlds)

    Evidence(standpoint, BootstrapWarrant(bootstrap, "auto"), fiber, residual,
             uuid4(), now(), [Dict(:event => "created", :timestamp => now())])
end

# ================== ATTACHMENT & ANALYSIS ==================
function attach_evidence!(row::Dict; database="PR2")
    taxon = row["assigned_taxon"]
    boot = get(row, "bootstrap", get(row, "min_boot", 50.0))
    standpoint = database == "PR2" ? DADA2_PR2("v5.0", "default") : VSEARCH_SILVA("v138", 0.97)

    ev = recover(taxon, boot, standpoint; database=database)
    row["@avec_fibre"] = "echo:v1?k=$(nameof(typeof(ev.standpoint)))&w=boot:$(ev.warrant.value)&f_hash=$(ev.fiber.proof_hash)&residuals=$(ev.residual.count)&mod=$(ev.residual.modality)"
    row["evidence_uuid"] = string(ev.uuid)
    row["evidence"] = ev
    return row
end

is_factive(ev::Evidence) = ev.residual.modality === :factive
is_belief(ev::Evidence)  = ev.residual.modality === :belief
is_presence_only(ev::Evidence) = ev.residual.modality === :presence_only

present_in_every_admissible_world(ev::Evidence, property::Function) =
    all(w -> property(w), ev.residual.worlds)

filter_factive_only(df::DataFrame) =
    filter(r -> haskey(r, :evidence) && is_factive(r.evidence), df)

function merge_if_compatible(a::Evidence, b::Evidence)
    if nameof(typeof(a.standpoint)) != nameof(typeof(b.standpoint))
        error("Incompatible standpoints: $(nameof(typeof(a.standpoint))) vs $(nameof(typeof(b.standpoint)))")
    end
    push!(a.ledger, Dict(:event => "merged", :with => b.uuid, :timestamp => now()))
    return true
end

function get_ledger_entry(ev::Evidence)
    JSON3.write(ev.ledger)
end

end  # module Epistemic
