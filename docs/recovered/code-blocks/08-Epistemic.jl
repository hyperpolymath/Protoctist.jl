# src/core/Epistemic.jl
module Epistemic
export Standpoint, Warrant, EchoFiber, ResidualEvidence, Evidence,
       recover, attach_evidence!, is_factive, present_in_every_admissible_world,
       filter_factive_only, merge_if_compatible, get_ledger

using DataFrames, UUIDs, Dates, SHA, JSON3

abstract type Standpoint end
struct DADA2_PR2 <: Standpoint; release::String; error_profile_hash::String; end
struct VSEARCH_SILVA <: Standpoint; release::String; identity_threshold::Float64; end

abstract type Warrant end
struct BootstrapWarrant <: Warrant; value::Float64; rank::String; end

struct EchoFiber
    admissible_candidates::Vector{String}
    proof_hash::String
end

abstract type Modality end
struct Factive <: Modality end
struct Belief <: Modality end
struct PresenceOnly <: Modality end

struct ResidualEvidence
    count::Int
    modality::Modality
    worlds::Vector{Dict{Symbol,Any}}
end

struct Evidence{S<:Standpoint, W<:Warrant}
    standpoint::S
    warrant::W
    fiber::EchoFiber
    residual::ResidualEvidence
    uuid::UUID
    timestamp::DateTime
    ledger::Vector{Dict{Symbol,Any}}
end

function recover(taxon::String, bootstrap::Float64, standpoint::Standpoint; database="PR2")
    candidates = bootstrap ≥ 90.0 ? [taxon] :
                 bootstrap ≥ 60.0 ? [taxon, taxon*" spp."] :
                 [taxon, taxon*" complex", "unclassified_"*taxon]
    proof_hash = bytes2hex(sha256(string(taxon, bootstrap, nameof(typeof(standpoint)), database)))[1:16]
    modality = bootstrap ≥ 90 ? Factive() : bootstrap ≥ 60 ? Belief() : PresenceOnly()
    worlds = [Dict(:taxon => c, :bootstrap => bootstrap) for c in candidates]
    fiber = EchoFiber(candidates, proof_hash)
    residual = ResidualEvidence(length(candidates), modality, worlds)
    Evidence(standpoint, BootstrapWarrant(bootstrap, "auto"), fiber, residual,
             uuid4(), now(), [Dict(:event => "created", :time => now())])
end

function attach_evidence!(row::Dict; database="PR2")
    taxon = row["assigned_taxon"]
    boot = get(row, "bootstrap", 50.0)
    standpoint = database == "PR2" ? DADA2_PR2("v5.0", "default") : VSEARCH_SILVA("v138", 0.97)
    ev = recover(taxon, boot, standpoint; database=database)
    row["@avec_fibre"] = "echo:v1?k=$(nameof(typeof(ev.standpoint)))&w=boot:$(ev.warrant.value)&f_hash=$(ev.fiber.proof_hash)&residuals=$(ev.residual.count)&mod=$(nameof(typeof(ev.residual.modality)))"
    row["evidence_uuid"] = string(ev.uuid)
    row["evidence"] = ev
    row
end

is_factive(ev::Evidence) = ev.residual.modality isa Factive
present_in_every_admissible_world(ev::Evidence, property::Function) = all(w -> property(w), ev.residual.worlds)
filter_factive_only(df::DataFrame) = filter(r -> haskey(r, :evidence) && is_factive(r.evidence), df)

function merge_if_compatible(a::Evidence, b::Evidence)
    nameof(typeof(a.standpoint)) == nameof(typeof(b.standpoint)) || error("Incompatible standpoints")
    push!(a.ledger, Dict(:event => "merged", :with => string(b.uuid), :time => now()))
    true
end

get_ledger(ev::Evidence) = JSON3.write(ev.ledger)

end
