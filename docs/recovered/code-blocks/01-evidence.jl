# src/core/evidence.jl
module Evidence

export Standpoint, Warrant, EchoFiber, ResidualSet, Evidence, attach_fibre!

using ..Types: Taxon, ASV, Bootstrap, IdentityScore
using JSON3, UUIDs, Dates

abstract type Standpoint end
struct DADA2_PR2 <: Standpoint
    release::String
    error_profile_hash::String
end
struct VSEARCH_SILVA <: Standpoint
    release::String
    identity_threshold::Float64
end
# Add more as needed (Nanopore, consensus, etc.)

abstract type Warrant end
struct BootstrapWarrant <: Warrant
    value::Float64          # 0–100
    rank::String
end
struct IdentityWarrant <: Warrant
    percent::Float64
end
struct CompositeWarrant <: Warrant
    bootstrap::Float64
    identity::Float64
end

struct EchoFiber
    admissible_candidates::Vector{Taxon}   # the "cloud"
    proof_hash::String                     # cryptographic witness
    standpoint::Standpoint
end

struct ResidualSet
    count::Int
    modality::Symbol   # :factive, :belief, :collapsed_residue, :sans_fibre
    admissible_worlds::Vector{Dict{Symbol,Any}}
end

struct Evidence
    standpoint::Standpoint
    warrant::Warrant
    fibre::EchoFiber
    residuals::ResidualSet
    timestamp::DateTime
    uuid::UUID
end

# Core recovery function (mirrors your Agda Recover formula)
function recover_evidence(taxon::Taxon, bootstrap::Float64, standpoint::Standpoint)
    candidates = compute_admissible_candidates(taxon, bootstrap, standpoint)  # stub — implement with real logic
    fiber = EchoFiber(candidates, generate_proof_hash(taxon, candidates), standpoint)
    
    modality = if bootstrap ≥ 90
        :factive
    elseif bootstrap ≥ 60
        :belief
    else
        :collapsed_residue
    end
    
    residuals = ResidualSet(length(candidates), modality, [Dict(:taxon => c) for c in candidates])
    
    Evidence(standpoint, BootstrapWarrant(bootstrap, taxon.rank), fiber, residuals, now(), uuid4())
end

function generate_proof_hash(taxon, candidates)
    # Simple deterministic hash for prototype — replace with proper cryptographic witness
    string(hash((taxon.name, length(candidates), nameof(typeof(taxon)))))
end

function compute_admissible_candidates(taxon::Taxon, bootstrap::Float64, standpoint::Standpoint)
    # Placeholder — replace with real lookup from your reference DB + residual logic
    # This is where your mathematician’s ResidualEvidence.agda spec is turned into code
    if bootstrap > 85
        return [taxon]                     # singleton — full identification
    else
        return [taxon, Taxon("related_spp", taxon.rank)]  # example cloud
    end
end

# Attach to a merged table row
function attach_fibre!(row::Dict)
    taxon = Taxon(row["assigned_taxon"], row["rank"])
    boot = get(row, "bootstrap", 50.0)
    standpoint = DADA2_PR2("PR2_v5.0", "errorhash123")
    
    ev = recover_evidence(taxon, boot, standpoint)
    row["@avec_fibre"] = "echo:v1?k=$(nameof(typeof(ev.standpoint)))&w=boot:$(ev.warrant.value)&f_hash=$(ev.fibre.proof_hash)&residuals=$(ev.residuals.count)&modality=$(ev.residuals.modality)"
    row["evidence_uuid"] = string(ev.uuid)
    return row
end

end  # module
