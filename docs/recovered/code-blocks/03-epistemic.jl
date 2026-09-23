# src/core/epistemic.jl
module Epistemic

export Standpoint, Warrant, EchoFiber, ResidualEvidence, Evidence, 
       attach_evidence!, is_factive, filter_factive_only, 
       present_in_every_admissible_world

abstract type Standpoint end
struct DADA2_PR2   <: Standpoint; release::String; end
struct VSEARCH_SILVA <: Standpoint; release::String; end

abstract type Warrant end
struct Bootstrap <: Warrant; value::Float64; rank::String; end
struct Identity  <: Warrant; value::Float64; end

struct EchoFiber
    admissible::Vector{String}   # the "cloud" — list of possible taxa
    proof_hash::String
end

struct ResidualEvidence
    count::Int
    modality::Symbol          # :factive, :belief, :presence_only, :sans_fibre
    worlds::Vector{Dict}
end

struct Evidence{S<:Standpoint, W<:Warrant}
    standpoint::S
    warrant::W
    fiber::EchoFiber
    residual::ResidualEvidence
end

# Core recovery function — mirrors your Agda Recover formula
function recover(taxon::String, bootstrap::Float64, standpoint::Standpoint)
    candidates = (bootstrap ≥ 85) ? [taxon] : [taxon, taxon*" spp."]
    
    fiber = EchoFiber(candidates, string(hash((taxon, bootstrap, standpoint))))
    
    modality = bootstrap ≥ 90 ? :factive : 
               bootstrap ≥ 60 ? :belief : :presence_only
    
    ResidualEvidence(length(candidates), modality, 
                    [Dict(:taxon => c, :bootstrap => bootstrap) for c in candidates])
end

function attach_evidence!(row::Dict)
    taxon = row["assigned_taxon"]
    boot  = get(row, "bootstrap", 50.0)
    stand = DADA2_PR2("PR2_v5")
    
    res = recover(taxon, boot, stand)
    ev  = Evidence(stand, Bootstrap(boot, row["rank"]), res.fiber, res)
    
    row["@avec_fibre"] = "echo:v1?k=$(typeof(stand))&w=boot:$(boot)&residuals=$(res.count)&mod=$(res.modality)&hash=$(res.fiber.proof_hash)"
    row["evidence"] = ev
    return row
end

# Useful predicates for analysis
is_factive(ev::Evidence) = ev.residual.modality === :factive
present_in_every_admissible_world(ev::Evidence, property::Function) = 
    all(w -> property(w), ev.residual.worlds)

filter_factive_only(df::DataFrame) = filter(r -> is_factive(r.evidence), df)

end
