module Epistemic

using JSON3, SHA, Dates

struct Standpoint
  tool::String; tool_ver::String; db_key::String; db_release::String
  params_hash::String; normalisation::String
end
struct Warrant
  boot_min::Union{Nothing,Float64}
  id::Union{Nothing,Float64}
  cov::Union{Nothing,Float64}
  depth::Union{Nothing,Int}
  neg_leak::Union{Nothing,Float64}
  chimera::Union{Nothing,Bool}
end
struct ProjectionY
  taxon::String; rank::String; feature_id::String; sample_id::String
end
struct Receipt
  v::String; k::Standpoint; w::Warrant; y::ProjectionY; sig::String
end

b64url(x) = replace(base64encode(x), ('+'=>'-','/'=>'_','='=>""))
canon(k,w,y) = JSON3.write((k=k,w=w,y=y); allow_inf=false)
make_receipt(k::Standpoint, w::Warrant, y::ProjectionY) = begin
  sig = "sha256:" * b64url(SHA.sha256(canon(k,w,y)))
  Receipt("echo:v1", k,w,y, sig)
end
verify_receipt(r::Receipt) = ("sha256:" * b64url(SHA.sha256(canon(r.k,r.w,r.y)))) == r.sig

# Example thresholds (configurable)
const DADA2_THRESH = Dict("Species"=>0.98, "Genus"=>0.88, "Class"=>0.75, "Division"=>0.70)
const VSEARCH_THRESH = Dict("Species"=>(id=0.99,cov=0.95), "Genus"=>(id=0.95,cov=0.90), "Division"=>(id=0.85,cov=0.80))

function highest_warranted_rank(k::Standpoint, w::Warrant)
  if k.tool == "dada2" && w.boot_min !== nothing
    for r in ("Species","Genus","Class","Division")
      if w.boot_min ≥ get(DADA2_THRESH, r, 2.0); return r; end
    end
  elseif k.tool == "vsearch" && w.id !== nothing && w.cov !== nothing
    for r in ("Species","Genus","Division")
      thr = get(VSEARCH_THRESH, r, nothing)
      if thr !== nothing && w.id ≥ thr.id && w.cov ≥ thr.cov; return r; end
    end
  end
  return nothing
end

# Gates (E): tune per study/run
passes_gates(w::Warrant; depth_min=100, leak_max=0.005) =
  (w.depth !== nothing && w.depth ≥ depth_min) &&
  (w.neg_leak === nothing || w.neg_leak ≤ leak_max) &&
  (w.chimera === nothing || w.chimera == false)

function epi_status(r::Receipt; depth_min=100, leak_max=0.005)
  if !verify_receipt(r); return :SansFibre; end
  if !passes_gates(r.w; depth_min, leak_max); return :Collapsed; end
  wr = highest_warranted_rank(r.k, r.w)
  return wr === nothing ? :Belief : :Factive
end

# Residual summary (fast; prefer vsearch top-hits if present)
struct ResidualInfo
  warranted_rank::Union{Nothing,String}
  status::Symbol
  residual_named_count::Int
  residual_has_novel::Bool
  contam_flag::Bool
end

function residual_info(r::Receipt; named_candidates_at_species::Int=0, gates_ok_for_novel::Bool=false)
  st = epi_status(r)
  wr = highest_warranted_rank(r.k, r.w)
  # novel only if: parent rank warranted (e.g., Genus), species fails threshold, and gates_ok
  has_novel = (!isnothing(wr) && wr != "Species" && named_candidates_at_species == 0 && gates_ok_for_novel)
  ResidualInfo(wr, st, named_candidates_at_species, has_novel, false)
end

end # module
