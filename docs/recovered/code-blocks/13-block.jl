using JSON3

struct NodeAnn
  f::Float64; b::Float64; res::Float64; novel::Int; contam::Int; counts::Int
end

# Suppose you have a Dict{String,NodeAnn} ann keyed by node name
function annotate_newick(newick::String, ann::Dict{String,NodeAnn})
  # naive approach: replace 'Name:' with 'Name[&...]:'
  return replace(newick, r"([A-Za-z0-9_]+):" => s -> begin
    name = match(r"([A-Za-z0-9_]+):", s.match).captures[1]
    if haskey(ann, name)
      a = ann[name]
      "$(name)[&f=$(a.f),b=$(a.b),res=$(a.res),novel=$(a.novel),contam=$(a.contam),counts=$(a.counts)]:"
    else
      s.match
    end
  end)
end
