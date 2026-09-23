include("core/Epistemic.jl")
using .Epistemic
for row in eachrow(merged_table)
    attach_evidence!(row; database=config.database)
end
