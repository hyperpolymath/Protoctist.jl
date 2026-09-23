include("core/Epistemic.jl")
using .Epistemic

function merge_taxa_with_evidence(table::DataFrame)
    for row in eachrow(table)
        attach_evidence!(row)
    end
    # Write @avec_fibre and evidence_uuid to DuckDB
end
