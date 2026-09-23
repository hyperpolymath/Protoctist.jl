include("core/evidence.jl")
using .Evidence

function merge_taxa_with_evidence(merged_table::DataFrame)
    for row in eachrow(merged_table)
        attach_fibre!(row)          # adds @avec_fibre and evidence_uuid
    end
    # Write to DuckDB with new columns
    DBInterface.execute(db, "ALTER TABLE merged ADD COLUMN IF NOT EXISTS avec_fibre TEXT")
    # ... insert updated rows
end
