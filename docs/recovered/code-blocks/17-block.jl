# In the forked MetaManifold
using EchoTypes          # Your mathematician's formal shadow
using Cladistics         # The cladistic engine
using Epistemic          # Our higher-level wrapper (to be written)

# Example usage in CladeCumulus
function render_cladistic(table::DataFrame)
    evidence = parse_fibre.(table.avec_fibre)        # from EchoTypes
    enriched = attach_evidence!.(evidence)           # from Epistemic
    tree = build_cumulative_tree(enriched)           # from Cladistics.jl
    render_with_clouds(tree)                         # uses residual counts and modalities
end
