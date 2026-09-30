# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# EPA-ng -> jplace -> one confident placement per query.
# Run from this directory:  julia --project=../.. run.jl

using Protoctist
const PIO = Protoctist.IO

b = PIO.read_jplace("epa_result.jplace")
println(length(b.placements), " candidate placements")

best = PIO.best_placements(b)                    # highest like_weight_ratio per query
for p in best
    println(rpad(p.name, 6), " edge ", p.edge_num, "  lwr ", p.like_weight_ratio)
end

confident = filter(p -> p.like_weight_ratio >= 0.9, best)
PIO.write_jplace(PIO.PlacementBatch(b.tree, confident, b.fields, b.version), "confident.jplace")
println("wrote confident.jplace with ", length(confident), " placements")
