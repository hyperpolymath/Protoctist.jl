# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# DADA2 + PR2 -> Protoctist -> iTOL.
# Run from this directory:  julia --project=../.. run.jl
# lineages.csv is what lineages.R writes; it is committed so this runs without R.

using Protoctist

run = load_protist_run("lineages.csv")
println(run)
for (taxon, reads) in clade_cumulus(run)[1:5]
    println(rpad(taxon, 16), reads)
end
for f in export_iqtree_annotated(run, "out")
    println("wrote ", f)
end
