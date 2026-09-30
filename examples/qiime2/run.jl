# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# QIIME 2 + SILVA -> Protoctist.
# Run from this directory:  julia --project=../.. run.jl
#
# taxonomy.tsv:      qiime tools export --input-path taxonomy.qza --output-path .
# feature-table.tsv: biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

using Protoctist

# "d__Bacteria; p__Proteobacteria" -> "Bacteria;Proteobacteria".
# An empty rank ("g__") becomes a blank, which Protoctist reads as unfilled.
strip_prefixes(taxon) =
    join((replace(strip(l), r"^[a-z]__" => "") for l in split(taxon, ';')), ';')

lineage = Dict{String,String}()
for line in Iterators.drop(eachline("taxonomy.tsv"), 1)
    id, taxon, _ = split(line, '\t')
    lineage[id] = strip_prefixes(taxon)
end

# One row per feature, reads summed over samples. biom writes counts as floats.
open("lineages.tsv", "w") do io
    println(io, "feature\tlineage\tcount")
    for line in eachline("feature-table.tsv")
        startswith(line, '#') && continue
        id, counts... = split(line, '\t')
        println(io, id, '\t', lineage[id], '\t', round(Int, sum(parse.(Float64, counts))))
    end
end

run = load_protist_run("lineages.tsv"; sep = '\t', ladder = :silva)
println(run)
for (taxon, reads) in clade_cumulus(run)
    println(rpad(taxon, 22), reads)
end
