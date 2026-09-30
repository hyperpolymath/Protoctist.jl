# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell <j.d.a.jewell@open.ac.uk>
#
# DADA2 + PR2: turn assignTaxonomy output into the lineage table Protoctist reads.
# Run from this directory:  Rscript lineages.R   (writes lineages.csv)
#
# The point of the recipe is the NA handling: DADA2 leaves unresolved ranks
# as NA, and they must become empty fields. Pasted as-is, they would become
# taxa literally named "NA".

# Stand-in for dada2::assignTaxonomy(seqs, "pr2_version_5.0.0_SSU_dada2.fasta.gz")
# and a seqtab.nochim: same shapes, three ASVs, two samples.
taxa <- matrix(c(
  "Eukaryota","TSAR","Alveolata","Ciliophora","Spirotrichea",NA,NA,NA,NA,
  "Eukaryota","TSAR","Rhizaria","Cercozoa",NA,"Glissomonadida",NA,NA,NA,
  "Eukaryota","Obazoa","Opisthokonta","Fungi",NA,NA,NA,NA,NA),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("ASV1","ASV2","ASV3"),
    c("Domain","Supergroup","Division","Subdivision","Class","Order","Family","Genus","Species")))
seqtab <- matrix(c(120, 40, 0, 35, 60, 5), nrow = 2, byrow = TRUE,
                 dimnames = list(c("S01","S02"), rownames(taxa)))

# --- the recipe ---
lineage <- apply(taxa, 1, function(r) paste(ifelse(is.na(r), "", r), collapse = ";"))
out <- data.frame(feature = rownames(taxa), lineage = lineage,
                  count = colSums(seqtab)[rownames(taxa)])
write.csv(out, "lineages.csv", row.names = FALSE)
