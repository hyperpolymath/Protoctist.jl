using Protoctist, EpistemicTypes, Cladistics

# 1) Load tables (with avec_fibre) and verify receipts
df = Protoctist.IO.read_ecsv("merged.ecsv")           # adds verified flag
pol = EpistemicTypes.ThresholdPolicy()                # or load from policy.yml

# 2) Compute epi_status + residuals per row (warrants + gates)
df.status = [EpistemicTypes.epi_status(r, pol) for r in df.receipt]
df.resid  = [EpistemicTypes.residual_info(r; named_candidates_at_species=df.named_cands[i],
                                          gates_ok_for_novel=df.novel_ok[i], pol=pol) for (i,r) in enumerate(df.receipt)]

# 3) Build taxonomy tree and roll up evidence
tree = Protoctist.Tree.build_taxonomy_tree(df.tax_path)  # PR2/SILVA ranks
annot = Protoctist.Tree.rollup_counts(tree, df, pol)     # per-node counts, f/b/res, novel/contam

# 4) Export annotated Newick + (if used) jplace + iTOL datasets
newick_str = Protoctist.Tree.annotate_newick(tree, annot)
write("tree.annotated.tree", newick_str)
Protoctist.IO.write_jplace(placements, "placements.jplace")       # if from EPA-ng/pplacer
Protoctist.IO.export_itol_bundle(tree, annot, "itol_bundle/")     # optional
