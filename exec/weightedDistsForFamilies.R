require(paranomeKsR)

message("USAGE: Rscript path/2/paranomeKsR/exec/weightedDistsForFamilies.R start_ind stop_ind out_dir")

input.args <- commandArgs(trailingOnly = TRUE)

start.ind <- as.integer(input.args[[1]])
stop.ind <- as.integer(input.args[[2]])
out.dir <- file.path(input.args[[3]])

families <- unique(chi.paralogous.fams[, 1])[start.ind:stop.ind]
dists.tbl <- data.frame(family = c(), tree.node = c(), weighted.distances = c(), 
    stringsAsFactors = FALSE)
for (fam.name in families) {
    fam.dists <- weightedDistsForFamily(fam.name)
    dists.tbl <- rbind(dists.tbl, data.frame(family = fam.name, tree.node = names(fam.dists$weighted.distances), 
        weighted.distances = fam.dists$weighted.distances, stringsAsFactors = FALSE))
    write.tree(fam.dists$cluster, file.path(out.dir, paste(fam.name, "_cluster.newick", 
        sep = "")))
}

write.table(dists.tbl, file.path(out.dir, paste("weighted_Ks_table_families_", start.ind, 
    "_", stop.ind, ".tsv", sep = "")), col.names = FALSE, row.names = FALSE, quote = FALSE, 
    sep = "\t")

message("DONE") 
