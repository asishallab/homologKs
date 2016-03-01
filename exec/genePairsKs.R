require(paranomeKsR)

message("USAGE: Rscript path/2/paranomeKsR/exec/genePairsKs.R blastp+_format6_result_table.tsv start_row stop_row coding_sequences.fasta redis_url out_path [intermediate_files_directory] [filter_MSA (TRUE|FALSE default: TRUE)] [KaKs_Calculator's method (default: YN)]")

#' Read and process input arguments:
input.args <- commandArgs(trailingOnly = TRUE)

start.row <- as.integer(input.args[[2]])
stop.row <- as.integer(input.args[[3]])
blastp.res.tbl <- read.table(input.args[[1]], stringsAsFactors = FALSE)[start.row:stop.row, 
    ]
cds <- readDNAStringSet(input.args[[4]])
redisConnect(input.args[[5]])

t.d <- if (length(input.args) > 6) {
    file.path(input.args[[7]])
} else tempdir()

if (length(input.args) > 7) {
    options(paranomeKsR.filter.MSA = as.logical(input.args[[8]])[[1]])
}

if (length(input.args) > 8) {
    options(paranomeKsR.ks.method = input.args[[9]])
}

#' Start the computation:
b.tbl <- blastp.res.tbl[with(blastp.res.tbl, which(V1 != V2 & V3 >= 30 & V4 >= 150)), 
    c("V1", "V2")]
b.tbl$V3 <- as.numeric(NA)
for (i in 1:nrow(b.tbl)) {
    x <- unlist(b.tbl[i, 1:2])
    pair.ks <- computeKsPipeline(x, cds, t.d)
    if (!is.null(pair.ks) && !is.na(pair.ks)) 
        b.tbl[which(b.tbl$V1 == x[[1]] & b.tbl$V2 == x[[2]]), "V3"] <- pair.ks
}

# Write results:
o.path <- file.path(input.args[[6]])
write.table(b.tbl, o.path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

message("DONE") 
