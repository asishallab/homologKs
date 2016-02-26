#' Checks wether in the current Redis server there is an entry for gene pair
#' 'gene.pair'
#'
#' @param Character vector of length two holding members A and B of the gene
#' pair.
#'
#' @return Boolean. TRUE of and only if the alphabetically sorted gene
#' accessions 'A_B' appear in the Redis store.
#' @export
genePairInRedis <- function(gene.pair) {
    g.pair <- paste(sort(gene.pair), collapse = "_")
    !is.null(redisGet(g.pair))
}

#' Puts the integer value 1 into the Redis Server using the alphabetically
#' sorted gene accessions in 'gene.pair' as key.
#'
#' @param Character vector of length two holding members A and B of the gene
#' pair.
#'
#' @return The output of calling rredis::redisSet(...)
#' @export
putGenePairIntoRedis <- function(gene.pair) {
    g.pair <- paste(sort(gene.pair), collapse = "_")
    redisSet(g.pair, 1)
}

#' Runs the pipeline to compute the Ks values of gene pair 'x' with coding
#' sequences held in 'cds'. Temporary files are stored in 't.d'.
#'
#' @param x character vector of gene IDs
#' @param cds an instance of Biostrings::XStringSet holding the unaligned
#' coding sequences of the genes in 'x'
#' @param t.d the directory to store the files in, default is tempdir()
#'
#' @return A numeric being the result of invoking 'KaKs_Calculator' on the
#' aligned coding sequences, or NA if an error occurres.
#' @export
computeKsPipeline <- function(x, cds, t.d = tempdir()) {
    p.n <- paste(sort(x), collapse = "_")
    tryCatch({
        if (!genePairInRedis(x)) {
            putGenePairIntoRedis(x)
            pair.cds.msa <- alignCodingSequencesPipeline(cds[x], t.d, p.n)
            axt.path <- file.path(t.d, paste(p.n, ".axt", sep = ""))
            writeLines(unlist(c(p.n, lapply(pair.cds.msa, toString))), con = axt.path)
            ks.out <- sub(".axt", "_KaKs_results.txt", axt.path, fixed = TRUE)
            system(paste("KaKs_Calculator -i", axt.path, "-o", ks.out))
            read.table(ks.out, header = TRUE, stringsAsFactors = FALSE)[[1, "Ks"]]
        }
    }, error = function(e) {
        warning("Computing Ks value of gene pair '", p.n, "' caused an error:\n", 
            e)
        NA
    })
}

#' Generates a distance matrix for genes in family 'fam.name' with distance
#' values (Ks) obtained from 'ks.tbl'.
#'
#' @param fam.name The name of the gene family to generate a distance matrix
#' for
#' @param fam.tbl An instance of base::data.frame with two columns. Columns one
#' must hold the family names and column two the gene members (IDs). Default is
#' 'chi.paralogous.fams'.
#' @param ks.tbl An instance of base::data.frame with three columns. Column one
#' and two hold gene identifier, and column three the measured distance (Ks)
#' value. Default is 'chi.paranome.ks'.
#' @param na.as.ks The value to use, if for a given gene pair no distance value
#' can be found in 'ks.tbl'. Default is 1e6
#'
#' @return An instance of base::matrix(numeric) holding the pairwise distances
#' (Ks) of the genes in family 'fam.name'. Obviously the matrix' diagonal will
#' be zero.
#' @export
familyDist <- function(fam.name, fam.tbl = chi.paralogous.fams, ks.tbl = chi.paranome.ks, 
    na.as.ks = 1e+06) {
    genes <- sort(unique(fam.tbl[which(fam.tbl[, 1] == fam.name), 2]))
    x <- ks.tbl[with(ks.tbl, which(V1 %in% genes & V2 %in% genes)), 3:4]
    genesKs <- setNames(x$V3, x$V4)
    .Call("familyDistCpp", genes, genesKs, na.as.ks)
}

#' Checks wether node 'tr.nd' is a tip or not in tree 'phylo.tr'.
#'
#' @param phylo.tr an instance of ape::phylo representing the tree in which to
#' look up the node.
#' @param tr.nd an integer identifying the node for which to determine wether
#' it is as tip or not.
#'
#' @return TRUE if and only if 'tr.nd' is a tip of tree 'phylo.tr', FALSE
#' otherwise.
#' @export
isTreeTip <- function(phylo.tr, tr.nd) {
    tr.tips <- Descendants(phylo.tr, getRoot(phylo.tr), type = "tips")
    tr.nd %in% tr.tips
}

#' Extracts the tip labels of the tips found in the subtree spanned by 'tr.nd'
#' in the tree 'phylo.tr'. If the argument node is a tip its label is returned.
#'
#' @param phylo.tr an instance of ape::phylo representing the tree
#' @param tr.nd an integer identifying the node of interest
#'
#' @return A character vector of tip labels 
#' @export
getDescendantTipsOrSelf <- function(phylo.tr, tr.nd) {
    i <- if (isTreeTip(phylo.tr, tr.nd)) 
        tr.nd else Descendants(phylo.tr, tr.nd, type = "tips")[[1]]
    phylo.tr$tip.label[i]
}

#' Identifies all given distances between relevant tips. Relevant means that
#' all gene pairs are considered between the left descending subtree and the
#' right descending one. This means (1) 'fam.tr' is expected to be a binary
#' tree in which each node has only two descending branches, and (2) the
#' weighted distance for a given node is computed from all pairs generated
#' between left and right descending subtrees. For more details see
#' Maere, Steven, Stefanie De Bodt, Jeroen Raes, Tineke Casneuf, Marc Van
#' Montagu, Martin Kuiper, and Yves Van de Peer. “Modeling Gene and Genome
#' Duplications in Eukaryotes.” Proceedings of the National Academy of Sciences
#' of the United States of America 102, no. 15 (April 12, 2005): 5454–59.
#' doi:10.1073/pnas.0501102102.
#'
#' @param fam.tr an instance of ape::phylo representing the tree of the
#' currently investigated gene family. The tree MUST be a binary tree in which
#' each inner node has exactly two descending branches.
#' @param tr.nd The inner tree node for which to compute the weighted distances
#' (Ks).
#' @param ks.tbl an instance of base::data.frame with three columns, in which
#' the first two hold gene identifier and the third holds the measured distance
#' (Ks) values. Default is 'chi.paranome.ks'.
#'
#' @return A numeric value computed as the mean of all distances of 'relevant'
#' gene pairs.
#' @export
weightedDistsForNode <- function(fam.tr, tr.nd, ks.tbl = chi.paranome.ks) {
    desc.nds <- Descendants(fam.tr, tr.nd, type = "children")
    gene.pairs <- expand.grid(getDescendantTipsOrSelf(fam.tr, desc.nds[[1]]), getDescendantTipsOrSelf(fam.tr, 
        desc.nds[[2]]), stringsAsFactors = FALSE)
    ks.vals <- ks.tbl[with(ks.tbl, with(gene.pairs, which(V1 %in% Var1 & V2 %in% 
        Var2 | V1 %in% Var2 & V2 %in% Var1))), 3]
    sum(ks.vals)/length(ks.vals)
}

#' Generates a named numeric vector in which for each inner node of 'fam.tr'
#' the weighted distance is held.
#'
#' @param fam.tr an instance of ape::phylo representing the families tree
#' @param ks.tbl an instance of base::data.frame with three columns, in which
#' the first two hold gene identifier and the third holds the measured distance
#' (Ks) values. Default is 'chi.paranome.ks'.
#'
#' @return A named numeric vector in which names are 'node_X' and values are
#' each inner nodes weighted distances.
#' @export
weightedDistsForTree <- function(fam.tr, ks.tbl = chi.paranome.ks) {
    tr.inner.nodes <- sort(setdiff(unlist(fam.tr$edge), Descendants(fam.tr, getRoot(fam.tr), 
        type = "tips")[[1]]))
    nd.names <- paste("node", tr.inner.nodes, sep = "_")
    setNames(as.numeric(unlist(lapply(tr.inner.nodes, function(tr.nd) weightedDistsForNode(fam.tr, 
        tr.nd, ks.tbl)))), nd.names)
}

#' Test function using RUnit to verify the mean distance values for inner tree
#' nodes.
#'
#' @param NONE
#'
#' @return TRUE if and only if no validation has failed.
#' @export
testWeightedDistsForTree <- function() {
    fam.name <- "fam300"
    fam.dists <- familyDist(fam.name)
    fam.clust <- hclust(as.dist(fam.dists), method = "single")
    fam.tr <- as.phylo(fam.clust)
    fam.ks <- weightedDistsForTree(fam.tr)
    checkEquals(fam.tr$Nnode, length(fam.ks))
    x <- getDescendantTipsOrSelf(fam.tr, 19)
    checkEquals(chi.paranome.ks[[with(chi.paranome.ks, which(V1 == x[[1]] & V2 == 
        x[[2]] | V1 == x[[2]] & V2 == x[[1]])), 3]], fam.ks[["node_19"]])
    x.tips <- getDescendantTipsOrSelf(fam.tr, 18)
    x.comp <- "CARHR048310.1"
    checkEquals(as.numeric(mean(chi.paranome.ks[with(chi.paranome.ks, which(V1 %in% 
        x.tips & V2 == x.comp | V1 == x.comp & V2 %in% x.tips)), 3])), as.numeric(fam.ks[["node_17"]]))
    # If we get to here, everything is fine:
    TRUE
} 
