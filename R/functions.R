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
#' @param cds an instance of Biostrings::DNAStringSet holding the unaligned
#' coding sequences of the genes in 'x'. Default is 'paranomeKsR::chi.cds'.
#' @param aas an instance of Biostrings::AAStringSet holding the amino acid
#' sequences for the genes in 'cds'. You can use Biostrings::translate or the
#' program MACSE to generate them from the coding sequences. Default is
#' 'paranomeKsR::chi.aas'.
#' @param t.d the directory to store the files in, default is tempdir()
#' @param codeml.call the call passed to system(...) in order to start codeml.
#' Default is 'codeml', use option 'paranomeKsR.codeml.call' to set another
#' default.
#'
#' @return A numeric being the result of invoking PAML's codeml on the aligned
#' coding sequences, or NA if an error occurres.
#' @references Yang, Z. and Nielsen, R. (2000) Mol. Biol. Evol., 17, 32-43.
#' @export
computeKsPipeline <- function(x, cds = chi.cds, aas = chi.aas, t.d = tempdir(), codeml.call = getOption("paranomeKsR.codeml.call", 
    "codeml")) {
    p.n <- paste(sort(x), collapse = "_")
    tryCatch({
        if (!genePairInRedis(x)) {
            putGenePairIntoRedis(x)
            pair.cds.msa.path <- alignCodingSequencesPipeline(cds[x], aas[x], t.d, 
                p.n)
            pair.codeml.in.out <- renderCodemlInputs(x, pair.cds.msa.path)
            orig.dir <- getwd()
            setwd(t.d)
            system(paste(codeml.call, pair.codeml.in.out[["in"]]))
            setwd(orig.dir)
            as.numeric(system(paste("tail -1", pair.codeml.in.out[["out"]], "| awk -F \"dS = \" '{print $2}'"), 
                intern = TRUE))
        }
    }, error = function(e) {
        message("ERROR: Computing Ks value of gene pair '", p.n, "' caused an error:\n", 
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
#' @param ks.tbl An instance of base::data.frame with four columns. Column one
#' and two hold gene identifier, and column three the measured distance (Ks)
#' value. Default is 'chi.paranome.ks'. Column four holds the gene pair's
#' identifiers concatonated in alphabetical order and separated by '_'.
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
    familyDistCpp(genes, genesKs, na.as.ks)
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
#' @param weight.func The function used to compute the weighted distance at the
#' family tree's node 'tr.nd'. Default is stats::median. Other functions can be
#' set as an option, use
#' 'options('ksParanomeR.dist.weight.func'=yourFunction)'.
#'
#' @return A numeric value computed as the mean of all distances of 'relevant'
#' gene pairs.
#' @export
weightedDistsForNode <- function(fam.tr, tr.nd, ks.tbl = chi.paranome.ks, weight.func = getOption("ksParanomeR.dist.weight.func", 
    stats::median)) {
    desc.nds <- Descendants(fam.tr, tr.nd, type = "children")
    gene.pairs <- expand.grid(getDescendantTipsOrSelf(fam.tr, desc.nds[[1]]), getDescendantTipsOrSelf(fam.tr, 
        desc.nds[[2]]), stringsAsFactors = FALSE)
    ks.vals <- ks.tbl[with(ks.tbl, with(gene.pairs, which(V1 %in% Var1 & V2 %in% 
        Var2 | V1 %in% Var2 & V2 %in% Var1))), 3]
    weight.func(ks.vals)
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

#' Pipeline to compute the weighted distances (Ks) for a family of genes.
#' Weighting is done based on each bipartition of subtrees at a given internal
#' node of the generated family tree, such that for each node the mean of the
#' distances is computed for all gene pairs, with one member from each subtree.
#'
#' @param fam.name the name of the family as it appears in 'fam.tbl'
#' @param fam.tbl a two column instance of base::data.frame where the first
#' column holds the family name and the second the identifier of gene members.
#' Default is 'chi.paralogous.fams'
#' @param dist.tbl an instance of base::data.frame with four columns. Column
#' one and two hold the gene identifier, column three the distance 'Ks' value,
#' and column four the concatonated and alphabetically sorted gene IDs joined
#' by '_'. Default is 'chi.paranome.ks'
#' @param na.dist.as.num A large value to be used for gene pairs where, because
#' of no significant similarity, 'dist.tbl' has no entry.
#'
#' @return A list with two named entries: 'cluster' an instance of ape::phylo
#' representing the families binary distance tree, and 'weighted.distances' a
#' named numeric vector with names 'node_i' and values the weighted distance
#' computed for the corresponding node.
#' @export
weightedDistsForFamily <- function(fam.name, fam.tbl = chi.paralogous.fams, dist.tbl = chi.paranome.ks, 
    na.dist.as.num = 1e+06) {
    fam.dists <- as.dist(familyDist(fam.name, fam.tbl, dist.tbl, na.dist.as.num))
    fam.tr <- as.phylo(hclust(fam.dists, method = "single"))
    fam.ks <- weightedDistsForTree(fam.tr, dist.tbl)
    list(cluster = fam.tr, weighted.distances = fam.ks)
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

#' Generates the tree and control file to run PAML's codeml in order to compute
#' the Ks of 'gene.pair'
#'
#' @param gene.pair a character of length 2 holding the identifiers of the
#' pair's genes.
#' @param cds.msa.path the file path to the gene pair's codon sequences
#' alignment
#' @param ks.out.path the file path to the gene pair's codeml output file,
#' defaut is generated from 'cds.msa.path' substituting '_CDS_MSA.fasta' with
#' '_codeml_out.txt'
#' @param cds.codeml.path the file path to the gene pair's codeml input file,
#' defaut is generated from 'cds.msa.path' substituting '_CDS_MSA.fasta' with
#' '_codeml_in.cnt'
#' @param tree.path the file path to the gene pair's tree file, defaut is
#' generated from 'cds.msa.path' substituting '_CDS_MSA.fasta' with
#' '_tree.newick'. The tree will be the trivial one: '(A,B);'.
#'
#' @return A list with two string entries: 'in' the path to the control file
#' used as input to PAML's codeml, and 'out' path to the output file generated
#' by PAML's codeml.
#' @export
renderCodemlInputs <- function(gene.pair, cds.msa.path, ks.out.path = sub("_CDS_MSA.fasta", 
    "_codeml_out.txt", cds.msa.path, fixed = TRUE), cds.codeml.path = sub("_CDS_MSA.fasta", 
    "_codeml_in.cnt", cds.msa.path, fixed = TRUE), tree.path = sub("_CDS_MSA.fasta", 
    "_tree.newick", cds.msa.path, fixed = TRUE)) {
    writeLines(paste("(", gene.pair[[1]], ",", gene.pair[[2]], ");", sep = ""), con = tree.path)
    brew(text = codeml.tmpl, output = cds.codeml.path)
    list(`in` = cds.codeml.path, out = ks.out.path)
}

#' Body of inline Rcpp function 'familyDistCpp'
#' @export
fun.bd <- "CharacterVector genes(sGenes);\nNumericVector genePairKs(sGenePairsKs);\nCharacterVector genePairs = genePairKs.names();\nint n( genes.size() );\nNumericMatrix distMtrx(n,n);\nstd::vector<std::string> pair(2);\nstd::string sPair = \"\";\nNumericVector defDist(sDefDist);\nfor ( int i=1; i<n; ++i ) {\n  for ( int j=0; j<i; ++j ) {\n    pair[0] = genes(i);\n    pair[1] = genes(j);\n    std::sort( pair.begin(), pair.end() );\n    sPair = pair[0] + \"_\" + pair[1];\n    bool present =  std::find(genePairs.begin(), genePairs.end(), sPair.c_str()) != genePairs.end();\n    if ( present ) {\n      distMtrx(i, j) = as<double>( genePairKs( sPair ) );\n    } else {\n      distMtrx(i, j) = defDist(0);\n    }\n  }\n}\n\nrownames(distMtrx) = genes;\ncolnames(distMtrx) = genes;\nreturn( wrap( distMtrx ) );"

#' Template for brew to generate the control file used as input for PAML's
#' codeml.
#' @export
codeml.tmpl <- "seqfile = <%= cds.msa.path %>\n    treefile = <%= tree.path %>\n     outfile = <%= ks.out.path %>\n       noisy = 0\n     verbose = 0\n     runmode = -2\n   cleandata = 1\n     seqtype = 1\n   CodonFreq = 2\n       model = 2\n     NSsites = 0\n       icode = 0\n       Mgene = 0\n   fix_kappa = 0\n       kappa = 2\n   fix_omega = 0\n       omega = 1\n   fix_alpha = 1\n       alpha = .0\n      Malpha = 0\n       ncatG = 4\n       clock = 0\n       getSE = 0\nRateAncestor = 0\n      method = 0" 
