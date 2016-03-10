.onLoad <- function(libname = find.package("paranomeKsR"), pkgname = "paranomeKsR") {
    data("chiParalogousFamilies", package = "paranomeKsR")
    data("chiParanomeKs", package = "paranomeKsR")
    data("chiSequences", package = "paranomeKsR")
    data("chiFamilyBasedWeightedKs", package = "paranomeKsR")
    data("chiPmeiFamiliesAndGenes", package = "paranomeKsR")
    # Lamentably compilation with proper Rcpp files somehow failes, usage of the
    # inline package works, however:
    assign("familyDistCpp", cxxfunction(signature(sGenes = "character", sGenePairsKs = "numeric", 
        sDefDist = "numeric"), body = fun.bd, plugin = "Rcpp"), envir = globalenv())
} 
