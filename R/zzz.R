.onLoad <- function(libname = find.package("paranomeKsR"), pkgname = "paranomeKsR") {
    data("chiParalogousFamilies", package = "paranomeKsR")
    data("chiParanomeKs", package = "paranomeKsR")
    assign("familyDistCpp", cxxfunction(signature(sGenes = "character", sGenePairsKs = "numeric", 
        sDefDist = "numeric"), body = fun.bd, plugin = "Rcpp"), envir = globalenv())
} 
