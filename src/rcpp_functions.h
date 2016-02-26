//Includes/namespaces

#include <Rcpp.h>

using namespace Rcpp ;

//' @title
//' familyDistCpp
//' @description
//' Generates a distance matrix for each pair of genes generated from 'sGenes'
//' using the distances (Ks) values held in 'sGenePairsKs'. If a pair is not
//' found in 'sGenePairsKs' the default distance (Ks) 'sDefDist' is used instead.
//' 
//' @param sGenes a character vector of gene identifiers 
//' 
//' @param sGenePairsKs a List in which the keys are gene pairs and the values
//' are their pre-computed distances (Ks). Gene pairs are encoded as their
//' concatonated identifiers in alphabetical order separated bu "_", e.g. "A_B"
//' and not "B_A".
//'
//' @param sDefDist a numeric vector holding a single value which is to be used
//' as the default distance for pairs for which no pre-computed distance can be
//' found in 'genePairsKs'. This should be a very large value, because
//' pre-computed distances are lacking for pairs with no significant sequence
//' similarity.
//' 
//' @export
// [[Rcpp::export]]
RcppExport SEXP familyDistCpp( SEXP sGenes, SEXP sGenePairsKs, SEXP sDefDist );
