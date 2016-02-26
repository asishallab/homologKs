#include "rcpp_functions.h"
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
//' @param sGenePairsKs a NumericVector in which the keys are gene pairs and
//' the values are their pre-computed distances (Ks). Gene pairs are encoded as
//' their concatonated identifiers in alphabetical order separated bu "_", e.g.
//' "A_B" and not "B_A".
//'
//' @param sDefDist a numeric vector holding a single value which is to be used
//' as the default distance for pairs for which no pre-computed distance can be
//' found in 'genePairsKs'. This should be a very large value, because
//' pre-computed distances are lacking for pairs with no significant sequence
//' similarity.
//' 
//' @export
SEXP familyDistCpp( SEXP sGenes, SEXP sGenePairsKs, SEXP sDefDist ){
  BEGIN_RCPP

     CharacterVector genes(sGenes);
     NumericVector genePairKs(sGenePairsKs);
     CharacterVector genePairs = genePairKs.names();
     int n( genes.size() );
     NumericMatrix distMtrx(n,n);
     std::vector<std::string> pair(2);
     std::string sPair = "";
     NumericVector defDist(sDefDist);
     for ( int i=1; i<n; ++i ) {
       for ( int j=0; j<i; ++j ) {
         pair[0] = genes(i);
         pair[1] = genes(j);
         std::sort( pair.begin(), pair.end() );
         sPair = pair[0] + "_" + pair[1];
         Rcpp::Rcout << sPair << std::endl;
         bool present =  std::find(genePairs.begin(), genePairs.end(), sPair.c_str()) != genePairs.end();
         if ( present ) {
           distMtrx(i, j) = as<double>( genePairKs( sPair ) );
         } else {
           distMtrx(i, j) = defDist(0);
         }
       }
     }

     rownames(distMtrx) = genes;
     colnames(distMtrx) = genes;
     return( wrap( distMtrx ) );

  END_RCPP
}
