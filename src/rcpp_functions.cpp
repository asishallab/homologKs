#include "rcpp_functions.h"
using namespace Rcpp ;

SEXP familyDistCpp( SEXP sGenes, SEXP sGenePairsKs, SEXP sDefDist ){
  BEGIN_RCPP

     CharacterVector genes(sGenes);
     List genePairKs(sGenePairsKs);
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
