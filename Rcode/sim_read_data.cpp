// Script written by P. D. Blischak (2016)
// SNP genotyping and analysis of genetic diversity in polyploids
// Paul D. Blischak, Laura S. Kubatko, Andrea D. Wolfe
/**************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;

/* Script will simulate total and reference read counts, as well as
   error rates and then returns them to R as a List. */

// [[Rcpp::export]]
List simReads(IntegerMatrix g, double coverage, int ploidy) {
  
  IntegerMatrix totReadMat(g.nrow(), g.ncol());
  IntegerMatrix refReadMat(g.nrow(), g.ncol());
  NumericVector errorRates(g.ncol());
  double gEpsilon = 0.0;
  
  for(int l = 0; l < g.ncol(); l++){
    errorRates[l] = ::Rf_rbeta(1, 200);
  }
  
  for(int i = 0; i < g.nrow(); i++){
    for(int l = 0; l < g.ncol(); l++){
      
      if(g(i,l) != -9){
        
        totReadMat(i,l) = 0;
        while(totReadMat(i,l) == 0){
          totReadMat(i,l) = ::Rf_rpois(coverage);
        }
        
        gEpsilon = (g(i,l) / (double) ploidy) * (1.0 - errorRates[l]) + (1.0 - g(i,l) / (double) ploidy) * errorRates[l];
        
        if(gEpsilon == 0)
          Rcout << gEpsilon;
        
        refReadMat(i,l) = ::Rf_rbinom(totReadMat(i,l), gEpsilon);
        
      } else {
        totReadMat(i,l) = -9;
        refReadMat(i,l) = -9;
      }
      
    }
  }
  
  return List::create(Named("totalReads") = totReadMat,
                      Named("referenceReads") = refReadMat,
                      Named("error") = errorRates);
}
