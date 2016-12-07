#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simFreq(int L, double lowBound, double highBound, double a, double b){

  double tmpFreq;
  NumericVector res(L);

  for(int l = 0; l < L; l++){

    tmpFreq = -1.0;

    while(tmpFreq < lowBound || tmpFreq > highBound)
      tmpFreq = ::Rf_rbeta(a, b);

    res[l] = tmpFreq;

  }

  return res;

}
