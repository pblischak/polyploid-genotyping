#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simConst(NumericVector f, int nGen){
  
  NumericVector res(f.size());
  for(int l = 0; l < f.size(); l++)
    res[l] = f[l];
  
  int popSize = 125000;
  int tmp = 0;
  
  for(int g = 0; g < nGen; g++){
    for(int l = 0; l < f.size(); l++){
      tmp = ::Rf_rbinom(popSize, res[l]);
      res[l] = tmp / (double) popSize;
    }
  }
  
  return res;
  
}

// [[Rcpp::export]]
NumericVector simGrowth(NumericVector f, int nGen){
  
  NumericVector res(f.size());
  for(int l = 0; l < f.size(); l++)
    res[l] = f[l];
  
  int popSize = 500;
  int tmp = 0;
  
  for(int g = 0; g < nGen; g++){
    popSize += 5;
    for(int l = 0; l < f.size(); l++){
      tmp = ::Rf_rbinom(popSize, res[l]);
      res[l] = tmp / (double) popSize;
    }
  }
  
  return res;
  
}
