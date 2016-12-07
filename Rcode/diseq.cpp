#include <Rcpp.h>
using namespace Rcpp;

double diseqLogLik(NumericVector gg, double freq, double phi,
                   int ind, int loc, int ploidy);
double lnBetaBinomialPdf(int size, int x, double alpha, double beta);
IntegerMatrix simDiseqGenos(NumericVector pp, NumericVector ff, SEXP ind, SEXP pldy);
int nonunif_int(const Rcpp::IntegerVector &vals, const Rcpp::NumericVector &probs);
int rbetabinom_cpp(const double &alpha, const double &beta, const int &nn);

NumericMatrix diseqCpp(NumericVector freqs, NumericVector phi,
                       IntegerMatrix tot, IntegerMatrix ref,
                       NumericVector err, SEXP pldy){

  int pos_lia, ploidy = as<int>(pldy), loci = tot.ncol(), ind = tot.nrow();
  double gEpsilon = 0.0;
  NumericMatrix res(freqs.size() * phi.size(), 2 + loci);
  NumericVector gLiks(ind * loci * (ploidy+1), 0.0);

  for(int l = 0; l < loci; l++){
    for(int i = 0; i < ind; i++){
      for(int a = 0; a <= ploidy; a++){

        pos_lia = l*ind*(ploidy+1) + i*(ploidy+1) + a;

        gEpsilon = (a / (double) ploidy) * (1 - err[l]) +
                   (1 - a / (double) ploidy) * err[l];

        gLiks[pos_lia] = ::Rf_dbinom(ref(i,l), tot(i,l), gEpsilon, 1);

      }
    }
  }

  int row = 0;
  for(int f = 0; f < freqs.size(); f++){
    for(int p = 0; p < phi.size(); p++){

      res(row, 0) = freqs[f];
      res(row, 1) = phi[p];

      for(int l = 0; l < loci; l++){
        res(row, 2 + l) = diseqLogLik(gLiks, freqs[f], phi[p], ind, l, ploidy);
      }

      row++;

      Rcout << row << "\n";
    }
  }


  return res;

}

double diseqLogLik(NumericVector gg, double freq, double phi,
                   int ind, int loc, int ploidy){

  int pos_lia;
  double logLik = 0.0, indLikSum = 0.0;
  NumericVector indLikVec(ploidy+1, 0.0);

  for(int i = 0; i < ind; i++){
    for(int a = 0; a <= ploidy; a++){

      pos_lia = loc*ind*(ploidy+1) + i*(ploidy+1) + a;

      indLikVec[a] = exp(gg[pos_lia] + lnBetaBinomialPdf(ploidy, a, freq * phi, (1 - freq) * phi));

    }

    indLikSum = sum(indLikVec);
    logLik += log(indLikSum);

  }

  return logLik;

}

double lnBetaBinomialPdf(int size, int x, double alpha, double beta){

  double res = ::Rf_lchoose(size, x)
               + ::Rf_lbeta(x + alpha, size - x + beta)
               - ::Rf_lbeta(alpha, beta);

  return res;

}

// [[Rcpp::export]]
IntegerMatrix simDiseqGenos(NumericVector pp, NumericVector ff, SEXP ind, SEXP pldy){

  int ploidy = as<int>(pldy);
  int nind = as<int>(ind);
  int nloci = pp.length();
  NumericVector phi = 1.0/ff - 1.0;
  IntegerMatrix gg(nind, nloci);
  double aa = 0;
  double bb = 0;

  for(int i = 0; i < nind; i++){
    for(int l = 0; l < nloci; l++){

      if(pp[l] == 1){
        gg(i,l) = ploidy;
      } else if(pp[l] == 0){
        gg(i,l) = 0;
      } else {

        aa = pp[l]*phi[i];
        bb = (1.0 - pp[l])*phi[i];

        gg(i,l) = rbetabinom_cpp(aa, bb, ploidy);

      }
    }
  }

  return gg;

}

int nonunif_int(const Rcpp::IntegerVector &vals, const Rcpp::NumericVector &probs){

  if(vals.size() != probs.size()){
    throw std::range_error("Vectors are not of the same length.");
  }

  double prob_sum = sum(probs);
  Rcpp::NumericVector normd_probs = probs/prob_sum;
  Rcpp::NumericVector cummulative_probs = Rcpp::cumsum(normd_probs);

  // Get random number between 0 and 1
  double samp = ::Rf_runif(0,1);

  int res = -999;

  if(samp < cummulative_probs[0]){
    res = vals[0];
  } else {

    for(int i = 1; i < vals.size(); i++){

      if(samp >= cummulative_probs[i-1] && samp < cummulative_probs[i]){
        res = vals[i];
        break;
      }

    }
  }

  if(res == -999){
    Rcpp::stop("Invalid random integer generated");
  }

  return res;

}

int rbetabinom_cpp(const double &alpha, const double &beta, const int &nn){

  Rcpp::IntegerVector k(nn + 1);
  Rcpp::NumericVector p(nn + 1);
  double tmp_prob;
  int res = -999;

  for(int i = 0; i <= nn; i++){
    k[i] = i;
    tmp_prob = ::Rf_lchoose(nn, i) + ::Rf_lbeta(i + alpha, nn - i + beta) - ::Rf_lbeta(alpha, beta);
    p[i] = exp(tmp_prob);
  }

  res = nonunif_int(k, p);

  if(res == -999){
    Rcpp::stop("Invalid rbetabinom integer generated.");
  }

  return res;

}
