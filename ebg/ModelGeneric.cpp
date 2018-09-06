#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "ModelGeneric.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

// Declare static member variables that are needed to pass member functions to Brent's method code
std::vector<double> ModelGeneric::_freqs, ModelGeneric::_errRates, ModelGeneric::_perSiteLogLik,
                    ModelGeneric::_theta0Loci, ModelGeneric::_theta1Loci, ModelGeneric::_theta2Loci, ModelGeneric::_thetaPrimeLoci,
                    ModelGeneric::_vLoci, ModelGeneric::_rLoci, ModelGeneric::_prevFreqs;

std::vector< std::vector< std::vector<double> > > ModelGeneric::_gExp, ModelGeneric::_gLiks;
std::vector< std::vector< std::vector< std::vector<double> > > > ModelGeneric::_gExp4D;
//std::vector<int> ModelGeneric::_totReads, ModelGeneric::_refReads;
std::vector< std::vector<int> > ModelGeneric::_totReads, ModelGeneric::_refReads;
double ModelGeneric::_tol = 1.0e-10, ModelGeneric::_currLogLik = 0.0, ModelGeneric::_prevLogLik = 0.0, ModelGeneric::_stopVal = 1.0e-5,
       ModelGeneric::_alpha = -1.0;
int ModelGeneric::_nInd = -999, ModelGeneric::_nLoci = -999, ModelGeneric::_ploidy = -999, ModelGeneric::_maxIters = 100, ModelGeneric::_currLoc = 0;
std::string ModelGeneric::_totalReadsFile = "none", ModelGeneric::_refReadsFile = "none", ModelGeneric::_errorRatesFile = "none", ModelGeneric::_prefix;
bool ModelGeneric::_quiet = 0, ModelGeneric::_brent = 0;

void ModelGeneric::getData(){

  _totReads.resize(_nInd);
  _refReads.resize(_nInd);
  for(int i = 0; i < _nInd; i++){
    _totReads[i].resize(_nLoci);
    _refReads[i].resize(_nLoci);
  }
  
  int readVal = -999, rowCounter = 0, colCounter = 0, count = 0;
  double errorVal = 0.0;

  std::ifstream _totalReadsStream(_totalReadsFile);
  std::ifstream _refReadsStream(_refReadsFile);
  std::ifstream _errorRatesStream(_errorRatesFile);

  if(_totalReadsStream.is_open()){
    for(int i = 0; i < _nInd; ++i){
      for(int l = 0; l < _nLoci; ++l){
        if(_totalReadsStream >> readVal){
          if(readVal == 0){
            std::cerr << "\n  Values of 0 in the total reads matrix are not allowed.\n"
                      << "  Missing data should be marked as -9. (" << i + 1 << ","  << l + 1
                      << ")\n\n";
            //exit(EXIT_FAILURE);
          } else {
            _totReads[i][l] = readVal;
            //std::cerr << readVal << std::endl;
          }
        } else {
          std::cerr << "Couldn't be read: " << i << "," << l << std::endl;
        }
      }
    }
  } else {
    std::cerr << "Could not open total reads file: " << _totalReadsFile << std::endl;
    exit(EXIT_FAILURE);
  }
  
  //colCounter = 0, rowCounter = 0;
  if(_refReadsStream.is_open()){
    for(int i = 0; i < _nInd; ++i){
      for(int l = 0; l < _nLoci; ++l){
        _refReadsStream >> _refReads[i][l];
        if(_refReads[i][l] == -9 && _totReads[i][l] != -9){
          _totReads[i][l] = -9;
        }
      }
    }
    /*while(_refReadsStream >> readVal){
      _refReads[rowCounter][colCounter] = readVal;
      colCounter++;
      
      if(colCounter == _nLoci){
        rowCounter++;
        //std::cerr << rowCounter << "," << colCounter << std::endl;
        colCounter = 0;
      }
    }*/
  } else {
    std::cerr << "Could not open ALT allele reads file: " << _refReadsFile << std::endl;
    exit(EXIT_FAILURE);
  }

  if(_errorRatesStream.is_open()){
    while(_errorRatesStream >> errorVal){
      count++;
      _errRates.push_back(errorVal);
      //std::cerr << count << std::endl;
    }
  } else {
    std::cerr << "Could not open error rates file: " << _errorRatesFile << std::endl;
    exit(EXIT_FAILURE);
  }


}

void ModelGeneric::checkInput(){

  if(!_quiet)
    std::cerr << "Checking input data...    ";

  double total = (double) _nInd * _nLoci;
  int missing = 0;

  if(_totReads.size() * _totReads[0].size() != (size_t) (_nInd * _nLoci)){
    std::cerr << "The size of the total reads matrix is incorrect (" << _totReads.size() << " not equal to num-ind * num-loci)." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(_refReads.size() * _refReads[0].size() != (size_t) (_nInd * _nLoci)){
    std::cerr << "The size of the ALT allele reads matrix is incorrect (" <<  _refReads.size() << " not equal to num-ind * num-loci)." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(_errRates.size() != (size_t) _nLoci){
    std::cerr << "The size of the error rates vector is incorrect (" << _errRates.size() << " not equal to num-loci)." << std::endl;
    exit(EXIT_FAILURE);
  }

  for(int i = 0; i < _nInd; i++){
    for(int l = 0; l < _nLoci; l++){
      if(_totReads[i][l] == MISSING)
        missing++;
    }
  }

  if(!_quiet)
    std::cerr << "            Good (" << missing/total * 100.0 <<" % missing data)"<< std::endl;

}

/* Code for Brent's method.*/
double ModelGeneric::local_min(double a, double b, double t, double (&f) (double x), double &x){
  double c;
  double d = 0.0;
  double e;
  double eps;
  double fu;
  double fv;
  double fw;
  double fx;
  double m;
  double p;
  double q;
  double r;
  double sa;
  double sb;
  double t2;
  double tol;
  double u;
  double v;
  double w;
//
//  C is the square of the inverse of the golden ratio.
//
  c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );

  eps = sqrt ( r8_epsilon ( ) );

  sa = a;
  sb = b;
  x = sa + c * ( b - a );
  w = x;
  v = w;
  e = 0.0;
  fx = f(x);
  fw = fx;
  fv = fw;

  for ( ; ; ){
    m = 0.5 * ( sa + sb ) ;
    tol = eps * r8_abs ( x ) + t;
    t2 = 2.0 * tol;
//
//  Check the stopping criterion.
//
    if ( r8_abs ( x - m ) <= t2 - 0.5 * ( sb - sa ) ){
      break;
    }
//
//  Fit a parabola.
//
    r = 0.0;
    q = r;
    p = q;

    if ( tol < r8_abs ( e ) ){
      r = ( x - w ) * ( fx - fv );
      q = ( x - v ) * ( fx - fw );
      p = ( x - v ) * q - ( x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q ){
        p = - p;
      }
      q = r8_abs ( q );
      r = e;
      e = d;
    }

    if ( r8_abs ( p ) < r8_abs ( 0.5 * q * r ) &&
         q * ( sa - x ) < p &&
         p < q * ( sb - x ) ){
//
//  Take the parabolic interpolation step.
//
      d = p / q;
      u = x + d;
//
//  F must not be evaluated too close to A or B.
//
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 ){
        if ( x < m ){
          d = tol;
        } else {
          d = - tol;
        }
      }
    }
//
//  A golden-section step.
//
    else {
      if ( x < m ) {
        e = sb - x;
      } else {
        e = sa - x;
      }
      d = c * e;
    }
//
//  F must not be evaluated too close to X.
//
    if ( tol <= r8_abs ( d ) ) {
      u = x + d;
    } else if ( 0.0 < d ){
      u = x + tol;
    } else {
      u = x - tol;
    }

    fu = f(u);
//
//  Update A, B, V, W, and X.
//
    if ( fu <= fx ){
      if ( u < x ){
        sb = x;
      } else {
        sa = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    } else {
      if ( u < x ) {
        sa = u;
      } else {
        sb = u;
      }

      if ( fu <= fw || w == x ){
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else if ( fu <= fv || v == x || v == w ){
        v = u;
        fv = fu;
      }
    }
  }
  return fx;
}

double ModelGeneric::norm(std::vector<double> &v){

  double tmp = 0.0;

  for(size_t i = 0; i < v.size(); i++){
    tmp += pow(v[i],2);
  }

  return sqrt(tmp);

}
