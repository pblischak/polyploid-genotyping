#ifndef MODEL_GENERIC_HPP
#define MODEL_GENERIC_HPP

class ModelGeneric {

public:

  ModelGeneric(){};
  ~ModelGeneric(){};

protected:

  static std::vector<double> _freqs, _prevFreqs, _errRates, _perSiteLogLik,
                             _theta0Loci, _theta1Loci, _theta2Loci,
                             _rLoci, _vLoci, _thetaPrimeLoci;
  static std::vector< std::vector< std::vector<double> > > _gLiks, _gExp;
  static std::vector< std::vector< std::vector< std::vector<double> > > > _gExp4D;
  // static std::vector<int> _totReads, _refReads;
  static std::vector< std::vector<int> > _totReads, _refReads;
  static double _tol, _currLogLik, _prevLogLik, _stopVal, _alpha;
  static int _nInd, _nLoci, _ploidy, _maxIters, _currLoc;
  static std::string _totalReadsFile, _refReadsFile, _errorRatesFile, _prefix;
  static bool _quiet, _brent;

  // Code for Brent's method modified from John Burkardt's code (brent.hpp, brent.cpp).
  double local_min(double a, double b, double t, double (&f) (double x), double &x);

  // Member functions for data parsing
  void getData();
  void checkInput();

  // Member utility functions
  double f_epsilon(int x, int y, double e);
  double r8_abs(double x); // Also from John Burkardt's BRENT code
  double r8_epsilon(); // Also from John Burkardt's BRENT code
  int gMax(std::vector<double> &v); // return the index of maximum entry in a vector
  double norm(std::vector<double> &v);

};

inline double ModelGeneric::f_epsilon(int x, int y, double e){
  return (x / (double) y) * (1.0 - e) + (1.0 - (x / (double) y)) * e;
}

inline double ModelGeneric::r8_epsilon(){
  const double value = 2.220446049250313E-016;
  return value;
}

inline double ModelGeneric::r8_abs(double x){
  double value;

  if ( 0.0 <= x ){
    value = x;
  } else {
    value = - x;
  }

  return value;
}

inline int ModelGeneric::gMax(std::vector<double> &v){

  if(v.size() == 0)
    return -1;

  double max_val = v[0];
  int max_index = 0;

  for(int i = 1; i < (int) v.size(); i++){
    if(v[i] > max_val){
      max_val = v[i];
      max_index = i;
    }
  }

  return max_index;

}

#endif //MODEL_GENERIC_HPP
