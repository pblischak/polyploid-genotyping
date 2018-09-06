#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
//#include <dlib/optimization.h>

#include "ModelGeneric.hpp"
#include "ModelDiseq.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

// Declare static member variable
std::vector<double> ModelDiseq::_phi, ModelDiseq::_perIndLogLik;
int ModelDiseq::_currInd;
std::vector<bool> ModelDiseq::_convergedInd, ModelDiseq::_convergedLoci;

//using namespace dlib;

ModelDiseq::ModelDiseq(int ac, char* av[]){
  _prefix = "diseq";
  for(int i = 2; i < ac; i++){
    if(strcmp(av[i],"--num-ind") == 0 || strcmp(av[i], "-n") == 0){
      _nInd = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--num-loci") == 0 || strcmp(av[i], "-l") == 0){
      _nLoci = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--ploidy") == 0 || strcmp(av[i], "-p") == 0){
      _ploidy = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--total-reads") == 0 || strcmp(av[i], "-t") == 0){
      _totalReadsFile = av[i + 1];
    } else if(strcmp(av[i],"--alt-reads") == 0 || strcmp(av[i], "-a") == 0){
      _refReadsFile = av[i + 1];
    } else if(strcmp(av[i],"--prefix") == 0){
      _prefix = av[i + 1];
    } else if(strcmp(av[i],"--error-rates") == 0 || strcmp(av[i], "-e") == 0) {
      _errorRatesFile = av[i + 1];
    } else if(strcmp(av[i],"--tol") == 0){
      _tol = atof(av[i + 1]);
    } else if(strcmp(av[i],"--iters") == 0){
      _maxIters = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--stop") == 0){
      _stopVal = atof(av[i+1]);
    } else if(strcmp(av[i],"--quiet") == 0 || strcmp(av[i], "-q") == 0){
      _quiet = 1;
    }
  }
  // Check the input for necessary arguments
  checkCommandLine();
  // Read in the data
  getData();
  // Check data
  checkInput();
  // Initialize the parameters for the ECM algorithm
  initParams();
}

void ModelDiseq::ecm(){
  double currLogLik = 0.0, prevLogLik = 0.0;
  bool convergence = 0;
  for(int j = 0; j < _maxIters; j++){
    prevLogLik = currLogLik;
    if(!_quiet)
      std::cerr << "Step: " << j + 1 << "\t";

    eStep();
    mStep();
    currLogLik = calcLogLik();
    if(!_quiet)
      std::cerr << std::setprecision(16) << "logLik: " << currLogLik << "\t" << "Diff: " << currLogLik - prevLogLik << std::endl;

    checkConvergence(convergence);
    if(convergence|| fabs(currLogLik - prevLogLik) < 1e-8)
      break;
  }
}

void ModelDiseq::checkConvergence(bool &con){
  int convergedLociCount = 0, convergedIndCount = 0;
  for(int l = 0; l < _nLoci; l++){
    if(_convergedLoci[l])
      convergedLociCount++;
  }
  for(int i = 0; i < _nInd; i++){
    if(_convergedInd[i])
      convergedIndCount++;
  }
  if(convergedLociCount == _nLoci && convergedIndCount == _nInd)
    con = 1;
}

void ModelDiseq::printOutput(){
  std::cerr << "Writing output files...";
  std::string genosFile = _prefix + "-genos.txt";
  std::string freqsFile = _prefix + "-freqs.txt";
  std::string phiFile = _prefix + "-F.txt";
  std::string plFile = _prefix + "-PL.txt";
  std::ofstream genosStream;
  genosStream.open(genosFile, std::ios::out);
  std::ofstream freqsStream;
  freqsStream.open(freqsFile, std::ios::out);
  std::ofstream phiStream;
  phiStream.open(phiFile, std::ios::out);
  std::ofstream plStream;
  plStream.open(plFile, std::ios::out);
  std::vector<double> tmp_val(_ploidy+1, 0.0), g_post_prob(_ploidy+1, 0.0);
  double tmp_val_sum = 0.0, phiConverted = 0.0, tmp_PL = 0.0;
  std::vector< std::vector<int> > genotypes;
  genotypes.resize(_nInd);
  for(int i = 0; i < _nInd; i++){
    genotypes[i].resize(_nLoci, -9);
  }
  int i3d = 0, i2d = 0;
  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;
      // Catch missing data and skip
      if(_totReads[i][l] == missing)
        continue;

      tmp_val_sum = 0.0;
      phiConverted = 1.0 / _phi[i] - 1.0;
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        tmp_val[a] = _gLiks[l][i][a] + r->lnBetaBinomPdf(_ploidy, a, _freqs[l] * phiConverted, (1 - _freqs[l]) * phiConverted);
        tmp_val_sum += exp(tmp_val[a]);
      }
      for(int a = 0; a <= _ploidy; a++){
        g_post_prob[a] = exp(tmp_val[a])/tmp_val_sum;
        tmp_PL = -10 * log10(g_post_prob[a]);
        if(tmp_PL == -0){
          plStream << -1 * tmp_PL;
        } else {
          plStream << tmp_PL;
        }
        if(a >= 0 && a < _ploidy){
          plStream << ",";
        }
      }
      plStream << "\t";
      genotypes[i][l] = gMax(g_post_prob);
    }
    plStream << std::endl;
  }
  for(int i = 0; i < _nInd; i++){
    for(int l = 0; l < _nLoci; l++){
      i2d = i * _nLoci + l;
      genosStream << genotypes[i][l] << "\t";
    }
    genosStream << std::endl;
  }
  for(int l = 0; l < _nLoci; l++)
    freqsStream << _freqs[l] << std::endl;

  for(int i = 0; i < _nInd; i++)
    phiStream << _phi[i] << std::endl;

  std::cerr << "Done." << std::endl;
}

/*******************************
 ** Private member functions  **
 *******************************/

void ModelDiseq::checkCommandLine(){
  if(!_quiet)
    std::cerr << "\nChecking command line input...        ";

  int errorCaught = 0;
  if(_nInd < 0){
    std::cerr << "\nMissing or invalid option for -n [--num-ind].\n";
    errorCaught++;
  }
  if(_nLoci < 0){
    std::cerr << "\nMissing or invalid option for -l [--num-loci].\n";
    errorCaught++;
  }
  if(_ploidy < 0){
    std::cerr << "\nMissing or invalid option for -p [--ploidy].\n";
    errorCaught++;
  }
  if(strcmp(_totalReadsFile.c_str(), "none") == 0){
    std::cerr << "\nMissing or invalid option for -t [--total-reads].\n";
    errorCaught++;
  }
  if(strcmp(_refReadsFile.c_str(), "none") == 0){
    std::cerr << "\nMissing or invalid option for -a [--alt-reads].\n";
    errorCaught++;
  }
  if(strcmp(_errorRatesFile.c_str(), "none") == 0) {
    std::cerr << "\nMissing or invalid option for -e [--error-rates].\n";
    errorCaught++;
  }
  if(errorCaught > 0){
    diseqUsage();
    exit(EXIT_FAILURE);
  }
  if(!_quiet)
    std::cerr << "Good" << std::endl;
}

void ModelDiseq::initParams(){
  _gExp.resize(_nLoci);
  _gLiks.resize(_nLoci);
  for(int l = 0; l < _nLoci; l++){
    _gExp[l].resize(_nInd);
    _gLiks[l].resize(_nInd);
    for(int i = 0; i < _nInd; i++){
      _gExp[l][i].resize(_ploidy + 1);
      _gLiks[l][i].resize(_ploidy + 1);
    }
  }
  _freqs.resize(_nLoci);
  _perSiteLogLik.resize(_nLoci);
  _perIndLogLik.resize(_nInd);
  _phi.resize(_nInd);
  _convergedInd.resize(_nInd);
  _convergedLoci.resize(_nLoci);
  int i2d = 0, i3d = 0; // indices for entries in 2D and 3D arrays stored as vectors
  double gEpsilon = 0.0; // error corrected genotype
  for(int l  = 0; l < _nLoci; l++){
    _freqs[l] = r->uniformRv();
    _convergedLoci[l] = 0;
  }
  for(int i = 0; i < _nInd; i++){
    _phi[i] = 0.1 * r->uniformRv() + 0.01;
    _convergedInd[i] = 0;
  }
  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        if(_totReads[i][l] == missing){
          _gLiks[l][i][a] = bad_lik;
        } else {
          gEpsilon = f_epsilon(a, _ploidy, _errRates[l]);
          _gLiks[l][i][a] = r->lnBinomPdf(_totReads[i][l], _refReads[i][l], gEpsilon);
          //std::cerr << gEpsilon << "    " << _gLiks[i3d] << std::endl;
        }
      }
    }
  }
}

double ModelDiseq::calcFreqLogLik(double x){
  int i3d = 0;
  double logLik = 0.0, alpha = 0.0, beta = 0.0, phiConverted = 0;
  for(int i = 0; i < _nInd; i++){
    // Catch missing data and skip
    if(_totReads[i][_currLoc] == missing)
      continue;

    phiConverted = 1.0 / _phi[i] - 1.0;
    if(phiConverted > 1000.0)
      phiConverted = 1000.0;

    alpha = x * phiConverted;
    beta  = (1 - x) * phiConverted;
    for(int a = 0; a <= _ploidy; a++){
      i3d = _currLoc * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
      logLik += _gExp[_currLoc][i][a] * (_gLiks[_currLoc][i][a] + r->lnBetaBinomPdf(_ploidy, a, alpha, beta));
      //std::cerr << _currLoc+1 << "\t" << i+1 << "\t" << a << "\t" << _gExp[i3d] << "\t" << r->lnBetaBinomPdf(_ploidy, a, alpha, beta) << std::endl;
    }
  }
  //if(_currLoc == 120)
  //  std::cerr << _currLoc+1 << "\t" << logLik << "\t" << x << std::endl;
  return -logLik;
}

double ModelDiseq::calcPhiLogLik(double x){
  int i3d = 0;
  double logLik = 0.0, alpha = 0.0, beta = 0.0, phiConverted = 1.0 / x - 1.0;
  for(int l = 0; l < _nLoci; l++){
    // Catch missing data and skip
    if(_totReads[_currInd][l] == missing)
      continue;

    alpha = _freqs[l] * phiConverted;
    beta  = (1 - _freqs[l]) * phiConverted;
    for(int a = 0; a <= _ploidy; a++){
      i3d = l * _nInd * (_ploidy + 1) + _currInd * (_ploidy + 1) + a;
      logLik += _gExp[l][_currInd][a] * (_gLiks[l][_currInd][a] + r->lnBetaBinomPdf(_ploidy, a, alpha, beta));
    }
  }
  return -logLik;
}

void ModelDiseq::eStep(){
  std::vector<double> tmp_exp(_ploidy + 1, 0.0);
  double tmp_exp_sum = 0.0, phiConverted = 0.0;
  int i3d = 0; // index for 3D array stored as vector
  bool print = 0;
  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      // Catch missing data and skip
      if(_totReads[i][l] == missing)
        continue;

      tmp_exp_sum = 0.0;
      phiConverted = 1.0 / _phi[i] - 1.0;
      if(phiConverted > 1000)
        phiConverted = 1000.0;

      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        tmp_exp[a] = _gLiks[l][i][a] + r->lnBetaBinomPdf(_ploidy, a, _freqs[l] * phiConverted, (1.0 - _freqs[l]) * phiConverted);
        //std::cerr << tmp_exp[a] << "    " << _gLiks[i3d] << "    " << r->lnBetaBinomPdf(_ploidy, a, _freqs[l] * phiConverted, (1 - _freqs[l]) * phiConverted) << std::endl;
        tmp_exp_sum += exp(tmp_exp[a]);
      }
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        _gExp[l][i][a] = exp(tmp_exp[a]) / tmp_exp_sum;
        if(print){
        std::cout << i+1 << "    " << l+1 << "    "
                  << a << "    " << r->betaBinomPdf(_ploidy, a, _freqs[l] * phiConverted, (1.0 - _freqs[l]) * phiConverted) << "    "
                  << tmp_exp[a] << "    " <<  _gExp[l][i][a] << "    "
                  << phiConverted << "    " << _freqs[l] << std::endl;
        }
      }
    }
  }
}

void ModelDiseq::eStepTwo(){
  std::vector<double> tmp_exp(_ploidy + 1, 0.0);
  double tmp_exp_sum = 0.0, phiConverted = 0.0, max_val = 0.0;
  int i3d = 0; // index for 3D array stored as vector
  bool print = 0;

  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){

      // Catch missing data and skip
      if(_totReads[i][l] == missing)
        continue;

      tmp_exp_sum = 0.0;
      phiConverted = 1.0 / _phi[i] - 1.0;
      /* Prevent phi from blowing up. This results in
         and inbreeding coefficient that is really small. */
      if(phiConverted > 1000)
        phiConverted = 1000.0;

      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        tmp_exp[a] = log(_gLiks[l][i][a]) + r->lnBetaBinomPdf(_ploidy, a, _freqs[l] * phiConverted, (1.0 - _freqs[l]) * phiConverted);
        //std::cerr << tmp_exp[a] << "    " << _gLiks[i3d] << "    " << r->betaBinomPdf(_ploidy, a, _freqs[l] * _phi[i], (1 - _freqs[l]) * _phi[i]) << std::endl;
      }

      max_val = *std::max_element(tmp_exp.begin(), tmp_exp.end());

      for(int a = 0; a <= _ploidy; a++)
        tmp_exp_sum += exp(tmp_exp[a] - max_val);

      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        _gExp[l][i][a] = exp(tmp_exp[a] - max_val - tmp_exp_sum);

        if(print){
        std::cout << i+1 << "    " << l+1 << "    "
                  << a << "    " << r->betaBinomPdf(_ploidy, a, _freqs[l] * phiConverted, (1.0 - _freqs[l]) * phiConverted) << "    "
                  << tmp_exp[a] << "    " <<  _gExp[l][i][a] << "    "
                  << phiConverted << "    " << _freqs[l] << std::endl;
        }

      }

    }
  }
}

void ModelDiseq::mStep(){
  double prev = 0.0;
  for(int l = 0; l < _nLoci; l++){
    if(_convergedLoci[l]){
      continue;
    } else {
      _currLoc = l;
      prev = _freqs[l];
      _perSiteLogLik[l] = local_min(0.0, 1.0, _tol, calcFreqLogLik, _freqs[l]);
    }
    if(sqrt(pow(_freqs[l] - prev, 2)) < _stopVal)
      _convergedLoci[l] = 1;
  }
  for(int i = 0; i < _nInd; i++){
    if(_convergedInd[i]){
      continue;
    } else {
      _currInd = i;
      prev = _phi[i];
      _perIndLogLik[i] = local_min(0.0, 1.0, _tol, calcPhiLogLik, _phi[i]);
    }
    if(sqrt(pow(_phi[i] - prev, 2)) < _stopVal)
      _convergedInd[i] = 1;
  }
}

/*
void ModelDiseq::mStepTwo(){
  double prev = 0.0;
  matrix<double, 1, 1> tmp_param_val;
  for(int l = 0; l < _nLoci; l++){
    if(_convergedLoci[l]){
      continue;
    } else {
      _currLoc = l;
      prev = _freqs[l];
      tmp_param_val(0,0) = _freqs[l];
      find_min_box_constrained(bfgs_search_strategy(),
                               objective_delta_stop_strategy(1e-7),
                               calcFreqLogLik, dlib::derivative(calcFreqLogLik), tmp_param_val, 1.0e-100, 1.0 - 1.0e-100);
      _freqs[l] = tmp_param_val(0,0);
      if(sqrt(pow(_freqs[l] - prev, 2)) < _stopVal)
        _convergedLoci[l] = 1;
    }
  }
  for(int i = 0; i < _nInd; i++){
    if(_convergedInd[i]){
      continue;
    } else {
      _currInd = i;
      prev = _phi[i];
      tmp_param_val(0,0) = _phi[i];
      find_min_box_constrained(bfgs_search_strategy(),
                               objective_delta_stop_strategy(1e-7),
                               calcPhiLogLik, dlib::derivative(calcPhiLogLik), tmp_param_val, 1.0e-100, 1.0 - 1.0e-100);
      _phi[i] = tmp_param_val(0,0);
      if(sqrt(pow(_phi[i] - prev, 2)) < _stopVal)
        _convergedInd[i] = 1;
    }
  }
}
*/

double ModelDiseq::calcLogLik(){
  int i3d = 0;
  double logLik = 0.0, tmpLik = 0.0, alpha = 0.0, beta = 0.0, phiConverted = 0.0;
  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      // Catch missing data and skip
      if(_totReads[i][l] == missing)
        continue;

      phiConverted = 1.0 / _phi[i] - 1.0;
      /* Prevent phi from blowing up. This results in
         and inbreeding coefficient that is really small. */
      if(phiConverted > 1000.0)
        phiConverted = 1000.0;

      alpha = _freqs[l] * phiConverted;
      beta  = (1 - _freqs[l]) * phiConverted;
      tmpLik = 0.0;
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        tmpLik += exp(_gLiks[l][i][a] + r->lnBetaBinomPdf(_ploidy, a, alpha, beta));
      }
      logLik += log(tmpLik);
    }
  }
  return logLik;
}
