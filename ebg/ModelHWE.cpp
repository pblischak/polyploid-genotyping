#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

#include "ModelGeneric.hpp"
#include "ModelHWE.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

std::vector<bool> ModelHWE::_convergedLoci;

ModelHWE::ModelHWE(int ac, char* av[]){

  _prefix = "hwe";

  // Parse command line options
  for(int i = 2; i < ac; i++){

    if(strcmp(av[i],"--num-ind") == 0 || strcmp(av[i], "-n") == 0){
      _nInd = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--num-loci") == 0 || strcmp(av[i], "-l") == 0){
      _nLoci = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--ploidy") == 0 || strcmp(av[i], "-p") == 0){
      _ploidy = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--total-reads") == 0 || strcmp(av[i], "-t") == 0){
      _totalReadsFile = av[i + 1];
    } else if(strcmp(av[i],"--ref-reads") == 0 || strcmp(av[i], "-r") == 0){
      _refReadsFile = av[i + 1];
    } else if(strcmp(av[i],"--error-rates") == 0 || strcmp(av[i], "-e") == 0) {
      _errorRatesFile = av[i + 1];
    } else if(strcmp(av[i],"--prefix") == 0){
      _prefix = av[i + 1];
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

void ModelHWE::em(){

  double currLogLik = 0.0, prevLogLik = 0.0;
  bool convergence = 0;

  for(int j = 0; j < _maxIters; j++){

    prevLogLik = currLogLik;

    if(!_quiet)
      std::cerr << "Step: " << j + 1 << "        ";

    eStep();
    mStep();

    currLogLik = calcLogLik();
    if(!_quiet)
      std::cerr << std::setprecision(16) << "logLik: " << currLogLik << "        " << "Diff: " << currLogLik - prevLogLik << std::endl;

    checkConvergence(convergence);
    if(convergence || fabs(currLogLik - prevLogLik) < 1e-8)
      break;

  }

}

/*****************************/
/* Private member functions */
/***************************/

/* */

void ModelHWE::checkCommandLine(){

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
    std::cerr << "\nMissing or invalid option for -r [--ref-reads].\n";
    errorCaught++;
  }

  if(strcmp(_errorRatesFile.c_str(), "none") == 0) {
    std::cerr << "\nMissing or invalid option for -e [--error-rates].\n";
    errorCaught++;
  }

  if(errorCaught > 0){
    hweUsage();
    exit(EXIT_FAILURE);
  }

  if(!_quiet)
    std::cerr << "Good" << std::endl;

}

void ModelHWE::checkConvergence(bool &con){

  int convergedLociCount = 0;
  for(int l = 0; l < _nLoci; l++){

    if(_convergedLoci[l])
      convergedLociCount++;

  }

  //std::cerr << "Percent convereged: " << convergedLociCount / (double) _nLoci << std::endl;
  if(convergedLociCount == _nLoci)
    con = 1;

}

void ModelHWE::initParams(){

  _gExp.resize(_nLoci * _nInd * (_ploidy + 1));
  _gLiks.resize(_nLoci * _nInd * (_ploidy + 1));
  _freqs.resize(_nLoci);
  _perSiteLogLik.resize(_nLoci);
  _convergedLoci.resize(_nLoci);

  int i2d = 0, i3d = 0; // indices for entries in 2D and 3D arrays stored as vectors
  double gEpsilon = 0.0; // error corrected genotype

  for(int l = 0; l < _nLoci; l++){
    _freqs[l] = r->uniformRv();
    _convergedLoci[l] = 0;
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;

        if(_totReads[i2d] == MISSING){
          _gLiks[i3d] = BADLIK;
        } else {
          gEpsilon = f_epsilon(a, _ploidy, _errRates[l]);
          _gLiks[i3d] = r->binomPdf(_totReads[i2d], _refReads[i2d], gEpsilon);
          // std::cout << gEpsilon << "    " << _gLiks[i3d] << std::endl;
        }

      }
    }
  }

}

void ModelHWE::eStep(){

  std::vector<double> tmp_exp(_ploidy + 1, 0.0);
  double tmp_exp_sum = 0.0;
  int i3d = 0; // index for 3D array stored as vector

  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){

      // Catch missing data and skip
      if(_totReads[i * _nLoci + l] == MISSING)
        continue;

      tmp_exp_sum = 0.0;

      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        tmp_exp[a] = _gLiks[i3d] * r->binomPdf(_ploidy, a, _freqs[l]);

        tmp_exp_sum += tmp_exp[a];

      }


      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        _gExp[i3d] = tmp_exp[a] / tmp_exp_sum;
      }

    }
  }

}

void ModelHWE::mStep(){

  int i3d = 0; // index for 3D array stored as vector
  double nChromosomes = 0, numerator = 0.0, prev = 0.0;

  for(int l = 0; l < _nLoci; l++){
    nChromosomes = 0;

    if(_convergedLoci[l]){
      continue;
    } else {
      numerator = 0.0;
      prev = _freqs[l];
      for(int i = 0; i < _nInd; i++){

        // Catch missing data and skip
        if(_totReads[i * _nLoci + l] == MISSING)
          continue;

        nChromosomes += _ploidy;
        for(int a = 1; a <= _ploidy; a++){
          i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
          numerator += (double) (a * _gExp[i3d]);
        }
      }

      _freqs[l] = (double) numerator / nChromosomes;

      if(sqrt(pow(_freqs[l] - prev, 2)) < _stopVal)
        _convergedLoci[l] = 1;

    }

  }
}

double ModelHWE::calcSiteLogLik(double x){

  int i3d = 0;
  double logLik = 0.0;

  for(int i = 0; i < _nInd; i++){

    // Catch missing data and skip
    if(_totReads[i * _nLoci + _currLoc] == MISSING)
      continue;

    for(int a = 0; a <= _ploidy; a++){
      i3d = _currLoc * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
      logLik += _gExp[i3d] * (log(_gLiks[i3d]) + r->lnBinomPdf(_ploidy, a, x));
    }
  }

  return -logLik;

}

void ModelHWE::mStepBrent(){

  for(int l = 0; l < _nLoci; l++){
    _currLoc = l;
    _perSiteLogLik[l] = local_min(0.0, 1.0, _tol, calcSiteLogLik, _freqs[l]);
  }

}

void ModelHWE::printOutput(){

  std::string genosFile = _prefix + "-genos.txt";
  std::string freqsFile = _prefix + "-freqs.txt";

  std::ofstream genosStream;
  genosStream.open(genosFile, std::ios::out);


  std::ofstream freqsStream;
  freqsStream.open(freqsFile, std::ios::out);

  std::vector<double> tmp_val(_ploidy+1, 0.0), g_post_prob(_ploidy+1, 0.0);
  double tmp_val_sum = 0.0;
  std::vector<int> genotypes(_nInd * _nLoci, -9);
  int i3d = 0, i2d = 0;

  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;
      tmp_val_sum = 0.0;

      // Catch missing data and skip
      if(_totReads[i2d] == MISSING)
        continue;

      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        tmp_val[a] = _gLiks[i3d] * r->binomPdf(_ploidy, a, _freqs[l]);
        tmp_val_sum += tmp_val[a];
      }

      for(int a = 0; a <= _ploidy; a++){
        g_post_prob[a] = tmp_val[a]/tmp_val_sum;
      }

      genotypes[i2d] = gMax(g_post_prob);

    }
  }

  for(int i = 0; i < _nInd; i++){
    for(int l = 0; l < _nLoci; l++){
      i2d = i * _nLoci + l;
      genosStream << genotypes[i2d] << "\t";
    }
    genosStream << std::endl;
  }

  for(int l = 0; l < _nLoci; l++)
    freqsStream << _freqs[l] << std::endl;

}

double ModelHWE::calcLogLik(){

  int i3d = 0;
  double logLik = 0.0, tmpLik = 0.0;

  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){

      // Catch missing data and skip
      if(_totReads[i * _nLoci + l] == MISSING)
        continue;

      tmpLik = 0.0;

      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        //std::cerr << _gExp[i3d] << "    " << log(_gLiks[i3d]) << "    " << r->binomPdf(_ploidy, a, _freqs[l]) << std::endl;

        tmpLik += _gLiks[i3d] * r->binomPdf(_ploidy, a, _freqs[l]);
        //std::cerr << log(tmpLik) << std::endl;

      }

      logLik += log(tmpLik);
      //std::cerr << logLik << std::endl;

    }
  }

  //_gExp[i3d] *
  return logLik;

}
