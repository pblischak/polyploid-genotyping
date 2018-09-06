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
    } else if(strcmp(av[i],"--alt-reads") == 0 || strcmp(av[i], "-a") == 0){
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
    std::cerr << "\nMissing or invalid option for -a [--alt-reads].\n";
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
  //_gExp.resize(_nLoci * _nInd * (_ploidy + 1));
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
  //_gLiks.resize(_nLoci * _nInd * (_ploidy + 1));
  _freqs.resize(_nLoci);
  _perSiteLogLik.resize(_nLoci);
  _convergedLoci.resize(_nLoci);
  int i2d = 0, i3d = 0; // indices for entries in 2D and 3D arrays stored as vectors
  double gEpsilon = 0.0; // error corrected genotype
  for(int l = 0; l < _nLoci; l++){
    _freqs[l] = r->uniformRv();
    while(_freqs[l] < min_freq || _freqs[l] > max_freq){
      _freqs[l] = r->uniformRv();
    }
    _convergedLoci[l] = 0;
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        if(_totReads[i][l] == MISSING){
          _gLiks[l][i][a] = BADLIK;
        } else {
          gEpsilon = f_epsilon(a, _ploidy, _errRates[l]);
          //_gLiks[l][i][a] = r->binomPdf(_totReads[i][l], _refReads[i][l], gEpsilon);
          _gLiks[l][i][a] = exp(r->lnBinomPdf(_totReads[i][l], _refReads[i][l], gEpsilon));
          // std::cout << gEpsilon << "    " << _gLiks[i3d] << std::endl;
          if(_gLiks[l][i][a] < max_small || isnan(_gLiks[l][i][a])){
            _gLiks[l][i][a] = max_small;
          }
        }
      }
    }
  }
}

void ModelHWE::eStep(){
  std::vector<double> tmp_exp(_ploidy + 1, 0.0);
  double tmp_exp_sum = 0.0, binom_prob = 0.0;
  int i3d = 0; // index for 3D array stored as vector
  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      // Catch missing data and skip
      if(_totReads[i][l] == MISSING)
        continue;

      tmp_exp_sum = 0.0;
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        binom_prob = r->binomPdf(_ploidy, a, _freqs[l]);
        if(isnan(binom_prob)){
          binom_prob = max_small;
        }
        tmp_exp[a] = _gLiks[l][i][a] * binom_prob;
        if(tmp_exp[a] < max_small){
          tmp_exp[a] = max_small;
        }
        tmp_exp_sum += tmp_exp[a];
      }
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        _gExp[l][i][a] = tmp_exp[a] / tmp_exp_sum;
        if(_gExp[l][i][a] < max_small || isnan(_gExp[l][i][a])){
          _gExp[l][i][a] = max_small;
        }
        
        if(_gExp[l][i][a] > almost_one){
          _gExp[l][i][a] = almost_one;
        }
      }
    }
  }
}

void ModelHWE::mStep(){
  int i3d = 0; // index for 3D array stored as vector
  double nChromosomes = 0, numerator = 0.0, prev = 0.0;
  bool check = 1;
  for(int l = 0; l < _nLoci; l++){
    nChromosomes = 0;
    if(_convergedLoci[l]){
      continue;
    } else {
      numerator = 0.0;
      prev = _freqs[l];
      for(int i = 0; i < _nInd; i++){
        // Catch missing data and skip
        if(_totReads[i][l] == MISSING)
          continue;

        nChromosomes += _ploidy;
        for(int a = 1; a <= _ploidy; a++){
          i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
          numerator += (double) (a * _gExp[l][i][a]);
        }
      }
      _freqs[l] = (double) numerator / nChromosomes;
      if(_freqs[l] < min_freq || isnan(_freqs[l])){
        _freqs[l] = min_freq;
        check = 0;
      }
      
      if(_freqs[l] > max_freq){
        _freqs[l] = max_freq;
        check = 0;
      }
      if(sqrt(pow(_freqs[l] - prev, 2)) < _stopVal && check)
        _convergedLoci[l] = 1;
    }
  }
}

double ModelHWE::calcSiteLogLik(double x){
  int i3d = 0;
  double logLik = 0.0;
  for(int i = 0; i < _nInd; i++){
    // Catch missing data and skip
    if(_totReads[i][_currLoc] == MISSING)
      continue;

    for(int a = 0; a <= _ploidy; a++){
      i3d = _currLoc * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
      if(_gLiks[_currLoc][i][a] > max_small){
        logLik += _gExp[_currLoc][i][a] * (log(_gLiks[_currLoc][i][a]) + r->lnBinomPdf(_ploidy, a, x));
      } else {
        logLik += _gExp[_currLoc][i][a] * (log_max_small + r->lnBinomPdf(_ploidy, a, x));
      }
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
  std::cerr << "Writing output files...";
  std::string genosFile = _prefix + "-genos.txt";
  std::string freqsFile = _prefix + "-freqs.txt";
  std::string plFile = _prefix + "-PL.txt";
  std::ofstream genosStream;
  genosStream.open(genosFile, std::ios::out);
  std::ofstream freqsStream;
  freqsStream.open(freqsFile, std::ios::out);
  std::ofstream plStream;
  plStream.open(plFile, std::ios::out);
  std::vector<double> tmp_val(_ploidy+1, 0.0), g_post_prob(_ploidy+1, 0.0);
  double tmp_val_sum = 0.0, tmp_PL = 0.0;
  std::vector< std::vector<int> > genotypes;
  genotypes.resize(_nInd);
  for(int i = 0; i < _nInd; i++){
    genotypes[i].resize(_nLoci, -9);
  }
  int i3d = 0, i2d = 0;
  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;
      tmp_val_sum = 0.0;
      // Catch missing data and skip
      if(_totReads[i][l] == MISSING)
        continue;

      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        tmp_val[a] = _gLiks[l][i][a] * r->binomPdf(_ploidy, a, _freqs[l]);
        tmp_val_sum += tmp_val[a];
      }
      for(int a = 0; a <= _ploidy; a++){
        g_post_prob[a] = tmp_val[a]/tmp_val_sum;
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

  std::cerr << "Done." << std::endl;
}

double ModelHWE::calcLogLik(){
  int i3d = 0;
  double logLik = 0.0, tmpLik = 0.0;
  double freqLik = 0.0;
  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      // Catch missing data and skip
      if(_totReads[i][l] == MISSING)
        continue;
      tmpLik = 0.0;
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
        freqLik = r->binomPdf(_ploidy, a, _freqs[l]);
        if(freqLik < max_small){
          freqLik = max_small;
        }
        if(_totReads[i][l] > 800)
          /*std::cerr << _totReads[i][l] << "," << _refReads[i][l] << "," << _freqs[l] << "    " 
                    << _gExp[l][i][a] << "    " << _gLiks[l][i][a] << "    " << log(_gLiks[l][i][a]) << "    " 
                    << r->binomPdf(_ploidy, a, _freqs[l]) << "    " 
                    << r->lnBinomPdf(_ploidy, a, _freqs[l]) << std::endl;*/
        
        tmpLik += _gLiks[l][i][a] * r->binomPdf(_ploidy, a, _freqs[l]);
        //std::cerr << log(tmpLik) << std::endl;
      }
      logLik += log(tmpLik);
      //std::cerr << logLik << std::endl;
    }
  }
  //_gExp[i3d] *
  return logLik;
}
