#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

#include "ModelGeneric.hpp"
#include "ModelAlloSNP.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

std::vector<double> ModelAlloSNP::_freqs1, ModelAlloSNP::_freqs2, ModelAlloSNP::_perSiteLogLikFreqs2;
int ModelAlloSNP::_ploidy1 = -999, ModelAlloSNP::_ploidy2 = -999;
std::string ModelAlloSNP::_refFreqsFile = "none";
std::vector<bool> ModelAlloSNP::_convergedLoci;

ModelAlloSNP::ModelAlloSNP(int ac, char* av[]){

  _prefix = "alloSNP";

  for(int i = 2; i < ac; i++){

    if(strcmp(av[i],"--freqs-file") == 0 || strcmp(av[i], "-f") == 0){
      _refFreqsFile = av[i + 1];
    } else if(strcmp(av[i],"--num-ind") == 0 || strcmp(av[i], "-n") == 0){
      _nInd = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--num-loci") == 0 || strcmp(av[i], "-l") == 0){
      _nLoci = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--ploidy1") == 0 || strcmp(av[i], "-p1") == 0){
      _ploidy1 = atoi(av[i + 1]);
    } else if(strcmp(av[i],"--ploidy2") == 0 || strcmp(av[i], "-p2") == 0){
      _ploidy2 = atoi(av[i + 1]);
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
      _stopVal = atof(av[i + 1]);
    } else if(strcmp(av[i],"--brent") == 0){
      _brent = 1;
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

  // Read in reference frequencies
  getRefFreqs();

}

void ModelAlloSNP::em(){

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
      std::cerr << std::setprecision(16) << "logLik: " << currLogLik << "\tDiff: " << currLogLik - prevLogLik << std::endl;

    checkConvergence(convergence);

    if(convergence)
      break;

  }

  if(!convergence && _brent){
    prevLogLik = currLogLik;

    if(!_quiet){
      std::cerr << "EM algorithm did not converge in " << _maxIters << " iterations. Switching to Brent's optimization." << std::endl;
      std::cerr << "Brent step:\t";
    }

    mStepBrent();
    currLogLik = calcLogLik();

    if(!_quiet)
      std::cerr << std::setprecision(10) << "-logLik: " << currLogLik << "\tDiff: " << currLogLik - prevLogLik << std::endl;

  }

}

void ModelAlloSNP::printOutput(){
  std::cerr << "Writing output files...";
  std::string gOneFile = _prefix + "-g1.txt";
  std::string gTwoFile = _prefix + "-g2.txt";
  std::string freqsTwoFile = _prefix + "-freqs2.txt";
  std::string plFile = _prefix + "-PL.txt";

  std::ofstream gOneStream;
  gOneStream.open(gOneFile, std::ios::out);
  std::ofstream gTwoStream;
  gTwoStream.open(gTwoFile, std::ios::out);
  std::ofstream freqsTwoStream;
  freqsTwoStream.open(freqsTwoFile, std::ios::out);
  std::ofstream plStream;
  plStream.open(plFile, std::ios::out);

  std::vector<double> tmp_val((_ploidy1 + 1)*(_ploidy2 + 1), 0.0), g_post_prob((_ploidy1 + 1)*(_ploidy2 + 1), 0.0);
  double tmp_val_sum = 0.0, tmp_PL = 0.0;
  std::vector<int> res(2);
  std::vector< std::vector<int> > genotypes1, genotypes2;
  genotypes1.resize(_nInd), genotypes2.resize(_nInd);
  for(int i = 0; i < _nInd; i++){
    genotypes1[i].resize(_nLoci, -9);
    genotypes2[i].resize(_nLoci, -9);
  }
  int i3d = 0, i2d = 0, i2d2 = 0;

  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;

      // Catch missing data and skip
      if(_totReads[i][l] == MISSING)
        continue;

      tmp_val_sum = 0.0;

      for(int a1 = 0; a1 <= _ploidy1; a1++){
        for(int a2 = 0; a2 <= _ploidy2; a2++){
          i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + (a1 + a2);
          i2d2 = a1 * (_ploidy2 + 1) + a2;

          tmp_val[i2d2] = _gLiks[l][i][a1+a2] * (r->binomPdf(_ploidy1, a1, _freqs1[l])
                                              *  r->binomPdf(_ploidy2, a2, _freqs2[l]));

          tmp_val_sum += tmp_val[i2d2];

        }
      }

      for(int a1 = 0; a1 <= _ploidy1; a1++){
        for(int a2 = 0; a2 <= _ploidy2; a2++){
          i2d2 = a1 * (_ploidy2 + 1) + a2;
          g_post_prob[i2d2] = tmp_val[i2d2] / tmp_val_sum;
          tmp_PL = -10 * log10(g_post_prob[i2d2]);
          if(tmp_PL == -0){
            plStream << -1 * tmp_PL;
          } else {
            plStream << tmp_PL;
          }
          if(((a1+1) * (a2+1)) < ((_ploidy1+1) * (_ploidy2+1))){
            plStream << ",";
          }
        }
      }
      plStream << "\t";
      gMax2(g_post_prob, _ploidy1 + 1, _ploidy2 + 1, res);
      genotypes1[i][l] = res[0];
      genotypes2[i][l] = res[1];
    }
    plStream << std::endl;
  }

  for(int i = 0; i < _nInd; i++){
    for(int l = 0; l < _nLoci; l++){
      i2d = i * _nLoci + l;
      gOneStream << genotypes1[i][l] << "\t";
      gTwoStream << genotypes2[i][l] << "\t";
    }
    gOneStream << std::endl;
    gTwoStream << std::endl;
  }

  for(int l = 0; l < _nLoci; l++)
    freqsTwoStream << _freqs2[l] << std::endl;

  std::cerr << "Done." << std::endl;
}

  /****************************/
 /* Private member functions */
/****************************/

void ModelAlloSNP::checkConvergence(bool &con){

  int convergedLociCount = 0;

  for(int l = 0; l < _nLoci; l++){

    if(_convergedLoci[l])
      convergedLociCount++;

  }

  //std::cerr << "Percent convereged: " << convergedLociCount / (double) _nLoci << std::endl;

  if(convergedLociCount == _nLoci)
    con = 1;

}

void ModelAlloSNP::checkCommandLine(){

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

  if(_ploidy1 < 0){
    std::cerr << "\nMissing or invalid option for -p1 [--ploidy1].\n";
    errorCaught++;
  }

  if(_ploidy2 < 0){
    std::cerr << "\nMissing or invalid option for -p2 [--ploidy2].\n";
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

  if(strcmp(_refFreqsFile.c_str(), "none") == 0) {
    std::cerr << "\nMissing or invalid option for -f [--freqs-file].\n";
    errorCaught++;
  }

  if(errorCaught > 0){
    alloSNPusage();
    exit(EXIT_FAILURE);
  }

  if(!_quiet)
    std::cerr << "Good" << std::endl;

}

void ModelAlloSNP::initParams(){

  _ploidy = _ploidy1 + _ploidy2;
  _gExp4D.resize(_nLoci);
  for(int l = 0; l < _nLoci; l++){
    _gExp4D[l].resize(_nInd);
    for(int i = 0; i < _nInd; i++){
      _gExp4D[l][i].resize(_ploidy1 + 1);
      for(int a1 = 0; a1 <= _ploidy1; a1++){
        _gExp4D[l][i][a1].resize(_ploidy2 + 1);
      }
    }
  }
  //_gExp.resize(_nLoci * _nInd * (_ploidy1 + 1) * (_ploidy2 + 1));
  _gLiks.resize(_nLoci);
  for(int l = 0; l < _nLoci; l++){
    _gLiks[l].resize(_nInd);
    for(int i = 0; i < _nInd; i++){
      _gLiks[l][i].resize(_ploidy + 1);
    }
  }
  //* _nInd * (_ploidy + 1));
  _freqs2.resize(_nLoci);

  _perSiteLogLikFreqs2.resize(_nLoci);
  _convergedLoci.resize(_nLoci);

  int i2d = 0, i3d = 0; // indices for entries in 2D and 3D arrays stored as vectors
  double gEpsilon = 0.0; // error corrected genotype

  for(int l = 0; l < _nLoci; l++){
    _freqs2[l] = r->uniformRv();
    _convergedLoci[l] = 0;
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;
      for(int a = 0; a <= _ploidy; a++){
        i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;

        if(_totReads[i][l] == MISSING){
          _gLiks[l][i][a] = BADLIK;
        } else {
          gEpsilon = f_epsilon(a, _ploidy, _errRates[l]);
          _gLiks[l][i][a] = r->binomPdf(_totReads[i][l], _refReads[i][l], gEpsilon);
          //std::cerr << gEpsilon << "    " << _gLiks[i3d] << std::endl;
        }

      }
    }
  }

}

void ModelAlloSNP::getRefFreqs(){

  std::ifstream _refFreqsStream(_refFreqsFile);

  double freqVal = 0.0;

  if(_refFreqsStream.is_open()){

    while(_refFreqsStream >> freqVal)
      _freqs1.push_back(freqVal);

  } else {
    std::cerr << "Could not open reference frequencies file: " << _refFreqsFile << std::endl;
    exit(EXIT_FAILURE);
  }

  if(_freqs1.size() != (size_t) _nLoci){
    std::cerr << "Number of loci in reference frequencies file is incorrect: " << std::endl
              << "Num. loci specified: " << _nLoci << std::endl
              << "Num. in reference frequencies: " << _freqs1.size() << std::endl << std::endl;

    exit(EXIT_FAILURE);
  }


}

void ModelAlloSNP::eStep(){

  std::vector<double> tmp_exp((_ploidy1 + 1) * (_ploidy2 + 1), 0.0);
  double tmp_exp_sum = 0.0;
  int i4d = 0, i3d = 0, i2d = 0; // index for 4D, 3D, and 2D array stored as vector

  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      tmp_exp_sum = 0.0;

      if(_totReads[i][l] == MISSING)
        continue;

      for(int a1 = 0; a1 <= _ploidy1; a1++){
        for(int a2 = 0; a2 <= _ploidy2; a2++){
          i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + (a1 + a2);
          i2d = a1 * (_ploidy2 + 1) + a2;

          tmp_exp[i2d] = _gLiks[l][i][a1+a2] * r->binomPdf(_ploidy1, a1, _freqs1[l])
                                             * r->binomPdf(_ploidy2, a2, _freqs2[l]);

          tmp_exp_sum += tmp_exp[i2d];

          if(std::isnan(tmp_exp[i2d]))
            std::cerr << l+1 << "\t" << _freqs2[l] <<  "\t" << r->binomPdf(_ploidy1, a1, _freqs1[l]) << "\t" << r->binomPdf(_ploidy2, a2, _freqs2[l]) << "\t" << tmp_exp[i2d] << std::endl;

        }
      }

      for(int a1 = 0; a1 <= _ploidy1; a1++){
        for(int a2 = 0; a2 <= _ploidy2; a2++){
          i4d = l * _nInd * (_ploidy1 + 1) * (_ploidy2 + 1) + i * (_ploidy1 + 1) * (_ploidy2 + 1) + a1 * (_ploidy2 + 1) + a2;
          i2d = a1 * (_ploidy2 + 1) + a2;

          _gExp4D[l][i][a1][a2] = tmp_exp[i2d] / tmp_exp_sum;
        }
      }

    }
  }

}

void ModelAlloSNP::mStep(){

  int i4d = 0;
  double nChromosomes = 0, numerator = 0.0, prev = 0.0;

  for(int l = 0; l < _nLoci; l++){

    nChromosomes = 0;

    if(_convergedLoci[l]){
      continue;
    } else {
      numerator = 0.0;
      prev = _freqs2[l];
      for(int i = 0; i < _nInd; i++){

        // Catch missing data and skip
        if(_totReads[i][l] == MISSING){
          //std::cerr << "Skip " << i+1 << "," << l+1 << std::endl;
          continue;
        }

        nChromosomes += _ploidy2;
        for(int a1 = 0; a1 <= _ploidy1; a1++){
          for(int a2 = 0; a2 <= _ploidy2; a2++){
            i4d = l * _nInd * (_ploidy1 + 1) * (_ploidy2 + 1) + i * (_ploidy1 + 1) * (_ploidy2 + 1) + a1 * (_ploidy2 + 1) + a2;
            numerator += (double) a2 * _gExp4D[l][i][a1][a2];
          }
        }

      }

      //std::cerr << nChromosomes << "\t" << _nInd * _ploidy2 << std::endl;

      _freqs2[l] = numerator / nChromosomes;

      if(sqrt(pow(_freqs2[l] - prev, 2)) < _stopVal)
        _convergedLoci[l] = 1;

      //if(std::isnan(numerator))
      //  std::cerr << l+1 << "\t" << numerator << "\t" << _freqs2[l] << std::endl;
    }
  }

}

void ModelAlloSNP::mStepBrent(){

  for(int l = 0; l < _nLoci; l++){
    _currLoc = l;

    if(_convergedLoci[l]){
      continue;
    } else {
      _perSiteLogLikFreqs2[l] = local_min(0.0, 1.0, _tol, calcFreqs2LogLik, _freqs2[l]);
    }

  }

}

double ModelAlloSNP::calcFreqs2LogLik(double x){

  int i3d = 0;
  double tmpLik = 0.0, logLik = 0.0;

  for(int i = 0; i < _nInd; i++){
    tmpLik = 0.0;

    // Catch missing data and skip
    if(_totReads[i][_currLoc] == MISSING)
      continue;


    for(int a1 = 0; a1 <= _ploidy1; a1++){
      for(int a2 = 0; a2 <= _ploidy2; a2++){
        i3d = _currLoc * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + (a1 + a2);
        tmpLik += _gLiks[_currLoc][i][a1+a2] * r->binomPdf(_ploidy1, a1, _freqs1[_currLoc])
                                             * r->binomPdf(_ploidy2, a2, x);
      }
    }

    logLik += log(tmpLik);

  }

  return -logLik;

}

double ModelAlloSNP::calcLogLik(){
  int i3d = 0;
  double logLik = 0.0, tmpLik = 0.0;

  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){

      // Catch missing data and skip
      if(_totReads[i][l] == MISSING)
        continue;

      for(int a1 = 0; a1 <= _ploidy1; a1++){
        for(int a2 = 0; a2 <= _ploidy2; a2++){
          i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + (a1 + a2);
          tmpLik += _gLiks[l][i][a1+a2] * (r->binomPdf(_ploidy1, a1, _freqs1[l]) * r->binomPdf(_ploidy2, a2, _freqs2[l]));
        }
      }

      logLik += log(tmpLik);
    }
  }

  return logLik;

}
