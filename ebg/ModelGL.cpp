#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

#include "ModelGeneric.hpp"
#include "ModelGL.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

ModelGL::ModelGL(int ac, char* av[]){

  _prefix = "gl";

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
    } else if(strcmp(av[i],"--prefix") == 0){
      _prefix = av[i + 1];
    } else if(strcmp(av[i],"--error-rates") == 0 || strcmp(av[i], "-e") == 0) {
      _errorRatesFile = av[i + 1];
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

}

void ModelGL::checkCommandLine(){

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
    gatkUsage();
    exit(EXIT_FAILURE);
  }

  if(!_quiet)
    std::cerr << "Good" << std::endl;

}

void ModelGL::printOutput(){

  std::string genosFile   = _prefix + "-genos.txt";
  std::string logDiffFile = _prefix + "-logDiff.txt";

  std::ofstream genosStream;
  genosStream.open(genosFile, std::ios::out);

  std::ofstream logDiffStream;
  logDiffStream.open(logDiffFile, std::ios::out);

  double gEpsilon = 0.0, curr_max_prob = 0;
  int i2d = 0, curr_max = 0;
  std::vector<double> tmp_g(_ploidy + 1, 0);

  for(int i = 0; i < _nInd; i++){
    for(int l = 0; l < _nLoci; l++){
      i2d = i * _nLoci + l;

      if(_totReads[i2d] == MISSING){
        genosStream   << "-9\t";
        logDiffStream << "-9\t";
        continue;
      }
      for(int a = 0; a <= _ploidy; a++){
        gEpsilon = f_epsilon(a, _ploidy, _errRates[l]);
        if(a == 0){
          tmp_g[a] = r->lnBinomPdf(_totReads[i2d], _refReads[i2d], gEpsilon);
          curr_max = 0;
          curr_max_prob = tmp_g[a];
        } else {
          tmp_g[a] = r->lnBinomPdf(_totReads[i2d], _refReads[i2d], gEpsilon);
          if(tmp_g[a] > curr_max_prob){
            curr_max = a;
            curr_max_prob = tmp_g[a];
          }
        }

      }

      genosStream << curr_max << "\t";
      std::sort(tmp_g.begin(), tmp_g.end(), std::greater<double>());
      logDiffStream << pow(10, (tmp_g[0] / log(10)) - (tmp_g[1] / log(10))) << "\t";

    }
    genosStream   << std::endl;
    logDiffStream << std::endl;
  }

}
