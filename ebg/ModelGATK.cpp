#include <vector>
#include <algorithm>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

#include "ModelGeneric.hpp"
#include "ModelGATK.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

ModelGATK::ModelGATK(int ac, char* av[]){

  _prefix = "gatk";

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
    } else if(strcmp(av[i],"--print-probs") == 0){
      _printProbs = 1;
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

  // Initialize parameters
  initParams();

}

void ModelGATK::checkCommandLine(){

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
    gatkUsage();
    exit(EXIT_FAILURE);
  }

  if(!_quiet)
    std::cerr << "Good" << std::endl;

}

void ModelGATK::initParams(){

  _gLiks.resize(_nLoci * _nInd * (_ploidy + 1));
  int i2d = 0, i3d = 0; // indices for entries in 2D and 3D arrays stored as vectors
  double gEpsilon = 0.0; // error corrected genotype

  for(int l = 0; l < _nLoci; l++){
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

void ModelGATK::printOutput(){

  std::string genosFile = _prefix + "-genos.txt";
  std::string probsFile = _prefix + "-probs.txt";
  std::string ratioFile = _prefix + "-ratio.txt";

  std::ofstream genosStream;
  genosStream.open(genosFile, std::ios::out);

  std::ofstream probsStream;
  if(_printProbs)
    probsStream.open(probsFile, std::ios::out);

  std::ofstream ratioStream;
  ratioStream.open(ratioFile, std::ios::out);

  std::vector<double> tmp_val(_ploidy+1, 0.0), g_post_prob(_ploidy+1, 0.0), ratio(_nInd * _nLoci, -9);
  double tmp_val_sum = 0.0;
  std::vector<int> genotypes(_nInd * _nLoci, -9);
  int i3d = 0, i2d = 0;

  for(int l = 0; l < _nLoci; l++){
    for(int i = 0; i < _nInd; i++){
      i2d = i * _nLoci + l;
      tmp_val_sum = 0.0;
      if(_printProbs)
        probsStream << l+1 << "\t" << i+1 << "\t";

      // Catch missing data and skip
      if(_totReads[i2d] == MISSING)
        continue;

        for(int a = 0; a <= _ploidy; a++){
          i3d = l * _nInd * (_ploidy + 1) + i * (_ploidy + 1) + a;
          tmp_val[a] = _gLiks[i3d];
          tmp_val_sum += tmp_val[a];
        }

        for(int a = 0; a <= _ploidy; a++){
          g_post_prob[a] = tmp_val[a]/tmp_val_sum;
          if(_printProbs)
            probsStream << g_post_prob[a] << "\t";
        }
        if(_printProbs)
          probsStream << std::endl;

        genotypes[i2d] = gMax(g_post_prob);
        std::sort(g_post_prob.begin(), g_post_prob.end(), std::greater<double>());
        ratio[i2d] = g_post_prob[0] / g_post_prob[1];

      }
    }

    for(int i = 0; i < _nInd; i++){
      for(int l = 0; l < _nLoci; l++){
        i2d = i * _nLoci + l;
        genosStream << genotypes[i2d] << "\t";
        ratioStream << ratio[i2d] << "\t";
      }
      genosStream << std::endl;
      ratioStream << std::endl;
    }

}
