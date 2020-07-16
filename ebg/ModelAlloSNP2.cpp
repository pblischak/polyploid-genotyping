#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>

#include "ModelGeneric.hpp"
#include "ModelAlloSNP2.hpp"
#include "MbRandom.hpp"
#include "main.hpp"

std::vector<double> ModelAlloSNP2::_freqs1, ModelAlloSNP2::_freqs2;
int ModelAlloSNP2::_ploidy1 = -999, ModelAlloSNP2::_ploidy2 = -999;
std::string ModelAlloSNP2::_refFreqsFile1 = "none";
std::string ModelAlloSNP2::_refFreqsFile2 = "none";

ModelAlloSNP2::ModelAlloSNP2(int ac, char* av[]){

  _prefix = "alloSNP2";

  for(int i = 2; i < ac; i++){

    if(strcmp(av[i],"--freqs-file1") == 0 || strcmp(av[i], "-f1") == 0){
      _refFreqsFile1 = av[i + 1];
    } else if(strcmp(av[i],"--freqs-file2") == 0 || strcmp(av[i], "-f2") == 0){
      _refFreqsFile2 = av[i + 1];
    }else if(strcmp(av[i],"--num-ind") == 0 || strcmp(av[i], "-n") == 0){
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

void ModelAlloSNP2::printOutput(){
  std::cerr << "Writing output files...";
  std::string gOneFile = _prefix + "-g1.txt";
  std::string gTwoFile = _prefix + "-g2.txt";
  std::string plFile = _prefix + "-PL.txt";

  std::ofstream gOneStream;
  gOneStream.open(gOneFile, std::ios::out);
  std::ofstream gTwoStream;
  gTwoStream.open(gTwoFile, std::ios::out);
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

  std::cerr << "Done." << std::endl;
}

  /****************************/
 /* Private member functions */
/****************************/

void ModelAlloSNP2::checkCommandLine(){

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

  if(strcmp(_refFreqsFile1.c_str(), "none") == 0) {
    std::cerr << "\nMissing or invalid option for -f1 [--freqs-file1].\n";
    errorCaught++;
  }
  
  if(strcmp(_refFreqsFile2.c_str(), "none") == 0) {
    std::cerr << "\nMissing or invalid option for -f2 [--freqs-file2].\n";
    errorCaught++;
  }

  if(errorCaught > 0){
    alloSNP2usage();
    exit(EXIT_FAILURE);
  }

  if(!_quiet)
    std::cerr << "Good" << std::endl;

}

void ModelAlloSNP2::initParams(){

  _ploidy = _ploidy1 + _ploidy2;
  _gLiks.resize(_nLoci);
  for(int l = 0; l < _nLoci; l++){
    _gLiks[l].resize(_nInd);
    for(int i = 0; i < _nInd; i++){
      _gLiks[l][i].resize(_ploidy + 1);
    }
  }

  int i2d = 0, i3d = 0; // indices for entries in 2D and 3D arrays stored as vectors
  double gEpsilon = 0.0; // error corrected genotype

  for(int l = 0; l < _nLoci; l++){
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

void ModelAlloSNP2::getRefFreqs(){

  std::ifstream _refFreqsStream1(_refFreqsFile1);
  std::ifstream _refFreqsStream2(_refFreqsFile2);

  double freqVal = 0.0;

  if(_refFreqsStream1.is_open()){

    while(_refFreqsStream1 >> freqVal)
      _freqs1.push_back(freqVal);

  } else {
    std::cerr << "Could not open reference frequencies file: " << _refFreqsFile1 << std::endl;
    exit(EXIT_FAILURE);
  }

  if(_freqs1.size() != (size_t) _nLoci){
    std::cerr << "Number of loci in reference frequencies file 1 is incorrect: " << std::endl
              << "Num. loci specified: " << _nLoci << std::endl
              << "Num. in reference frequencies: " << _freqs1.size() << std::endl << std::endl;

    exit(EXIT_FAILURE);
  }
  
  if(_refFreqsStream2.is_open()){

    while(_refFreqsStream2 >> freqVal)
      _freqs2.push_back(freqVal);

  } else {
    std::cerr << "Could not open reference frequencies file: " << _refFreqsFile2 << std::endl;
    exit(EXIT_FAILURE);
  }

  if(_freqs2.size() != (size_t) _nLoci){
    std::cerr << "Number of loci in reference frequencies file 2 is incorrect: " << std::endl
              << "Num. loci specified: " << _nLoci << std::endl
              << "Num. in reference frequencies: " << _freqs2.size() << std::endl << std::endl;

    exit(EXIT_FAILURE);
  }

}
