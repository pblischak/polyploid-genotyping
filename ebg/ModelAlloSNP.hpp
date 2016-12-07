#ifndef MODEL_ALLO_SNP_HPP
#define MODEL_ALLO_SNP_HPP

/* Inherits common parameters from ModelGeneric class. */

class ModelAlloSNP : public ModelGeneric {

public:
  ModelAlloSNP(int ac, char* av[]);
  void em();
  void printOutput();

private:
  //Member variables
  static std::vector<double> _freqs1, _freqs2, _perSiteLogLikFreqs2;
  static int _ploidy1, _ploidy2;
  static std::string _refFreqsFile;
  static std::vector<bool> _convergedLoci;

  // Member functions for data parsing
  void getRefFreqs();
  void checkCommandLine();
  void initParams();
  double calcLogLik();

  //Member functions
  static double calcFreqs2LogLik(double x);
  void eStep();
  void mStep();
  void mStepBrent();
  void gMax2(std::vector<double> &v, int row, int col, std::vector<int> &result);
  void checkConvergence(bool &con);

};

inline void ModelAlloSNP::gMax2(std::vector<double> &v, int row, int col, std::vector<int> &result){

  int currMaxRow = 0, currMaxCol = 0;
  double currMax = 0.0;

  if(result.size() != 2){
    std::cout << "You done messed up...\n" << std::endl;
    exit(EXIT_FAILURE);
  }

  for(int r = 0; r < row; r++){
    for(int c = 0; c < col; c++){

      if(v[r * col + c] > currMax){
        currMax = v[r * col + c];
        currMaxRow = r;
        currMaxCol = c;
      }

    }
  }

  result[0] = currMaxRow;
  result[1] = currMaxCol;

}

#endif //MODEL_ALLO_SNP_HPP
