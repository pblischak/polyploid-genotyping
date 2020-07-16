#ifndef MODEL_ALLO_SNP2_HPP
#define MODEL_ALLO_SNP2_HPP

/* Inherits common parameters from ModelGeneric class. */

class ModelAlloSNP2 : public ModelGeneric {

public:
  ModelAlloSNP2(int ac, char* av[]);
  void printOutput();

private:
  //Member variables
  static std::vector<double> _freqs1, _freqs2;
  static int _ploidy1, _ploidy2;
  static std::string _refFreqsFile1, _refFreqsFile2;

  // Member functions for data parsing
  void getRefFreqs();
  void checkCommandLine();
  void initParams();

  //Member functions
  void gMax2(std::vector<double> &v, int row, int col, std::vector<int> &result);

};

inline void ModelAlloSNP2::gMax2(std::vector<double> &v, int row, int col, std::vector<int> &result){

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

#endif //MODEL_ALLO_SNP2_HPP
