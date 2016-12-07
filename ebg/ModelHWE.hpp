#ifndef MODEL_HWE_HPP
#define MODEL_HWE_HPP

class ModelHWE : public ModelGeneric {

public:
  ModelHWE(int ac, char* av[]); /* Constructor for setting up model from command line arguments. */
  void em(); /* Wrapper function for EM algorithm that is called in main. */
  void printOutput();

private:
  // Member variables
  static std::vector<bool> _convergedLoci;

  // Member functions for EM algorithm
  void eStep();
  void mStep();
  void mStepBrent();
  static double calcLogLik();

  static double calcSiteLogLik(double x);

  // Member functions for data parsing
  void checkCommandLine();
  void initParams();
  void checkConvergence(bool &con);

};

#endif //MODEL_HWE_HPP
