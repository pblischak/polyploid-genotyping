#ifndef MODEL_DISEQ_HPP
#define MODEL_DISEQ_HPP

/* Inherits common parameters from ModelGeneric class. */

class ModelDiseq : public ModelGeneric {

public:
  ModelDiseq(int ac, char* av[]); /* Constructor for setting up model. */
  void ecm(); /* Wrapper function for ECM that is called in main. */
  void printOutput();

private:
  // Member variables
  static std::vector<double> _phi, _prevPhi, _perIndLogLik, _theta0Ind,
                             _theta1Ind, _theta2Ind,
                             _rInd, _vInd, _thetaPrimeInd;
  static int _currInd;
  static std::vector<bool> _convergedLoci, _convergedInd;

  // Member functions for data parsing
  void checkCommandLine();
  void initParams();

  // Member functions
  static double calcFreqLogLik(double x);
  static double calcPhiLogLik(double x);
  static double calcLogLik();
  void eStep();
  void eStepTwo();
  void mStep();
  //void mStepAccel();
  void checkConvergence(bool &con);

};

#endif //MODEL_DISEQ_HPP
