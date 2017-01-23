#ifndef MODEL_GATK_HPP
#define MODEL_GATK_HPP

class ModelGATK : public ModelGeneric {

public:
  ModelGATK(int ac, char* av[]); /* Constructor for setting up model. */
  void printOutput(); /* Calculate posterior probs and print. */

private:
  void checkCommandLine();
  void initParams();
  bool _printProbs = 0;

};

#endif // MODEL_GATK_HPP
