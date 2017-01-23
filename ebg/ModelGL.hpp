#ifndef MODEL_GL_HPP
#define MODEL_GL_HPP


class ModelGL : public ModelGeneric {

public:
  ModelGL(int ac, char* av[]); /* Constructor for setting up model. */
  void printOutput(); /* Calculate posterior probs and print. */

private:
  void checkCommandLine();

};

#endif // MODEL_GL_HPP
