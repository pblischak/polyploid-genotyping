#include <iostream>
#include <string>
#include <cstring>
#ifdef _OPENMP
  #include <omp.h>
#endif

#include "MbRandom.hpp"
#include "ModelGeneric.hpp"
#include "ModelHWE.hpp"
#include "ModelDiseq.hpp"
#include "ModelAlloSNP.hpp"
#include "ModelAlloSNP2.hpp"
#include "ModelGATK.hpp"
#include "ModelGL.hpp"
#include "main.hpp"

MbRandom *r = new MbRandom;

int main(int argc, char* argv[]){
  std::string selectedModel = "none";
  if(argc < 2){
    std::cerr << "\nNot enough arguments.\n" << std::endl;
    mainUsage();
    exit(EXIT_FAILURE);
  }
  if(strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0){
    mainUsage();
    exit(EXIT_SUCCESS);
  }
  if(strcmp(argv[1],"-V") == 0 || strcmp(argv[1],"--version") == 0){
    std::cerr << "\nThis is ebg version " << VERSION << " (" << VERSIONDATE << ").\n" << std::endl;
    exit(EXIT_SUCCESS);
  }
  if(strcmp(argv[1], "hwe") == 0 && (strcmp(argv[2], "--help") == 0 || strcmp(argv[2], "-h") == 0)){
    hweUsage();
    exit(EXIT_SUCCESS);
  }
  if(strcmp(argv[1], "diseq") == 0 && (strcmp(argv[2], "--help") == 0 || strcmp(argv[2], "-h") == 0)){
    diseqUsage();
    exit(EXIT_SUCCESS);
  }
  if(strcmp(argv[1], "alloSNP") == 0 && (strcmp(argv[2], "--help") == 0 || strcmp(argv[2], "-h") == 0)){
    alloSNPusage();
    exit(EXIT_SUCCESS);
  }
  if(strcmp(argv[1], "alloSNP2") == 0 && (strcmp(argv[2], "--help") == 0 || strcmp(argv[2], "-h") == 0)){
    alloSNP2usage();
    exit(EXIT_SUCCESS);
  }
  if(strcmp(argv[1], "gatk") == 0 && (strcmp(argv[2], "--help") == 0 || strcmp(argv[2], "-h") == 0)){
    gatkUsage();
    exit(EXIT_SUCCESS);
  }
  if(strcmp(argv[1], "gl") == 0 && (strcmp(argv[2], "--help") == 0 || strcmp(argv[2], "-h") == 0)){
    glUsage();
    exit(EXIT_SUCCESS);
  }
  if(strcmp(argv[1], "hwe") != 0 && strcmp(argv[1], "diseq") != 0 && strcmp(argv[1], "alloSNP") != 0 && strcmp(argv[1], "alloSNP2") != 0 && strcmp(argv[1], "gatk") != 0 && strcmp(argv[1], "gl") != 0){
    std::cerr << "\nA valid model was not selected. The options are: hwe, diseq, alloSNP, alloSNP2, gatk, and gl.\n" << std::endl;
    mainUsage();
    exit(EXIT_FAILURE);
  } else {
    selectedModel = argv[1];
  }
  //if(strcmp(selectedModel.c_str(), "none") == 0){
  //  std::cerr << "\nA model was not selected. The options are: freqs, diseq, and alloSNP.\n" << std::endl;
  //} else
  if(argc <= 2){
    std::cerr << "Too few arguments.\n" << std::endl;
    mainUsage();
    exit(EXIT_FAILURE);
  }

  if(strcmp(selectedModel.c_str(), "hwe") == 0){
    ModelHWE mHWE(argc, argv);
    mHWE.em();
    mHWE.printOutput();
  } else if(strcmp(selectedModel.c_str(), "diseq") == 0){
    ModelDiseq mDiseq(argc, argv);
    mDiseq.ecm();
    mDiseq.printOutput();
  } else if(strcmp(selectedModel.c_str(), "alloSNP") == 0){
    ModelAlloSNP mAlloSNP(argc, argv);
    mAlloSNP.em();
    mAlloSNP.printOutput();
  } else if(strcmp(selectedModel.c_str(), "alloSNP2") == 0){
    ModelAlloSNP2 mAlloSNP2(argc,argv);
    mAlloSNP2.printOutput();
  }else if(strcmp(selectedModel.c_str(), "gatk") == 0){
    ModelGATK mGATK(argc, argv);
    mGATK.printOutput();
  } else if(strcmp(selectedModel.c_str(), "gl") == 0){
    ModelGL mGL(argc, argv);
    mGL.printOutput();
  } else {
    std::cerr << "\nThe selected model (" << selectedModel << ") is not a valid option.\n"
              << "Please choose one of the following: hwe, diseq, alloSNP, gatk, gl.\n" << std::endl;
  }
  return 0;
}

void mainUsage(){
  std::cerr << "\nUsage: ebg <model> [model options]\n" << std::endl;
  std::cerr << "Main options:" << std::endl;
  std::cerr << "  -h [--help]        Prints a help message" << std::endl;
  std::cerr << "  -V [--version]     Print version information" << std::endl;
  std::cerr << "\nAvailable models:" << std::endl;
  std::cerr << "  hwe                Estimate genotypes under Hardy-Weinberg" << std::endl;
  std::cerr << "  diseq              Estimate genotypes assuming H-W disequilibrium" << std::endl;
  std::cerr << "  alloSNP            Estimate genotypes within the subgenomes of an allopolyploid using one parent's allele frequencies" << std::endl;
  std::cerr << "  alloSNP2           Estimate genotypes within the subgenomes of an allopolyploid using both parent's allele frequencies" << std::endl;
  std::cerr << "  gatk               Estimate genotypes using the GATK model (flat genotype prior)." << std::endl;
  std::cerr << "  gl                 Estimate genotypes directly using genotype likelihoods." << std::endl;
  std::cerr << "\nOptions for each model can be found by typing: ./ebg <model> --help\n" << std::endl;
}

void hweUsage(){
  std::cerr << "\nUsage: ebg hwe [hwe options]\n" << std::endl;
  std::cerr << "Required options:" << std::endl
            << "  -n [--num-ind]      <int>        Number of individuals" << std::endl
            << "  -l [--num-loci]     <int>        Number of loci" << std::endl
            << "  -p [--ploidy]       <int>        Ploidy level" << std::endl
            << "  -t [--total-reads]  <string>     Total read counts file" << std::endl
            << "  -a [--alt-reads]    <string>     ALT allele read counts file" << std::endl
            << "  -e [--error-rates]  <string>     Sequencing error rates file" << std::endl
            << "\nAdditional options:" << std::endl
            << "  --tol               <double>     Tolerance for Brent's method (default = 1e-10)" << std::endl
            << "  --iters             <int>        Number of iterations to run EM algorithm (default = 100)" << std::endl
            << "  --stop              <double>     Stop value for EM algorithm parameter updates (default = 1e-5)" << std::endl
            << "  --prefix            <string>     Prefix for output files (default = hwe)" << std::endl
            << "  --quiet                          Flag to suppress output to stdout\n" << std::endl;
}

void diseqUsage(){
  std::cerr << "\nUsage: ebg diseq [diseq options]\n" << std::endl;
  std::cerr << "Required options:" << std::endl
            << "  -n [--num-ind]      <int>        Number of individuals" << std::endl
            << "  -l [--num-loci]     <int>        Number of loci" << std::endl
            << "  -p [--ploidy]       <int>        Ploidy level" << std::endl
            << "  -t [--total-reads]  <string>     Total read counts file" << std::endl
            << "  -a [--alt-reads]    <string>     ALT allele read counts file" << std::endl
            << "  -e [--error-rates]  <string>     Sequencing error rates file" << std::endl
            << "\nAdditional options:" << std::endl
            << "  --tol               <double>     Tolerance for Brent's method (default = 1e-10)" << std::endl
            << "  --iters             <int>        Number of iterations to run ECM algorithm (default = 100)" << std::endl
            << "  --stop              <double>     Stop value for ECM algorithm parameter updates (default = 1e-5)" << std::endl
            << "  --prefix            <string>     Prefix for output files (default = diseq)" << std::endl
            << "  --quiet                          Flag to suppress output to stdout\n" << std::endl;
}

void alloSNPusage(){
  std::cerr << "\nUsage: ebg alloSNP [alloSNP options]\n" << std::endl;
  std::cerr << "Required options:" << std::endl
            << "  -f [--freqs-file]   <string>     File with reference allele frequencies" << std::endl
            << "  -n [--num-ind]      <int>        Number of individuals" << std::endl
            << "  -l [--num-loci]     <int>        Number of loci" << std::endl
            << "  -p1 [--ploidy1]     <int>        Ploidy level of subgenome one" << std::endl
            << "  -p2 [--ploidy2]     <int>        Ploidy level of subgenome two" << std::endl
            << "  -t [--total-reads]  <string>     Total read counts file" << std::endl
            << "  -a [--alt-reads]    <string>     ALT allele read counts file" << std::endl
            << "  -e [--error-rates]  <string>     Sequencing error rates file" << std::endl
            << "\nAdditional options:" << std::endl
            << "  --tol               <double>     Tolerance for Brent's method (default = 1e-10)" << std::endl
            << "  --iters             <int>        Number of iterations to run EM algorithm (default = 100)" << std::endl
            << "  --stop              <double>     Stop value for EM algorithm parameter updates (default = 1e-5)" << std::endl
            << "  --prefix            <string>     Prefix for output files (default = alloSNP)" << std::endl
            << "  --brent                          Flag to use Brent's method if the EM algorithm doesn't converge" << std::endl
            << "  --quiet                          Flag to suppress output to stdout\n" << std::endl;
}

void alloSNP2usage(){
  std::cerr << "\nUsage: ebg alloSNP2 [alloSNP2 options]\n" << std::endl;
  std::cerr << "Required options:" << std::endl
            << "  -f1 [--freqs-file1] <string>     File with parent 1 reference allele frequencies" << std::endl
            << "  -f2 [--freqs-file2] <string>     File with parent 2 reference allele frequencies" << std::endl
            << "  -n [--num-ind]      <int>        Number of individuals" << std::endl
            << "  -l [--num-loci]     <int>        Number of loci" << std::endl
            << "  -p1 [--ploidy1]     <int>        Ploidy level of subgenome one" << std::endl
            << "  -p2 [--ploidy2]     <int>        Ploidy level of subgenome two" << std::endl
            << "  -t [--total-reads]  <string>     Total read counts file" << std::endl
            << "  -a [--alt-reads]    <string>     ALT allele read counts file" << std::endl
            << "  -e [--error-rates]  <string>     Sequencing error rates file" << std::endl
            << "\nAdditional options:" << std::endl
            << "  --prefix            <string>     Prefix for output files (default = alloSNP)" << std::endl
            << "  --quiet                          Flag to suppress output to stdout\n" << std::endl;
}

void gatkUsage(){
  std::cerr << "\nUsage: ebg gatk [hwe options]\n" << std::endl;
  std::cerr << "Required options:" << std::endl
            << "  -n [--num-ind]      <int>        Number of individuals" << std::endl
            << "  -l [--num-loci]     <int>        Number of loci" << std::endl
            << "  -p [--ploidy]       <int>        Ploidy level" << std::endl
            << "  -t [--total-reads]  <string>     Total read counts file" << std::endl
            << "  -a [--alt-reads]    <string>     ALT allele read counts file" << std::endl
            << "  -e [--error-rates]  <string>     Sequencing error rates file" << std::endl
            << "\nAdditional options:" << std::endl
            << "  --prefix            <string>     Prefix for output files (default = gatk)" << std::endl
            << "  --print-probs                    Flag to alow printing of posterior probabilities of genotypes" << std::endl
            << "  --quiet                          Flag to suppress output to stdout\n" << std::endl;
}

void glUsage(){
  std::cerr << "\nUsage: ebg gl [hwe options]\n" << std::endl;
  std::cerr << "Required options:" << std::endl
            << "  -n [--num-ind]      <int>        Number of individuals" << std::endl
            << "  -l [--num-loci]     <int>        Number of loci" << std::endl
            << "  -p [--ploidy]       <int>        Ploidy level" << std::endl
            << "  -t [--total-reads]  <string>     Total read counts file" << std::endl
            << "  -a [--alt-reads]    <string>     ALT allele read counts file" << std::endl
            << "  -e [--error-rates]  <string>     Sequencing error rates file" << std::endl
            << "\nAdditional options:" << std::endl
            << "  --prefix            <string>     Prefix for output files (default = gl)" << std::endl
            << "  --quiet                          Flag to suppress output to stdout\n" << std::endl;
}
