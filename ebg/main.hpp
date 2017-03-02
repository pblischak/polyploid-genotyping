#ifndef MAIN_HPP
#define MAIN_HPP

void mainUsage();
void hweUsage();
void diseqUsage();
void alloSNPusage();
void gatkUsage();
void glUsage();

#define VERSION "0.2.2-alpha"
#define VERSIONDATE "March 2017"

#define MISSING -9
#define BADLIK  -9999.999

const double max_small = 1.0e-100; // smallest allowable value for EM algorithms.
const double bad_lik   = -9999.999;
const int    missing   = -9;

extern MbRandom *r;

#endif //MAIN_HPP
