#ifndef MAIN_HPP
#define MAIN_HPP

void mainUsage();
void hweUsage();
void diseqUsage();
void alloSNPusage();
void alloSNP2usage();
void gatkUsage();
void glUsage();

#define VERSION "0.4.0-alpha"
#define VERSIONDATE "July 2020"

#define MISSING -9
#define BADLIK  -9999.999

const double max_small  = 1.0e-300; // smallest allowable value for EM algorithms.
const double log_max_small = -300;
const double big_positive = 10e300;
const double big_negative = -10e300;
const double min_freq = 1.0e-20;
const double max_freq = 1.0 - min_freq;
const double almost_one = 1 - max_small;
const double bad_lik    = -9999.999;
const int    missing    = -9;

extern MbRandom *r;

#endif //MAIN_HPP
