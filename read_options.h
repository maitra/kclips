#ifndef READOPTIONS_H
#define READOPTIONS_H
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <ctype.h>
#include <math.h>
#include "string.h"

void usage(char *s);

short read_options(int agc, char **agv, char **infil, char **oufil, int *n, int *m, int *kmin, int *kmax, int *niter, int *ntries, int *shortiters, double *alpha, int *robust, int *initmeth, int *scale, int *tW, int *no_chg); 

#endif /*READOPTIONS_H*/

