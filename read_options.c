#include <getopt.h>
#include <ctype.h>
#include "constants.h"
#include <stdio.h>
#include <stdlib.h>
#include "read_options.h"

void usage(char *s)
{
	fprintf(stderr, "Usage: %s [-v] (flag for verbose, default no) \n -i input file \n -o output files suffix \n -n number of records/observations \n -m number of features/dimensions \n -k minimum number of groups \n [-K] maximum number of groups (default: same as k) \n [-N] maximum number of iterations (default: 100) \n [-s] number of short iterations (default: 10) \n [-d] should deterministic initialization be done in addition to random initialization? Default: no. Note: Only random initialization is done for cases wih m >= 9 because of the computational complexity involved in the initialization algorithm.\n [-I] number of random initalizations (for the random initializer), default: n*k*m\n [-r] robust scale estimate (flag, default no)\n [-a] alpha value (default, 0.05) \n [-S] should coordinates be scaled or not (flag, default none)\n [-W] should determinant of W be used in the calculation of the objective function (flag, default no: trace of W is used)\n", s);
	//"vi:o:n:m:k:K:I:t:s:r:a:W:S"
	fprintf(stderr, "Note: [.] indicates optional arguments\n\tfor details, see Maitra and Ramler (2009).\n");
}

short read_options(int agc, char **agv, char **infil, char **oufil, int *n, int *m, int *kmin, int *kmax, int *niter, int *ntries, int *shortiters, double *alpha, int *robust, int *initmeth, int *scale, int *tW, int *no_chg) {
	char c;
	short verbose = 0;

	*tW = 1;
	*no_chg = 1;
	*scale = 0;
	*niter = 1000;
	*shortiters = 10;
        *ntries = 0;
	*alpha = 0.05;
	*initmeth = 0;
	*robust = 0;
/*	double value, alpha = 0.05, finalW, **mu, w4sb[]={3., 4.5, 6., 7.5, 9.};*/
 //	nw = 5;
//	max_eliminull = -1;
//	min_eliminull = 0;

	/* Note that the colon after the character argument means that it has to have an argument, not that the character itself is optional */
	
	while ((c = getopt(agc, agv, "vi:o:n:m:k:K:dN:s:ra:Sc:WI:"))!= EOF) {
		switch (c) {
		case 'I':
			*ntries = atoi(optarg);
			break;
		case 'd':
			*initmeth = 1;
			break;
		case 'W':
			*tW = 0;
			break;
		case 'a':
			*alpha = atof(optarg);
			break;
		case 'r':
			*robust = 1;
			break;
		case 'v':
			verbose = 1;
			break;
		case 'i':
			*infil = optarg;
			break;
		case 'o':
			*oufil = optarg;
			break;
		case 'k':
			*kmin = atoi(optarg);
			break;
		case 'K':
			*kmax = atoi(optarg);
			break;
		case 'm':
			*m = atoi(optarg);
			break;
		case 'n':
			*n = atoi(optarg);
			break;
		case 'N':
			*niter = atoi(optarg);
			break;
		case 'c':
			*no_chg = (size_t) atoi(optarg);
			break;
		case 's':
			*shortiters = atoi(optarg);
			break;
		case 'S':
			*scale = 1;
			break;
		case '?':
			if (optopt == 'i' || optopt == 'o' || optopt == 'm' || optopt == 'k' || optopt == 'n' || optopt == 'N' || optopt == 'I' || optopt == 's' || optopt == 'a')
				fprintf (stderr, "Option -%c requires an argument.\n", optopt);
			else
				if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
			usage(agv[0]);
			exit(1);
		default:
			fprintf(stderr, "-i -o -m -n -k needed, -v -N -I -K -c -s -a -S -W -d optional\n");
			usage(agv[0]);
			exit(1);
		}
	}
	if (!*infil || !*oufil || !*kmin || !*m || !*n) {
		printf("Missing mandatory arguments\n");
		usage(agv[0]);
		exit(1);
	}
	if (!*kmax)
		*kmax = *kmin;
	if (verbose) {
		printf ("input file = %s, output class file = %s\n", *infil, *oufil);
		printf ("output means file = %s, mink = %i, maxk = %i max iterations = %i, no change iterations = %i\n", *oufil, *kmin, *kmax, *niter, *no_chg);
		printf ("scaling by range done for each coordinate = %d\n", *scale);
		printf ("number of random initializations = %d\n", *ntries);

	}
	return verbose;
}

