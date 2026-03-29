#include<stdio.h>
#include<stdlib.h>
#include "array.h"
#include "random.h"
#include <math.h>
#include "scatter_group.h"
#include "kclips.h"
#include <time.h>
#include "read_options.h"

#define Inf 1e+140

/*void preprocess_data(double **data, int n, int m);*/

int main(int argc, char **argv)
{
	time_t t1; /*,t2,t3;*/
	unsigned int SEED1; /* SEED2;*/

	short verb;
	int i, j, m = 0, nn = 0, robust, iterations, shortiters, no_chg, scale;
	int max_eliminull, min_eliminull, nw = 5;
	int minK = 0, maxK = 0, trW, optK, ntries, initmeth;
	double value, alpha, finalW, w4sb[]={3., 4.5, 6., 7.5, 9.};
	FILE *finp;
	char *infil = NULL, *otxt = NULL, *otext, *mtext, *ptext, *jtext;

	time(&t1);
	SEED1 = t1;
/*	time(&t2);
	time(&t3);
	SEED2 = t2*t3;
	set_seed(SEED1,SEED2); */

/*	printf("seed = %zu %ld %ld\n", SEED1, random(), random());*/
	setseed(&SEED1);
/*  	printf("seed = %zu %ld\n", SEED1, random());*/

	max_eliminull = -1;
	min_eliminull = 0;

	verb = read_options(argc, argv, &infil, &otxt, &nn, &m, &minK, &maxK, &iterations, &ntries, &shortiters, &alpha, &robust, &initmeth, &scale, &trW, &no_chg);
	
	if (verb)
		printf ("input file = %s, output file prefix = %s\n", infil, otxt);

	if (verb)
		printf("%d %d %d %d %d\n", nn, m, minK, maxK, ntries);

	otext = malloc((strlen(otxt) + strlen("-class.out") + 10) * sizeof(char));
	mtext = malloc((strlen(otxt) + strlen("-means.out") + 10) * sizeof(char));
	ptext = malloc((strlen(otxt) + strlen("-robustscalepar.out") + 10) * sizeof(char));
	jtext = malloc((strlen(otxt) + strlen("-obj.out") + 10) * sizeof(char));

	if (!otext || !mtext || !ptext || !jtext) {
		fprintf(stderr, "Could not allocate memory for output filenames\n");
		free(otext);
		free(mtext);
		free(ptext);
		exit(1);
	} else {
		int *class, *nc;
		double **data, **mu;
		FILE *fout;
		
		sprintf(otext, "%s-class.out", otxt);
		sprintf(mtext, "%s-means.out", otxt);
		sprintf(ptext, "%s-robustscalepar.out", otxt);
		sprintf(jtext, "%s-obj.out", otxt);

		if (verb)
			printf("Begin reading in data\n");
		
		MAKE_MATRIX(data, nn, m);
		MAKE_VECTOR(class, nn);
		if (verb)
			printf("%s\n", infil);
		if ((finp=fopen(infil, "r"))==NULL) {
			printf("Error opening input file: is it correctly located?\n");
			printf("Exiting...\n");
			exit(1);
		}
		for (i = 0; i < nn; i++) 
			for (j = 0; j < m; j++) 
				fscanf(finp,"%lf",&data[i][j]);
		fclose(finp);
		
		if (verb)
			printf("Done reading in data\n");
		
/*		if (preprocess) 
			preprocess_data(data, nn, m); available in preprocess.c*/
		
		MAKE_MATRIX(mu, maxK, m);
		MAKE_VECTOR(nc, maxK);
		
		value = kclips(data, nn, m, minK, maxK, alpha, nc, class, mu, w4sb, 
			       shortiters, iterations, initmeth, min_eliminull, 
			       max_eliminull, robust, &finalW, nw, trW, ntries, &optK);
		
/*
  w4sb: the weights used in the robust bi-weight estimator of standard deviation
  (w in page 342, last line of Maitra and Ramler (2009).
  nw: the number of candidate w's in w4sb (5 by default: {3.,4.5,6.,7.5,9.}
  traceW: if trace rather than determinant of W is used.


*/
		printf("optK = %d value of objective function = %f\n", optK, value);
		fout = fopen(jtext, "w");
		fprintf(fout,"%d %f", optK, value);
		fclose(fout);

		fout = fopen(otext, "w");
		for (i = 0; i < nn; i++) 
			fprintf(fout, "%d ", class[i]);
		fclose(fout);
	
		fout = fopen(mtext, "w");
		for (i = 0; i < optK; i++) {
			for (j = 0; j < m;j++) 
				fprintf(fout, "%f ", mu[i][j]);
			fprintf(fout, "\n");
		}  
		fclose(fout);
		
		fout = fopen(ptext, "w");
		if (robust)
			fprintf(fout,"%f ", finalW);
		fclose(fout);

		FREE_VECTOR(nc);
		FREE_MATRIX(mu);
		FREE_VECTOR(class);  
		FREE_MATRIX(data);
		free(otext);
		free(mtext);
		free(ptext);
		free(jtext);
	}
	

	return EXIT_SUCCESS;
}
// ./run_kclips_estk -i /home/maitra/Datasets/att-faces/temp.txt -o /home/maitra/Datasets/att-faces/kclips-results/radonDWT -n 400 -m 280 -k 1 -K 100
