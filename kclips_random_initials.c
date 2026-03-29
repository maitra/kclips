/*
  This function provides a randomly initialized kclips run, with an evaluation
  of the objective function */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "array.h"
#include "constants.h"
#define PI 3.141593
#define Inf 1e+140
#define CEIL(a) ((int)(a) + (((a) > ((int)(a))) ? 1 : 0))

#include "scatter_cluster.h"
#include "rmath-stalone/rfns.h"

void kclips_rnd_initials(double **x, SIZE_T n, SIZE_T p, SIZE_T nclass, 
			 double **center, double prob, SIZE_T numnn);

void copy_optimal_details(int n, int p, int K, int *nc, int *optnc, 
			  int *class, int *optclass, int *ikeep, int *optikeep,
			  double **mu, double **optmu, double finW,
			  double *optfinW);

double kclips_rndinit_one_run(double **x, int n, int p, int nclass, double **Mu,
			      int *nc, int *class, double alpha, int *ikeep, 
			      int robust, int iters,
			      int maxelim, int minelim, double w4sb[], 
			      double *finalW, int nw, int traceW)
/*
  w4sb: the weights used in the robust bi-weight estimator of standard deviation
  (w in page 342, last line of Maitra and Ramler (2009).
  nw: the number of candidate w's in w4sb (5 by default: {3.,4.5,6.,7.5,9.}
  traceW: if trace rather than determinant of W is used. */
{
	double radk;
	SIZE_T flag = 0;
	double detW = Inf;
	int temp_elim;
	int *currkeep;

	MAKE_VECTOR(currkeep, n);

	if (maxelim < 0) 
		temp_elim = CEIL(n/(1.*p*nclass)); 
	else 
		temp_elim=maxelim;
//	printf("nclass = %d minelim = %d maxelim = %d temp_elim = %d\n", nclass, minelim, maxelim, temp_elim);
	
	while (!flag) {
		double currfinalW, currdetW;
		int nstar = 0;
		
		temp_elim = MAX(temp_elim, 2); /* guarantee at least two points
						  in the neighborhood to be a 
						  valid kclips initial candidate
					       */
		kclips_rnd_initials(x, n, p, nclass, Mu, 0.25, temp_elim);
		/* obtain initial candidate for the kclips algorithm*/

		radk = scatter_cluster(x, n, p, Mu, nclass, class, 
				       nc, iters, alpha, currkeep, robust, 
				       w4sb, &currfinalW, nw);
		for(int i = 0; i < n; i++)  
			nstar += currkeep[i]; 

		if (!nstar) 
			flag = 0; /* no points inside the cores -- not sure if
				     this can happen -- so we go again*/
		else {
			if (traceW)
				currdetW = -trace_crit(x, n, p, nclass, Mu, 
						       radk, currkeep, 
						       class,  
						       SQ(radk)/chisq_quantile(alpha, 
									       (double)n, 0, 0));
			else 
				currdetW = -det_crit(x, n, p, nclass, Mu, 
						     radk, currkeep, class);
			if (isfinite(currdetW)) {
					detW = currdetW;
					ikeep = currkeep;
					flag = 1;
				}
		}
	}
	FREE_VECTOR(currkeep);
	return detW;
}

	
	
			

double kclips_rndinit_one_run_orig(double **x, int n, int p, int nclass, double **Mu, 
			      int *nc, int *class, double alpha, int *ikeep, 
			      int robust, int iters,
			      int maxelim, int minelim, double w4sb[], 
			      double *finalW, int nw, int traceW)
/*
  
  POINTLESS and NOT USED any more.  Not sure what I was trying to do here earlier.

  w4sb: the weights used in the robust bi-weight estimator of standard deviation
  (w in page 342, last line of Maitra and Ramler (2009).
  nw: the number of candidate w's in w4sb (5 by default: {3.,4.5,6.,7.5,9.}
  traceW: if trace rather than determinant of W is used. */
{
	double radk;
	SIZE_T penalty = 1,  flag = 0;
	double **currMu, detW = Inf;
	int *currnc, *currclass, *currkeep, temp_elim, cnt = 0;
	
	MAKE_MATRIX(currMu, nclass, p);
	MAKE_VECTOR(currnc, nclass);
	MAKE_VECTOR(currclass, n);
	MAKE_VECTOR(currkeep, n);


	if (maxelim < 0) 
		temp_elim = CEIL(n/(1.*p*nclass)); 
	else 
		temp_elim=maxelim;
//	printf("nclass = %d minelim = %d maxelim = %d temp_elim = %d\n", nclass, minelim, maxelim, temp_elim);
	
	while (!flag && (temp_elim >= minelim)) {
		double currfinalW, currdetW;
		int nstar = 0;
		temp_elim = MAX(temp_elim, 2); /* guarantee at least two points
						  in the neighborhood to be a 
						  valid kclips initial candidate
					       */
		kclips_rnd_initials(x, n, p, nclass, currMu, 0.25, temp_elim);
		/* obtain initial candidate for the kclips algorithm*/


		radk = scatter_cluster(x, n, p, currMu, nclass, currclass, 
				       currnc, iters, alpha, currkeep, robust, 
				       w4sb, &currfinalW, nw);
		penalty++;
		for(int i = 0; i < n; i++)  
			nstar += currkeep[i]; 

		if (!nstar) 
			flag = 1; /* no points inside the cores -- not sure if
				     this can happen */
		else {
			if (traceW)
				currdetW = -trace_crit(x, n, p, nclass, currMu, 
						       radk, currkeep, 
						       currclass,  
						       SQ(radk)/chisq_quantile(alpha, 
									       (double)n, 0, 0));
			else 
				currdetW = -det_crit(x, n, p, nclass, currMu, 
						     radk, currkeep, currclass);
			if (currdetW < detW) {
				detW = currdetW;
				copy_optimal_details(n, p, nclass, currnc, nc, 
						     currclass, class, currkeep,
						     ikeep, currMu, Mu, 
						     currfinalW, finalW);
				cnt = 0;
			}
			else cnt++;
		}
		temp_elim -= CEIL(exp(log((double)CEIL(n/(1.*p*nclass)))-0.15*penalty));
//		printf("temp_elim = %d, detW = %f nstar = %d, penalty = %d\n", temp_elim, detW,nstar, penalty);
		if (cnt >= 10) /*no improvement in 10 consecutive runs*/
			flag = 1;
	}
	
	/* Now determine the remaining area. */
	FREE_MATRIX(currMu);
	FREE_VECTOR(currnc);
	FREE_VECTOR(currclass);
	FREE_VECTOR(currkeep);
	return detW;

}
	
double kclips_rndinit_initials(double **x, int n, int p, int nclass, 
			       double **Mu, int *nc, int *class, double alpha, 
			       int *ikeep, int robust, int shortiters, 
			       int longiter, int maxelim, int minelim, 
			       double w4sb[], double *finalW, int nw, 
			       int traceW, int ntries)
/*
  runs the kclips_rndinit_one_run ntries number of times (each to lax 
  convergence (measured by setting of shortiters, number of iterations) and
  then the best result up to longiter many times (to finer convergence).
  Returns the GOF measure according to if trace or determinant criterion is
  used.
*/
{
	double radk, **currMu, detW = Inf, currdetW, currfinalW;
	SIZE_T i, nstar = 0;
	int *currnc, *currclass, *currkeep;
	
	MAKE_MATRIX(currMu, nclass, p);
	MAKE_VECTOR(currnc, nclass);
	MAKE_VECTOR(currclass, n);
	MAKE_VECTOR(currkeep, n);

	for (i = 0; i < ntries; i++) {
		currdetW =  kclips_rndinit_one_run(x, n, p, nclass, currMu, 
						   currnc, currclass, alpha, 
						   currkeep, robust, 
						   shortiters, maxelim, 
						   minelim, w4sb,
						   &currfinalW, nw, traceW);
//		printf("currdetW = %f\n", currdetW);

		if (currdetW < detW) {
			detW = currdetW;
			copy_optimal_details(n, p, nclass, currnc, nc, 
					     currclass, class, currkeep, ikeep, 
					     currMu, Mu, currfinalW, finalW);
		}
	}

	FREE_MATRIX(currMu);
	FREE_VECTOR(currnc);
	FREE_VECTOR(currclass);
	FREE_VECTOR(currkeep);

	radk = scatter_cluster(x, n, p, Mu, nclass, class, nc, longiter, 
			       alpha, ikeep, robust, w4sb, finalW, nw);
	
	for (int i = 0; i < n; i++)  
		nstar += ikeep[i]; 
//	printf("nstar = %zu, radk = %f traceW = %d\n", nstar, radk, traceW);
//	if (!nstar) 
//	flag = 1; /* no points inside the cores -- not sure if this can happen */
//	else {
	if (nstar) {
		if (traceW)
			detW = -trace_crit(x, n, p, nclass, Mu, radk, ikeep, 
					   class, SQ(radk)/chisq_quantile(alpha, 
									  (double)n, 0, 0));
		else 
			detW = -det_crit(x, n, p, nclass, Mu, radk, ikeep, class);
	}
//	printf("detW = %f\n", detW);
	return detW;
}
