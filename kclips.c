#include<stdio.h>
#include<stdlib.h>
#include "array.h"
#include "random.h"
#include <math.h>
#include "scatter_group.h"
#include <time.h>

#define Inf 1e+140

void copy_optimal_cluster_info(int n, int p, int K, int *nc, int *optnc, 
			       int *class, int *optclass, double **mu, 
			       double **optmu, double finW, double *optfinW);

void copy_optimal_details(int n, int p, int K, int *nc, int *optnc, 
			  int *class, int *optclass, int *ikeep, int *optikeep,
			  double **mu, double **optmu, double finW,
			  double *optfinW);

double kclips_rndinit_initials(double **x, int n, int p, int nclass, 
			       double **Mu, int *nc, int *class, double alpha, 
			       int *ikeep, int robust, int shortiters, 
			       int longiter, int maxelim, int minelim, 
			       double w4sb[], double *finalW, int nw, 
			       int traceW, int ntries);

double  kclips_fixedK(double **x, int n, int p, int k, double alpha, int *nc,
		      int *classids, double **mu, double w4sb[], int shortiters,
		      int iters, int initmeth, int minelim, int maxelim, 
		      int robust, double *finalW, int nw, int trW, int ntries)
{
	double temp;
	int *keep, nstar = 0, i;

	if (!ntries)
		ntries = n * k * p;
	
	MAKE_VECTOR(keep, n);
	
	if ((initmeth) && (p < 9))  { /*both initmethods*/
		double **medoids;
		double currfinalW, **currmu, currtemp;
		int *currclassids, *currkeep, currnstar = 0, *currnc;
		
		MAKE_MATRIX(medoids, k, p);
	
		temp = scatter_group(x, n, p, k, mu, nc, classids, alpha,
				     keep, robust, iters, medoids, maxelim,
				     minelim, w4sb, finalW, nw, trW);
		
		for (i = 0; i < n; i++) 
			nstar += keep[i];
//		printf("nstar = %d, temp = %f\n", nstar, temp);

//		temp += (1.418 - 0.5*log(1.*nstar))*nstar*p/(1.*n) + 1.0*log(k+1.0/p);
//		printf("temp = %f\n", temp);

		if(temp <=(-.1*Inf)) 
			temp=Inf;

		FREE_MATRIX(medoids);
		
		MAKE_MATRIX(currmu, k, p);
		MAKE_VECTOR(currkeep, n);
		MAKE_VECTOR(currnc, k);
		MAKE_VECTOR(currclassids, n);
		
		currtemp = kclips_rndinit_initials(x, n, p, k, currmu, 
						   currnc, currclassids,
						   alpha, currkeep,
						   robust, shortiters,
						   iters, maxelim,
						   minelim, w4sb,
						   &currfinalW, nw, 
						   trW, ntries);
		for (i = 0; i < n; i++)
				currnstar += currkeep[i];
		
//		currtemp += (1.418 - 0.5*log(1.*currnstar))*currnstar*p/(1.*n) + 1.0*log(k+1.0/p);
			
		if(currtemp <=(-.1*Inf)) 
			currtemp=Inf;
		
		if (currtemp < temp) {
			temp = currtemp;
			
			copy_optimal_details(n, p, k, currnc, nc, 
					     currclassids, classids, 
					     currkeep, keep, currmu, 
					     mu, currfinalW, finalW);
		}				
		
		FREE_MATRIX(currmu);
		FREE_VECTOR(currkeep);
		FREE_VECTOR(currnc);
		FREE_VECTOR(currclassids);
	}
	else {
		temp = kclips_rndinit_initials(x, n, p, k, mu, nc, classids,
					       alpha, keep, robust, shortiters,
					       iters, maxelim, minelim, w4sb, 
					       finalW, nw, trW, ntries);
		for (i = 0; i < n; i++)
			nstar += keep[i];
			
//		temp += (1.418 - 0.5*log(1.*nstar))*nstar*p/(1.*n) + 1.0*log(k+1.0/p);
		if(temp <= (-.1*Inf)) 
			temp=Inf;
	}

	for (i = 0; i < n; i++) {
		classids[i]++;
		classids[i] *= keep[i];
	}

	FREE_VECTOR(keep);
	return temp;
}

double  kclips(double **x, int n, int p, int mink, int maxk, double alpha, 
	       int *nc, int *classids, double **mu, double w4sb[], 
	       int shortiters, int iters, int initmeth, int minelim, 
	       int maxelim, int robust, double *finalW, int nw, int trW, 
	       int ntries, int *optK)
{
	double temp, **currmu;
	int *currnc, *currclassids;

	*optK = mink;

//	MAKE_MATRIX(mu, mink, p);
//	MAKE_VECTOR(nc, mink);
	temp = kclips_fixedK(x, n, p, mink, alpha, nc, classids, mu, w4sb, 
			     shortiters, iters, initmeth, minelim, 
			     maxelim, robust, finalW, nw, trW, ntries);
	

	MAKE_MATRIX(currmu, maxk, p);
	MAKE_VECTOR(currnc, maxk);
	MAKE_VECTOR(currclassids, n);	

	for (int k = mink + 1; k <= maxk; k++)	{
		double currtemp, currfinalW;

		currtemp = kclips_fixedK(x, n, p, k, alpha, currnc, 
					 currclassids, currmu, w4sb, shortiters,
					 iters, initmeth, minelim, maxelim, 
					 robust, &currfinalW, nw, trW, ntries);
//		printf("k = %i, temp = %f\n", k, currtemp);
		if (currtemp < temp) {
			temp = currtemp;
			(*optK) = k;
			copy_optimal_cluster_info(n, p, k, currnc, nc, 
						  currclassids, classids,
						  currmu, mu, currfinalW, 
						  finalW); 
		}
	}
	FREE_VECTOR(currclassids);
	FREE_VECTOR(currnc);
	FREE_MATRIX(currmu);
	
	return temp;
}
