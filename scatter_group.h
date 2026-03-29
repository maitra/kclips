#ifndef _SCATTER_GROUP_H
#define _SCATTER_GROUP_H
#include "mat_vec.h"

#define CEIL(a) ((int)(a) + (((a) > ((int)(a))) ? 1 : 0))

void vassign_closest(double **a, int m, int n, double **c, int k, int *ic1, 
		     int *nc);
int starts(int n,int m,double **x,int nclus,int *ningrp,
		   int *grpids,double **medoids, int elim,int iro);
void initials(double **x,int n,int p,int nclass,int *nc,double **Mu,
	      double **LTSigma,int *class);
void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);
double **reducedata(double **x, int n, int p, int k, int *ikeep,int *class, 
		    int *rnc, int nstar, int *rclass);
void get_ranges(int n, int p, double **X, double **ranges);
double get_volume(int n, int p, double **x, double radius, int nstar, int K,
		  double **mu);

double trace_crit(double **data,int n, int p, int K,double **mu,double rad,
		  int *ikeep,int *idns,double sig2);

double scatter_cluster(double **a, int m, int n, double **c, int k, int *ic1,
		       int *nc,int iter, double alpha, int *ikeep, int robust, 
		       /*int istarts,*/
		       double w4sb[], double *finalW, int nw);

double scatter_group(double **x, int n, int p, int nclass, double **Mu, int *nc,
		     int *class, double alpha, int *ikeep, int robust, 
		     int iters, double **medoids, int maxelim, int minelim,
		     double w4sb[], double *finalW, int nw, int traceW);

#endif /* _SCATTER_GROUP_H */
