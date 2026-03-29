#ifndef SCATTER_CLUSTER_H
#define SCATTER_CLUSTER_H
#define PI 3.141593
#include "mat_vec.h"

double sbi(int n, double *x, double c, int mcent);

double det_crit(double **data,int n, int p, int K,double **mu,double rad,
		int *ikeep,int *idns);
int initials(double **x,int n,int p,int nclass,int *nc,
	     double **Mu,double **LTSigma,int *class);
double determinant(double *LTSigma,int n);
void get_ranges(int n, int p, double **X, double **ranges);
double get_volume(int n,int p,double **x,double radius,int nstar,int K,
		  /*int exact,*/ double **mu/*, int nsim*/);
double trace_crit(double **data,int n, int p, int K,double **mu,double rad,
		  int *ikeep,int *idns,double sig2);
double scatter_cluster_sbw(double **a, int m, int n, double **c, int k, 
			   int *ic1,int *nc, int iter, double alpha, 
			   int *ikeep, int robust, /*int istarts,*/
			   double *w4sb, double *finalW, int nw);
double scatter_cluster_sd(double **a, int m, int n, double **c, int k, 
			  int *ic1,int *nc, int iter, double alpha, 
			  int *ikeep, int robust, /*int istarts,*/
			  double *w4sb, double *finalW);
double scatter_cluster(double **a, int m, int n, double **c, int k, int *ic1,
		       int *nc, int iter, double alpha, int *ikeep, int robust, 
		       /*int istarts,*/
		       double w4sb[], double *finalW, int nw);

void vassign_closest(double **a, int m, int n, double **c, int k, int *ic1, 
		     int *nc);
#endif /*SCATTER_CLUSTER_H*/
