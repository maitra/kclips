#ifndef _KCLIPS_H
#define _KCLIPS_H

double  kclips(double **x, int n, int p, int mink, int maxk, double alpha, 
	       int *nc, int *classids, double **mu, double w4sb[], 
	       int shortiters, int iters, int initmeth, int minelim, 
	       int maxelim, int robust, double *finalW, int nw, int trW, 
	       int ntries, int *optK);


double  kclips_fixedK(double **x, int n, int p, int k, double alpha, int *nc,
		      int *classids, double **mu, double w4sb[], int shortiters,
		      int iters, int initmeth, int minelim, int maxelim, 
		      int robust, double *finalW, int nw, int trW, int ntries);

#endif /* _KCLIPS_H */
