#include "array.h"
#include "scatter_cluster.h"
#include "rmath-stalone/rfns.h"

#define BIG 1e+40

	   
double sampd(double **a, int m, int n,  int k, double **c, int*ic1, int *nc, 
	     int irobust, double cnstnt);
double reduced_cluster_means(double **a, int m, int n, double **c, int k, 
			     int *ic1, int *nc,int *ikeep, double chival, 
			     int rob, double cnst);
void vassign_closest(double **a, int m, int n, double **c, int k, int *ic1, 
		     int *nc);

double scatter_cluster_starts(double **a, int m, int n, double **c, int k, 
			      int *ic1,int *nc, int iter, double alpha, 
			      int *ikeep, int robust, double constant)
{
/* This function calls the assign_closest function and reduced cluster means 
   to iteratively reassign the points. */
  
	int i=0;
	double chival=0.,eps=0.00001,diff=BIG,curcrit,oldcrit=-1000000000.,
		newrad=0.;
	chival = chisq_quantile(alpha,(double)(n),0,0);   /* Calls Rmath lib */

	vassign_closest(a, m, n, c, k, ic1, nc);
/*  printf("Running CIPS\n");*/

	while((i<iter)&&(diff>eps)){
		newrad=reduced_cluster_means(a, m, n, c, k, ic1, nc, ikeep, chival, 
					     robust, constant);
		vassign_closest(a, m, n, c, k, ic1, nc);

		curcrit=newrad;
		diff=fabs((curcrit-oldcrit)/curcrit);
		oldcrit=curcrit;
		i++;
	}

	return newrad;
}


