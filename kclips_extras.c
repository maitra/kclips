#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "constants.h"
#include "random.h"
#include "array.h"

#define Inf 1e+140
#define SQ(x) ((x) * (x))
#define MAX(a, b) ((a) > (b) ? (a):(b))

SIZE_T assign_closest(double *X, SIZE_T p, SIZE_T nclass, double **Mu);

void get_ranges(int n, int p, double **X, double **ranges){

/* This function finds the limits of the ranges of the data.
   n: the # of obs in X
   p: the dimension of X
   X: the data
   ranges: a 2 x p matrix of minimums and maximums
*/

	int i,j;

	for(i=0;i<p;i++){
		ranges[0][i] = X[1][i];
		ranges[1][i] = X[1][i];

		for(j=1; j<n; j++){
			if(ranges[0][i] > X[j][i]){
				ranges[0][i] = X[j][i];
			}
			if(ranges[1][i] < X[j][i]){
				ranges[1][i] = X[j][i];
			}
		}
	}

	return;
}


double get_volume(int n, int p, double **x, double radius, int nstar, int K,
		  /*int exact, //never used in Ivan's code -- modified to 
		   be according to situation*/ 
		  double **mu 
		  /* , int nsim //approximation number of simulations, changed
		     to be SQ(p)*10000, internal to code */)
{
/* This function calculates the remaining p-dimensional "volume" after existing clusters have been found. */
	int i, k=0, texact /*taken out by RM, brought back by RM*/;
	double **ranges, uni=1.0, vol = 0., prop = 1.;

/* First get total rectanglish volume of data. */
	MAKE_MATRIX(ranges,2,p);
	
	get_ranges(n, p, x, ranges);
  
	for(i=0;i<p;i++)
		uni = uni * (ranges[1][i] - ranges[0][i]);
//	printf("uni = %f\n", uni); /*this is the rectangle volume */

	/* Now need to check if any spheres are outside of box. */
	texact=1; //RM dropped brought back by RM*/
	for(k = 0; k < K; k++){ //already dropped  by IR, brought back by RM
		for(i = 0; (i < p && texact); i++)
			if (MAX(ranges[1][i] - mu[k][i], mu[k][i] - ranges[0][i]) < radius)
				texact=0;
	}
	if (texact) { //RM took out and brought back
		double ptemp= 0.0, gam = 1.0, sphere = 0.0; //RM moved up */
		ptemp = (p/2.0+1.0);
		while(ptemp > 1.0){
			ptemp -= 1.0;
			gam *= ptemp;
		}
		if(ptemp != 1.0) 
			gam *= 0.5 * sqrt(PI);
		
		sphere = pow(radius,1.*p) * pow(PI, (p/2)) / gam;
//		printf("sphere = %f\n", sphere);
		vol = (uni-(double)K*sphere) * prop;
	}
	else{
		//    printf("Obtaining Approximation of Volume\n");
		SIZE_T j, count=0, flag, nsim = SQ(p)*10000;
		double  *unif,  *tdist;

		MAKE_VECTOR(unif, p);
		MAKE_VECTOR(tdist, K);
		for (i = 0; i < nsim; i++) {
			flag = 0;
			for (k = 0; k < K; k++) 
				tdist[k] = 0.0;
			for (j = 0; j < p; j++){
				unif[j] = runir(ranges[0][j], ranges[1][j]);
				for (k = 0; k < K; k++)
					tdist[k] += SQ(unif[j] - mu[k][j]);
			}
			for(k = 0; ((k < K) && !flag); k++)
				if (sqrt(tdist[k]) < radius) 
					flag=1;
			if (!flag) 
				count++;
		}
		FREE_VECTOR(unif);
		vol = uni*(double)(count)/(double)(nsim);

		FREE_VECTOR(tdist);

	}

	FREE_MATRIX(ranges);
	
	/*  printf("Volume = %e\n",vol);*/
	return vol;
}


/* This function copies the important information of a clustering to 
   the "optimal clustering".  */
void copy_optimal_cluster(int n,int p, int K, int *nc, int *optnc, int *class,
			  int *optclass, int *ikeep, int *optikeep,
			  double **mu, double **optmu, double **starting_mu,
			  double **optstarting_mu,double finW,double *optfinW)
{
	int i,j,k;
	for(i=0;i<n;i++){
		optclass[i]=class[i];
		optikeep[i]=ikeep[i];
	}
	for(k=0;k<K;k++){
		optnc[k]=nc[k];
		for(j=0;j<p;j++){
			optmu[k][j]=mu[k][j];
			optstarting_mu[k][j]=starting_mu[k][j];
		}
	}
	(*optfinW)=finW;
	return;
}

/* This function copies the important information of a clustering to 
   the "optimal clustering".  */
/* similar to copy_optimal_cluster but no keeping starting mu's. Not sure what 
   they do in the other routine */

void copy_optimal_details(int n, int p, int K, int *nc, int *optnc, 
			  int *class, int *optclass, int *ikeep, int *optikeep,
			  double **mu, double **optmu, double finW,
			  double *optfinW) 
{
	int i,j,k;
	for(i=0;i<n;i++){
		optclass[i]=class[i];
		optikeep[i]=ikeep[i];
	}
	for(k=0;k<K;k++){
		optnc[k]=nc[k];
		for(j=0;j<p;j++)
			optmu[k][j]=mu[k][j];
	}
	(*optfinW)=finW;
	return;
}


/* This function copies the important information of a clustering to 
   the "optimal clustering".  no ikeep here */

void copy_optimal_cluster_info(int n, int p, int K, int *nc, int *optnc, 
			       int *class, int *optclass, double **mu, 
			       double **optmu, double finW, double *optfinW) 
{
	int i,j,k;
	for (i = 0; i < n; i++) 
		optclass[i] = class[i];
	for (k = 0; k < K; k++) {
		optnc[k] = nc[k];
		for(j = 0; j < p; j++)
			optmu[k][j] = mu[k][j];
	}
	(*optfinW) = finW;
	return;
}


double **reducedata(double **x, int n, int p, int k, int *ikeep,int *class, 
		    int *rnc, int nstar, int *rclass)
{
/*  This function takes in a data set x and a vector ikeep containing
    the information of whether of not the point is scattered and returns a matrix
    that contains only the points that are kept.  It also updates the "cluster
    sizes" in rnc (a vector that needs to be passed into the function) and the 
    classification.
*/
	int i,j,istar=0;
	double **reduced;
	MAKE_MATRIX(reduced,nstar,p);

	for(i=0;i<n;i++){
		if(ikeep[i]){
			rclass[istar] = class[i];
			for(j=0;j<p;j++){
				reduced[istar][j]=x[i][j];
			}
			istar++;
		}
	}
	return reduced;
}

void vassign_closest(double **x, SIZE_T n, SIZE_T p, double **Mu, 
		     SIZE_T nclass, SIZE_T *ind, SIZE_T *nc)
{
/* This function assigns the observations in the matrix a to the closest 
   cluster based on the centers in the matrix Mu. rewritten by RM*/
	SIZE_T i;
	for (i = 0; i < nclass; i++) 
		nc[i] = 0;
	for(i = 0; i < n; i++) {
		ind[i] = assign_closest(x[i], p, nclass, Mu);
		nc[ind[i]]++;
	}
}



void zero_ivect(int n, int *vec)
{
	int i;
	for(i=0;i<n;i++) 
		vec[i]=0;
}

void zero_dvect(int n, double *vec)
{
	int i;
	for(i=0;i<n;i++) 
		vec[i]=0.0;
}

void zero_imat(int n,int p, int **mat)
{
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++)
			mat[i][j]=0;
	}
}
void zero_dmat(int n,int p, double **mat)
{
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<p;j++)
			mat[i][j]=0.0;
	}
}

double quickvolume(int n,int p,double **x,double radius,int nstar,int K,
		   double **mu,double rect)
{
	double vol=0.0,ptemp,gam,sphere;
	gam=1.;
	ptemp=(p/2.+1.);

	while(ptemp>1.){
		ptemp-=1.;
		gam*= ptemp;
	}

	if(ptemp<0.6) gam*=0.5*sqrt(PI);
  
	sphere=pow(radius,1.*p)*pow(PI,(p/2.))/gam;
  
	vol = rect-(double)K*sphere;

	return vol;
}

