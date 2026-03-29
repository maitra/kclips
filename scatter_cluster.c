#include "array.h"
#include "scatter_cluster.h"
#include "quantile.h"
#include "rmath-stalone/rfns.h"

#define BIG 1e+40
#define Inf 1e+140

SIZE_T stabilize_fa(double *w, SIZE_T p);

double sampsd(double **a, int m, int n,  int k, double **c, int *ic1, int *nc, 
	      int irobust, double cnstnt)
{
/* This computes the sample standard deviation for a k cluster problem.  
   It assumes 
   homogeneous variance across clusters. 
   m = # of data points
   n = dimension of data
   k = # of clusters
   c = cluster means
   ic1 = current cluster id
   nc = vector of current cluster sizes
   irobust = 0 if usual sd calculated, 1 for bi-weight estimator
   cnstnt = constant to be passed through to biweight function (w of Maitra and Ramler (2009)
*/
	int i,j,l,ii;
	double sum=0.0;

	if(irobust){  /* Use this loop for biweight estimator of variance */
		double **tempa, **tempaT, *vecA, **medians;
		MAKE_MATRIX(medians,k,n);/*1*/
		MAKE_MATRIX(tempa, m, n);/*2*/ /* more memory than ever needed */
		MAKE_MATRIX(tempaT, n, m);/*3*/ /* more memory than needed 
						 but saving computer time*/

		for(l = 0; l < k; l++){ /* for calculate medians for each
					   coordinate and cluster */
			i = 0;
			for(ii = 0; ii < m; ii++){
				if(ic1[ii]==l){
					for(j=0;j<n;j++)
						tempa[i][j]=a[ii][j];
					i++;
				}
			}
			matrpose(tempa, nc[l], n, tempaT);

			for(j=0;j<n;j++){
				if (nc[l]) 
					medians[l][j]=median(tempaT[j], nc[l]);
				else
					medians[l][j] = 0;
			}
		}
		FREE_MATRIX(tempa);
		FREE_MATRIX(tempaT);/*3*/
		
		MAKE_VECTOR(vecA,(m*n));/*4*/
		for(i=0;i<m;i++){
			for(j=0;j<n;j++)
				vecA[n*i+j]=a[i][j]-medians[ic1[i]][j];
		}

		FREE_MATRIX(medians);/*1*/
   
		sum=sbi((m*n), vecA, cnstnt, 1);

		FREE_VECTOR(vecA);/*4*/
	}
	else {  /* Use this loop for "regular" (sd) estimate of variance */
		for(i = 0; i<m;i++){
			for(j=0;j<n;j++){
				sum+=(a[i][j]-c[ic1[i]][j])*(a[i][j]-c[ic1[i]][j]);
			}
		}
		sum=sqrt(sum/(m*n));
	}
	return sum;
}



double reduced_cluster_means(double **a, int m, int n, double **c, int k, 
			     int *ic1, int *nc,int *ikeep, double chival, 
			     int rob, double cnst)
/* This funtion recalculates the cluster centers based a only a portion of 
   the data
   from a given-% confidence region.
   m = # of observations
   n = dimensionality
   conf = the confidence level to use to create the region
   ikeep = length n vector that tells whether or not a point of inside or
   outside of a confidence region
   = 1 if in cluster, 0 otherwise
   rob = 1 if using robust estimator for variance 0 otherwise
   cnst = constant needed for robust estimator
*/
{
	int i,j,l, *tempnc;
	double sqdist = 0;
	double sd = sampsd(a, m, n, k, c, ic1, nc, rob, cnst);
	double radsq = SQ(sd)*chival;

        /* Find the square of the radius of the "confidence region" */
	/* Now decide which points in each cluster are in the confidence region. */
	MAKE_VECTOR(tempnc,k);/*2*/
	for (l=0;l<k;l++) 
		tempnc[l] = 0;
	
	for(i=0;i<m;i++){
		sqdist = 0.;
		for(j=0;j<n;j++)
			sqdist += SQ((a[i][j]-c[ic1[i]][j]));
		
		if(sqdist <= radsq)
			ikeep[i]=1;
		else 
			ikeep[i]=0;
		tempnc[ic1[i]] += ikeep[i];
	}
	
	/* Reset cluster centers for those still containing points. */
	for(l=0;l<k;l++)
		if (tempnc[l])
			for(j=0;j<n;j++)
				c[l][j]=0.0;
	
	for (i=0;i<m;++i) 
		if((ikeep[i]) && (tempnc[ic1[i]])) 
			for (j=0;j<n;++j) 
				c[ic1[i]][j] += a[i][j];
		
	for(l=0;l<k;l++)
		if(tempnc[l]>0)
			for (j=0;j<n;j++) 
				c[l][j] /= tempnc[l];

	FREE_VECTOR(tempnc);
	return sqrt(radsq);
}


double scatter_cluster(double **a, int m, int n, double **c, int k, int *ic1,
		       int *nc, int iter, double alpha, int *ikeep, 
		       int robust, double *w4sb, double *finalW, int nw)
/*		       int istarts, // never used by Ivan, always set to 1*/
{
	/* 
	   RM: This is just a wrapper to the correct function depending on 
	   whether robust is null or not, i.e.  whether the standard deviation
	   used to compute the cores is the usual sd or the bi-weight 
	   estimator.  
	   It returns the radius of the spheres of the clusters.
	   istarts decides if we determine convergence based on relative change
	   in radius (seems to be hard-coded in by default) or change in 
	   determinant of the W matrix (never used since calling istarts is 
	   always set at 1). So this has been dropped from Ivan's program.
	*/ 
	if (robust) 
		return scatter_cluster_sbw(a,m,n,c,k,ic1,nc,iter,alpha,ikeep,
					   robust, /*istarts,*/ w4sb, finalW, nw);
	else 
		return scatter_cluster_sd(a,m,n,c,k,ic1,nc,iter,alpha,ikeep, 
					  robust, /*istarts,*/w4sb, finalW);
}

double scatter_cluster_sbw(double **a, int m, int n, double **c, int k, 
			   int *ic1, int *nc, int iter, double alpha, 
			   int *ikeep, int robust, /*int istarts,*/
			   double w4sb[], double *finalW, int n_w)
{
/* This function calls the assign_closest function and reduced cluster means 
   to iteratively reassign the points. */
  
	int i=0, w, i1, i2;
	double chival, newrad=0., optw=3.0,optval=BIG, eps = 1e-5, diff=0.,
		curcrit,oldcrit=0.,tempW,**c_init;

	chival = chisq_quantile(alpha,(double)(n),0,0);  
	/* Calls local rmath lib */

	MAKE_MATRIX(c_init,k,n);
	/* Need to provide same starting values for all runs of w4sb, hence
	 keep these aside */
	for(i1=0;i1<k;i1++)
		for(i2=0;i2<n;i2++)
			c_init[i1][i2]=c[i1][i2];

	for(w = 0; w < n_w; w++){
		/* reset "kill" parameters  */
		diff=1.*BIG;
		curcrit=-10.*BIG;
		oldcrit=-1.*BIG;
		newrad=0.;
		tempW=w4sb[w]; 

		for(i1 = 0; i1 < k; i1++){
			for(i2 = 0; i2 < n; i2++)
				c[i1][i2] = c_init[i1][i2];
		}

		if (k == 1) {
			for (i = 0; i < m; i++)
				ic1[i] = 0;
			nc[0] = m;
		}
		else
			vassign_closest(a, m, n, c, k, ic1, nc);

		for (i = 0; ((i < iter) && (diff > eps)); i++) {
			newrad = reduced_cluster_means(a, m, n, c, k, ic1, nc,
						     ikeep, chival, robust, 
						     tempW);
			if (k == 1) {
				for (i = 0; i < m; i++)
					ic1[i] = 0;
				nc[0] = m;
			}

			if (k > 1)
				vassign_closest(a, m, n, c, k, ic1, nc);
			curcrit=newrad;
			diff=fabs((curcrit-oldcrit)/curcrit);
			oldcrit=curcrit;
		}
		
		if(curcrit<optval){
			optval=curcrit;
			optw=tempW;
		}
	}
	/* reset "kill" parameters  */
	eps=0.00001; 
	diff=1.*BIG;
	curcrit=-10.*BIG;
	oldcrit=-1.*BIG;
	newrad=0.;
	(*finalW)=optw;
	for(i1=0;i1<k;i1++)
		for(i2=0;i2<n;i2++)
			c[i1][i2]=c_init[i1][i2];
	vassign_closest(a, m, n, c, k, ic1, nc);
	
	for (i = 0; ((i < iter) && (diff > eps)); i++) {
		newrad=reduced_cluster_means(a, m, n, c, k, ic1, nc, ikeep, 
					     chival, robust, optw);
		vassign_closest(a, m, n, c, k, ic1, nc);
		curcrit=newrad;
		diff=fabs((curcrit-oldcrit)/curcrit);
		oldcrit=curcrit;
	}

	FREE_MATRIX(c_init);

	return newrad;
} /* end scatter_cluster_sbw */

double scatter_cluster_sd(double **a, int m, int n, double **c, int k, 
			  int *ic1, int *nc, int iter, double alpha, 
			  int *ikeep,int robust, /*int istarts,*/double w4sb[], 
			  double *finalW)
{
/* This function calls the assign_closest function and reduced cluster means 
   to iteratively reassign the points. */
  	int i=0;
	double chival = 0., eps=0.00001,diff=BIG,curcrit,oldcrit=-1000000000.,
		newrad = 0;

	chival=chisq_quantile(alpha,(double)(n),0,0);   /* Calls Rmath lib */
	vassign_closest(a, m, n, c, k, ic1, nc);


	for (i = 0; ((i < iter) && (diff > eps)); i++) {
		newrad=reduced_cluster_means(a, m, n, c, k, ic1, nc, ikeep, 
					     chival, robust, 9.0);
		curcrit=newrad;
		diff=fabs((curcrit-oldcrit)/curcrit);
		oldcrit=curcrit;
	}
	return newrad;
} /* end scatter_cluster_sd */

double det_crit(double **data, int n, int p, int K,double **mu,double rad,
		int *ikeep, int *idns)
{
	double val,**fmu,**fLTSigma,*w,detW,vol;
	int i, j, nstar=0, *fnc, *iddc;

	MAKE_VECTOR(fnc,K+1);

	for(i = 0; i < (K + 1);i++) 
		fnc[i]=0;
	MAKE_VECTOR(iddc,n);
	for(i=0;i<n;i++){
		nstar+=ikeep[i];
		iddc[i]=ikeep[i]*(idns[i]+1);
		fnc[iddc[i]]++;
	}

	/* First get W matrix  */
	MAKE_MATRIX(fmu,K+1,p);
	MAKE_MATRIX(fLTSigma,K+1,p*(p+1)/2);

/*	printf("nstar = %d\n", nstar);

	print_ivector(fnc, K+1, " %d ");*/
	initials(data,n,p,K+1,fnc,fmu,fLTSigma,iddc);
	
	FREE_VECTOR(iddc);
	MAKE_VECTOR(w,p*(p+1)/2);
	for(j=0;j<p*(p+1)/2;j++) 
		w[j]=0.0;

	for(i=1;i<(K+1);i++) {
		for (j=0;j<(p*(p+1)/2);j++) 
			w[j]+=(fnc[i]-1)*fLTSigma[i][j];
	}

	i = logLTdeterminant(w, p, &detW);

	if (i) {
		if (!stabilize_fa(w, p)) {
			i = logLTdeterminant(w, p, &detW);
		}
		else
			detW=Inf;
	}

/*		detW = Inf;*/
	
	/* Now get remaining area (i.e. volume of region)  */

	vol = get_volume(n, p, data, rad, nstar, K, mu); /* changed to reflect
							    simplifications in
							    get_volume*/


	val = nstar/(1.*n) * detW/2 + (1-nstar/(1.*n))*log(vol);

	FREE_VECTOR(w);
	FREE_MATRIX(fLTSigma);
	FREE_MATRIX(fmu);
	FREE_VECTOR(fnc);

	return val;

}


double trace_crit(double **data,int n, int p, int K, double **mu, double rad,
		  int *ikeep, int *idns, double sig2)
{

	int i, j, *idtc, nstar = 0;
	double traceW=0., vol, val;

//	printf("sig2 = %f\n", sig2);
	
	MAKE_VECTOR(idtc, n);
	for (i = 0; i < n; i++){
		nstar += ikeep[i];
		idtc[i] = ikeep[i] * (idns[i] + 1);
		if (ikeep[i])
			for(j = 0; j < p; j++)
				traceW += SQ(data[i][j] - mu[idtc[i]-1][j]);
	}

	vol = get_volume(n, p, data, rad, nstar, K, mu); /* changed to reflect
							    changes in get_vol*/
//	printf("K = %d traceW = %f vol = %f\n", K, traceW, vol);
     
	val=nstar*p/2.0*log(2*PI+sig2)+1./(2.*sig2)*traceW+(n-nstar)*log(vol);  
	FREE_VECTOR(idtc);
	return val;
}

