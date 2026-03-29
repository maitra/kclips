/*
  This program has the utilities which does initialization using random sampling
 */

#include <stdio.h>
#include <stdlib.h>
#include "array.h"
#include "quantile.h"
#include "knnp.h"
#include "random.h"
#include "mat_vec.h"

#define SQ(x) ((x) * (x))

#define Inf 1e+140

void remove_duplicates(fdata **ma, SIZE_T *na, SIZE_T nc, SIZE_T dc)
{
	fdata *sa = NULL;
	KEY **key;
	SIZE_T *index = malloc(*na * sizeof index), i, k = 0;
	
	for (int i = 0; i < *na; i++) 
		index[i] = 0;

	key = knnp_bridge(ma, ma, sa, *na, *na, nc, 2, dc);

	for (i = 0; i < *na; i++)
		if (!index[i] && (key[i][1].d < 1e-10)) {
			if (index[key[i][1].k])
					index[i] = 1; /* this means that the 
							 "closest neighbor" 
							 has been identified as 
							 a duplicate, so the
							 current index should 
							 be tagged as the
							 duplicate. */
			else 
				index[key[i][1].k] = 1;
			/*current index is the first representative of the 
			  duplicate so should be included but not the duplicate
			  itself */
		}
	
	if (sa != NULL) 
		free1a(sa, *na, sizeof(fdata));
	
	for (i = 0; i < *na; i++)
		if (!index[i] && k!=i) {
			for (int j = 0; j < nc; j++)
				ma[k][j] = ma[i][j];
			k++;
		}

	free2a(key, *na, 2, sizeof(KEY));

	*na = k;

//	printf("k = %zu \n", k);

	free(index);
}



/*------------------------------------------------------------------------------------*/

KEY **knnp_bridge(fdata **ma, fdata **mb, fdata *sa, SIZE_T na, SIZE_T nb,
		  SIZE_T nc, SIZE_T nk, int dc)
{
	SIZE_T np = get_nb_cores(), cs = get_cache_size(); 
	/* see knnp.h for arguments here */
/*		SIZE_T np = 4, cs = 4096*1024; */
	return knnp(ma, mb, sa, na, nb, nc, nk, np, cs, dc, 0);
}

ThinSample thinned_sample(fdata **x, SIZE_T n, SIZE_T p, fdata prob, SIZE_T nk)
{
/*  nk = number of neighbors included inside the core of a cluster */
	fdata **ma, *sa = NULL, thrsh;
	SIZE_T dc = 3, i, j, k = 0;
	KEY **key;
	ThinSample Y;

//	printf("Entering thinned sample\n");
//	printf("nk = %zu\n", nk);

	MAKE_MATRIX(ma, n, p);

	for (i = 0; i < n; i++)
		for (j = 0; j < p; j++)
			ma[i][j] = x[i][j];

	remove_duplicates(ma, &n, p, dc);

	key = knnp_bridge(ma, ma, sa, n, n, p, nk + 1, dc);

	if (sa != NULL) 
		free1a(sa, n, sizeof(fdata));

	sa = malloc(n*sizeof(fdata));
	for (i = 0; i < n; i++) {
		sa[i] = 0;
//		printf("n = %zu crashing after\n", i);
		for (j = 1; j <= nk; j++) 
			sa[i] += key[i][j].d; /* total distance of all nk nearest neighbors of the ith observation*/
		sa[i] /= nk; /*average distance of all nk NN's of ith obs.*/
//		printf("sa = %f nk = %zu n = %zu\n", sa[i], nk, n);
	}
//	printf("crashing after here\n");
	free2a(key, n, nk + 1, sizeof(KEY));

	thrsh = select_kth_smallest(sa, n, prob*n + 1); 
/*kth smallest average distance among all observations: here k is prob*n+1.*/


	Y.p = p;  
	MAKE_VECTOR(Y.pi, n);
	MAKE_MATRIX(Y.x, n, Y.p);

	for (i = 0; i < n; i++)
		if (sa[i] <= thrsh) {
			for (j = 0; j < p; j++)
				Y.x[k][j] = ma[i][j];
			Y.pi[k] = sa[i]; /*average squared distance*/
			k++;
		}
	free(sa);
	FREE_MATRIX(ma);
	Y.n = k;

	Y.sum = 0;
	for (i = 0; i < Y.n; i++)
		Y.sum += 1/Y.pi[i];

//	printf("Exiting thinned sample\n");

	return Y;
}


double distance_from_closest_center(double *X, SIZE_T p, SIZE_T nclass, 
				    double **Mu)  
{
	SIZE_T j, l;
	double temp, dum1;

	temp = Inf;
	for (l = 0; l < nclass; l++) {
		dum1 = 0.;
		for (j = 0; j < p; j++)
			dum1 += SQ(X[j] - Mu[l][j]);
		if (dum1 < temp) 
			temp = dum1;
	}
	return temp;
}



void kclips_rnd_initials(double **x, SIZE_T n, SIZE_T p, SIZE_T nclass, 
			 double **center, double prob, SIZE_T numnn)
{
/* This function provides an initial candidate for the kclips algorithm taking 
   a random set of (valid) observations. */

	SIZE_T i, j, k;
	double u = runi(), sum = 0, cumsum = 0, *dist;
	ThinSample Y;

	Y = thinned_sample(x, n, p, prob, numnn);

	/* get the first seed w.p. inversely proportional to the average
	   distance to the numnn-nearest neighbors. */
	u *= Y.sum;
	
	for (i = 0; (i < Y.n) && (cumsum < u); i++)
		cumsum += 1/Y.pi[i];

	for (j = 0; j < p; j++)
		center[0][j] = Y.x[i-1][j]; /* selected the first center */
	
	/* now for all the other centers */
	
	MAKE_VECTOR(dist, Y.n);
	for (k = 1; k < nclass; k++)  {
		for (i = 0; i < Y.n; i++)
			dist[i] = distance_from_closest_center(Y.x[i], Y.p,
							       k, center);
		/* next center now obtained w.p. proportional to squared 
		   distance from closest center.*/
		sum = 0;
		for (i = 0; i < Y.n; i++)
			sum += dist[i];
		u = runi() * sum;
		cumsum = 0;
		for (i = 0; (i < Y.n) && (cumsum < u); i++)
			cumsum += dist[i];
		for (j = 0; j < p; j++) 
			center[k][j] = Y.x[i-1][j];
	}
	FREE_VECTOR(Y.pi);
	FREE_MATRIX(Y.x);
	FREE_VECTOR(dist);
}
