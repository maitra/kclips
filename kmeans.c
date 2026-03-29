/**
 * @file kmeans.c
 *
 * Implement Hartigan and Wong's Algorithm AS 136 published in Applied 
 * Statistics (1979) 28(1):100
 *
 * Code based on Applied Statistics algorithms (C) Royal Statistical Society
 * 1979. Adapted for C by Ranjan Maitra, Baltimore, 07/12/02
 * identical to kmns.c except for indx which is a number here.
 * Further adjusted by K. Dorman, Ames, IA, 9/12.
 */

#include <float.h>
#include <string.h>
#include "array.h"
#include "kmeans.h"
#include "util.h"



const char *KMEANS_ERROR_STRING[KMEANS_NUMBER_FAULTS] = {
	"Success",				/* first 4 are #kmeans() */
	"Inferred null cluster",		/* return codes */
	"Requested K=1 or K>n, more clusters than observations",
	"Unused error",
	"Initialization method failed",
	"Number of errors",
	"Exceeded maximum iterations"
};


void optra(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, SIZE_T *ic1,
	SIZE_T *ic2, SIZE_T *nc, double *an1, double *an2, SIZE_T *live,
	SIZE_T *cd, double *d, SIZE_T *indx);

void qtran(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, SIZE_T *ic1,
	SIZE_T *ic2, SIZE_T *nc, double *an1, double *an2, SIZE_T *live,
	double *d, SIZE_T *indx);

/**
 * Hartigan and Wong's algorithm.
 *
 * @param a mXn data of m observations, n dimensions
 * @param m number of observations
 * @param n dimension of data
 * @param c matrix of
 * @param k number of clusters
 * @param icl cluster assignment of each observation
 * @param nc number of members in each cluster
 * @param iter maximum number of iterations
 * @param wss k-means criterion minimized, within-sums-of-squares
 * @param ifault error code
 * @return error status
 */
void
kmeans(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, SIZE_T *ic1,
	    SIZE_T *nc, SIZE_T iter, double *wss, int *ifault)
{
	SIZE_T i, j, l;		/* indices */
	SIZE_T indx = 0;	/* no. observations evaluated without swap */
	SIZE_T *live;		/* live set information */
	SIZE_T *cd;		/* recompute distance */
	SIZE_T *ic2;		/* second closest centers */
	double db;		/* to calculate distances */
	double *d;		/* penalized distance to closest center for each obs. */
	double dt[2];		/* distance to closest & next closest centers */
	double *an1;		/* for each cluster nk/(nk+1) */
	double *an2;		/* for each cluster nk/(nk-1) */
	
	*ifault = KMEANS_NO_ERROR;

	/* invalid or trivial input */
	if (k <= 1 || k >= m) {

		/* user requests 1 cluster */
		if (k == 1) {
			nc[0] = m;
			/* compute centers [PARALLEL] */
			for (j = 0; j < n; j++) 
				c[0][j] = 0;	/* initialize */
			for (i = 0; i < m; i++) {
				ic1[i] = 0;	/* assign obs. i to cluster 1 */
				for (j = 0; j < n; j++) 
					c[0][j] += a[i][j];
			}
			for (j = 0; j < n; j++) 
				c[0][j] /= m;	/* mean */
			/* compute within-sum-of-squares [PARALLEL] */
			wss[0] = 0;
			for (i = 0; i < m; i++)
				for (j = 0; j < n; j++) 
					wss[0] += SQ(a[i][j] - c[0][j]);

		/* user requests m clusters */
		} else if (k == m)
			for (i = 0; i < m; i++) {
				nc[i] = 1;	/* 1 member per cluster */
				memcpy(c[i], a[i], n * sizeof **c);
				ic1[i] = i;
				wss[i] = 0;
			}

		/* user requests too few or too many clusters */
		else 
			*ifault = KMEANS_CALLER_INPUT_ERROR;
	} /* handle k<=1 or k>= m */
	else { /* now handle the more interesting case of 1 < k < m */
	/* allocate second-closest center vector */
		MAKE_VECTOR(ic2, m);
		
		/* find two closest centres, ic1[i] and ic2[i], for obs. i [PARALLEL] */
		for (i = 0; i < m; i++) {
			
		/* initially guess first and second center */
			ic1[i] = 0;
			ic2[i] = 1;
			
			/* compute distances to first and second */
			for (l = 0; l < 2; ++l) {
				dt[l] = 0.;
				for (j = 0; j < n; ++j)
					dt[l] += SQ(a[i][j] - c[l][j]);
			}
			
			/* swap first and second, if second closer */
			if (dt[0] > dt[1]) {
				ic1[i] = 1;
				ic2[i] = 0;
				db = dt[0];
				dt[0] = dt[1];
				dt[1] = db;
			}
			
			/* for remaining centers */
			for (l = 2; l < k; ++l) {
				
				/* compute distance */
				db = 0.;
				for (j = 0; j < n; ++j) {
					db += SQ(a[i][j] - c[l][j]);
					if (db > dt[1])	/* losing already */
						break;
				}
				
				/* it is closer */
				if (db < dt[1]) {
					
					/* second closest so far */
					if (db >= dt[0]) {
						dt[1] = db;
						ic2[i] = l;
						/* closest so far */
					} else {
						dt[1] = dt[0];
						ic2[i] = ic1[i];
						dt[0] = db;
						ic1[i] = l;
					}
				}
			}
		} /* find closest centers */
		
		/* compute new cluster centroids [PARALLEL] */
		for (l = 0; l < k; ++l) {	/* initialize */
			nc[l] = 0;
			for (j = 0; j < n; ++j)
				c[l][j] = 0.;
		}
		for (i = 0; i < m; ++i) {	/* compute */
			nc[ic1[i]]++;
			for (j = 0; j < n; ++j)
				c[ic1[i]][j] += a[i][j];
		}
		
		/* allocate data */
		MAKE_VECTOR(an1, k);
		MAKE_VECTOR(an2, k);
		MAKE_VECTOR(d, m);
		MAKE_VECTOR(live, k);
		MAKE_VECTOR(cd, k);
		
		/* initialize pre-computed values needed by algorithm */
		for (l = 0; l < k; ++l) {
		
			/* error if initialization routine produces empty cluster */
			if (nc[l] == 0) {
				*ifault = KMEANS_NULL_CLUSTER_ERROR;
				goto FINISHED;
			}
			
			/* finish centroid calculation */
			for (j = 0; j < n; ++j)
				c[l][j] /= nc[l];
			
			/* dis(i,l) * ans2[l] = cost of add obs. i to cluster l */
			an2[l] = nc[l] / (nc[l] + 1.0);
			
			/* dis(i,1) * ans1[l] = cost of keeping obs. i in cluster l */
			an1[l] = DBL_MAX;	/* infinite cost if l singlet obs. i */
			if (nc[l] > 1)
				an1[l] = nc[l] / (nc[l]-1.0);
			
			/* live[l] stores obs. last transferred to/from this cluster */
			live[l] = m + 1;	/* all clusters start live */
			cd[l] = m + 1;		/* all distances need computing */
		}
		
		/* iterate optimal and quick transfer stages until max. iterations */
		for (i = 0; i < iter; i++) {
			/* optimal transfer stage: each point reallocated,
			 * if appropriate to cluster which best reduces wss
			 */
			optra(a, m, n, c, k, ic1, ic2, nc, an1, an2, live, cd, d, &indx);
			
			/* m observations processed and no change, we're done */
			if (indx == m)
			break;
			
			/* quick transfer stage: consider moving each point i to ic2[i]
			 * and keep looping until no more change.
			 */
			qtran(a, m, n, c, k, ic1, ic2, nc, an1, an2, live, d, &indx);
			
			/* two clusters: no need for additional optimal transfer */
			if (k==2)
				break;
		} 
		
		/* maximum iterations exceeded */
		if ((indx != m) && (k != 2))
			*ifault = KMEANS_EXCEED_ITER_WARNING;
		
		/* compute centroids and wss for each cluster [PARALLEL] */
		/* KSD TODO: centroids already computed! */
		for (l = 0; l < k; ++l) {	/* initialize */
			wss[l] = 0.;
			for (j = 0; j < n; ++j)
			c[l][j] = 0.;
		}
		
		for (i = 0; i < m; ++i)		/* compute center */
			for (j = 0; j < n; ++j)
				c[ic1[i]][j] += a[i][j];
		
		for (j = 0; j < n; ++j) {	/* compute wss */
		for (l = 0; l < k; ++l)
			c[l][j] /= (double) nc[l];
		for (i = 0; i < m; ++i)
			wss[ic1[i]] += SQ(a[i][j] - c[ic1[i]][j]);
		}
		
	FINISHED:
	/* free memory */
		FREE_VECTOR(ic2);
		FREE_VECTOR(an1);
		FREE_VECTOR(an2);
		FREE_VECTOR(d);
		FREE_VECTOR(live);
		FREE_VECTOR(cd);
	}

} /* kmeans */

/**
 * Optimal-transfer stage of K-means algorithm.
 *
 * This is the optimal transfer stage.  Each point is re-allocated, if
 * necessary, to the cluster that will induce a maximum reduction in the
 * within-cluster sum of squares.  If cluster L was updated in the last
 * quick-transfer stage, it belongs to the live set throughout this stage.
 * Otherwise, at each step, it is not in the live set if it has not been
 * updated in the last $m$ optimal transfer steps.
 *
 * @param a data observations
 * @param m number of observations
 * @param n number of dimensions
 * @param c centers
 * @param k number of clusters
 * @param ic1 closest and
 * @param ic2 next closest center
 * @param nc number in each cluster
 * @param an1 factor for removing 1 from each cluster
 * @param an2 factor for adding 1 to each cluster
 * @param live set to i+1 if i transferred to/from this cluster
 * @param cd indicates if dis to closest center NOT needed
 * @param d decrement in wss upon removing sequence from its current centroid
 * @param indx first 0, number of consecutive observations not transferred
 */
void optra(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, SIZE_T *ic1,
	SIZE_T *ic2, SIZE_T *nc, double *an1, double *an2, SIZE_T *live,
	SIZE_T *cd, double *d, SIZE_T *indx)
{
	SIZE_T i, j, l, ll, i1;
	SIZE_T l1;	/* closest cluster */
	SIZE_T l2;	/* second closest cluster */
	double r2, de, rr, al1, al2, alt, alw;

	/* one loop through all observations */
	for (i = 0, i1=1; i < m; ++i, ++i1) {
		(*indx)++;
		l1 = ic1[i];

		/* no transfer if observation i is only member of its cluster */
		if (nc[l1] != 1) {
			l2 = ic2[i];
			ll = l2;

			/* centroid c[l1] updated: re-compute membership cost;
			 * qtran unsyncs cd and live by computing membership
			 * cost and leaving clusters lives */
			if (i < cd[l1]) { /* && i < live[l1]) { */
				de = 0.;
				for (j = 0; j < n; ++j)
					de += SQ(a[i][j] - c[l1][j]);
				d[i] = de * an1[l1];
			}

			/* find cluster with min. join cost, start with l2;
			 * it may also become next closest center */
			de = 0.;
			for (j = 0; j < n; ++j)
				de += SQ(a[i][j] - c[l2][j]);
			r2 = de * an2[l2];

			for (l = 0; l < k; ++l) {
				/* observations in live set are compared
				 * against all clusters; observations outside
				 * live set are compared only against live
				 * clusters 
				 */
				if ((i1 >= live[l1] && i1 >= live[l]) || l == l1 || l == ll)
					continue;

				rr = r2 / an2[l];
				de = 0.;
				j = 0;
				while (j < n) {
					de += SQ(a[i][j] - c[l][j]);
					if (de >= rr)	/* worse */
						break;
					j++;
				}
				if (j == n) {	/* better exchange */
					r2 = de * an2[l];
					l2 = l;
				}
			}

			/* no transfer improves wss, l2 is new ic2[i] */
			if (r2 >= d[i])
				ic2[i] = l2;

			/* transfer i from l1 to l2, l1 is new ic2[i] */
			else {
				*indx = 0;
				live[l1] = m + i1;
				live[l2] = m + i1;
				cd[l1] = live[l1];	/* now cd and live */
				cd[l2] = live[l2];	/* have same info */

				/* update centers [PARALLEL] */
				al1 = (double) nc[l1];
				alw = al1 - 1.0;
				al2 = (double) nc[l2];
				alt = al2 + 1.0;
				for (j = 0; j < n; ++j) {
					c[l1][j] = (c[l1][j] * al1 - a[i][j]) / alw;
					c[l2][j] = (c[l2][j] * al2 + a[i][j]) / alt;
				}

				/* update cluster counts & assignments */
				nc[l1]--;
				nc[l2]++;
				ic1[i] = l2;
				ic2[i] = l1;

				/* update cost multipliers */
				an2[l1] = alw / al1;
				an1[l1] = DBL_MAX;
				if (alw > 1.0)
					an1[l1] = alw / (alw - 1.0);
				an1[l2] = alt / al2;
				an2[l2] = alt / (alt + 1.0);
			}
		}

		/* no transfers; algorithm has converged, so nothing to do */
		if (*indx == m)
			return;
	}

	/* finished loop, but there was transfer; adjust last obs. transferred
	 * to/from l */
	for (l = 0; l < k; ++l) {
		live[l] -= m;
		cd[l] = 0;	/* qtran will do calculations based on live */
		/* in addition, when qtran finishes all distances to closest
		 * centers are computed; no need for recomputation until first
		 * optimal transfer */
	}

	return;
} /* optra */

/**
 * Quick transfer stage of Hartigan and Wong's K-means algorithm.
 *
 * This is the quick transfer stage. ic1[i] is the cluster observation i
 * belongs to.  ic2[i] is the next closest cluster, the one i is most likely to
 * transfer to. For each observation i, ic1[i] & ic2[i] are switched, if
 * necessary, to reduce within-cluster sum of squares.  The cluster centres are
 * updated after each step. In the optimal transfer stage, live[l] indicates
 * which observation was last transferred in/out of cluster l.  Later
 * observations have compared against the new center; earlier observations have
 * not yet.
 *
 * @param a data observations
 * @param m number of observations
 * @param n number of dimensions
 * @param c centers
 * @param k number of clusters
 * @param ic1 closest and
 * @param ic2 next closest center
 * @param nc number in each cluster
 * @param an1 factor for removing
 * @param an2 factor for adding
 * @param live set to i if i last observation transferred to/from cluster
 * @param d distance to assigned centroid for each point
 * @param indx first 0, number of consecutive observations not transferred
 */
void qtran(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, SIZE_T *ic1,
	SIZE_T *ic2, SIZE_T *nc, double *an1, double *an2, SIZE_T *live, 
	double *d, SIZE_T *indx)
{
	SIZE_T i, j;
	SIZE_T l1, l2;		/* closest and next closest centers */
	SIZE_T istep = 1;	/* current iteration */
	SIZE_T loop_cnt = 1;	/* keep looping as long as transfer */
	double r2, da, al1, al2, alt, alw;
  
	istep = 0;
	while (1) {
		for (i = 0; i < m; i++, istep++, loop_cnt++) {

			/* no transfer during last loop; we're done */
			if (loop_cnt > m) {
				/* there was at least one quick transfer 
				 * last update index needs to be on single loop scale
				 */
				if (!(*indx))
					for (j=0; j<k ;j++)
						/* affected clusters revived */
						if (live[j] > m)
							live[j] = m + 1;
				return;
			}

			l1 = ic1[i];
			l2 = ic2[i];
			/* no transfer if obs. i is sole member of cluster */
			if (nc[l1] != 1) {

				/* l1 changed; recompute membership cost */
				if (istep <= live[l1]) {
					da=0.;
					for (j = 0; j < n; ++j)
						da += SQ(a[i][j] - c[l1][j]);
					d[i] = da * an1[l1];
				}

				/* closest centers changed; transfer possible */
				if (istep < live[l1] || istep < live[l2]) {
					/* compute join cost to next closest */
					r2 = d[i] / an2[l2];
					da = 0.;
					j = 0;
					while (j < n) {
						da += SQ(a[i][j] - c[l2][j]);
						if (da >= r2)
							break;
						j++;
					}

					/* it's closer: transfer! */
					if (j == n) {
						loop_cnt = 0;
						*indx = 0;
						live[l1] = istep + m;
						live[l2] = istep + m;

						/* update centers */
						al1 = (double) nc[l1];
						alw = al1-1.0;
						al2 = (double) nc[l2];
						alt = al2+1.0;
						for (j = 0; j < n; ++j) {
							c[l1][j] = (c[l1][j] * al1 - a[i][j]) / alw;
							c[l2][j] = (c[l2][j] * al2 + a[i][j]) / alt;
						}

						/* update cluster counts & assignments */
						nc[l1]--;
						nc[l2]++;
						ic1[i] = l2;
						ic2[i] = l1;

						/* update cost multipliers */
						an2[l1] = alw / al1;
						an1[l1] = DBL_MAX;
						if (alw > 1.0)
							an1[l1] = alw / (alw - 1.0);
						an1[l2] = alt / al2;
						an2[l2] = alt / (alt + 1.0);
					}
				}
			}
		}
	}
} /* qtran */

/**
 * Returns human-friendly description of a k-means error.
 *
 * @param err integer representation of error; something returned by
 *            #kmeans() or some k-means initialiation routine.
 * @return descriptive string explaining error or NULL if invalid error
 */
const char *kmeans_error(int err) {
	if (err > 0 && err < KMEANS_NUMBER_ERRORS)
		return KMEANS_ERROR_STRING[err];
	else
		return NULL;
} /* kmeans_error */
