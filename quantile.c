/**
 * Compute a given set of quantiles for a one-dimensional arrays
 *
 * Ranjan Maitra, Ames, IA 50011
 * 2013/05/24
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "quantile.h"

void quantile(SIZE_T n, double *x, double *p, double *q, SIZE_T numqs, 
	      SIZE_T type)
/* Similar function to the Splus/R quantile function. Returns in q[i]
   the p[i]th quantile of the dataset in the n-dimensional vector x. The 
   quantiles are calculated using different types as per Hyndman and Fan (The  
   American Statistician, 1996). 
   Note that Hyndman and Fan (1996) recommend Type 8, so, by default if the
   types are not numbered 1 through 9, we use a type 8 calculation.
   Also, that R returns type
   7 quantiles by default: to get type 8 quantiles, use type=8 in R.
   For future development, if more type quantiles is desired, just change a,
   b to get the different types for continuous quantiles (for types 4-9),
   for the dicontinuous type (1-3) quantiles, make changes.
   Here, numqs represent the number of    quantiles desired to be returned */
{
	if ((numqs == 1) && (*p == 0.5))
		*q = median(x, n);
	else {
		SIZE_T *srtd_idx;
		
		srtd_idx = malloc(n * sizeof*(srtd_idx));
		
		HeapOrder(x, n, srtd_idx, n);
		
	
		for (SIZE_T jj = 0; jj < numqs; jj++) {
			double eta = p[jj] * n, b = 1/3., a = b, h = 0.5;
			SIZE_T j = eta;

			switch(type) {
			case 1:	
				if (eta > j)
					h = 1;
				else 
					h = 0;
				break;
			case 2: 
				if (eta > j)
						h = 1;
				else
					h = 0.5;
				break;
			case 3:
				if ((eta == j) && (!(j % 2)))
					h = 0;
				else
					h = 1;
				break;
			default: {
				switch (type) {
				case 4:
					a = 0;
					b = 1;
					break;
				case 5:
					a = 0.5;
					b = 0.5;
					break;
				case 6:
					a = 0;
					b = 0;
					break;
				case 7:
					a = 1;
					b = 1;
					break;
				case 9:
					a = 0.375;
					b = 0.375;
					break;
				}
				eta = a + p[jj] * (n + 1 - a - b);
				j = eta;
				h = eta - j;
			}
			}
			if (j < 1) 
				q[jj] = x[srtd_idx[0]];
			else 
				if (j >= n)
					q[jj] = x[srtd_idx[n-1]];
				else 
					q[jj] = (1 - h) * x[srtd_idx[j - 1]] + 	h * x[srtd_idx[j]];
		}
		free(srtd_idx);
	}
}

double trimmed_mean(SIZE_T n, double *x, double leftprop, double rightprop)
/* This function calculates the trimmed mean of x after dropping the lowest 
   "left" proportion of observations and the upper "right" proportion of 
   observations. The sum of "left" and "right, each assumed to be between 0 
   and 1, is further assumed to be between 0 and 1 and the code does not check
   for consistency.*/
{
	double sum = 0.0;
	
	if (!((SIZE_T)(leftprop * n + 0.5)) ||  (!((SIZE_T)(rightprop * n + 0.5)))) { 
		/* no need for truncating because all are included */
		for (SIZE_T i = 0; i < n; i++) 
			sum += x[i];
		sum /= n;
	}
	
	else {
		SIZE_T rightk, *index, nx = 0;
		
		index = malloc(n * sizeof index);
		
		rightk = rightprop*n;
		
		HeapOrder(x, n, index, rightk); /* don't need beyond the first 
						   rightk smallest elements */
		
		for (SIZE_T i = leftprop * n; i < rightk; i++) {
			sum += x[index[i]];
			nx++;
		}
		free(index);
		sum /= nx;
	}
	return sum;
}

