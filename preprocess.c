#include <stdlib.h>
#include <stdio.h>
#include "quantile.h"
#include "random.h"

/* This function calculates the minimum distance between any two values in a 
   coordinate. If the smallest distance between any two subsequent numbers is
   zero, we add a small random number which is uniform(0, 10^-3) times that of 
   the smallest positive difference between any two values in the coordinate.
*/

double preprocess_coordinate(double *x, int n)
{
	double minval, minposval, duf;

	HeapInSituSort(x, n, n);
	
	minval = minposval = x[1] - x[0];
	
	for (int i = 1; i < n - 1; i++) {
		duf = x[i+1] - x[i];
		if (duf < minval) {
			if (minval != 0)
				minval = duf;
		}
		if (duf > 0 ) {
			if (minposval == 0) 
				minposval = duf;
			else {
				if (duf < minposval)
					minposval = duf;
			}
		}
	}
	if ((minposval > 0) && (minval == minposval)) 
		return 0;
	else 
		return minposval;
}
	
void preprocess_data(double **data, int n, int m)
{
	int i, j;
	double x[n];

	for (j = 0; j < m; j++) {
		double minposval;
		for (i = 0; i < n; i++)
			x[i] = data[i][j];
		minposval = preprocess_coordinate(x, n);
		if (minposval > 0) {
			minposval *= 0.01;
			for (i = 0; i < n; i++)
				data[i][j] += minposval * runi();
		}
	}
	return;	
}
