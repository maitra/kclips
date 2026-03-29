/* These set of functions are used to obtain a robust estimate of the scale 
   parameter based on a bi-weight estimator (See Understanding Robust and 
   Exploratory Data Analysis by Hoaglin, Mosteller & Tukey for more details.)*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "array.h"
#include "quantile.h"

#define SQ(x) ((x) * (x))
#define CUBE(x) ((x) * SQ(x))

double MAD(int n, double *x, int mcent){
/* This function calculates the median absolute deviation from the median 
   of a vector x. 
   If mcent = 1 then it is assumed that the vector x has already been 
   centered around it's median
*/
	int i;
	double *absdev, mad;
	MAKE_VECTOR(absdev,n);
  
	if(mcent==0){
		double sampmed = median(x, n);
		for(i=0;i<n;i++){
			absdev[i]=fabs((x[i]-sampmed));
		}
	}
	else{
		for(i=0;i<n;i++){
			absdev[i]=fabs(x[i]);
		}    
	}
	mad=median(absdev, n);
	FREE_VECTOR(absdev);
	return mad;
}


double biweight(double x){
/* Calculates the biweight function. */
	double bi=0.0;
	if(fabs(x)<1.0) 
		bi = x * SQ(1 - SQ(x*x));
	return bi;
}

double sbi(int n, double *x, double c, int mcent){
/* This function calculates the biweight estimate of the scale (i.e. standard 
   deviation) for a vector x of length n.  Here, c represents the constant 
   used in determining the weight of the observations.  c = 9 or 6 is typical 
   as it gives zero weight to observations 6 (4) standard deviations away 
   from the median.
   mcent controls where the observations are median-centered or no.
*/

	int i;
	double sbi=0.0, M=0.0,mad=0.0, *u,num=0.0,dem=0.0;

  
	if (!mcent) 
		M = median(x, n);
  	mad = MAD(n, x, mcent);
	MAKE_VECTOR(u,n);
//	printf("MAD = %f\n", mad);

//	printf("c = %f\n", c);

	for(i=0;i<n;i++)
		u[i] = (x[i] - M) / (c * mad);
	  
	for(i=0;i<n;i++)
		if(SQ(u[i])<1.0){
			num+= SQ(x[i] - M) * SQ(SQ(1 - SQ(u[i])));
			dem+=(1 - SQ(u[i])) * (1 - 5 * SQ(u[i]));
		}
	
	FREE_VECTOR(u);
	num=sqrt((double)n*num);
	dem=fabs(dem);
	sbi=num/dem;
	return sbi;
}

