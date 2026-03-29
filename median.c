/*
  This function calculates the median of an array in an efficient way
  using an adaptation of the QuickSelect algorithm provided in
  "Numerical recipes in C", Second Edition, 
  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
  This is a modification of the code by Nicolas Devillard - 1998 for
  odd n.  The modification handles the case for n even and odd in one shot, 
  and hence is very efficient.

  A nice description of the QuickSelection algorithm is provided
  by Ryan Tibshirani (www.stat.cmu.edu/~ryantibs) in a manuscript written by
  him (on binmedian) in 2008. 

  I also wrote a new function called select_kth_smallest which goes through
  the same steps as median but provides the kth smallest value in the array x. 

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"

#define element_type double

#define SWAP_ELEMENTS(a,b) { register element_type t=(a);(a)=(b);(b)=t; }


element_type median_in_situ(element_type *arr, SIZE_T n) 
{
	/* This function returns the median. It is similar in scope to the 
	   previous function median, except that the array is returned back 
	   semi-sorted with all the entries below the median position smaller
	   than the median and all entries above the median position larger
	   than the median. */
	
	SIZE_T low = 0, high = n - 1;
	SIZE_T rightmed = n >> 1, leftmed = (n - 1) >> 1;
	SIZE_T middle, ll, hh;
	element_type median;

	for (;high > (low + 1);) {
		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) >> 1;
		if (arr[middle] > arr[high])  
			SWAP_ELEMENTS(arr[middle], arr[high]) ;
		if (arr[low] > arr[high])       
			SWAP_ELEMENTS(arr[low], arr[high]) ;
		if (arr[middle] > arr[low])     
			SWAP_ELEMENTS(arr[middle], arr[low]) ;

		/* Swap low item (now in position middle) into position 
		   (low+1) */
		SWAP_ELEMENTS(arr[middle], arr[low+1]) ;

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do 
				ll++; 
			while 
				(arr[low] > arr[ll]) ;
			do 
				hh--; 
			while 
				(arr[hh]  > arr[low]) ;

			if (hh < ll)
				break;

			SWAP_ELEMENTS(arr[ll], arr[hh]) ;
		}

		/* Swap middle item (in position low) back into correct position */
		SWAP_ELEMENTS(arr[low], arr[hh]) ;

		/* Re-set active partition */
		
		if (hh <= leftmed)
			low = ll;
		if (hh >= rightmed)
			high = hh - 1;
	}

	if ((high == low + 1) /* Two elements only */ && (arr[low] > arr[high]))
		SWAP_ELEMENTS(arr[low], arr[high]) ;
	/* else it is all one or two elements only */
	median = (arr[leftmed] + arr[rightmed])/2;

/*	for (SIZE_T i = 0; i < ((leftmed + rightmed) >> 1) - 1; i++)
		if ((arr[i] > median) || (arr[n - i - 1] < median))
		printf("in median_in_situ: %d %d %f %d %f median = %f\n", n, i, arr[i], n-i-1, arr[n-i-1], median );*/
	return median;
}


element_type select_kth_smallest(element_type *x, SIZE_T n, SIZE_T k) 
{
	SIZE_T low = 0, high = n - 1;
	SIZE_T middle, ll, hh;
	element_type *arr, kthval;

	arr = malloc(n * sizeof(element_type));

	memcpy(arr, x, n * sizeof(element_type));

	for (;high > (low + 1);) {
		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) >> 1;
		if (arr[middle] > arr[high])  
			SWAP_ELEMENTS(arr[middle], arr[high]) ;
		if (arr[low] > arr[high])       
			SWAP_ELEMENTS(arr[low], arr[high]) ;
		if (arr[middle] > arr[low])     
			SWAP_ELEMENTS(arr[middle], arr[low]) ;

		/* Swap low item (now in position middle) into position 
		   (low+1) */
		SWAP_ELEMENTS(arr[middle], arr[low+1]) ;

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
			do 
				ll++; 
			while (arr[low] > arr[ll]) ;
			do 
				hh--; 
			while (arr[hh]  > arr[low]) ;

			if (hh < ll)
				break;

			SWAP_ELEMENTS(arr[ll], arr[hh]) ;
		}

		/* Swap middle item (in position low) back into correct position */
		SWAP_ELEMENTS(arr[low], arr[hh]) ;

		/* Re-set active partition */
		
		if (hh <= k)
			low = ll;
		if (hh >= k)
			high = hh - 1;
	}

	if ((high == low + 1) /* Two elements only */ && (arr[low] > arr[high]))
		SWAP_ELEMENTS(arr[low], arr[high]) ;
	/* else it is all one or two elements only */
	kthval = arr[k];
	free(arr);
	return kthval;
}

element_type median(element_type *x, SIZE_T n) 
{
	element_type *arr, median;
	
	arr = malloc(n * sizeof(element_type));

	memcpy(arr, x, n * sizeof(element_type));
	
	median = median_in_situ(arr, n);
	
	free(arr);
	return median;
}
