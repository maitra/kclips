#include "array.h"
#include<stdio.h>
#include<values.h>
#include "quantile.h"

#define SWAP_ELEMENTS(a,b) {register double temp=(a);(a)=(b);(b)=temp;}
#define MIN(a, b) ((a) < (b) ? (a):(b))
/*Declaring heap globally so that we do not need to pass it as an argument every time*/
/* Heap used here is Min Heap */

/*Insert an element SIZE_To the heap */
void InsertInSitu(double *heap, SIZE_T *heapSize, double element)
{/* question if element needs to be a pointer (don't think so) */
        SIZE_T now = ++(*heapSize);
        heap[*heapSize - 1] = element; /*Insert in the last place*/
        /*Adjust its position*/
        while (((now >> 1) >= 1) && (heap[(now >> 1) - 1] > element))   {
		heap[now - 1] = heap[(now >> 1) - 1];
                now = now >> 1;
        }
        heap[now - 1] = element;
}

double DeleteMinInSitu(double *heap, SIZE_T *heapSize)
{
        /* heap[1] is the minimum element. So we remove heap[1]. Size of the 
	   heap is decreased. 
           Now heap[1] has to be filled. We put the last element in its place 
	   and see if it fits.
           If it does not fit, take minimum element among both its children and
	   replace parent with it.
           Again see if the last element fits in that place.*/
/*	This function returns the index of the first element of the heap: this
 is the index of the next smallest element.*/

	double lastElement, minElement = heap[0];
	SIZE_T child, now;

        lastElement = heap[--(*heapSize)];
        /* now refers to the index at which we are now */
        for(now = 1; now*2 <= *heapSize ; now = child)        {
                /* child is the index of the element which is minimum among 
		   both the children */ 
                /* Indexes of children are i*2 and i*2 + 1*/
                child = now*2;
                /*child!=heapSize beacuse heap[heapSize+1] does not exist, 
		  which means it has only one child */
                if(child != (*heapSize) && heap[child] < heap[child - 1] ) 
                        child++;
                /* To check if the last element fits or not it suffices to 
		   check if the last element is less than the minimum element 
		   among both the children*/
                if(lastElement > heap[child - 1])     {      
                        heap[now - 1] = heap[child - 1];
		}
		else /* It fits there */
			break;
	}
	heap[now - 1] = lastElement;
        return minElement;
}

void HeapInSituSort(double *array, SIZE_T size, SIZE_T first_k) 
/* provides the indices of the first k smallest values in the array */
{
	SIZE_T iter, heapSize = 1; 

        for(iter = 1; iter < size; iter++) 
                InsertInSitu(array, &heapSize, array[iter]);
        for(iter = 0; iter < first_k; iter++)
                array[size - iter - 1] = DeleteMinInSitu(array, &heapSize);
	for(iter = 0; iter < MIN(first_k, size >> 1); iter++)
		SWAP_ELEMENTS(array[iter], array[size - iter - 1]);
	
}

