#include<stdio.h>
#include<stdlib.h>
#include<values.h>
#include "quantile.h"

/*Insert an element SIZE_To the heap */
void Insert(double *heap, SIZE_T *heapSize, double element, SIZE_T *index, SIZE_T iter)
{
        SIZE_T now = ++(*heapSize);
        heap[*heapSize - 1] = element; /*Insert in the last place*/
	index[*heapSize - 1] = iter;
        /*Adjust its position*/
        while (((now >> 1) >= 1) && (heap[(now >> 1) - 1] > element))   {
                heap[now - 1] = heap[(now >> 1) - 1];
		index[now - 1] = index[(now >> 1) - 1];
                now = now >> 1;
        }
        heap[now - 1] = element;
	index[now - 1] = iter;
}

SIZE_T DeleteMin(double *heap, SIZE_T *index, SIZE_T *heapSize)
{
        /* heap[1] is the minimum element. So we remove heap[1]. Size of the 
	   heap is decreased. 
           Now heap[1] has to be filled. We put the last element in its place 
	   and see if it fits.
           If it does not fit, take minimum element among both its children and
	   replace parent with it.
           Again see if the last element fits in that place.
	   This function returns the index of the first element of the heap: 
	   this is the index of the next smallest element.*/

	double lastElement;
	SIZE_T child, now, lastindex, index_min = index[0];

        lastElement = heap[(*heapSize) - 1];
	lastindex = index[(*heapSize) - 1]; /* */
	(*heapSize)--;
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
			index[now - 1] = index[child - 1]; /* */
		}
		else /* It fits there */
			break;
	}
	heap[now - 1] = lastElement;
	index[now - 1] = lastindex; /* */
        return index_min;
}

void HeapOrder(double *array, SIZE_T size, SIZE_T *sorted_index, SIZE_T first_k) 
/* provides the indices of the first k smallest values in the array */
{
	double *heap;
	SIZE_T *index, iter = 0, heapSize = 1; 

	heap = malloc(size * sizeof *(heap));
	index = malloc(size * sizeof *(index));

	heap[0] = array[iter];/*	initialize heap*/
	index[0] = iter;

        for(iter = 1; iter < size; iter++) 
                Insert(heap, &heapSize, array[iter], index, iter);

        for(iter = 0; iter < first_k; iter++)
                sorted_index[iter] = DeleteMin(heap, index, &heapSize);
	
	free(index);
	free(heap);
}

