#ifndef QUANTILE_H
#define QUANTILE_H

#include "constants.h"

#define element_type double

element_type median(element_type *x, SIZE_T n);

element_type select_kth_smallest(element_type *x, SIZE_T n, SIZE_T k);

void HeapOrder(double *array, SIZE_T size, SIZE_T *sorted_index, SIZE_T first_k);

void HeapInSituSort(double *array, SIZE_T size, SIZE_T first_k);

void quantile(SIZE_T n, double *x, double *p, double *q, SIZE_T numqs, SIZE_T type);

double trimmed_mean(SIZE_T n, double *x, double leftprop, double rightprop);

#endif /* QUANTILE_H */
