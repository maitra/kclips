#include "array.h"
#include "mat_vec.h"
#include "cephes_eigens.h"

void ordered_eigens(double *a, double **evec, double *eval, SIZE_T p);

/*
##
## given LTmatrix that is (near-singular) use the principal components method
## of factor analysis to stabilize it. We use 99% of the variation explained 
## by the PCs.
##
*/

SIZE_T stabilize_fa(double *w, SIZE_T p)
{
	double *e, **ev, sum = 0, sums =0, W[p], alpha = 0.99;
	SIZE_T l, ll, flag=0;

	for (SIZE_T j = 0; j < p ; j++)
		W[j] = w[j*(j+1)/2+j];

	MAKE_VECTOR(e, p);
	MAKE_MATRIX(ev, p, p);
   
	ordered_eigens(w, ev, e, p);
	/* note ordered_eigens changes the w */

	for (l=0; (l<p) && (e[l]>0) ; l++)
		sum+=e[l];

	for(ll = 0; ((ll < l) && (sums < alpha*sum) && (e[ll] > 0)); ll++) 
		sums += e[ll];

	if (ll) {
	
		for (SIZE_T j = 0; j < p; j++) {
			for (SIZE_T m = 0; m <= j; m++) {
				double sum1 = 0;
				for (SIZE_T k = 0; k < ll; k++) 
					sum1 += e[k]*ev[j][k]*ev[m][k];
				w[j*(j+1)/2 + m] = sum1;
			}
		}

		for (l = 0; l < p; l++)
			w[l*(l+1)/2+l] = W[l];
	}
	else
		flag = 1;
	FREE_VECTOR(e);
	FREE_MATRIX(ev);
	return flag;
}
