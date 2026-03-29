#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "array.h"

void confusion(int n, int k1, int k2, int *clid1, int *clid2, int **conf
	       ,int start);

/* This functions calculates both the Rand and Adjusted Rand statistic
for comparing two clusterings on the same dataset.
It returns the value fo the adjusted Rand and gives the Rand for free.
 */


double arand(int N, int k1, int k2, int *clid1, int *clid2,int start, double *rand)
{
		 
  int i,j,**n,sumk1sq, *sumk1, sumk2sq, *sumk2,sumsq,sumk1k2sq;
  double discordant, Rand, adjrand,Eindex;


  MAKE_MATRIX(n,k1,k2);

  confusion(N,k1,k2,clid1,clid2,n,start);

  MAKE_VECTOR(sumk1,k1);
  sumk1sq=0;
  for(i=0;i<k1;i++){
    sumk1[i]=0;
    for(j=0;j<k2;j++){
      sumk1[i]+=n[i][j];
    }
    sumk1sq+=sumk1[i]*sumk1[i];
  }

  MAKE_VECTOR(sumk2,k2);
  sumk2sq=0;
  for(j=0;j<k2;j++){
    sumk2[j]=0;
    for(i=0;i<k1;i++){
      sumk2[j]+=n[i][j];
    }
    sumk2sq+=sumk2[j]*sumk2[j];
  }


  sumsq=0;
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      sumsq+=n[i][j]*n[i][j];
    }
  }


  FREE_MATRIX(n);


  sumk1k2sq=0;
  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      sumk1k2sq+=sumk1[i]*sumk1[i]*sumk2[j]*sumk2[j];
    }
  }

  FREE_VECTOR(sumk1);
  FREE_VECTOR(sumk2);
  
  Eindex=sumk1k2sq/(1.0*N*(N-1)+N*1.0*N/(1.0*N-1.0)) -1.0*(sumk2sq+sumk1sq)/(1.0*N-1.0);
  Eindex*=2.;
  Eindex/=1.0*N*(1.0*N-1.0);

  discordant=0.5*(sumk1sq+sumk2sq)-1.0*sumsq;

  Rand=1.0-discordant/(1.0*N*(N-1)/2.0);

  (*rand)=Rand;

  adjrand=(Rand-Eindex)/(1.0-Eindex);

  return adjrand;
}
