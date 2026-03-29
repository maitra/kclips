#include <stdio.h>
#include <stdlib.h>
#include "array.h"
#include "mat_vec.h"

/* This function creates the confusion (also known as a contigency table or 
 discrepancy matrix based on two potentially different classification, clid1 
and clid2.  It stores these in the pre-allocated matrix conf. */


void confusion(int n, int k1, int k2, int *clid1, int *clid2, int **conf
,int start){
  /* n: number of points clustered
  k1: number of clusters in first clustering
  k2: number of clusters in second clustering
  clid1: cluster id's for first clustering
  clid2: cluster id's for second clustering
  conf: k1xk2 preallocated matrix
  start: indicates how first cluster is labeled (i.e. starts with a 1 or 0)  */

  int i,j;

  for(i=0;i<k1;i++){
    for(j=0;j<k2;j++){
      conf[i][j]=0;
    }
  }

  for(i=0;i<n;i++){
    conf[(clid1[i]-start)][(clid2[i]-start)]++;
  }
  return;
}

