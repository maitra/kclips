#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mat_vec.h"
#include "array.h"
#include "quantile.h"
#define Inf 1e+140

void adiss(int n, int p, double **X, int *clid, int *nclus, double *avgd)
{
/* This function will create the average dissimilarity vector
 for all points in a cluster for the entire dataset.
*/

  int i1,i2,j;
  double temp;

  for(i1=0;i1<n;i1++){
    avgd[i1]=0.0;
    for(i2=0;i2<n;i2++){
      if(clid[i1]==clid[i2] && i1!=i2){
        temp=0.0;
        for(j=0; j<p;j++){
          temp+= SQ((X[i1][j]-X[i2][j]));
        }
        avgd[i1]+=sqrt(temp);
      }
    }
    avgd[i1]=avgd[i1]/nclus[clid[i1]];   /* ??? might need clid[i1]-1 */
  }
  return;
}


void kmedoid(int n, int m, double **data, int *clustid, int *nc, int K, 
	     double **medoids)
{
/*  This function perform the k-medoids algorithm on a dataset for a given
    number of clusters K.  It assumes some type of initial classification and
    only gives the medoids for each cluster.
*/
  int i,j,k, *minid;
  double *average,tempmin;
  /* First need to find the medoids for each cluster */

  MAKE_VECTOR(average,n);
  adiss(n,m,data,clustid,nc,average);
  MAKE_VECTOR(minid,K);
  for(k=0;k<K;k++){
    minid[k]=0;
    tempmin=Inf;
    for(i=0;i<n;i++){
      if(clustid[i]==k){
        if(tempmin>average[i]){
          tempmin=average[i];
          minid[k]=i;
        }
      }
    }
    for(j=0;j<m;j++) {medoids[k][j]=data[minid[k]][j];}
  }
  FREE_VECTOR(average);
  FREE_VECTOR(minid);

  return;

}



