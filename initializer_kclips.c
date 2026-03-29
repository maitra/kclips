#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mat_vec.h"
#include "array.h"
#include "constants.h"
#include "quantile.h"
#define Inf 1e+140

double scatter_cluster_starts(double **a, SIZE_T m, SIZE_T n, double **c, 
			      SIZE_T k, SIZE_T *ic1, SIZE_T *nc, 
			      SIZE_T iter, double alpha, SIZE_T *ikeep, 
			      SIZE_T robust, double constant);

void kmedoid(SIZE_T n, SIZE_T m, double **data, SIZE_T *clustid, SIZE_T *nc, 
	     SIZE_T K, double **medoids);

void hclassify(SIZE_T n, SIZE_T m, double **x, SIZE_T hcrit, SIZE_T nclass, 
	       SIZE_T *class);

void kmeans(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, SIZE_T *ic1,
	    SIZE_T *nc, SIZE_T iter, double *wss,SIZE_T *ifault);
SIZE_T initials(double **x,SIZE_T n,SIZE_T p,SIZE_T nclass,SIZE_T *nc,
		double **Mu,double **LTSigma,SIZE_T *class);
double det_crit(double **data, SIZE_T n, SIZE_T p, SIZE_T K, double **mu,
		double rad, SIZE_T *ikeep, SIZE_T *idns);
/*
double logdetw(SIZE_T nclass, SIZE_T p, double **LTSigma, SIZE_T *nc)
{
	double *temp,detW;
	SIZE_T i,j;
	MAKE_VECTOR(temp,p*(p+1)/2);
	for(j=0;j<(p*(p+1)/2);j++) temp[j]=0;
	for(i=0;i<nclass;i++) {
		for (j=0;j<(p*(p+1)/2);j++) temp[j]+=(nc[i]-1)*LTSigma[i][j];
	}
	i = LTdeterminant(temp,p, &detW); 
	FREE_VECTOR(temp);
	if (i)
		return Inf;
	else 
	return detW;
}
*/
/* write the digits of n in its broken down "base" into buf */
void break_down(SIZE_T n, SIZE_T *base, SIZE_T *buf, SIZE_T buflen)
{
        SIZE_T i;
        for (i=0; i<buflen; i++) {
		buf[i] = n%base[i];
		n /= base[i];
        }
}

SIZE_T assign_closest(double *X, SIZE_T p, SIZE_T nclass, double **Mu)
{
	SIZE_T j, l, class = 0;
	double temp, dum1;

	temp = Inf;
	for (l = 0; l < nclass; l++) {
		dum1 = 0.;
		for (j = 0; j < p; j++)
			dum1 += SQ(X[j] - Mu[l][j]);
		if (dum1 < temp) {
			temp = dum1;
			class = l;
		}
	}
	return class;
}

/*
void eliminulls(double **x,SIZE_T n,SIZE_T p,SIZE_T *nclass,double **Mu, 
		    double **Centers, SIZE_T kk)
{ 
//  This routine eliminates all those centers which have representation less 	  than or equal to kk and returns the remaining centers 

	SIZE_T i,j,*nc,k=0,newtotcl=(*nclass), *ind;

	CMAKE_VECTOR(nc,(*nclass));
	MAKE_VECTOR(ind, n);
	
	for(i=0;i<n;i++) {
		ind[i] = assign_closest(x[i], p, (*nclass), Mu);
		nc[ind[i]]++;
	}
	for(i=0;i<(*nclass);i++) 
		if(nc[i]<=kk) 
			newtotcl--;

	for(i=0;i<(*nclass);i++){
		if (nc[i]>kk) {
			for(j=0;j<p;j++) Centers[k][j]=Mu[i][j];
			k++;
		}
	}
	(*nclass)=newtotcl;
	FREE_VECTOR(ind);
	FREE_VECTOR(nc);
}
*/


void eliminulls(double **x, SIZE_T n, SIZE_T p, SIZE_T *nclass, double ***Mu,
		  SIZE_T kk) 
{	/*This function eliminates all those centers which have   
	  representation less than or equal to kk and returns the remaining 
	  centers in the same pointer.*/

	SIZE_T i, j, k = 0, *ind, *nc, newtotcl = (*nclass);

	MAKE_VECTOR(ind, n);
	CMAKE_VECTOR(nc, (*nclass));	/* calloc version */

	for(i = 0; i < n; i++) {
		ind[i] = assign_closest(x[i], p, (*nclass), *Mu);
		nc[ind[i]]++;
	}
	for(i = 0; i < (*nclass); i++)
		if (!nc[i])
			newtotcl--;
	if (newtotcl < *nclass) { // some classes have no representation
		for (i = 0; i < (*nclass); i++)
			if (nc[i] > kk) {
				for (j = 0; j < p; j++)
					(*Mu)[k][j] = (*Mu)[i][j];
				for (j = 0; j < n; j++) 
					if (ind[j] == i)
						ind[j] = k;
				k++; 
			}
		(*nclass) = newtotcl;
	}

	FREE_VECTOR(ind);
	FREE_VECTOR(nc);
}



void partialfill(double **full, double **partial, SIZE_T nf, SIZE_T np, SIZE_T m){
	SIZE_T i,j;
	for(i=0;i<np;i++){
		for(j=0;j<m;j++){
			full[i][j] = partial[i][j];
		}
	}
	return;
}

SIZE_T starts_in_original_domain(SIZE_T n,SIZE_T m,double **y,SIZE_T nclus,SIZE_T *ningrp,
			      SIZE_T *grpids, double **meds, SIZE_T elim, SIZE_T irob)
{    
	double **cent,*dum1,*probs,*qtl,*cctr,**dumy,**dumcctr,
		*wss1,**mu,**tempcent;
	SIZE_T i,j,k,*counts,*ncl,totcl=1,sumcl=0,*dum2,*buf,*ning,kopt,
		*grclass,*cum1,maxiter=100000,ind=1,complete=0,index=0,partn;   

	MAKE_VECTOR(ncl,m);/*1*/
	i=(int) ceil(pow(n*m,1.0/(m+1)));
	for(j=0;j<m;j++)   ncl[j]=i;  /*Ivan's Change*/
	for(j=0;j<m;j++) {
		totcl*=ncl[j];
		sumcl+=ncl[j];
	}

	SIZE_T box;
	box=n*((int)ceil(sqrt(m)));

	MAKE_VECTOR(cctr,m);/*3*/
	MAKE_VECTOR(dum1,m*(m+1)/2);/*4*/
  
	FREE_VECTOR(dum1);/*4*/
	FREE_VECTOR(cctr);/*3*/

	MAKE_VECTOR(dum1,n);/*5*/
	MAKE_VECTOR(cctr,sumcl);/*6*/
	MAKE_MATRIX(dumy,n,1);/*7*/
	MAKE_VECTOR(dum2,n);/*8*/

	k=0;
	for(i=0;i<m;i++) {
		MAKE_VECTOR(ning,ncl[i]);/*9*/
		MAKE_VECTOR(qtl,ncl[i]);/*10*/
		MAKE_MATRIX(dumcctr,ncl[i],1);/*11*/
		MAKE_VECTOR(wss1,ncl[i]);/*12*/
		for(j=0;j<n;j++) dum1[j]=y[j][i];
		MAKE_VECTOR(probs,ncl[i]);/*13*/
		for(j=0;j<ncl[i];j++) 
			probs[j]=(2*j+1)/(2.*ncl[i]);
		quantile(n, dum1, probs, qtl, ncl[i], 8);
		for(j=0;j<n;j++) dumy[j][0]=y[j][i];
		for(j=0;j<ncl[i];j++) dumcctr[j][0]=qtl[j];
		kmeans(dumy,n,1,dumcctr,ncl[i],dum2,ning,maxiter,wss1,&kopt);
		for(j=0;j<ncl[i];j++) cctr[k++]=dumcctr[j][0];
		FREE_VECTOR(probs);/*13*/
		FREE_VECTOR(wss1);/*12*/
		FREE_VECTOR(ning);/*9*/
		FREE_VECTOR(qtl);/*10*/
		FREE_MATRIX(dumcctr);/*11*/
	}
	FREE_VECTOR(dum2);/*8*/
	FREE_MATRIX(dumy);/*7*/
	FREE_VECTOR(dum1);/*5*/


	MAKE_VECTOR(cum1,m);/*14*/
	cum1[0]=0;
	for(i=1;i<m;i++) cum1[i]=cum1[i-1]+ncl[i-1];

	MAKE_MATRIX(cent,(box),m);/*15*/
	MAKE_VECTOR(buf,m);/*16*/

	while(complete <=totcl){
		for (i=index; i<(box); i++) {
			complete++;
			break_down(i, ncl, buf, m);
			for(j=0;j<m;j++) cent[i][j]=cctr[cum1[j]+buf[j]];
		}
		partn=(box);
		MAKE_MATRIX(tempcent, (box), m);
		eliminulls(y,n,m,&partn,&cent,elim);/*17*/
//		partialfill(ccent,tempcent,(box),partn,m);
		FREE_MATRIX(tempcent);/*18*/
		index=partn;
	}

	/*printf("Grid Finished\n");*/

	FREE_VECTOR(cum1);/*14*/
	FREE_VECTOR(cctr);/*6*/
	FREE_VECTOR(buf);/*16*/
	FREE_VECTOR(ncl);/*1*/

	eliminulls(y,n,m,&partn,&cent,elim);/*18*/
	totcl=partn;

	if(totcl>=nclus) {
		double std,val,optval=Inf,atop=0.05,wtop=6.,**smu;
		SIZE_T *ikp,i11,i12,*medoidid,*nmed;

		ind=0;

		/* Ivan's change: change from kmeans to scatter_cluster 04-17-06 */

		MAKE_VECTOR(ikp,n);/*19*/
		MAKE_VECTOR(counts,totcl);/*20*/
		MAKE_MATRIX(smu,totcl,m);

		cpy(cent,totcl,m,smu);

		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts, 
					   maxiter, 0.01, ikp, irob, 3.);
		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   1, 0.05, ikp, irob, 6.);
		val=det_crit(y,n,m,totcl,cent,std,ikp,grpids);
		if(val<optval) {
			atop=0.01;
			wtop=3.;
			optval=val;
		}
		cpy(cent,totcl,m,smu);



		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   maxiter, 0.01, ikp, irob, 9.);
		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   1, 0.05, ikp, irob, 6.);
		val=det_crit(y,n,m,totcl,cent,std,ikp,grpids);
		if(val<optval) {
			atop=0.01;
			wtop=9.;
			optval=val;
		}
		cpy(cent,totcl,m,smu);

		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   maxiter, 0.05, ikp, irob, 3.);
		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   1, 0.05, ikp, irob, 6.);
		val=det_crit(y,n,m,totcl,cent,std,ikp,grpids);
		if(val<optval) {
			atop=0.05;
			wtop=3.;
			optval=val;
		}
		cpy(cent,totcl,m,smu);



		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   maxiter, 0.05, ikp, irob, 9.);
		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   1, 0.05, ikp, irob, 6.);
		val=det_crit(y,n,m,totcl,cent,std,ikp,grpids);
		if(val<optval) {
			atop=0.05;
			wtop=9.;
			optval=val;
		}
		cpy(cent,totcl,m,smu);


		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   maxiter, 0.1, ikp, irob, 3.);
		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   1, 0.05, ikp, irob, 6.);
		val=det_crit(y,n,m,totcl,cent,std,ikp,grpids);
		if(val<optval) {
			atop=0.1;
			wtop=3.;
			optval=val;
		}
		cpy(cent,totcl,m,smu);



		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   maxiter, 0.1, ikp, irob, 9.);
		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   1, 0.05, ikp, irob, 6.);
		val=det_crit(y,n,m,totcl,cent,std,ikp,grpids);
		if(val<optval) {
			atop=0.1;
			wtop=9.;
			optval=val;
		}
		cpy(cent,totcl,m,smu);


		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   maxiter, 0.2, ikp, irob, 3.);
		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   1, 0.05, ikp, irob, 6.);
		val=det_crit(y,n,m,totcl,cent,std,ikp,grpids);
		if(val<optval) {
			atop=0.2;
			wtop=3.;
			optval=val;
		}
		cpy(cent,totcl,m,smu);



		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   maxiter, 0.2, ikp, irob, 9.);
		std=scatter_cluster_starts(y, n, m, smu, totcl, grpids, counts,
					   1, 0.05, ikp, irob, 6.);
		val=det_crit(y,n,m,totcl,cent,std,ikp,grpids);
		if(val<optval) {
			atop=0.2;
			wtop=9.;
			optval=val;
		}
		cpy(cent,totcl,m,smu);


/*    printf("alpha=%f : w=%f\n",atop,wtop);*/
		std=scatter_cluster_starts(y, n, m, cent, totcl, grpids, counts,
					   maxiter, atop, ikp, irob, wtop);


		FREE_MATRIX(smu);
		FREE_VECTOR(counts);/*20*/
		FREE_VECTOR(ikp);/*19*/

		MAKE_VECTOR(grclass,totcl);/*21*/
		hclassify(totcl,m,cent,2,nclus,grclass);

		MAKE_VECTOR(nmed,nclus);/*22*/
		MAKE_VECTOR(medoidid,n);/*23*/
		for(i11=0;i11<nclus;i11++) nmed[i11]=0;
		for(i11=0;i11<n;i11++){
			medoidid[i11]=grclass[grpids[i11]];
			nmed[medoidid[i11]]++;
		}

		MAKE_MATRIX(mu,nclus,m);/*24*/
		kmedoid(n, m, y, medoidid, nmed, nclus, mu);
		FREE_VECTOR(medoidid);
		FREE_VECTOR(nmed);


		for(i11=0;i11<nclus;i11++){
			for(i12=0;i12<m;i12++){
				meds[i11][i12]=mu[i11][i12];
			}
		}

		FREE_VECTOR(grclass);/*21*/
		FREE_MATRIX(cent);/*18*/
		for(i=0;i<n;i++) grpids[i]=assign_closest(y[i],m,nclus,mu);
		FREE_MATRIX(mu);/*24*/

	}
	else    FREE_MATRIX(cent);
 
	return ind;
}

SIZE_T starts(SIZE_T n,SIZE_T m,double **x,SIZE_T nclus,SIZE_T *ningrp,SIZE_T *grpids,
	   double **meds,SIZE_T elim,SIZE_T iro)
{
	SIZE_T ind,i,j;
	double **y;
	MAKE_MATRIX(y,n,m);
	for(i=0;i<n;i++){
		for(j=0;j<m;j++) y[i][j]=x[i][j];
	}
	ind=starts_in_original_domain(n,m,y,nclus,ningrp,grpids,meds,elim,iro);

	FREE_MATRIX(y);
	return ind;
}
