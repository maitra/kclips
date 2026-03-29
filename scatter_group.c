/*
  This routine clusters a dataset into groups using the scatter_cluster 
  algorithm with 
  start points provided by the starts_via_svd routine. The function returns 
  the classification ids as well as the estimated centers.
*/
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "array.h"
#define PI 3.141593
#define Inf 1e+140
#include "scatter_group.h"
#include "rmath-stalone/rfns.h"

double scatter_group(double **x,int n,int p,int nclass,double **Mu,int *nc,
		     int *class,double alpha, int *ikeep, int robust, 
		     int iters, double **medoids, int maxelim, 
		     int minelim, double w4sb[], double *finalW, int nw, int traceW)
/*
  w4sb: the weights used in the robust bi-weight estimator of standard deviation
  (w in page 342, last line of Maitra and Ramler (2009).
  nw: the number of candidate w's in w4sb (5 by default: {3.,4.5,6.,7.5,9.}
  traceW: if trace rather than determinant of W is used.
*/
{
	int j, flag=0;
	double **LTSigma, ret=3., radk;

	MAKE_MATRIX(LTSigma,nclass,p*(p+1)/2);
	if (nclass == 1) {
		nc[0] = n;
		for (j = 0; j < n;j++) 
			class[j]=0;
    
		meandispersion(x, n, p, Mu[0], LTSigma[0]); 
/*		initials(x, n, p, nclass, nc, Mu, LTSigma, class);
	
		vassign_closest(x, n, p, Mu, nclass, class, nc); */

/*		print_dvector(Mu[0], p, "%f ");
		print_dvector(LTSigma[0], p, "%f "); */

		radk = scatter_cluster(x, n, p, Mu, nclass, class, nc,
				       iters, alpha, ikeep, robust, 
				       w4sb, finalW, nw);

//		printf("radk %f\n", radk);

		/* end of loop for given K  when nclass = 1*/
		/* begin case estimated K, but nclass = 1*/
		if (!traceW) {
			/* First reduce data to get W matrix. */
			int nstar=0, *rnc, tint=0, *id,i;
			double **tempdata, *w, detW, volume;
			int *rclass;
				
			MAKE_VECTOR(rnc,nclass);/*1*/
			for (i = 0; i < nclass; i++)
				rnc[i]=0;
			MAKE_VECTOR(id,n);
			for (i = 0; i < n; i++){
				nstar += ikeep[i];
				tint = (rnc[class[i]]+1*ikeep[i]);
				rnc[class[i]] = tint;
				id[i] = ikeep[i] * (class[i]+1);
			} 

			MAKE_VECTOR(rclass,nstar);
			tempdata=reducedata(x, n, p, nclass, ikeep, 
					    class, rnc, nstar, rclass);
			initials(tempdata, nstar, p, nclass, rnc, Mu, 
				 LTSigma,rclass);

//			print_dmatrix(Mu, nclass, p, "%f ");

//			print_dmatrix(LTSigma, nclass, p*(p+1)/2, "%f ");

			FREE_MATRIX(tempdata);
			FREE_VECTOR(rclass); 
			MAKE_VECTOR(w,p*(p+1)/2);
			for(j=0;j<p*(p+1)/2;j++) 
				w[j]=0.0;
				
			for(i=0;i<nclass;i++) {
				for (j=0;j<(p*(p+1)/2);j++) 
					w[j]+=(rnc[i]-1)*LTSigma[i][j];
			}
			
//			print_dvector(w, p*(p+1)/2, "%f ");
			
			i = logLTdeterminant(w, p, &detW);

//			printf("%d \n", i);

			if (i) 
				detW = Inf;
			FREE_VECTOR(rnc);
			FREE_VECTOR(w);
			FREE_VECTOR(id);
			/* Now determine the remaining area. */
				
			volume = get_volume(n, p, x, radk, nstar, nclass, Mu);
				
			if(volume <= 0.0) 
				volume = 1.0;
//			printf("nstar/n = %f\n",(1.*nstar)/(1.*n));
			ret = (1.0*nstar) / (1.0*n) * (0.5*detW - log(volume)) + log(volume) + log(nclass + 1./p);

//			printf("%f\n", ret);
		} /* end of case for nclass = 1 and detW calculations */
		else {
			/* This uses the trace (i.e. objective function)
			   to determine the numberof clusters. */
			/* First reduce data to get W matrix. */
			int nstar=0, i;
			for (i = 0; i < n; i++)   
				nstar += ikeep[i];
				
			double detw;
			detw = trace_crit(x, n, p,nclass,Mu,radk, ikeep,class,  
					  SQ(radk)/chisq_quantile(alpha,(double)(n),0,0));
				
			ret=detw;
		} /*end of case for nclass = 1 and traceW calculations */
	} /*end of case for nclass = 1*/
	else {  /* start with general case of classes > 1,  needs initialization */
		int temp_elim, penalty;
		double radk;
/*		double dum = p*nclass;
		fprintf(stderr, "p = %d\n", p);
		fprintf(stderr, "dum = %f\n", dum);*/
		
		flag=1;
		
		if (maxelim < 0) 
/*			fprintf(stderr, "here %d\n", n/(p*nclass) + ((n%(p*nclass))?1:0));
			fprintf(stderr, "here %d\n", CEIL(n/(1.0*p*nclass))); */
			temp_elim = CEIL(n/(1.*p*nclass)); 
		else 
			temp_elim=maxelim;
		
		penalty = 1;
		
		while(flag && (temp_elim>=minelim)){
			j = starts(n,p,x,nclass,nc,class,medoids,temp_elim,robust);
			if (!j) {
				vassign_closest(x, n, p, medoids, nclass, class, nc);
				flag=0;
			}
			else{
/*				subt = MAX(1, temp_elim-CEIL(exp(log((double)CEIL(n/(1.*p*nclass)))-0.15*penalty)));
				temp_elim -= subt; */
				
				temp_elim -= CEIL(exp(log((double)CEIL(n/(1.*p*nclass)))-0.15*penalty));
				temp_elim = MAX(temp_elim, 1);

//				printf("temp_elim = %d\n", temp_elim);
				
				penalty++;
				/*printf("Note: singular starting point - Reducing value of eliminulls to %d\n",temp_elim);*/
			}
		}

		if (flag) {
			printf("Note: singular starting point -- Cannot reduce eliminull further -- Program Aborted\n");
			ret=Inf;
		}
		else {
			int i11, i12;
			for(i11=0;i11<nclass;i11++){
				for(i12=0;i12<p;i12++)
					Mu[i11][i12]=medoids[i11][i12];
			}
      
			vassign_closest(x, n, p, Mu, nclass, class, nc);
			/*      for(i11=0;i11<nclass;i11++) printf("%d\n",nc[i11]);*/

			radk=scatter_cluster(x, n, p, Mu, nclass, class,nc, iters, alpha,  
					     ikeep, robust, w4sb,finalW,nw);
      
			/*estimating the objective function for no. of clusters*/
			if (traceW) {
				/* This uses the trace (i.e. objective function) to determine the numberof clusters. */
				int nstar=0, i;
				double detw;
				initials(x,n,p,nclass,nc,Mu,LTSigma,class);
				radk=scatter_cluster(x, n, p, Mu, nclass, class,nc,iters, alpha, 
						     ikeep, robust, w4sb,finalW,nw);
				/* First reduce data to get W matrix. */
				for(i=0;i<n;i++)  
					nstar+=ikeep[i];
				detw = trace_crit(x, n, p, nclass, Mu, radk, ikeep,class,  
						  radk*radk/chisq_quantile(alpha,(double)(n),0,0));
					
				/* Now determine the remaining area. */
				ret=detw;
			}
			else {
				int nstar=0,*id,i;
				/* First get W matrix  */
				double **fmu,**fLTSigma, *w, detW, volume;
				int *fnc;

				MAKE_VECTOR(fnc,nclass+1);
				for(i=0;i<nclass+1;i++) 
					fnc[i]=0;
				MAKE_VECTOR(id,n);
				for(i=0;i<n;i++){
					nstar+=ikeep[i];
					id[i]=ikeep[i]*(class[i]+1);
					fnc[id[i]]++;
				}
				MAKE_MATRIX(fmu,nclass+1,p);
				MAKE_MATRIX(fLTSigma,nclass+1,p*(p+1)/2);
				initials(x,nstar,p,nclass+1,fnc,fmu,fLTSigma,id);
				MAKE_VECTOR(w,p*(p+1)/2);
				for(j=0;j<p*(p+1)/2;j++) 
					w[j]=0.0;

				for(i=1;i<(nclass+1);i++) {
					for (j=0;j<(p*(p+1)/2);j++) 
						w[j]+=(fnc[i]-1)*fLTSigma[i][j];
				}

				i = logLTdeterminant(w,p, &detW);
					
				if (i) 
					detW = Inf;
				FREE_VECTOR(w);
				FREE_MATRIX(fmu);
				FREE_MATRIX(fLTSigma);
				FREE_VECTOR(fnc);
				/* Now determine the remaining area. */
				FREE_VECTOR(id);
				volume=get_volume(n, p, x, radk, nstar, nclass, Mu);
				if (volume <= 0.0) 
					volume = 1.0;
//				fprintf(stderr, "volume = %f nstar/n = %f\n", volume, (1.*nstar)/(1.*n));
				ret=nstar/(1.*n)*(detW/2 - log(volume)) + log(volume) + log(nclass + 1./p);
			}
		}
	}
	FREE_MATRIX(LTSigma);
//	printf("objective function value = %f\n", ret);
	return ret;
}


