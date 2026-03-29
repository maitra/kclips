/*------------------------------------------------------------------------------------*/
/*                                                                                    */
/*    Copyright (C) 2008, 2009 and 2010, Georges Quénot, LIG-CNRS.                    */
/*    Version 1.02 Last revision: March 31, 2010.                                     */
/*                                                                                    */
/*    This file is part of KNNLSB (K Nearest Neighbors Linear Scan Baseline).         */
/*                                                                                    */
/*    KNNLSB is free software: you can redistribute it and/or modify                  */
/*    it under the terms of the GNU Lesser General Public License as                  */
/*    published by the Free Software Foundation, either version 3 of                  */
/*    the License, or (at your option) any later version.                             */
/*                                                                                    */
/*    KNNLSB is distributed in the hope that it will be useful,                       */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of                  */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   */
/*    GNU Lesser General Public License for more details.                             */
/*                                                                                    */
/*    You should have received a copy of the GNU Lesser General Public                */
/*    License along with KNNLSB.  If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                                    */
/*------------------------------------------------------------------------------------*/

#include "knnp.h"

/*------------------------------------------------------------------------------------*/

#define BS 4096
#define RB_MAX 1024

/*------------------------------------------------------------------------------------*/

SIZE_T UsedMem = 0;
SIZE_T MaxUsedMem = 0;

SIZE_T UsedMemory()
{
	return(UsedMem);
}

SIZE_T MaxUsedMemory()
{
	return(MaxUsedMem);
}

void *malloc1a(SIZE_T n0, SIZE_T size)
{
	void *v1;
	if (n0 == 0) {
		fprintf(stderr,"Error in malloc1a: n0 = 0\n");
		return(NULL);
	}
	if (!(v1 = (void *) malloc(n0*size))) {
		fprintf(stderr,"Error in malloc1a: malloc() failed\n");
		return(NULL);
	}
	UsedMem += n0*size;
	if (UsedMem > MaxUsedMem) MaxUsedMem = UsedMem;
	return(v1);
}

void free1a(void *v1, SIZE_T n0, SIZE_T size)
{
	free(v1);
	UsedMem -= n0*size;
}

void *malloc2a(SIZE_T n0, SIZE_T n1, SIZE_T size)
{
	SIZE_T i0;
	char *v1;
	void **v2;
	if (n0 == 0) {
		fprintf(stderr,"Error in malloc2a: n0 = 0\n");
		return(NULL);
	}
	if (n1 == 0) {
		fprintf(stderr,"Error in malloc2a: n1 = 0\n");
		return(NULL);
	}
	if (!(v2 = (void **) malloc(n0*sizeof(void *)))
	    || !(v1 = (char *) malloc(n0*n1*size))) {
		fprintf(stderr,"Error in malloc2a: malloc() failed\n");
		return(NULL);
	}
	for (i0 = 0; i0 < n0; i0++) v2[i0] = v1+i0*n1*size;
	UsedMem += n0*sizeof(void *)+n0*n1*size;
	if (UsedMem > MaxUsedMem) MaxUsedMem = UsedMem;
	return(v2);
}

void free2a(void *v2, SIZE_T n0, SIZE_T n1, SIZE_T size)
{
	free(*((char **) v2));
	free(v2);
	UsedMem -= n0*sizeof(void *)+n0*n1*size;
}

double get_time()
{
	struct timeval t;
	gettimeofday(&t,(struct timezone *) NULL);
	return(t.tv_sec+((double) 0.000001)*t.tv_usec);
}

/*------------------------------------------------------------------------------------*/

SIZE_T get_nb_cores()
{
	SIZE_T n = 0;
	char line[BS];
	FILE *fp;
	if ((fp = fopen("/proc/cpuinfo","r")) == NULL) return (1);
	while (fgets(line,BS,fp)) n += !strncmp(line,"processor",9);
	fclose(fp); /* brought in by RM */
	return(n);
}

SIZE_T get_cache_size()
{
	SIZE_T n = 512;
	char line[BS],c;
	FILE *fp;
	if ((fp = fopen("/proc/cpuinfo","r")) == NULL) return (1);
	while (fgets(line,BS,fp))
		if (!strncmp(line,"cache size",10)) sscanf(line+11,"%c %d",&c,&n);
	fclose(fp); /* brought in by RM */
	return(n*1024);
}

/*------------------------------------------------------------------------------------*/

fdata *getnorm(fdata **m, SIZE_T n, SIZE_T nc, int dc)
{
	SIZE_T i,ic;
	double si;
	fdata *s,*v;
	if ((s = malloc1a(n,sizeof(fdata))) == NULL) return(NULL);
	for (v = *m, i = 0; i < n; s[i++] = si, v  += nc) {
		for (si = 0, ic = 0; ic < nc; ic++) si += v[ic]*v[ic];
		if (dc) si = si ? 1/sqrt(si) : 1;
	}
	return(s);
}

int keyCompare(const void *key1, const void *key2)
{
	if (((KEY *) key1)->d < ((KEY *) key2)->d) return(-1);
	if (((KEY *) key1)->d > ((KEY *) key2)->d) return(+1);
	return(0);
}

/*------------------------------------------------------------------------------------*/

void *knns(void *par)
{
	SIZE_T ic,ik,ia,ib,ik0,ik1;
	int cont;
	double d,sai,dv,vac,vbc;
	fdata *va,*vb;
	KNNPARAMS *kpar = (KNNPARAMS *) par;
	KEY **key = kpar->key;
	fdata **ma = kpar->ma;
	fdata **mb = kpar->mb;
	fdata *sa = kpar->sa;
	fdata *sb = kpar->sb;
	fdata *dm = kpar->dm;
	SIZE_T la = kpar->la;
	SIZE_T ha = kpar->ha;
	SIZE_T lb = kpar->lb;
	SIZE_T hb = kpar->hb;
	SIZE_T nc = kpar->nc;
	SIZE_T nk = kpar->nk;
	int dc = kpar->dc;
	for (ib = lb; ib < hb; ib++) {
		dm[ib] = FDATA_MAX;
		for (ik = 0; ik < nk; ik++) {
			key[ib][ik].k = -1;
			key[ib][ik].d = FDATA_MAX;
		}
	}
	for (ia = la; ia < ha; ia++) {
		va = ma[ia];
		if (dc == 3) {
			for (ib = lb; ib < hb; ib++) {
				vb = mb[ib];
				for (d = 0, ic = 0; ic < nc; ic++) {
					dv = vb[ic]-va[ic];
					d += dv*dv;
				}
				if (d < dm[ib]) {
					ik0 = 0; ik1 = nk-1; ik = (ik0+ik1)/2; cont = 1;
					while (cont) {
						if (d > key[ib][ik].d) ik0 = ik;
						else ik1 = ik;
						cont = (ik0+1 < ik1);
						ik = (ik0+ik1)/2;
					}
					if (d < key[ib][ik0].d) ik = ik0; else ik = ik1;
					if (ik < nk-1) memmove(key[ib]+ik+1,key[ib]+ik,(nk-1-ik)*sizeof(KEY));	
					key[ib][ik].k = ia;
					key[ib][ik].d = d;
					dm[ib] = key[ib][nk-1].d;
				}
			}
		} else if (dc == 2) {
			for (ib = lb; ib < hb; ib++) {
				vb = mb[ib];
				for (d = 0, ic = 0; ic < nc; ic++) {
					dv = (vbc = vb[ic])-(vac = va[ic]);
					d += dv*dv/(1.0e-30+vbc+vac);
				}
				if (d < dm[ib]) {
					ik0 = 0; ik1 = nk-1; ik = (ik0+ik1)/2; cont = 1;
					while (cont) {
						if (d > key[ib][ik].d) ik0 = ik;
						else ik1 = ik;
						cont = (ik0+1 < ik1);
						ik = (ik0+ik1)/2;
					}
					if (d < key[ib][ik0].d) ik = ik0; else ik = ik1;
					if (ik < nk-1) memmove(key[ib]+ik+1,key[ib]+ik,(nk-1-ik)*sizeof(KEY));	
					key[ib][ik].k = ia;
					key[ib][ik].d = d;
					dm[ib] = key[ib][nk-1].d;
				}
			}
		} else {
			if (sa) {
				sai = sa[ia];
			} else {
				for (sai = 0, ic = 0; ic < nc; ic++) sai += va[ic]*va[ic];
				if (dc) sai = sai ? 1/sqrt(sai) : 1;
			}
			for (ib = lb; ib < hb; ib++) {
				vb = mb[ib];
				for (d = 0, ic = 0; ic < nc; ic++) d += vb[ic]*va[ic];
				if (dc) d = -d*sb[ib]*sai; else d = sb[ib]+sai-2*d;
				if (d < dm[ib]) {
					ik0 = 0; ik1 = nk-1; ik = (ik0+ik1)/2; cont = 1;
					while (cont) {
						if (d > key[ib][ik].d) ik0 = ik;
						else ik1 = ik;
						cont = (ik0+1 < ik1);
						ik = (ik0+ik1)/2;
					}
					if (d < key[ib][ik0].d) ik = ik0; else ik = ik1;
					if (ik < nk-1) memmove(key[ib]+ik+1,key[ib]+ik,(nk-1-ik)*sizeof(KEY));	
					key[ib][ik].k = ia;
					key[ib][ik].d = d;
					dm[ib] = key[ib][nk-1].d;
				}
			}
		}
	}
	return(NULL);
}

KEY **knnps(fdata **ma, fdata **mb, fdata *sa, SIZE_T na, SIZE_T nb, SIZE_T nc,
	    SIZE_T nk, SIZE_T np, SIZE_T cs, int dc, int sta)
{
	SIZE_T ib,ik,ip,lb,hb,rb,iq,nq,*km,kr;
	double t0 = 0,t;
	fdata *sb = NULL,*dm,dr;
	pthread_t *tid;
	KEY **key,**keyp;
	KNNPARAMS *par;
	if (na < nk) return(NULL);
	if (sta) t0 = get_time();
	if (np == 0) np = get_nb_cores();
	if (cs == 0) cs = get_cache_size();
	if (np > na/nk) np = na/nk;
	if ((tid = malloc1a(np,sizeof(pthread_t))) == NULL) return(NULL);
	if ((par = malloc1a(np,sizeof(KNNPARAMS))) == NULL) return(NULL);
	if ((km = malloc1a(np,sizeof(SIZE_T))) == NULL) return(NULL);
	if ((dc < 2) && ((sb = getnorm(mb,nb,nc,dc)) == NULL)) return(NULL);
	if ((key = malloc2a(nb,nk,sizeof(KEY))) == NULL) return(NULL);
	if ((rb = cs/(2*(nc*sizeof(fdata)+np*nk*sizeof(KEY)))) == 0) rb = 1;
	if (rb > nb) rb = nb;
	if (rb > RB_MAX) rb = RB_MAX;
	nq = nb/rb;
	nq += (nq*rb < nb);
	if (sta > 1) printf("%d blocks of %d vectors\n",nq,rb);
	if ((keyp = malloc2a(np*rb,nk,sizeof(KEY))) == NULL) return(NULL);
	if ((dm = malloc1a(np*rb,sizeof(fdata))) == NULL);
	for (ip = 0; ip < np; ip++) {
		par[ip].key = keyp+ip*rb;
		par[ip].ma = ma;
		par[ip].mb = mb;
		par[ip].sa = sa;
		par[ip].sb = sb;
		par[ip].dm = dm+ip*rb;
		par[ip].la = (ip*na)/np;
		par[ip].ha = (ip+1 == np) ? na : ((ip+1)*na)/np;
		par[ip].nc = nc;
		par[ip].nk = nk;
		par[ip].np = np;
		par[ip].dc = dc;
	}
	for (iq = 0; iq < nq; iq++) {
		lb  = iq*rb;
		hb  = (iq+1 == nq) ? nb : (iq+1)*rb;    
		for (ip = 0; ip < np; ip++) {
			par[ip].lb = lb;
			par[ip].hb = hb;
		}
		for (ip = 1; ip < np; ip++) if (pthread_create(tid+ip,NULL,knns,par+ip)) {
				printf("error creating thread.");
				abort();
			}
		knns(par);
		for (ip = 1; ip < np; ip++) if (pthread_join(tid[ip],NULL)) {
				printf("error joining thread.");
				abort();
			}
		for (ib = lb; ib < hb; ib++) {
			for (ip = 0; ip < np; ip++) km[ip] = 0;
			for (ik = 0; ik < nk; ik++) {
				kr = 0;
				dr = keyp[ib-lb][km[0]].d;
				for (ip = 1; ip < np; ip++) if (dr > keyp[ip*rb+ib-lb][km[ip]].d) {
						kr = ip;
						dr = keyp[ip*rb+ib-lb][km[ip]].d;
					}
				key[ib][ik].k = keyp[kr*rb+ib-lb][(km[kr]++)].k;
				key[ib][ik].d = dr;
			}
		}
		for (ip = 0; ip < np; ip++) {
			par[ip].key -= rb;
			par[ip].dm -= rb;
		}
	}
	for (ib = 0; ib < nb; ib++) {
		for (ik = 0; ik < nk; ik++) if (dc == 1) {
				if (key[ib][ik].d < -1) key[ib][ik].d = -1;
				if (key[ib][ik].d > +1) key[ib][ik].d = +1;
				key[ib][ik].d = acos(-key[ib][ik].d);
			} else {
				if (key[ib][ik].d < 0) key[ib][ik].d = 0;
				key[ib][ik].d = sqrt(key[ib][ik].d);
			}
	}
	free1a(tid,np,sizeof(pthread_t));
	free1a(par,np,sizeof(KNNPARAMS));
	free1a(km,np,sizeof(SIZE_T));
	if (dc < 2) free1a(sb,nb,sizeof(fdata));
	free1a(dm,np*rb,sizeof(fdata));
	free2a(keyp,np*rb,nk,sizeof(KEY));
	if (sta > 0) {
		printf("%d-thread: %10.6f seconds, ",np,t = get_time()-t0);
		if (dc <= 1) printf("%5.2f Gflops, ",2*0.000000001*na*nb*nc/t);
		if (dc == 2) printf("%5.2f Gflops, ",5*0.000000001*na*nb*nc/t);
		if (dc == 3) printf("%5.2f Gflops, ",3*0.000000001*na*nb*nc/t);
		printf("%6.3f GBytes/s, ",0.000000001*na*nq*nc*sizeof(fdata)/t);
		printf("%7.3f GBytes/s\n",0.000000001*na*nb*nc*sizeof(fdata)/(t*np));
	}
	return(key);
}

/*------------------------------------------------------------------------------------*/

void *knnl(void *par)
{
	SIZE_T ic,ia,ib;
	double d,sai,dv,vac,vbc;
	fdata *va,*vb;
	KNNPARAMS *kpar = (KNNPARAMS *) par;
	KEY **key = kpar->key;
	fdata **ma = kpar->ma;
	fdata **mb = kpar->mb;
	fdata *sa = kpar->sa;
	fdata *sb = kpar->sb;
	SIZE_T la = kpar->la;
	SIZE_T ha = kpar->ha;
	SIZE_T lb = kpar->lb;
	SIZE_T hb = kpar->hb;
	SIZE_T nc = kpar->nc;
	int dc = kpar->dc;
	for (ia = la; ia < ha; ia++) {
		va = ma[ia];
		if (dc == 3) {
			for (ib = lb; ib < hb; ib++) {
				vb = mb[ib];
				for (d = 0, ic = 0; ic < nc; ic++) {
					dv = vb[ic]-va[ic];
					d += dv*dv;
				}
				key[ib][ia].k = ia;
				key[ib][ia].d = d;
			}
		} else if (dc == 2) {
			for (ib = lb; ib < hb; ib++) {
				vb = mb[ib];
				for (d = 0, ic = 0; ic < nc; ic++) {
					dv = (vbc = vb[ic])-(vac = va[ic]);
					d += dv*dv/(1.0e-30+vbc+vac);
				}
				key[ib][ia].k = ia;
				key[ib][ia].d = d;
			}
		} else {
			if (sa) {
				sai = sa[ia];
			} else {
				for (sai = 0, ic = 0; ic < nc; ic++) sai += va[ic]*va[ic];
				if (dc) sai = sai ? 1/sqrt(sai) : 1;
			}
			for (ib = lb; ib < hb; ib++) {
				vb = mb[ib];
				for (d = 0, ic = 0; ic < nc; ic++) d += vb[ic]*va[ic];
				if (dc) d = -d*sb[ib]*sai; else d = sb[ib]+sai-2*d;
				key[ib][ia].k = ia;
				key[ib][ia].d = d;
			}
		}
	}
	for (ib = lb; ib < hb; ib++) {
		qsort(key[ib]+la,ha-la,sizeof(KEY),keyCompare);
	}
	return(NULL);
}

KEY **knnpl(fdata **ma, fdata **mb, fdata *sa, SIZE_T na, SIZE_T nb, SIZE_T nc,
	    SIZE_T nk, SIZE_T np, SIZE_T cs, int dc, int sta)
{
	SIZE_T ib,ik,ip,lb,hb,rb,iq,nq,*km,kr;
	double t0=0,t;
	fdata *sb = NULL,dr;
	pthread_t *tid;
	KEY **key,**keyp;
	KNNPARAMS *par;
	if (na < nk) return(NULL);
	if (sta) t0 = get_time();
	if (np == 0) np = get_nb_cores();
	if (cs == 0) cs = get_cache_size();
	if ((tid = malloc1a(np,sizeof(pthread_t))) == NULL) return(NULL);
	if ((par = malloc1a(np,sizeof(KNNPARAMS))) == NULL) return(NULL);
	if ((km = malloc1a(np,sizeof(SIZE_T))) == NULL) return(NULL);
	if ((dc < 2) && ((sb = getnorm(mb,nb,nc,dc)) == NULL)) return(NULL);
	if ((key = malloc2a(nb,nk,sizeof(KEY))) == NULL) return(NULL);
	if ((rb = cs/(2*(nc*sizeof(fdata)))) == 0) rb = 1;
	if (rb > nb) rb = nb;
	if (rb > RB_MAX) rb = RB_MAX;
	nq = nb/rb;
	nq += (nq*rb < nb);
	if (sta > 1) printf("%d blocks of %d vectors\n",nq,rb);
	if ((keyp = malloc2a(rb,na,sizeof(KEY))) == NULL) return(NULL);
	for (ip = 0; ip < np; ip++) {
		par[ip].key = keyp;
		par[ip].ma = ma;
		par[ip].mb = mb;
		par[ip].sa = sa;
		par[ip].sb = sb;
		par[ip].la = (ip*na)/np;
		par[ip].ha = (ip+1 == np) ? na : ((ip+1)*na)/np;
		par[ip].nc = nc;
		par[ip].nk = nk;
		par[ip].np = np;
		par[ip].dc = dc;
	}
	for (iq = 0; iq < nq; iq++) {
		lb  = iq*rb;
		hb  = (iq+1 == nq) ? nb : (iq+1)*rb;    
		for (ip = 0; ip < np; ip++) {
			par[ip].lb = lb;
			par[ip].hb = hb;
		}
		for (ip = 1; ip < np; ip++) if (pthread_create(tid+ip,NULL,knnl,par+ip)) {
				printf("error creating thread.");
				abort();
			}
		knnl(par);
		for (ip = 1; ip < np; ip++) if (pthread_join(tid[ip],NULL)) {
				printf("error joining thread.");
				abort();
			}
		for (ib = lb; ib < hb; ib++) {
			for (ip = 0; ip < np; ip++) km[ip] = par[ip].la;
			for (ik = 0; ik < nk; ik++) {
				kr = -1;
				dr = FDATA_MAX;
				for (ip = 0; ip < np; ip++) {
					if ((km[ip] < par[ip].ha) &&(dr > keyp[ib-lb][km[ip]].d)) {
						kr = ip;
						dr = keyp[ib-lb][km[ip]].d;
					}
				}
				key[ib][ik].k = keyp[ib-lb][km[kr]++].k;
				key[ib][ik].d = dr;
			}
		}
		for (ip = 0; ip < np; ip++) {
			par[ip].key -= rb;
		}
	}
	for (ib = 0; ib < nb; ib++) {
		for (ik = 0; ik < nk; ik++) if (dc == 1) {
				if (key[ib][ik].d < -1) key[ib][ik].d = -1;
				if (key[ib][ik].d > +1) key[ib][ik].d = +1;
				key[ib][ik].d = acos(-key[ib][ik].d);
			} else {
				if (key[ib][ik].d < 0) key[ib][ik].d = 0;
				key[ib][ik].d = sqrt(key[ib][ik].d);
			}
	}
	free1a(tid,np,sizeof(pthread_t));
	free1a(par,np,sizeof(KNNPARAMS));
	free1a(km,np,sizeof(SIZE_T));
	if (dc < 2) free1a(sb,nb,sizeof(fdata));
	free2a(keyp,rb,na,sizeof(KEY));
	if (sta > 0) {
		printf("%d-thread: %10.6f seconds, ",np,t = get_time()-t0);
		if (dc <= 1) printf("%5.2f Gflops, ",2*0.000000001*na*nb*nc/t);
		if (dc == 2) printf("%5.2f Gflops, ",5*0.000000001*na*nb*nc/t);
		if (dc == 3) printf("%5.2f Gflops, ",3*0.000000001*na*nb*nc/t);
		printf("%6.3f GBytes/s, ",0.000000001*na*nq*nc*sizeof(fdata)/t);
		printf("%7.3f GBytes/s\n",0.000000001*na*nb*nc*sizeof(fdata)/(t*np));
	}
	return(key);
}

/*------------------------------------------------------------------------------------*/

KEY **knnp(fdata **ma, fdata **mb, fdata *sa, SIZE_T na, SIZE_T nb, SIZE_T nc,
	   SIZE_T nk, SIZE_T np, SIZE_T cs, int dc, int sta)
{
	KEY **key;
	key = (nk < 1024) ? knnps(ma,mb,sa,na,nb,nc,nk,np,cs,dc,sta)
		: knnpl(ma,mb,sa,na,nb,nc,nk,np,cs,dc,sta);
	return(key);
}

/*------------------------------------------------------------------------------------*/

KEY **knnc(fdata **ma, fdata **mb, SIZE_T na, SIZE_T nb, SIZE_T nc, SIZE_T nk, int dc,
	   int sta)
{
	SIZE_T ia,ib,ic,ik;
	double t0 = 0,t,dd,da,db,ds;
	fdata *va,*vb,dv;
	KEY **key,*keyp;
	if (sta) t0 = get_time();
	if ((key = malloc2a(nb,nk,sizeof(KEY))) == NULL) return(NULL);
	if ((keyp = malloc1a(na,sizeof(KEY))) == NULL) return(NULL);
	for (ib = 0; ib < nb; ib++) {
		vb = mb[ib];
		for (ia = 0; ia < na; ia++) {
			va = ma[ia];
			if (dc == 2) {
				for (dd = 0, ic = 0; ic < nc; ic++) {
					dv = va[ic]-vb[ic];
					dd += dv*dv/(1.0e-30+va[ic]+vb[ic]);
				}
			} else if (dc == 1) {
				for (da = db = ds = 0, ic = 0; ic < nc; ic++) {
					da += va[ic]*va[ic];
					db += vb[ic]*vb[ic];
					ds += va[ic]*vb[ic];
				}
			} else {
				for (dd = 0, ic = 0; ic < nc; ic++) {
					dv = va[ic]-vb[ic];
					dd += dv*dv;
				}
			}
			keyp[ia].k = ia;
			keyp[ia].d = (dc == 2) ? sqrt(dd) :
				(dc == 1) ? acos(ds/sqrt(da*db ? da*db : 1)) : sqrt(dd);
		}
		qsort(keyp,na,sizeof(KEY),keyCompare);
		for (ik = 0; ik < nk; ik++) {
			key[ib][ik].k = keyp[ik].k;
			key[ib][ik].d = keyp[ik].d;
		}
	}
	free1a(keyp,na,sizeof(KEY));
	if (sta) {
		printf("1-thread: %10.6f seconds, ",t = get_time()-t0);
		printf("%5.2f Gflops, ",0.000000002*na*nb*nc*(1+(dc == 2))/t);
		printf("%6.3f GBytes/s, ",0.000000001*na*nb*nc*sizeof(fdata)/t);
		printf("%7.3f GBytes/s\n",0.000000001*na*nb*nc*sizeof(fdata)/t);
	}
	return(key);
}

void knncheck(KEY **key0, KEY **key1, SIZE_T nb, SIZE_T nk)
{
	SIZE_T ib,ik;
	fdata d,dm = 0;
	for (ib = 0; ib < nb; ib++) {
		for (ik = 0; ik < nk ; ik++) {
			if (key0[ib][ik].k != key1[ib][ik].k) {
				printf("key0[%d][%d].k = %d, ",ib,ik,key0[ib][ik].k);
				printf("key0[%d][%d].d = %f\n",ib,ik,key0[ib][ik].d);
				printf("key1[%d][%d].k = %d, ",ib,ik,key1[ib][ik].k);
				printf("key1[%d][%d].d = %f\n",ib,ik,key1[ib][ik].d);
			}
			if (dm < (d = fabs(key0[ib][ik].d- key1[ib][ik].d))) dm = d;
		}
	}
	printf("Ckeck: maximum distance difference = %.15f\n\n",dm);
}

/*------------------------------------------------------------------------------------*/
