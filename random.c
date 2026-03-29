/* change to Atkinson's method for simulating from Poisson variables for large
   lambda 
   modified Ranjan Maitra, Ames, IA 50014. 2005/10/31.

   change to use the generic C function random. 
   Further modification by Ranjan Maitra, Ames, IA 50014, 2005/11/16.

   TODO: should change the uniform generator to use the AS algorithm, perhaps 
   using the inbuilt random function to initialize.

*/
#include <stdio.h>
#include <stdlib.h>	
#include <math.h>	
#include "random.h"

/*************************************/
void setseed(unsigned int *s)
{
	/*srand48(s)*/
	srandom(s);
}

/*************************************/
long genseed(void)
{
	return((long)random());
}
/*************************************/
double runi(void)
/*  old version ...
    {
    return ((double)rand() + 0.5)/((double)RAND_MAX + 1.0);
    }
*/
{
	/*double drand48();*/
	/*   static int idnum = -1;*/
	/*   return drand48(&idnum);*/
	/*   return drand48();*/

	return ((double)random() + 0.5)/((double)RAND_MAX + 1.0);
}
/*************************************/
double runir(double a,double b)
{
	return (b-a)*runi()+a;
}
/*************************************/
int runii(int na,int nb)
{
	return (int)(runi()*(nb-na+1.))+na;
}
/*************************************/
double rnor(double mu,double sd)
{
	double e,v1,v2,w;
	do {
		v1=2*runi()-1.;
		v2=2*runi()-1.;
		w=v1*v1+v2*v2;      
	}
	while (w>1.);
	e=sqrt((-2.*log(w))/w);
	if(runi() < 0.5)
		return v2*e*sd+mu;
	else
		return v1*e*sd+mu;
}

/*************************************/
int rpois(double mu)
{ 
	double c,b,a,p=1.,k,ck1,ck2,u1,u2,x,y;
	int n=0;
	if(mu<=0)return 0;
	if (mu<30) {
		mu=exp(-mu);
		do{
			p=p*runi();
			n++;
		}while(p>=mu);
		return n-1;
	}
	else {
		c=.767-3.36/mu;
		b=PI/sqrt(3.*mu);
		a=b*mu;
		if (c<=0) {
			printf("Error in Poisson deviate generation.\n");
			exit(1);
		}
		k=log(c)-mu-log(b);
		ck1=0.;
		do {
			ck2=0.;
			do {
				u1=runi();
				x=(a-log(.1e-18+(1.-u1)/u1))/b;
				if(x>-.5) ck2=1.;
			} while (ck2<0.5);
			n = (int)(x+.5);
			u2=runi();
			y=1+exp(a-b*x);
			ck1= a-b*x+log(.1e-18+u2/(y*y));
			ck2= k+n*log(.1e-18+mu)-lgamma(n+1.);
			if(ck1<=ck2) ck1=1.;
		} while (ck1<0.5);
		return n;
	}
}

/*************************************/
double rexp(double lambda)
{
	return -log(runi())/lambda;
}
/*************************************/
double rcauchy(double loc,double scale)
{
	return tan((runi()-0.5)*PI)*scale+loc;
}
/*************************************/
double rncauchy(double loc,double scale)
{
	double nscale;
	nscale = scale*.67;
	return tan((runi()-0.5)*PI)*nscale+loc;
}
/*************************************/
double rgamma(double alpha)
{
	double r1,r2,aa,x,w,c1,c2,c3,c4,c5;
	if(alpha<=0.)return 0.;
	if(alpha == 1.)return rexp(1.);
	if(alpha<1){
		aa=(alpha+EE)/EE;
		do{
			r1=runi();
			r2=runi();
			if(r1>1./aa){
				x = -log(aa*(1.-r1)/alpha);
				if(r2<pow(x,(alpha-1.)))return x;
			}
			else{
				x = pow((aa*r1),(1./alpha));
				if(r2<exp(-x))return x;
			}
		}while(r2<2);
	}
	else{
		c1=alpha-1;
		c2=(alpha-1./(6.*alpha))/c1;
		c3=2./c1;
		c4=c3+2.;
		c5=1./sqrt(alpha);
		do{
			do{
				r1=runi();
				r2=runi();
				if(alpha>2.5)r1=r2+c5*(1.-1.86*r1);
			}while(r1<=0 || r1 >= 1);
			w=c2*r2/r1;
			if(c3*r1+w+1/w <= c4)return c1*w;
			if(c3*log(r1)-log(w)+w<1)return c1*w;
		}while(r2<2);
	}
	exit(1);
}
/*************************************/
double rbeta(double alpha,double beta)
{
	double r1;
	if(alpha <=0. || beta <= 0.) return 0.;
	r1=rgamma(alpha);
	return r1/(r1+rgamma(beta));
}
/*************************************/
double rlogistic(double loc,double scale)
{
	return 1./(1.-runi())-1;
}
/*************************************/
double rlnorm(double norm,double norsd)
{
	return exp(rnor(norm,norsd));
}
/*************************************/
int rbin(int n,double p)
{
	int j=0,k;
	for(k=0;k<n;k++)if(runi()<p)j++;
	return j;
}
/*************************************/
double rweibull(double gamma)
{
	if(gamma<=0)return 0.;
	return pow(rexp(1.),(1/gamma));
}
/*************************************/
double rchisq(double t)
{
	return rgamma(t/2)*2.;
}
/*************************************/
double rf(double t,double u)
{
	if(t<=0 || u <= 0)return 0.;
	return rchisq(t)*u/(t*rchisq(u));
}
/*************************************/
double rstudent(double t)
{
	if(t<=0.)return 0.;
	return rnor(0.,1.)/sqrt(rchisq(t)/t);
}
/*************************************/

