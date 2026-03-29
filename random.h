#ifndef _RANDOM_H
#define _RANDOM_H
/* header for random.c */

#define PI 3.141592653589793
#define EE 2.718281828459046

long int random(void);   /* declare in order to avoid undef with ISO */
void srandom(unsigned int *seed); /* declare in order to avoid undef with ISO */

double   runi(void);
double runir(double a,double b);
int runii(int na,int nb);
double rnor(double mu,double sd);
int rpois(double mu);
double rexp(double lambda);
double rcauchy(double loc,double scale);
double rncauchy(double loc,double scale);
double rgamma(double alpha);
double rbeta(double alpha,double beta);
double rlogistic(double loc,double scale);
double rlnorm(double norm,double norsd);
int rbin(int n,double p);
double rweibull(double gamma);
double rchisq(double t);
double rf(double t,double u);
double rstudent(double t);
void  setseed(unsigned int *s);
long  genseed(void);

#endif /* _RANDOM_H*/
