#ifndef _RFNS_H
#define _RFNS_H
double normal_quantile(double p, double mu, double sigma, int lower_tail, int log_p);
double normal_cdf(double x, double mu, double sigma, int lower_tail, int log_p);
void normal_cdf_both(double x, double *cum, double *ccum, int i_tail, int log_p);
double normal_pdf(double x, double mu, double sigma, int give_log);


double stirling_error(double n);
double bd0_stalone(double x, double np);
double poisson_pdf_raw(double x, double lambda, int give_log);
double poisson_pdf(double x, double lambda, int give_log);
double poisson_cdf(double x, double lambda, int lower_tail, int log_p);

double gamma_cdf_raw (double x, double alph, int lower_tail, int log_p);
double gamma_cdf(double x, double alph, double scale, int lower_tail, int log_p);
double gamma_pdf(double x, double shape, double scale, int give_log);
double chisq_quantile(double p, double df, int lower_tail, int log_p);
double gamma_quantile(double p, double alpha, double scale, int lower_tail, 
		      int log_p);


double ftruncate(double x);
#endif /*_RFNS_H*/
