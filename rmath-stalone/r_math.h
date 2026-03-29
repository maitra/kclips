/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000--2007  R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#ifndef __R_MATH_H__
#define __R_MATH_H__

#include <math.h>
#include <float.h>
#include "util.h"

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif


#define ME_NONE         0
/*      no error */
#define ME_DOMAIN       1
/*      argument out of domain */
#define ME_RANGE        2
/*      value out of range */
#define ME_NOCONV       4
/*      process did not converge */
#define ME_PRECISION    8
/*      does not have "full" precision */
#define ME_UNDERFLOW    16
/*      and underflow occured (important for IEEE)*/

#define ML_NEGINF       ((-1.0) / 0.0)
#define ML_POSINF       (1.0 / 0.0)


#define ML_ERR_return_NAN { ML_ERROR(ME_DOMAIN, ""); return NAN; }
#define ISNAN(x) (isnan(x)!=0)


/* For a long time prior to R 2.3.0 ML_ERROR did nothing.
   We don't report ME_DOMAIN errors as the callers collect ML_NANs into
   a single warning.
 */
#define ML_ERROR(x, s) {                                                       \
	if(x > ME_DOMAIN) {                                                    \
		char *msg = "";                                                \
		switch(x) {                                                    \
			case ME_DOMAIN:                                        \
				msg = "argument out of domain in '%s'\n";      \
				break;                                         \
			case ME_RANGE:                                         \
				msg = "value out of range in '%s'\n";          \
				break;                                         \
			case ME_NOCONV:                                        \
				msg = "convergence failed in '%s'\n";          \
				break;                                         \
			case ME_PRECISION:                                     \
				msg = "full precision may not have been "      \
					"achieved in '%s'\n";                  \
				break;                                         \
			case ME_UNDERFLOW:                                     \
				msg = "underflow occurred in '%s'\n";          \
				break;                                         \
		}                                                              \
		error(x, NO_EXIT, __FILE__, __func__, __LINE__, msg, s);       \
	}                                                                      \
}


/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)                 \
	if (log_p) {                                           \
		if(p > 0)                                      \
			ML_ERR_return_NAN;                     \
		if(p == 0) /* upper bound*/                    \
			return lower_tail ? _RIGHT_ : _LEFT_;  \
		if(p == ML_NEGINF)                             \
			return lower_tail ? _LEFT_ : _RIGHT_;  \
	}                                                      \
	else { /* !log_p */                                    \
		if(p < 0 || p > 1)                             \
			ML_ERR_return_NAN;                     \
		if(p == 0)                                     \
			return lower_tail ? _LEFT_ : _RIGHT_;  \
		if(p == 1)                                     \
			return lower_tail ? _RIGHT_ : _LEFT_;  \
	}

#define R_D_Lval(p) (lower_tail ? (p) : (0.5 - (p) + 0.5))  /*  p  */
#define R_D_Cval(p) (lower_tail ? (0.5 - (p) + 0.5) : (p))  /*  1 - p */

#define R_DT_qIv(p) (log_p ? (lower_tail ? exp(p) : - expm1(p)) \
                               : R_D_Lval(p))

#define R_DT_CIv(p) (log_p ? (lower_tail ? -expm1(p) : exp(p)) \
                               : R_D_Cval(p))




double expm1(double x);
double log1p(double x);
double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p);
int attribute_hidden chebyshev_init(double *dos, int nos, double eta);
double attribute_hidden chebyshev_eval(double x, const double *a, const int n);

#endif
