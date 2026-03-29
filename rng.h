/**
 * @file rng.h
 * @author Karin S. Dorman
 * @date Last updated
 *   in $Revision: 4002 $
 *   on $Date: 2013-03-08 11:26:39 -0600 (Fri, 08 Mar 2013) $
 *   by $Author: kdorman $
 *
 * In house random number generator.  To use, call allocate_rng() specifying
 * one of the types in the random number generator enum.  Do not forget to
 * free_rng() when done.
 */

#ifndef __RNG_H__
#define __RNG_H__

#include <stdlib.h>
#include <malloc.h>
#include <limits.h>
#include <stdarg.h>

/* You must define MATHLIB_STANDALONE before including R library headers.  See
 * <http://cran.r-project.org/doc/manuals/R-exts.html#Standalone-Mathlib> for
 * more information.
 */
#ifdef __USE_R_LIBRARY__
#define MATHLIB_STANDALONE 1
#include <Rmath.h>
#endif

#define SEED_SRAND srand(hash_time())
#define STD_UNIF rng_unif(global_rng)
#define STD_NORM rng_norm(global_rng, 0, 1)

/**
 * Enum defining random number generator types.
 */
enum {DEFAULT_RNG,	/*!< default RNG is AS183 */
	AS183,		/*!< Applied Statistics Algorithm #183 */
	MARSAGLIA,	/*!< Marsaglia random number generator */
	RRNG,		/*!< Use R's random number generator */
	RAND		/*!< Use C's rand() */
};

typedef struct _rng rng;

extern rng *global_rng;

/**
 * Random number generator.  This struct identifies a type of random number
 * (see rng enum) generator and its seeds.  Different types of random number
 * generators use different numbers of seeds.
 */
struct _rng {
	short type;		/*!< type of random number generator */
	size_t nseeds;		/*!< number of seeds */
	unsigned int *seeds;	/*!< list of seeds */
};

/* function pointer to universal rng generators of different types */
double (*rng_unif)(rng *);
double (*rng_norm)(rng *, double, double);

/* function pointer used to seed rng */
void (*seed_rng)(rng *, ...);

/* function pointer extracts current seeds; most rng's do not need it */
void (*get_seeds)(rng *);

/* function pointer to rng free function */
void (*free_rng)(rng *);

/* You can imagine other functions you can add, like
double (*rexp)(rng *, double);
double (*rnorm)(rng *, double, double);
double (*rgamma)(rng *, double, double);
double (*ar_generator)(rng *, double (*)(double, va_list), ...);
...
You can have multiple implementations for each, allowing the sophisticated user
to choose her algorithm, but choosing sensible defaults for naive users.
 */


/* function pre-declarations */
rng *allocate_rng(short);
int hash_time(void);
int write_seeds(rng *, FILE *);

#endif
