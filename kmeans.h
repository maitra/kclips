#ifndef KMEANS_H
#define KMEANS_H

#include "constants.h"
#include "util.h"

/**
 * Include macros for allocating vectors and matrices.  Also define base
 * allocators to handle errors according to the program's style.  There are two
 * versions, one that immediately returns with the MEMORY_ALLOCATION_ERROR and
 * another that that simply issues an error.  The caller is responsible for
 * checking failed allocation.  By default, the base macros MAKE_1ARRAY and
 * CMAKE_1ARRAY are set to the return options.  Redefine wherever you do not
 * wish immediate return, for example when you want to perform extra clean-up
 * before leaving a function that has failed.
 */
#include "array.h"
#define MAKE_1ARRAY_RETURN(a,n) do {                                           \
	(a) = malloc((n) * sizeof *(a));                                       \
	if ((a) == NULL)                                                       \
		return error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__,       \
			__func__, __LINE__, NULL);                             \
} while (0)
#define MAKE_1ARRAY_ERROR(a,n) do {                                            \
	(a) = malloc((n) * sizeof *(a));                                       \
	if ((a) == NULL)                                                       \
		error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__, __func__,    \
			__LINE__, NULL);                                       \
} while (0)
#define CMAKE_1ARRAY_RETURN(a,n) do {                                          \
	(a) = calloc((n), sizeof *(a));                                        \
	if ((a) == NULL)                                                       \
		return error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__,       \
			__func__, __LINE__, NULL);                             \
} while (0)
#define CMAKE_1ARRAY_ERROR(a,n) do {                                           \
	(a) = calloc((n), sizeof *(a));                                        \
	if ((a) == NULL)                                                       \
		error(MEMORY_ALLOCATION_ERROR, NO_EXIT, __FILE__, __func__,    \
			__LINE__, NULL);                                       \
} while (0)

/* avoid gcc warnings by undef'ing */
#undef MAKE_1ARRAY
#undef CMAKE_1ARRAY
/* before reseting to rewritten macros */
#define MAKE_1ARRAY MAKE_1ARRAY_ERROR
#define CMAKE_1ARRAY CMAKE_1ARRAY_ERROR

/**
 * Default maximum number of K-means iterations.
 * The user can specify any number of maximum iterations when calling 
 * #kmeans(), but for consistency across the code, this define can be 
 * used.
 */
#define MAX_KMEANS_ITERATIONS	1000

/**
 * Default buffer size for buffered reads.
 */
#define DEFAULT_BUFFER_SIZE	1000

/**
 * Verbosity levels.
 */
enum {
	SILENT,		/*!< no output to screen */
	TERSE,		/*!< minimal output */
	VERBOSE,	/*!< most verbose */
	NUM_VERBOSITY_LEVELS
};

/**
 * How to standardize/scale the data.
 */
#define NO_TRANSFORM	0	/* interpretable bit field */
#define CENTER_DATA	1	/* 1 to center; 2 to scale; */
#define SCALE_DATA	2	/* 3 to standardize data */
#define DATA_TO_PC	4	/* convert to principal compmonents */
#define SCALE_PC	8	/* scale principal compmonents! */
enum {
	NO_SCALING,
	RANGE_SCALING,		/*!< divide by range */
	STANDARD_ERROR_SCALING,	/*!< divide by std. dev. */
};

/**
 * Errors produced by kmeans().
 */
enum {
	KMEANS_NO_ERROR,		/*!< error-free condition */
	KMEANS_NULL_CLUSTER_ERROR,	/*!< K-means produced null cluster */
	KMEANS_CALLER_INPUT_ERROR,	/*!< invalid caller input */
	KMEANS_UNUSED_ERROR,		/*!< not used */
	KMEANS_NO_SEEDS,		/*!< initialization gives 0 seeds */
	KMEANS_NUMBER_ERRORS,		/*!< number of K-means errors */
	KMEANS_EXCEED_ITER_WARNING,	/*!< K-means exceeded max iterations */
	KMEANS_NUMBER_FAULTS		/*!< total number of fault codes */
};

/** 
 * Human-friendly character strings describing each K-means error.
 */
extern const char *KMEANS_ERROR_STRING[];

typedef struct _options options;
typedef struct _data data;
typedef struct _model model;

/**
 * Options structure.
 */
struct _options {
	/* data i/o */
	const char *path;		/*!< path to files */
	const char *xname;		/*!< data filename */
	const char *idname;		/*!< classification output filename */
	const char *wssname;		/*!< wss output filename */
	const char *muname;		/*!< filename with mu */
	const char *auto_rfile;		/*!< rng seed filename */
	const char *rfile;		/*!< rng seed file for repeatability */
	const char *tfile;		/*!< runtime filename */
	const char *iter_file;		/*!< iteration count filename */
	const char *init_tfile;		/*!< initialization time filename */
	const char *efile;		/*!< error filename */
	const char *cfile;		/*!< estimated centroid filename */
	const char *sfile;		/*!< generated seeds filename */
	char *full_filename;		/*!< data file name */
	FILE *finp;			/*!< FILE pointer to data file */
	char *mu_full_filename;		/*!< full path mu filename */
	FILE *fmup;			/*!< FILE pointer to true clusters */
	int write_seeds;		/*!< write seeds to seed file */

	/* run settings */
	SIZE_T K;			/*!< number of clusters [TODO: more] */
	SIZE_T n_datasets;		/*!< number of datasets in data file */
	int conserve_memory;		/*!< computationally intensive, low
					 *   memory version */
	SIZE_T n_max_iter;		/*!< maximum number of iterations */
	double run_seconds;		/*!< no. secs to run stoch. methods */
	double *run_seconds_per;	/*!< as above but per dataset */
	int seed_from_file;		/*!< indicate if seed(s) from file */
	unsigned int seed1;		/*!< random number seeds */
	unsigned int seed2;		/*!< random number seeds */
	int verbosity;			/*!< verbosity level */

	/* initialization method */
	int init_method;		/*!< initialization method; see enums */
	const char *method_name;	/*!< human-friendly name of method */

	/* initialization method parameters */
	int transform_data;		/*!< how to transform data */
	int scaling_method;		/*!< RANGE or STANDARD_ERROR _SCALING */
	double as70_subsample_prop;	/*!< subsample size; prop. of data::n */
	double as70_sphere_radius;	/*!< compute density in spheres */
	double as70_separation_radius;	/*!< pick seeds separated by this */
	int ha75_dim;			/*!< pick single dim for partition */
	int hh93a_dim;			/*!< pick single dim for splitting */
	int hh93a_cluster_choice;	/*!< how to pick cluster to split */
	int bh67_steinley;		/*!< BH67 with Steinley's recommend. */
	double lw67_subsample_prop;	/*!< subsample size; prop. of data::n */
	SIZE_T lw67_subsample_size;	/*!< subsample size */
	SIZE_T lw67_min_subsample_size;	/*!< minimum subsample size: Kx */
	FILE *lw67_estimate_ss_size;	/*!< estimate subsample size */
	int lw67_scaling;		/*!< scaling method */
	double use_tg74;		/*!< use TG74 variant of BW67 */
	double dr96_lower_bound;	/*!< min. no. centers per kd-bucket */
	int ks69_transform_data;	/*!< data transformation to apply */
	int ks69_scaling;		/*!< scaling to apply */
	int ks69_first_seed;		/*!< how to pick first seed */
	int ks69_first_seed_default;	/*!< how to pick first seed */
	int ks69_next_seed;		/*!< how to pick first seed */
	int ks69_next_seed_default;	/*!< how to pick first seed */
	int hs85_use_centroid;		/*!< use centroid as first seed */
	SIZE_T kr90_subsample_size;	/*!< subsample for KR90 */
	double bf98_subsample_prop;	/*!< subsample size; prop. of data::n */
	int rh07_both_methods;		/*!< take best of RH07a and RH07b */
	int rh07_remove_outliers;	/*!< manage outliers */
	double rh07_outlier_proportion;	/*!< proportion of outliers */
	double osy12_alpha;		/*!< \in [1,2] */
	double osy12_tolerance;		/*!< precision to stop ICA iterations */
	int osy12_deflate;		/*!< use deflation: req'd if data::p > model::K */
	SIZE_T osy12_iterations;	/*!< maximum number of iterations in ICA */
	SIZE_T rh07_leaf_size;		/*!< minimum k-d tree bucket size */
	SIZE_T ltty08_dim;		/*!< precision of integer transform */
	SIZE_T hcsz09_mp;		/*!< no. nearest neighbors (mp param) */
	int av07_greedy;		/*!< choose greedy K-means++ */
	double hk05_alpha;		/*!< init. prob. to move observations */
	double hk05_beta;		/*!< prop. decrease in move prob. */
	SIZE_T hk05_Bs;			/*!< stop if Bs iters. w/ no change */
	SIZE_T hk05_Bl;			/*!< stop after Bl iterations */
	double ma09_alpha;		/*!< percent variation explained */
	SIZE_T ma09_num_projections;	/*!< max. no. of projections */
	double hh93b_prop_var;		/*!< prop. of variance explained */
	SIZE_T hh93b_max_pc;		/*!< max. number of projects */
	double orss06b_rho;		/*!< sqrt(epsilon) parameter of ORSS06b
					 * method, measures support for cluster
					 * k over cluster k+1 solution. */
	int orss06b_use_av07;		/*!< AV07 replace ORSS06b first stage */
	int rnd_init;			/*!< no. of iterations to run r+k-means;
					 * negative indicates to convergence */
	int dm13_smart;			/*!< smart seed selection */
	int dm13_greedy;		/*!< greedy seed selection */
	double dm13_precision;		/*!< precision on chi-squared calcs */
	int ka04_mitra_misuse;		/*!< use KA04 misuse of Mitra 2002 */
};

/**
 * Data structure.
 */
struct _data {
	SIZE_T n;		/*!< number of observations */
	SIZE_T p;		/*!< number of dimensions */
	double **x;		/*!< data matrix */
	double *dis;		/*!< stored pairwise distances */
};

/**
 * Model structure.  Variables defining model to fit and variables for k-means
 * algorithm.
 */
struct _model {
	SIZE_T K;		/*!< number of clusters */
	SIZE_T real_K;		/*!< realized K */
	int ifault;		/*!< error status of k-means */
	SIZE_T niter;		/*!< no. of iterations */

	double **means;		/*!< initialization seeds */

	SIZE_T *iclass;		/*!< inferred class assignments */
	SIZE_T *nc;		/*!< inferred class sizes */
	double *w;		/*!< cluster within-sums-of-squares */
	double wsum;		/*!< within-sum-of-squares */
	double init_time;	/*!< initialization time used */
	double kmeans_time;	/*!< time to run k-means */
	double time;		/*!< internal timer */
};


void kmeans(double **a, SIZE_T m, SIZE_T n, double **c, SIZE_T k, SIZE_T *ic1,
	SIZE_T *nc, SIZE_T iter, double *wss, int *ifault); 

const char *kmeans_error(int err);

#endif /* KMEANS_H */
