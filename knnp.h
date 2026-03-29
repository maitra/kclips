#ifndef KNNP_H
#define KNNP_H

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

/*#define fdata float
  #define FDATA_MAX FLT_MAX  */

#define fdata double   /* modified by RM */
#define FDATA_MAX DBL_MAX /* modified by RM*/


/*------------------------------------------------------------------------------------*/
#include <sched.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>   /* added by RM */
#include "constants.h"

/*------------------------------------------------------------------------------------*/

typedef struct {
	int k;
	fdata d;
} KEY;

typedef struct thinsample {
	int n;
	int p;
	double *pi;
	double **x;
	double sum; /* sum of the pi's*/
} ThinSample;


typedef struct {
	KEY **key;
	fdata **ma;
	fdata **mb;
	fdata *sa;
	fdata *sb;
	fdata *dm;
	SIZE_T la;
	SIZE_T ha;
	SIZE_T lb;
	SIZE_T hb;
	SIZE_T nc;
	SIZE_T nk;
	SIZE_T np;
	int dc;
} KNNPARAMS;

SIZE_T UsedMemory();
SIZE_T MaxUsedMemory();
void *malloc1a(SIZE_T n0, SIZE_T size);
void free1a(void *v1, SIZE_T n0, SIZE_T size);
void *malloc2a(SIZE_T n0, SIZE_T n1, SIZE_T size);
void free2a(void *v2, SIZE_T n0, SIZE_T n1, SIZE_T size);

SIZE_T get_nb_cores();
SIZE_T get_cache_size();

fdata *getnorm(fdata **m, SIZE_T n, SIZE_T nc, int dc);
KEY **knnp(fdata **ma, fdata **mb, fdata *sa, SIZE_T na, SIZE_T nb, SIZE_T nc,
	   SIZE_T nk, SIZE_T np, SIZE_T cs, int dc, int sta);

KEY **knnc(fdata **ma, fdata **mb, SIZE_T na, SIZE_T nb, SIZE_T nc, SIZE_T nk, int dc,
	   int sta);
void knncheck(KEY **key0, KEY **key1, SIZE_T nb, SIZE_T nk);

/*

  SYNOPSIS
  #include "knnp.h"

  fdata *getnorm(fdata **m, SIZE_T n, SIZE_T nc, int dc);
       
  KEY **knnp(fdata **ma, fdata **mb, fdata *sa, SIZE_T na, SIZE_T nb,
  SIZE_T nc, SIZE_T nk, SIZE_T np, SIZE_T cs, int dc, int sta);

  DESCRIPTION
  knnp() computes the nk nearest neighbors of the mb vectors
  within the ma vectors.

  na :  number of "document" vectors;
  nb :  number of "query" vectors;
  nc :  number of vector components;
  nk :  number of searched nearest neighbors;
  ma :  document vectors, two-dimensional (na*nc) array of fdata;
  mb :  query vectors, two-dimensional (nb*nc) array of fdata;
      
  sa :  pre-computed norm of document vectors using getnorm or NULL.
       
  np :  number of threads to be run in parallel or 0 to let the
  program set it to the number of available cores if
  get_nb_cores() works.
  cs :  specify the cache size on the system or 0 to let the program
  get it from the system if get_cache_size() works.
	    
  dc :  0 : use Euclidian distance with the scalar product trick;
  1 : use angle between vectors;
  2 : use chi square distance; 
  3 : use Euclidian distance without the scalar product trick.
       
  sta : print some statistcal information : knnp execution time,
  number of floating point operations per second, main
  memory bandwidth, cache bandwidth (per core).
       
  fdata is defined in knnp.h. It may be either float or double.
  FDATA_MAX should be defined accordingly as FLT_MAX or DBL_MAX.
  Default is double. // changed by RM
       
  getnorm() pre-computes the norms of a set of vectors. The square
  of the norm for the Euclidian distance and the inverse of the
  norm for the "cosine" distance is actually computed and returned.
       
  RETURN VALUE
  knnp() returns a two-dimensional (nb*nk) arrays of KEY structures:
       
  typedef struct {
  int k;
  fdata d;
  } KEY;

  There is one vector of KEY structures for each query vector.
       
  k : index of the jth nearest neighbor;
  d : distance to the jth nearest neighbor;
       
  getnorm() returns a one-dimensional array of fdata.
       
  On failure for any reason, knnp() and getnorm() return NULL.
  
*/

KEY **knnp_bridge(fdata **ma, fdata **mb, fdata *sa, SIZE_T na, SIZE_T nb,
		  SIZE_T nc, SIZE_T nk, int dc);

/* similar knnp (indeed is a wrapper to knnp) but has a reduced argument set. 
   All arguments of knnp which are also in knnp_bridge are the same*/


/*------------------------------------------------------------------------------------*/

void remove_duplicates(fdata **ma, SIZE_T *na, SIZE_T nc, SIZE_T dc);


#endif /* KNNP_H */

