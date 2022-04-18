/****************************************************************************
 * RealPaver v. 0.4                                                         *
 *--------------------------------------------------------------------------*
 * Author: Laurent Granvilliers                                             *
 * Copyright (c) 1999-2003 Institut de Recherche en Informatique de Nantes  *
 * Copyright (c) 2004      Laboratoire d'Informatique de Nantes Atlantique  *
 *--------------------------------------------------------------------------*
 * RealPaver is distributed WITHOUT ANY WARRANTY. Read the associated       *
 * COPYRIGHT file for more details.                                         *
 *--------------------------------------------------------------------------*
 * domain.h                                                                 *
 ****************************************************************************/

#ifndef __domain_h
#define __domain_h

#include "union_interval.h"

/* definitions */

typedef IBItv  IBDom;       /* interval domain for a variable */
typedef IBItv *IBDomains;   /* array of domains of variables */

#define IBDomV(d,var) d[var]


/* functions */

/* memory allocation */
IBDomains IBNewD     (int ndom);

/* memory desallocation */
void      IBFreeD    (IBDomains d);

/* dcopy := d, ndom is the number of domains */
void      IBCopyD    (IBDomains dcopy, IBDomains d, int ndom);

/* memory allocation of d' s.t. d' := d */
IBDomains IBNewCopyD (IBDomains d, int ndom);

/* d is empty ? */
int       IBEmptyD   (IBDomains d, int ndom);

/* domain of var := [x1,x2] */
void      IBSetD     (IBDomains d, int var, double x1, double x2);

/* d[i] := [x,x], i=0,...,ndom-1 */
void      IBSetDtoR  (IBDomains d, double x, int ndom);

/* output of d */
void      IBPrintD   (IBDomains d, int ndom, int digits);

#endif
