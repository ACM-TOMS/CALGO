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
 * constant.h                                                               *
 ****************************************************************************/

#ifndef __constant_h
#define __constant_h

#include "evaluator.h"


/*--- definitions ---*/

struct IBCo
{
  char *name;    /* name of constant */
  IBItv value;   /* value of constant */
};

struct IBConst
{
    struct IBCo *a;    /* array of constants */
    int N;             /* number of constants */
    int Nfree;         /* number of empty places in the array */
};
typedef struct IBConst *IBConstants;


#define IBConstAllocUnit 20


/*--- functions ---*/

/* memory management */
IBConstants IBNewConstants  ();
void        IBFreeConstants (IBConstants a);

/* adds in a a new constant 'name' with value eval(f) */
int         IBAddConstant   (IBConstants a, char *name, IBTree *f);

/* return the value of constant 'name' */
IBInterval *IBGetConstant   (IBConstants a, char *name);

#endif
