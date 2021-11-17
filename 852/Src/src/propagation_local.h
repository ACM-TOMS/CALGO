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
 * propagation_local.h                                                      *
 ****************************************************************************/

#ifndef __propagation_local_h
#define __propagation_local_h


#include "narrowing_newton.h"
#include "narrowing_hullbox.h"


/*------ Management of modified domains */
typedef struct
{
  int *dom;       /* indexes of modified domains */
  int N;          /* number of modified domains  */
} IBDmodified;

#define IBDMdom(d,i)  d->dom[i]
#define IBDMnb(d)     d->N

/* memory management */
IBDmodified *IBDMnew  (int n);
void         IBDMfree (IBDmodified *d);


/*------ Management of propagation over a set of projections */
typedef struct
{
  IBProjection **proj;   /* array of projections */
  int first;             /* active projections between first... */
  int end;               /* ...and end */
  int size;
  int N;
} IBPropagationList;

#define IBPLproj(l,i)    l->proj[i]
#define IBPLfirst(l)     l->first
#define IBPLend(l)       l->end
#define IBPLsize(l)      l->size
#define IBPLnbelem(l)    l->N


/* memory management */
IBPropagationList *IBPLnew  (int n);
void               IBPLfree (IBPropagationList *l);


/*------ Management of propagation given a set of modified domains */
typedef struct
{
 unsigned long *a;
 int N;
} IBPropagationGlobal;

#define IBPGlobalNb(l) l->N
#define IBPGlobalValue(l,i) l->a[i]

/* memory management */
IBPropagationGlobal *IBPGlobalNew  (int n);
void                 IBPGlobalFree (IBPropagationGlobal *l);


/*------ Management of propagation over a set of constraints */
typedef struct
{
  IBConstraint **ctr;
  int first;
  int end;
  int size;
  int N;
} IBPropagationListCtr;

#define IBPLCctr(l,i)  l->ctr[i]
#define IBPLCfirst(l)  l->first
#define IBPLCend(l)    l->end
#define IBPLCsize(l)   l->size
#define IBPLCnbelem(l) l->N

/* memory management */
IBPropagationListCtr *IBPLCnew  (int n);
void                  IBPLCfree (IBPropagationListCtr *l);



/*--- Propagation algorithm HC3 with decomposition ---*/
int IBFilteringHC3decomp (IBDomains d, IBDmodified *dmodified);


/*--- Propagation algorithm HC3 without decomposition ---*/
int IBFilteringHC3       (IBDomains d, IBDmodified *dmodified);


/*--- Propagation algorithm HC4 ---*/
int IBFilteringHC4       (IBDomains d, IBDmodified *dmodified);


/*--- Propagation algorithm HC4 and interval Newton ---*/
int IBFilteringHC4Newton (IBDomains d, IBDmodified *dmodified);


/*--- Propagation algorithm BC3 ---*/
int IBFilteringBC3       (IBDomains d, IBDmodified *dmodified);


/*--- Propagation algorithm BC3 and interval Newton ---*/
int IBFilteringBC3Newton (IBDomains d, IBDmodified *dmodified);


/*--- Propagation algorithm BC4 ---*/
int IBFilteringBC4       (IBDomains d, IBDmodified *dmodified);


/*--- Propagation algorithm BC5 ---*/
int IBFilteringBC5       (IBDomains d, IBDmodified *dmodified);


#endif
