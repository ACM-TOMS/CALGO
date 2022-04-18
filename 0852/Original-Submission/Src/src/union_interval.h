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
 * union_interval.h                                                         *
 ****************************************************************************/

#ifndef __union_interval_h
#define __union_interval_h

#include "interval.h"


/*------ definitions ------*/

typedef struct
{
  int N;       /* cardinal of the union */
  int Nfree;   /* number of unused intervals in the union */
  IBItv *t;    /* the ordered set of disjoint intervals */
}
IBUnion;

#define IBUnionN(u)     u->N
#define IBUnionNf(u)    u->Nfree
#define IBUnionI(u)     u->t      /* array of intervals */
#define IBUnionItv(u,j) u->t[j]   /* interval number j in u */

#define IBUnionAllocUnit 2             /* intervals created at each memory allocation */
#define IBIsEmptyU(u) IBUnionN(u)==0   /* empty union ? */


/*------ functions ------*/

/* allocation of a union of N intervals */
IBUnion *IBNewU        (int N);

/* allocation, desallocation */
IBUnion *IBNewEmptyU   ();
void     IBFreeU       (IBUnion *u);
IBUnion *IBNewCopyU    (IBUnion *u);
void     IBResetU      (IBUnion *u);
void     IBReallocU    (IBUnion *u);

/* i := hull(u) */
void     IBHullU       (IBUnion *u, IBItv i);

/* returns 1 if x is in u */
int IBDoubleInU        (IBUnion *u, double x);

/* returns the greatest number in u that is <= x
    hypothesis 1 : x not in u
    hypothesis 2 : such a number is supposed to exist, not (x < u) */
double IBPrevDoubleInU (IBUnion *u, double x);

/* returns the smallest number in u that is >= x
    hypothesis 1 : x not in u
    hypothesis 2 : such a number is supposed to exist, not (x > u) */
double IBNextDoubleInU (IBUnion *u, double x);

/* insertion of i in u[j] */
void     IBShiftRightU (IBUnion *u, int j, IBItv i);

/* deletion of interval u[j] */
void     IBShiftLeftU  (IBUnion *u, int j);

/* u := u union i, i can be modified */
void     IBUnionIU     (IBUnion *u, IBItv i);

/* u := u inter i, i can be modified
   returns 1 if the intersection is not empty */
int      IBInterIU     (IBUnion *u, IBItv i);

/* u := u inter ui, ui can be modified
   returns 1 if the intersection is not empty */
int      IBInterUU     (IBUnion *u, IBUnion *ui);

/* output of u */
void     IBWriteU      (FILE *out, IBUnion *u, int digits, int mode);



/* arithmetic operations and elementary functions */

/* Generic type for interval union functions */
typedef IBUnion * (* IBEvalFwdNode)(IBUnion *, IBUnion *);

IBUnion *IBAddUU       (IBUnion *u1, IBUnion *u2);
IBUnion *IBSubUU       (IBUnion *u1, IBUnion *u2);
IBUnion *IBNegU        (IBUnion *u1, IBUnion *useless);
IBUnion *IBMulUU       (IBUnion *u1, IBUnion *u2);
IBUnion *IBDivUU       (IBUnion *u1, IBUnion *u2);
IBUnion *IBSqrU        (IBUnion *u1, IBUnion *useless);
IBUnion *IBSqrtU       (IBUnion *u1, IBUnion *useless);
IBUnion *IBPowU        (IBUnion *u1, IBUnion *n);
IBUnion *IBDivRelUU    (IBUnion *u1, IBUnion *u2);
IBUnion *IBNthRootRelU (IBUnion *u1, IBUnion *n);
IBUnion *IBLogU        (IBUnion *u1, IBUnion *useless);
IBUnion *IBExpU        (IBUnion *u1, IBUnion *useless);

/* auxiliary functions */
IBUnion *IBDivRelOneII (IBItv i1, IBItv i2);
IBUnion *IBNthRootOneI (IBItv i1, IBItv n);
IBUnion *IBCoshRelOneI (IBItv i1);
IBUnion *IBSinRelOneI  (IBItv i1, IBItv dom);
IBUnion *IBTanRelOneI  (IBItv i1, IBItv dom);
IBUnion *IBAcosRelOneI  (IBItv i1);
IBUnion *IBAsinRelOneI  (IBItv i1);
IBUnion *IBAtanRelOneI  (IBItv i1);


/* Evaluation operations used in Algorithm HC4 */
void IBhc4AddUI        (IBUnion *result, IBUnion *u, IBItv itv);
void IBhc4SubUI        (IBUnion *result, IBUnion *u, IBItv itv);
void IBhc4SubIU        (IBUnion *result, IBItv itv, IBUnion *u);
void IBhc4NegU         (IBUnion *result, IBUnion *u);
void IBhc4MulUI        (IBUnion *result, IBUnion *u, IBItv itv);
void IBhc4DivRelUI     (IBUnion *result, IBUnion *u, IBItv itv);
void IBhc4DivRelIU     (IBUnion *result, IBItv itv, IBUnion *u);
void IBhc4NthRootRelUI (IBUnion *result, IBUnion *u, IBItv itv);
void IBhc4SqrUI        (IBUnion *result, IBUnion *u);
void IBhc4ExpUI        (IBUnion *result, IBUnion *u);
void IBhc4LogUI        (IBUnion *result, IBUnion *u);
void IBhc4MinUI        (IBUnion *result, IBUnion *u, IBItv itv1, IBItv itv2);
void IBhc4MaxUI        (IBUnion *result, IBUnion *u, IBItv itv1, IBItv itv2);
void IBhc4SinhRelUI    (IBUnion *result, IBUnion *u);
void IBhc4AsinhRelUI   (IBUnion *result, IBUnion *u);
void IBhc4CoshRelUI    (IBUnion *result, IBUnion *u);
void IBhc4AcoshRelUI   (IBUnion *result, IBUnion *u);
void IBhc4TanhRelUI    (IBUnion *result, IBUnion *u);
void IBhc4AtanhRelUI   (IBUnion *result, IBUnion *u);
void IBhc4CosRelUI     (IBUnion *result, IBUnion *u, IBItv dom);
void IBhc4SinRelUI     (IBUnion *result, IBUnion *u, IBItv dom);
void IBhc4TanRelUI     (IBUnion *result, IBUnion *u, IBItv dom);
void IBhc4AcosRelUI    (IBUnion *result, IBUnion *u);
void IBhc4AsinRelUI    (IBUnion *result, IBUnion *u);
void IBhc4AtanRelUI    (IBUnion *result, IBUnion *u);
#endif
