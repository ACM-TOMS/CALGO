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
 * narrowing_newton.h                                                       *
 ****************************************************************************/

#ifndef __narrowing_newton_h
#define __narrowing_newton_h


#include "constant.h"


/*--------------------------------------------------------------------------*
 *                     MATRICES (double, interval)                          *
 ---------------------------------------------------------------------------*/

/*-- Matrices of doubles --*/

typedef struct
{
   int nb;
   double **coeff;    /* square matrix nb*nb */
}
IBMDouble;

/* memory management and copy */
IBMDouble *IBMDoubleNew         (int n);
IBMDouble *IBMDoubleNewZero     (int n);
IBMDouble *IBMDoubleNewIdentity (int n);
void       IBMDoubleFree        (IBMDouble *mat);
void       IBMDoubleCopy        (IBMDouble *mcopy, IBMDouble *msource);


/*-- Matrices of intervals --*/

typedef struct
{
   int nb;
   IBItv **coeff;     /* square matrix nb*nb */
}
IBMInterval;


/* memory management and copy */
IBMInterval *IBMIntervalNew     (int n);
IBMInterval *IBMIntervalNewZero (int n);
void         IBMIntervalFree    (IBMInterval *mat);
void         IBMIntervalCopy    (IBMInterval *mcopy, IBMInterval *msource);


/*-- functions --*/
void IBMDImul         (IBMInterval *m, IBMDouble *m1, IBMInterval *m2);
void IBMDDmul         (IBDomains dnew, IBMDouble *m, IBDomains d);
int  IBMDoubleInverse (IBMDouble *inv, IBMDouble *m);




/*--------------------------------------------------------------------------*
 *                           INTERVAL METHODS                               *
 ---------------------------------------------------------------------------*/

/*-- Gauss-Seidel iteration for the linear system Av = b */
int IBGaussSeidelIteration(IBMInterval *A, IBDomains v, IBDomains b);


/*-- Narrowing operator using Interval Newton */
int IBNarrowIntervalNewton (IBDomains d);


/*-- Returns 1 if solution d is safe, 0 otherwise --*/
int IBSafeSolutionIntervalNewton(IBDomains d);

#endif
