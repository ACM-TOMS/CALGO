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
 * narrowing_newton.c                                                       *
 ****************************************************************************/

#include "narrowing_newton.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


extern IBVariables variables;      /* global array of constrained variables */
extern IBOperations operations;    /* global array of operations */
extern IBConstraints constraints;  /* global array of constraints */

extern IBMInterval *IBMIzero,        /* Zero matrix of intervals */
                   *IBMIjacobian,    /* Jacobian matrix */
                   *IBMIfinalsystem;
extern IBMDouble   *IBMDmidjacobian, /* Center of Jacobian matrix */
                   *IBMDzero,        /* Zero matrix of doubles */
                   *IBMDidentity,    /* Identity matrix of doubles */
                   *IBMDinverse;

extern IBDomains   IBdnwt1,
                   IBdnwt2,
                   IBdnwt3,
                   IBdnwt4;

extern int IBComputableIntervalNewton;        /* flag: 1 if Taylor can be computed */



IBMDouble *IBMDoubleNew(int n)
/***************************************************************************
*  Allocation of a square matrix n*n of doubles
*/
{
  IBMDouble *mat;
  int i;

  mat = (IBMDouble *)malloc(sizeof(IBMDouble));
  mat->nb = n;
  mat->coeff = (double **)malloc(n*sizeof(double *));
  for( i=0; i<n; i++ )
  {
    mat->coeff[i] = (double *)malloc(n*sizeof(double));
  }
  return( mat );
}


IBMDouble *IBMDoubleNewZero(int n)
/***************************************************************************
*  Allocation of a square matrix n*n of doubles
*/
{
  IBMDouble *mat;
  int i, j;
  mat = IBMDoubleNew(n);            /* allocation of matrix */

  for( i=0; i<n; i++ )
  {
    for( j=0; j<n; j++ )
    {
      mat->coeff[i][j] = 0.0;       /* initialization of intervals */
    }
  }
  return( mat );
}


IBMDouble *IBMDoubleNewIdentity(int n)
/***************************************************************************
*  Allocation of a square matrix n*n of doubles
*/
{
  IBMDouble *mat;
  int i, j;
  mat = IBMDoubleNew(n);            /* allocation of matrix */

  for( i=0; i<n; i++ )
  {
    for( j=0; j<n; j++ )
    {
      if( i==j ) mat->coeff[i][j] = 1.0;
      else       mat->coeff[i][j] = 0.0;
    }
  }
  return( mat );
}


void IBMDoubleFree(IBMDouble *mat)
/***************************************************************************
*  Desallocation of a square matrix of doubles
*/
{
  int i;
  for( i=0; i<mat->nb; i++ )
  {
    free(mat->coeff[i]);
  }
  free(mat->coeff);
  free(mat);
}


void IBMDoubleWrite(FILE *out, IBMDouble *m)
/***************************************************************************
*  To write m on out
*/
{
  int i,j;

  for( i=0; i<m->nb; i++ )
  {
    fprintf(out,"( ");

    for( j=0; j<m->nb; j++ )
    {
      fprintf(out,"%4g ",m->coeff[i][j]);
    }
    fprintf(out,")\n");
  }
}


void IBMDoubleCopy(IBMDouble *mcopy, IBMDouble *msource)
/***************************************************************************
*  Copy msource in mcopy
*/
{
  int i;
  for( i=0; i<msource->nb; i++ )
  {
    memcpy(mcopy->coeff[i],msource->coeff[i],msource->nb*sizeof(double));
  }
}


IBMInterval *IBMIntervalNew(int n)
/***************************************************************************
*  Allocation of a square matrix n*n of intervals
*/
{
  IBMInterval *mat;
  int i;

  mat = (IBMInterval *)malloc(sizeof(IBMInterval));
  mat->nb = n;
  mat->coeff = (IBItv **)malloc(n*sizeof(IBItv *));
  for( i=0; i<n; i++ )
  {
    mat->coeff[i] = (IBItv *)malloc(n*sizeof(IBItv));
  }
  return( mat );
}


IBMInterval *IBMIntervalNewZero(int n)
/***************************************************************************
*  Allocation of a square matrix n*n of intervals
*/
{
  IBMInterval *mat;
  int i, j;
  mat = IBMIntervalNew(n);                /* allocation of matrix */

  for( i=0; i<n; i++ )
  {
    for( j=0; j<n; j++ )
    {
      IBSetI(mat->coeff[i][j],0.0,0.0);   /* initialization of intervals */
    }
  }
  return( mat );
}


void IBMIntervalFree(IBMInterval *mat)
/***************************************************************************
*  Desallocation of a square matrix of intervals
*/
{
  int i;
  for( i=0; i<mat->nb; i++ )
  {
    free(mat->coeff[i]);
  }
  free(mat->coeff);
  free(mat);
}



void IBMIntervalWrite(FILE *out, IBMInterval *m)
/***************************************************************************
*  To write m on out, used for debugging
*/
{
  int i,j;

  for( i=0; i<m->nb; i++ )
  {
    fprintf(out,"( ");

    for( j=0; j<m->nb; j++ )
    {
      IBWriteI(stdout,m->coeff[i][j],9,IBPrintIntervalBounds);
      fprintf(out," ");
    }
    fprintf(out,")\n");
  }
}


void IBMIntervalCopy(IBMInterval *mcopy, IBMInterval *msource)
/***************************************************************************
*  Copy msource in mcopy
*/
{
  int i;
  for( i=0; i<msource->nb; i++ )
  {
    memcpy(mcopy->coeff[i],msource->coeff[i],msource->nb*sizeof(IBItv));
  }
}


void IBMDImul(IBMInterval *m, IBMDouble *m1, IBMInterval *m2)
/***************************************************************************
*  m := m1*m2  (square matrices)
*/
{
  int i, j, k;
  IBItv itv;
  IBMIntervalCopy(m,IBMIzero);    /* m := 0 */

  for( i=0; i<m->nb; i++ )
  {
    for( j=0; j<m->nb; j++ )
    {
      for( k=0; k<m->nb; k++ )
      {
        if( (!(IBIsZeroI(m2->coeff[k][j]))) && (m1->coeff[i][k]!=0.0) )
	{
          IBMulRIinternal(itv,m1->coeff[i][k],m2->coeff[k][j]);
          IBAddII(m->coeff[i][j],m->coeff[i][j],itv);
	}
      }
    }
  }
}


void IBMDDmul(IBDomains dnew, IBMDouble *m, IBDomains d)
/***************************************************************************
*  dnew := m*d
*/
{
  int i,j;
  IBItv itv;

  for( i=0; i<m->nb; i++ )
  {
    IBSetI(IBDomV(dnew,i),0.0,0.0);
    for( j=0; j<m->nb; j++ )
    {
      if( (!(IBIsZeroI(IBDomV(d,j)))) && (m->coeff[i][j]!=0.0) )
      {
        IBMulRIinternal(itv,m->coeff[i][j],IBDomV(d,j));
        IBAddII(IBDomV(dnew,i),IBDomV(dnew,i),itv);
      }
    }
  }
}


int IBMDoubleInverse(IBMDouble *inv, IBMDouble *m)
/***************************************************************************
*  inv := inverse of m
*/
{
  int i, j, k, l, n = m->nb;
  double max, max2, coeff, pivot;
  double **B = inv->coeff,
         **M = m->coeff;

  IBMDoubleCopy(inv,IBMDidentity);    /* inv := Identity */

  k = 0;
  while( k<n )                /* k from row 0 to row n-1 */
  {
    /* Search for the maximal pivot */
    max = fabs(M[k][k]);
    l = k;
    for( i=k+1; i<n; i++ )    /* i from row k+1 to row n */
    {
      max2 = fabs(M[i][k]);
      if( max<max2 )
      {
        max = max2; l = i;
      }
    }

    /* THE MATRIX IS SINGULAR then 0 is returned */
    if( max == 0.0 )
    {
      return( 0 );
      /* OPTIMISATION PROPOSED BY HANSEN : NOT USED HERE
        l = k;
        M[k][k] = 0.000001;
      */
    }

    if( l!=k )                /* Inverse the lines l and k in M and B */
    {
      for( i=k; i<n; i++ )    /* i from column k to column n-1 */
      {
        max = M[l][i];
        M[l][i] = M[k][i];
        M[k][i] = max;
      }
      for( i=0; i<n; i++ )
      {
        max = B[l][i];
        B[l][i] = B[k][i];
        B[k][i] = max;
      }
    }

    /* The pivot is now M[k][k] */
    pivot = M[k][k];

    /* Divides the lines k in M and B by the pivot */
    for( i=k; i<n; i++ ) M[k][i] /= pivot;
    for( i=0; i<n; i++ ) B[k][i] /= pivot;

    for( j=k+1; j<n; j++ )  /* lines under k */
    {
      coeff = M[j][k];
      for( i=k; i<n; i++ ) M[j][i] -= coeff*M[k][i];
      for( i=0; i<n; i++ ) B[j][i] -= coeff*B[k][i];
    }
    k++;
  }

  for( k=n-1; k>0; k-- )
  {
    for( j=k-1; j>=0; j-- )
    {
      coeff = M[j][k];
      M[j][k]  = 0;
      for( i=n-1; i>=0; i-- ) B[j][i] -= coeff*B[k][i];
    }
  }
  return( 1 );
}


int IBMIDiagonallyDominant(IBMInterval *m)
/***************************************************************************
*  Returns 1 if m is diagonally dominant
*/
{
  double mag, mig;
  int i, j;

  for( i=0; i<m->nb; i++ )
  {
    mig = IBMin(fabs(IBMinI(m->coeff[i][i])),fabs(IBMaxI(m->coeff[i][i])));
    mag = 0.0;
    for( j=0; j<m->nb; j++ )
    {
      if( j!=i ) mag += IBMax(fabs(IBMinI(m->coeff[i][j])),
                              fabs(IBMaxI(m->coeff[i][j])));
    }
    if( mig < mag ) return( 0 );
  }
  return( 1 );
}


int IBGaussSeidelChooseRow(IBMInterval *A, int i)
/***************************************************************************
*  Choose the best row for Gauss-Seidel iteration over variable i
*  Returns -1 if no line can be chosen (each has a zero on column i)
*/
{
  int j, line;

  line = i;  /* line to contract domain of variable i */

  if( IBIsZeroI(A->coeff[i][i]) )
  {
    line = -1; j = 0;
    while( (line==-1) && (j<IBVnb(variables)) )
    {
      if( !IBIsZeroI(A->coeff[j][i]) ) line = j;
      else j++;
    }
  }
  return( line );
}


int IBGaussSeidelIteration(IBMInterval *A, IBDomains v, IBDomains b)
/***************************************************************************
*  Gauss-Seidel iteration for system Av=b s.t. v is already initialized
*  returns 0 if failure, 0 otherwise
*/
{
  int i, j, l;
  IBItv itv, itv2;

  for( i=0; i<IBVnb(variables); i++ )
  {
    l = IBGaussSeidelChooseRow(A,i);

    if( l>=0 )
    {
      IBCopyI(itv2,IBDomV(b,l));                     /* itv2 := b[l] */
      for( j=0; j<i; j++ )
      {
        if( !IBIsZeroI(A->coeff[l][j]) )
	{
          IBMulII(itv,A->coeff[l][j],IBDomV(v,j));
          IBSubII(itv2,itv2,itv);                    /* itv2 := b[l] - ... - A[l,j]*x[j] - ... */
	}
      }
      for( j=i+1; j<IBVnb(variables); j++ )
      {
        if( !IBIsZeroI(A->coeff[l][j]) )
	{
          IBMulII(itv,A->coeff[l][j],IBDomV(v,j));
          IBSubII(itv2,itv2,itv);                    /* itv2 := b[l] - ... - A[l,j]*x[j] - ... */
	}
      }

      IBExtDivInterII(IBDomV(v,i),itv2,A->coeff[l][i]);   /* division by A[l,i] and intersection */

      if( IBEmptyI(IBDomV(v,i)) ) return( 0 );
    }
  }
  return( 1 );
}


int IBNarrowIntervalNewton(IBDomains d)
/***************************************************************************
*  Narrowing operator based on the Gauss-Seidel method applied on a system
*  obtained by preconditionning a first order Taylor expansion of the initial
*  system in the case it contains at least as much constraints as variables
*
*  if there are n variables and m>n constraints then the n first constraints
*  are considered
*/
{
  int i, j, isinvertible;
  IBConstraint *ctr;
  IBTree *f;
  IBDomains mid, fmid, finalvector, vector;
  IBItv itv;

#if SOFTWARE_PROFILE
  IBClockBegin(IBClockINwt);
#endif

  mid         = IBdnwt1;
  fmid        = IBdnwt2;
  finalvector = IBdnwt3;
  vector      = IBdnwt4;

  /* Creation of mid = midpoint(d) */
  IBRoundDown();
  for( i=0; i<IBVnb(variables); i++ )
  {
    IBMinI(IBDomV(mid,i)) = IBMidI(IBDomV(d,i));
  }
  IBRoundUp();
  for( i=0; i<IBVnb(variables); i++ )
  {
    IBMaxI(IBDomV(mid,i)) = IBMidI(IBDomV(d,i));
  }

  /* Initialization of the Jacobian matrix and its center matrix to ZERO */
  IBMIntervalCopy(IBMIjacobian,IBMIzero);
  IBMDoubleCopy(IBMDmidjacobian,IBMDzero);


  /* Computation of derivatives and creation of the Jacobian matrix
     and its center matrix */
  for( i=0; i<IBVnb(variables); i++ )
  {
    ctr = IBCCtr(constraints,i);

    IBTevalAll(IBCfunc(ctr),d);
    if( IBEmptyI(IBTfwd(IBCfunc(ctr))) )
    {
      return 0;
    }

    IBCderiv(ctr);

    IBRoundDown();
    for( j=0; j<IBCNbVar(ctr); j++ )
    {
      IBCopyI(IBMIjacobian->coeff[i][IBCVglobvar(ctr,j)],IBCVderiv(ctr,j));
      IBMDmidjacobian->coeff[i][IBCVglobvar(ctr,j)] = IBMidI(IBCVderiv(ctr,j));
    }
  }


  /* Creation of vector -f(m(X)) */
   for( i=0; i<IBVnb(variables); i++ )
  {
    f = IBCfunc(IBCCtr(constraints,i));
    IBTevalAll(f,mid);
    if( IBEmptyI(IBTfwd(f)) )
    {
      return 0;
    }

    IBNegI(IBDomV(fmid,i),IBTfwd(f),itv);
  }

  /* preconditionning only if the Jacobian is not diagonally dominant */
   /*  if( IBMIDiagonallyDominant(IBMIjacobian) )
  {
    isinvertible = 0;
  }
  else
  {
   */

   /* The diagonal dominance is not taken into account
      ex. x=+y=x-y=0, x,y in [-1,1]
          => Jacobian diagonally dominant
          => But need for preconditionning !
   */

    if( !IBMDoubleInverse(IBMDinverse,IBMDmidjacobian) )
    {
      isinvertible = 0;
    }
    else isinvertible = 1;
  /*
  }
  */


  if( isinvertible )
  {
    IBMDImul(IBMIfinalsystem,IBMDinverse,IBMIjacobian); /* inverse*jacobian */
    IBMDDmul(finalvector,IBMDinverse,fmid);             /* inverse*fmid */
  }
  else
  {
    IBMIntervalCopy(IBMIfinalsystem,IBMIjacobian);
    IBCopyD(finalvector,fmid,IBVnb(variables));
  }

  /* FINAL SYSTEM is now: IBMIfinalsystem*vector = finalvector */

  /* Creation of the initial vector = d - mid */
  for( i=0; i<IBVnb(variables); i++ )
  {
    IBSubII(IBDomV(vector,i),IBDomV(d,i),IBDomV(mid,i));
  }


  /* narrowing through Gauss-Seidel iterations */
  if( !IBGaussSeidelIteration(IBMIfinalsystem,vector,finalvector) )
  {
#if SOFTWARE_PROFILE
    IBClockEnd(IBClockINwt);
#endif

    return( 0 );
  }
  else
  {
    /* domain <- domain   inter   vector+midpoint(domain) */
    for( i=0; i<IBVnb(variables); i++ )
    {
      IBAddII(itv,IBDomV(vector,i),IBDomV(mid,i));
      IBInterII(IBDomV(d,i),IBDomV(d,i),itv);
      if( IBEmptyI(IBDomV(d,i)) )
      {
#if SOFTWARE_PROFILE
        IBClockEnd(IBClockINwt);
#endif

        return( 0 );
      }
    }
  }

#if SOFTWARE_PROFILE
  IBClockEnd(IBClockINwt);
#endif

  return( 1 );
}


int IBSafeSolutionIntervalNewton(IBDomains d)
/***************************************************************************
*  Returns 1 if solution d is safe, 0 otherwise
*/
{
  IBConstraint *ctr;
  int i, j, k, globvar, globvar2;
  IBItv itv, itv2, itv3, sum;
  IBDomains mid;

  if( IBComputableIntervalNewton )
  {
    /* mid := Center of domains */
    mid = IBdnwt1;
    IBRoundDown();
    for( i=0; i<IBVnb(variables); i++ )
    {
      IBMinI(IBDomV(mid,i)) = IBMidI(IBDomV(d,i));
    }
    IBRoundUp();
    for( i=0; i<IBVnb(variables); i++ )
    {
      IBMaxI(IBDomV(mid,i)) = IBMidI(IBDomV(d,i));
    }

    for( i=0; i<IBVnb(variables); i++ )
    {
      /* Test with constraint i and variable i */
      ctr = IBCCtr(constraints,i);
      globvar = i;

      j = IBCGlobvarToLocvar(ctr,i);
      if( j==-1 )   /* variable i not in constraint i */
      {
        return( 0 );
      }

      IBTevalAll(IBCfunc(ctr),d);   /* evaluation of f(d) */
      if( IBEmptyI(IBTfwd(IBCfunc(ctr))) )
      {
        return 0;
      }

      IBCderiv(ctr);                /* derivation */

      IBTevalAll(IBCfunc(ctr),mid); /* evaluation of f(midpoint(d)) */
      if( IBEmptyI(IBTfwd(IBCfunc(ctr))) )
      {
        return 0;
      }

      /* 0 in the gradient ? */
      if( IBDoubleInI(IBCVderiv(ctr,i),0.0) )
      {
        return( 0 );
      }

      /* Test with local variable j (i.e. variable i) of constraint i */
      globvar = IBCVglobvar(ctr,j);

      IBCopyI(sum,IBTfwd(IBCfunc(ctr)));

      for( k=0; k<j; k++ )   /* local variable k of c */
      {
        globvar2 = IBCVglobvar(ctr,k);
        IBSubII(itv2,IBDomV(d,globvar2),IBDomV(mid,globvar2));  /* itv2 := d[k] - midpoint[k] */
        IBMulII(itv3,IBCVderiv(ctr,k),itv2);      /* itv3 := df/dk (d) * (d[k] - midpoint[k]) */
        IBAddII(sum,sum,itv3);                  /* Sum_k ( df/dk (d) * (d[k] - midpoint[k]) ) */
      }
      for( k=j+1; k<IBCNbVar(ctr); k++ )        /* remainder of the sum */
      {
        globvar2 = IBCVglobvar(ctr,k);
        IBSubII(itv2,IBDomV(d,globvar2),IBDomV(mid,globvar2));
        IBMulII(itv3,IBCVderiv(ctr,k),itv2);
        IBAddII(sum,sum,itv3);
      }

      IBDivII(itv3,sum,IBCVderiv(ctr,j));      /* itv3 := Sum_k ( df/dk (d) * (d[k] - midpoint[k]) ) / df/dk[i] */
      IBSubII(itv,IBDomV(mid,globvar),itv3);   /* itv := midpoint[i] - Sum_k ( df/dk (d) * (d[k] - midpoint[k]) ) / df/dk[i] */

      if( !IBIncludedII(itv,IBDomV(d,globvar)) )  /* itv is not included in d[i] ? */
      {
        return( 0 );
      }
    }

    return( 1 );   /* all the tests have succeeded */
  }
  else return( 0 );
}
