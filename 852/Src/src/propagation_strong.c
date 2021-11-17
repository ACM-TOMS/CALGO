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
 * propagation_strong.c                                                     *
 ****************************************************************************/

#include "propagation_strong.h"

extern double      IBPragmaPrecision3B;  /* precision of 3B consistency */
extern IBVariables variables;            /* global array of variables */
extern double*     IBwidth3B;            /* used to adapt the precision of 3B consistency */
extern long        IBPragmaMaxTime;      /* stop after IBPragmaMaxTime milliseconds */


int IB3BReviseLeft(IBDomains d, int var, IBDmodified *dmodified,
                   IBLocalPropagation f2b, double w, double* out)
/***************************************************************************
*  Returns - 0 if d is inconsistent
*          - 1 if the left bound of d[var] is inconsistent; in this case,
*            *out is equal to the reduced left bound
*          - 2 if the left bound of d[var] is 3B(w) consistent
*
*  If the result is 1 and 2 then [min d[var], min d[var] + w] is consistent
*  If it is equal to 2, an optimization consists in saving the new computed
*  left bound in *out
*
*  The proof is done by using the 2B consistency algorithm f2b
*/
{
  double right,
         save = IBMinI(IBDomV(d,var));
  int allDomain = 0;

  IBRoundUp();
  right = save + w;

  IBDMnb(dmodified) = IBVnb(variables);

  if( right<IBMaxI(IBDomV(d,var)))
  {
    IBMaxI(IBDomV(d,var)) = right;    /* f2b is applied over the left bound of d[var]*/
  }
  else
  {
    allDomain = 1;                    /* f2b is applied over d */
  }

  /* left bound inconsistent ? */
  if( !f2b(d,dmodified) )
  {
    if( allDomain )
    {
      return 0;                      /* d is inconsistent */
    }
    else
    {
      *out = right;                  /* d[var] is consistent but its left bound is inconsistent */
      return 1;
    }
  }
  else
  {
    *out = IBMinI(IBDomV(d,var));    /* optimization: the new value of the left bound of d[var] */
    return 2;
  }
}


int IB3BReviseRight(IBDomains d, int var, IBDmodified *dmodified,
                    IBLocalPropagation f2b, double w, double* out)
/***************************************************************************
*  Returns - 0 if d is inconsistent
*          - 1 if the right bound of d[var] is inconsistent; in this case,
*            *out is equal to the reduced right bound
*          - 2 if the right bound of d[var] is 3B(w) consistent
*
*  If the result is 1 and 2 then [max d[var] - w, max d[var]] is consistent
*  If it is equal to 2, an optimization consists in saving the new computed
*  right bound in *out
*
*  The proof is done by using the 2B consistency algorithm f2b
*/
{
  double left,
         save = IBMaxI(IBDomV(d,var));
  int allDomain = 0;

  IBRoundDown();
  left = save - w;

  IBDMnb(dmodified) = IBVnb(variables);

  if( left>IBMinI(IBDomV(d,var)))
  {
    IBMinI(IBDomV(d,var)) = left;     /* f2b is applied over the right bound of d[var]*/
  }
  else
  {
    allDomain = 1;                    /* f2b is applied over d */
  }

  /* right bound inconsistent ? */
  if( !f2b(d,dmodified) )
  {
    if( allDomain )
    {
      return 0;                      /* d is inconsistent */
    }
    else
    {
      *out = left;                   /* d[var] is consistent but its right bound is inconsistent */
      return 1;
    }
  }
  else
  {
    *out = IBMaxI(IBDomV(d,var));    /* optimization: the new value of the right bound of d[var] */
    return 2;
  }
}


int IBFiltering3BOneStep(IBDomains d, IBDmodified *dmodified, IBLocalPropagation f2b,
                         IBDomains dcopy, double *w)
/***************************************************************************
*  One step of 3B consistency for each variable bound
*  Returns - 0 if d is inconsistent
*          - 1 if d is reduced and it is not 3B(w) consistent
*          - 2 if d is reduced and it is 3B(w) consistent
*/
{
  int i,
      n,
    result = 2;   /* a priori, d is 3B(w) consistent */
  double bound,
         wRatio = 5.0;

  /* Contraction of all domains */
  for( i=0; i<IBVnb(variables) && (IBClockObserve(IBClockSolve)<=IBPragmaMaxTime); i++ )
  {
    w[i] /= wRatio;
    if( w[i]<IBPragmaPrecision3B ) w[i] = IBPragmaPrecision3B;

    /* Copy of d in dcopy, d must not be modified by the call to IB3BReviseLeft */
    IBCopyD(dcopy,d,IBVnb(variables));
    if( !(n=IB3BReviseLeft(dcopy,i,dmodified,f2b,IBwidth3B[i],&bound)) )
    {
      return 0;       /* d is inconsistent */
    }
    else if( n==1 )
    {
      IBMinI(IBDomV(d,i)) = bound;
      result = 1;     /* modification, left bound not 3B(w) consistent */
    }
    else
    {
      IBMinI(IBDomV(d,i)) = bound;  /* optimization */
    }

    IBCopyD(dcopy,d,IBVnb(variables));
    if( !(n=IB3BReviseRight(dcopy,i,dmodified,f2b,IBwidth3B[i],&bound)) )
    {
      return 0;       /* d is inconsistent */
    }
    else if( n==1 )
    {
      IBMaxI(IBDomV(d,i)) = bound;
      result = 1;     /* modification, right bound not 3B(w) consistent */
    }
    else
    {
      IBMaxI(IBDomV(d,i)) = bound;  /* optimization */
    }
  }
  return( result );
}


int IBFilteringWeak3B(IBDomains d, IBDmodified *dmodified, IBLocalPropagation f2b)
/***************************************************************************
*  Weak 3B consistency, using the 2B consistency algorithm f2b
*
*  Returns - 0 if d is inconsistent
*          - 1 if d is consistent (and reduced or not)
*/
{
  int i;
  double bound,
         wRatio = 5.0;     /* 5.0 originates from the experiments */
  IBDomains dcopy;

  IBDMnb(dmodified) = IBVnb(variables);
  if( !f2b(d,dmodified) )  /* reduction of d */
  {
    return( 0 );
  }

  /* Initialization of widths used for the precision at bounds */
  for( i=0; i<IBVnb(variables); i++ )
  {
    IBwidth3B[i] = IBWidthI(IBDomV(d,i));
  }

  dcopy = IBNewD(IBVnb(variables));
  if( !IBFiltering3BOneStep(d,dmodified,f2b,dcopy,IBwidth3B) )
  {
    IBFreeD(dcopy);
    return( 0 );
  }
  else
  {
    IBFreeD(dcopy);
    return( 1 );
  }
}


int IBFiltering3B(IBDomains d, IBDmodified *dmodified, IBLocalPropagation f2b)
/***************************************************************************
*  3B consistency, using the 2B consistency algorithm f2b
*
*  Returns - 0 if d is inconsistent
*          - 1 if d is consistent (and reduced or not)
*/
{
  int i, n, fixedPoint;
  double bound,
         wRatio = 5.0;     /* 5.0 originates from the experiments */
  IBDomains dcopy;

  IBDMnb(dmodified) = IBVnb(variables);
  if( !f2b(d,dmodified) )  /* reduction of d */
  {
    return( 0 );
  }

  /* Initialization of widths used for the precision at bounds */
  for( i=0; i<IBVnb(variables); i++ )
  {
    IBwidth3B[i] = IBWidthI(IBDomV(d,i));
  }

  dcopy = IBNewD(IBVnb(variables));
  fixedPoint = 0;
  while( (!fixedPoint) && (IBClockObserve(IBClockSolve)<=IBPragmaMaxTime) )
  {
    /* Detection of fixed-point: d is consistent and all the parameters
       IBwidth3B[i] are small enough, i.e., <= IBPragmaPrecision3B */

    fixedPoint = 1;

    if( !(n=IBFiltering3BOneStep(d,dmodified,f2b,dcopy,IBwidth3B)) )
    {
      IBFreeD(dcopy);
      return( 0 );         /* inconsistency */
    }
    else if( n==1 )
    {
      fixedPoint = 0;
    }
    else                   /* n==2: no reduction of d, n==1: reduction of d */
    {
      i = 0;
      while( fixedPoint && (i<IBVnb(variables)) )
      {
        if( IBwidth3B[i]>IBPragmaPrecision3B )
	{
          fixedPoint = 0;
	}
	else i++;
      }
    }
  }
  IBFreeD(dcopy);
  return( 1 );
}
