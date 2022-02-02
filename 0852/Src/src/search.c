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
 * search.c                                                                 *
 ****************************************************************************/


#include "search.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern IBVariables   variables;           /* array of constrained variables */
extern IBConstraints constraints;         /* global array of constraints */

extern IBPropagation      IBfilter;       /* propagation algorithm */
extern IBLocalPropagation IBfilter2B;     /* 2B-based propagation algorithm used by IBfilter */
extern IBBisectVar        IBbisect;       /* choice function for the bisected variable */
extern IBBisectArity      IBsplit;        /* bisection step: generation of sub-domains */

extern IBMInterval *IBMIjacobian;         /* Jacobian matrix computed in Interval Newton,
                                             used for the MAX REDUCTION strategy */

extern double IBPragmaPrecision;          /* desired precision */
extern int    IBPragmaBisection;          /* bisection stretegy */
extern int    IBPragmaStyleInterval;      /* used for interval printing */
extern int    IBPragmaIntervalDigits;     /* number of digits for interval printing */
extern long   IBPragmaNbGeneratedDomains; /* number of domains examined */
extern unsigned long   IBPragmaMaxSolution; /* maximum number of solutions computed
                                             in the bisection process */
extern int    IBPragmaHullMode;           /* true the result is the hull of all the output boxes */
extern long   IBPragmaMaxTime;            /* stop after IBPragmaMaxTime milliseconds */
extern int    IBPragmaSubpaving;          /* 1 if a subpaving is computed */

extern int    IBPragmaNumberBisection;    /* arity of each bisection operation */


int IBBisectVariableLF(IBDomains d, IBDomains dold, int var)
/***************************************************************************
*  Stragegy LARGEST-FIRST
*
*  Returns the variable which domain will be bisected
*  var is the previous bisected variable in the round-robin strategy
*/
{
  double w;
  int i;

  w   = IBWidthI(IBDomV(d,0)); /* variable 0 is necessarily a user variable */
  var = 0;

  for( i=1; i<IBVnb(variables); i++ )
  {
    if( IBIsDomainBisectable(variables,i,d,w) )
    {
      w = IBWidthI(IBDomV(d,i));
      var = i;
    }
  }
  return( var );
}


int IBBisectVariableRR(IBDomains d, IBDomains dold, int var)
/***************************************************************************
*  Stragegy ROUND-ROBIN
*
*  Returns the variable which domain will be bisected
*  var is the previous bisected variable in the round-robin strategy
*/
{
  var = (var+1)%IBVnb(variables);
  while( !IBIsDomainBisectable(variables,var,d,IBPragmaPrecision) )
  {
    var = (var+1)%IBVnb(variables);
  }
  return( var );
}


int IBBisectVariableMN(IBDomains d, IBDomains dold, int var)
/***************************************************************************
*  Stragegy MAX REDUCTION
*
*  Returns the variable which domain will be bisected
*  var is the previous bisected variable in the round-robin strategy
*/
{
  double w, wvar, sum;
  int i, j;

  var = -1;
  wvar = 0.0;

  for( i=0; i<IBVnb(variables); i++ )
  {
    if( IBIsDomainBisectable(variables,i,d,IBPragmaPrecision) )
    {
      sum = 0.0;
      for( j=0; j<IBVnb(variables); j++ )
      {
        sum += IBMax(fabs(IBMinI(IBMIjacobian->coeff[i][j])),
                     fabs(IBMaxI(IBMIjacobian->coeff[i][j])));
      }
      w = IBWidthI(IBDomV(dold,i)) * sum;

      if( w>wvar )
      {
        wvar = w;
        var = i;
      }
    }
  }
  if( var==-1 )
  {
    return( IBBisectVariableRR(d,dold,var) );
  }
  else return( var );
}


int IBSolutionOuterBox(IBDomains d)
/***************************************************************************
*  d is a solution if no domain d[i] can be bisected !
*/
{
  int i;
  for( i=0; i<IBVnb(variables); i++ )
  {
    if( IBIsDomainBisectable(variables,i,d,IBPragmaPrecision) )
      return( 0 );
  }
  return( 1 );
}


int IBSolutionInnerBox(IBDomains d)
/***************************************************************************
*  Returns 1 if d represents an inner box of the user's model
*/
{
  int i;

  for( i=0; i<IBCNbCtr(constraints); i++ )
  {
    if( IBCPartOfModel(IBCCtr(constraints,i)) )
    {
      if( !IBCInnerBoxOfConstraint(IBCCtr(constraints,i),d) )
      {
        return( 0 );
      }
    }
  }

  return( 1 );
}


double IBPrecisionSolution(IBDomains d)
/***************************************************************************
*  Returns max(width(d[i]))
*/
{
  double size = 0.0,
         size2;
  int i;

  IBRoundUp();
  for( i=0; i<IBVnb(variables); i++ )
  {
    if( IBIsUserVar(variables,i) && ((size2=IBWidthI(IBDomV(d,i))) > size) )
      size = size2;
  }
  return( size );
}


void IBDListBisect2(IBDList *dlist, int var, long *nbdom)
/***************************************************************************
*  Bisection in 2 parts of the domain of var in dlist->first
*
*  *nbdom is the number of generated domains in the bisection process
*/
{
  IBDomains d1, d2;

  if( IBPragmaSubpaving )
  {
    /* Bisected domain: the first element in dlist
       => it is deplaced at the end of the list */
    d1 = IBDListAddLastSugar(dlist,var,IBDListGetFirstDomain(dlist));
    IBDListRemoveFirstSugar(dlist);
  }
  else
  {
    /* Bisected domain: the last element in dlist */
    d1 = IBDListGetLastDomain(dlist);
    IBDListSetLastVar(dlist,var);
  }

  d2 = IBDListAddLastSugar(dlist,var,d1);   /* copy of d1 */
  IBMinI(IBDomV(d1,var)) = IBMaxI(IBDomV(d2,var)) = IBMidI(IBDomV(d1,var));
  *nbdom += 2;
}


void IBDListBisect3(IBDList *dlist, int var, long *nbdom)
/***************************************************************************
*  Bisection in 3 parts of the domain of var in stack->first
*
*  *nbdom is the number of generated domains in the bisection process
*/
{
  double third, twothirds, x, y;
  IBDomains d1, d2, d3;

  if( IBPragmaSubpaving )
  {
    /* Bisected domain: the first element in dlist */
    d1 = IBDListGetFirstDomain(dlist);
  }
  else
  {
    /* Bisected domain: the last element in dlist */
    d1 = IBDListGetLastDomain(dlist);  
  }

  third     = IBMax(IBThirdI(IBDomV(d1,var)),IBMinI(IBDomV(d1,var)));
  twothirds = IBMin(IBTwoThirdsI(IBDomV(d1,var)),IBMaxI(IBDomV(d1,var)));

  if( (third==IBMinI(IBDomV(d1,var))) ||
      (twothirds==IBMaxI(IBDomV(d1,var))) ||
      (third==twothirds) ) {
    /* IBDomV(d1,var) contains at most three floating point-numbers */
    IBDListBisect2(dlist,var,nbdom);      /* bisection in two parts */
    return;
  }

  if( IBPragmaSubpaving )
  {
    /* Bisected domain: the first element in dlist
       => it is deplaced at the end of the list */
    d1 = IBDListAddLastSugar(dlist,var,IBDListGetFirstDomain(dlist));
    IBDListRemoveFirstSugar(dlist);
  }
  else
  {
    /* Bisected domain: the last element in dlist */
    d1 = IBDListGetLastDomain(dlist);
    IBDListSetLastVar(dlist,var);
  }

  d2 = IBDListAddLastSugar(dlist,var,d1);  /* copy of d1 */
  d3 = IBDListAddLastSugar(dlist,var,d1);  /* copy of d1 */
  IBMaxI(IBDomV(d1,var)) = IBMinI(IBDomV(d2,var)) = third;
  IBMaxI(IBDomV(d2,var)) = IBMinI(IBDomV(d3,var)) = twothirds;
  *nbdom += 3;
}


void IBComputeHull(IBDomains dhull, IBDomains d, int nvar, int nsol, double *prec, int isInner)
/***************************************************************************
*  
*/
{
  int i;
  double w;

  if( nsol==1 )  /* first computed output box */
  {
    IBCopyD(dhull,d,nvar);
    if( isInner )
    {
      *prec = 0.0;
    }
    else
    {
      *prec = IBPrecisionSolution(d);
    }
  }
  else
  {
    if( (!isInner) && ((w=IBPrecisionSolution(d))>*prec) )
    {
      *prec = w;
    }

    for( i=0; i<nvar; ++i )
    {
      IBMinI(IBDomV(dhull,i)) = IBMin(IBMinI(IBDomV(dhull,i)),IBMinI(IBDomV(d,i)));
      IBMaxI(IBDomV(dhull,i)) = IBMax(IBMaxI(IBDomV(dhull,i)),IBMaxI(IBDomV(d,i)));
    }
  }
}


int IBBisection(IBDomains d, int Nobisect, int* completeProcess)
/***************************************************************************
*  Bisection algorithm
*/
{
  IBDList* dlist;
  IBDmodified *dmodified;
  IBDomains df, dold, dhull;
  long nbsol = 0;
  char info[20];
  double maxWidthHull;
  int isInner;

  IBDListRemoveSugar removeDomain;
  IBDListGetDomain   getDomain;

  /* Initialization: all domains are declared modified */
  dmodified         = IBDMnew(IBVnb(variables));
  IBDMnb(dmodified) = IBVnb(variables);
  dold              = IBNewD(IBVnb(variables));  /* used in the strategy MAX REDUCTION */

  printf("INITIAL BOX\n");
  IBWriteVdom(stdout,d,variables,IBPragmaIntervalDigits,IBPragmaStyleInterval);

  /* Only a filtering is enforced according to the user's demand */
  if( Nobisect )
  {
    if( (* IBfilter)(d,dmodified,IBfilter2B) )
    /* Success of filtering */
    {
      if( IBSolutionInnerBox(d) )
      {
        printf("\nINNER BOX 1\n");
      }
      else
      {
        printf("\nOUTER BOX 1\n");
      }
      IBWriteVdom(stdout,d,variables,IBPragmaIntervalDigits,IBPragmaStyleInterval);
      printf("\n  precision: %.3g, ",IBPrecisionSolution(d));
      _IBprintlong(info,IBClockObserve(IBClockSolve),1);
      printf("elapsed time: %s ms\n",info);
      nbsol++;
    }
    IBDMfree(dmodified);
    IBFreeD(dold);
    *completeProcess = 1;
    return( nbsol );
  }

  if( IBPragmaHullMode )
  {
    dhull = IBNewD(IBVnb(variables));
  }

  /* The list of domains */
  dlist = IBDListNew();
  IBDListAddLastSugar(dlist,-1,d);

  /* Two bisection modes:
     - Subpaving: dlist is a queue (remove first, add last)
     - Points: dlist is a stack (remove last, add last)
  */
  if( IBPragmaSubpaving )
  {
    removeDomain = IBDListRemoveFirstSugar;
    getDomain    = IBDListGetFirstDomain;
  }
  else
  {
    removeDomain = IBDListRemoveLastSugar;
    getDomain    = IBDListGetLastDomain;
  }

  while( IBBisectIterate(IBPragmaSubpaving,
                         nbsol,
                         IBDListNbElements(dlist),
                         IBPragmaNumberBisection,
                         IBPragmaMaxSolution) &&
	(IBClockObserve(IBClockSolve)<=IBPragmaMaxTime) )
  {
    df = getDomain(dlist);  /* the domain is not removed from dlist */

    if( IBPragmaBisection==IBBisectMaxNarrow )
    {
      IBCopyD(dold,df,IBVnb(variables));
    }

    if( (* IBfilter)(df,dmodified,IBfilter2B) )
    {
      /* success of filtering */
      if( (isInner=IBSolutionInnerBox(df)) || IBSolutionOuterBox(df) )
      {
        nbsol++;
	if( IBPragmaHullMode )
	{
          IBComputeHull(dhull,df,IBVnb(variables),nbsol,&maxWidthHull,isInner);
	}
        else
	{
	  if( isInner )
	  {
            printf("\nINNER BOX %d\n",nbsol);            
	  }
          else
	  {
            if( IBSafeSolutionIntervalNewton(df) ) printf("\nSAFE ");
	    else printf("\n");
            printf("OUTER BOX %d\n",nbsol);
	  }

          IBWriteVdom(stdout,df,variables,IBPragmaIntervalDigits,IBPragmaStyleInterval);
          printf("\n  precision: %.3g, ",IBPrecisionSolution(df));
          _IBprintlong(info,IBClockObserve(IBClockSolve),1);
          printf("elapsed time: %s ms\n",info);
	}

	removeDomain(dlist); /* the domain is removed from dlist */
      }
      else
      {
       (* IBsplit)(dlist,
                   (* IBbisect)(df,dold,IBDListGetLastVar(dlist)),
                   &IBPragmaNbGeneratedDomains);
      }
    }
    else
    {
      /* failure of filtering */
      removeDomain(dlist); /* the domain is removed from dlist */
    }

    /* Initialization of the next propagation step */
    if( IBDListNbElements(dlist) )
    {
      /* Since the filtering is not idempotent then it is often more efficient
         to propagate for the whole system rather than the subpart depending on
         the bisected domain. Otherwise, propagation only wrt. the bisected domain:
         IBDMnb(dmodified) = 1;
         IBDMdom(dmodified,0) = IBStackDomFirstV(stack); */
      IBDMnb(dmodified) = IBVnb(variables);
    }
  }

  if( IBPragmaSubpaving )
  {
    while( IBDListNbElements(dlist) )
    {
      df = removeDomain(dlist);
      nbsol++;
      if( IBPragmaHullMode )
      {
        /* this is not an inner box */
        IBComputeHull(dhull,df,IBVnb(variables),nbsol,&maxWidthHull,IBSolutionInnerBox(df));
      }
      else
      {
	if( IBSolutionInnerBox(df) )
	{
          printf("\nINNER BOX %d\n",nbsol);            
	}
        else
	{
          if( IBSafeSolutionIntervalNewton(df) ) printf("\nSAFE ");
	  else printf("\n");
          printf("OUTER BOX %d\n",nbsol);
	}

        IBWriteVdom(stdout,df,variables,IBPragmaIntervalDigits,IBPragmaStyleInterval);
        printf("\n  precision: %.3g, ",IBPrecisionSolution(df));
        _IBprintlong(info,IBClockObserve(IBClockSolve),1);
        printf("elapsed time: %s ms\n",info);
      }
    }
  }


  if( IBDListNbElements(dlist) )
  {
    *completeProcess = 0;
  }
  else
  {
    *completeProcess = 1;
  }


  if( IBPragmaHullMode && (nbsol>0))
  {
    if( IBSolutionInnerBox(dhull) )
    {
      printf("\nINNER");            
    }
    else
    {
      printf("\nOUTER");
    }

    printf(" BOX: HULL of %d boxes\n",nbsol);
    IBWriteVdom(stdout,dhull,variables,IBPragmaIntervalDigits,IBPragmaStyleInterval);
    printf("\n  precision: %.3g, ",IBPrecisionSolution(dhull));
    _IBprintlong(info,IBClockObserve(IBClockSolve),1);
    printf("elapsed time: %s ms\n",info);
    printf("  precision of largest box: %.3g\n",maxWidthHull);

    IBFreeD(dhull);
  }

  IBDMfree(dmodified);
  IBDListFree(dlist);
  IBFreeD(dold);
  return( nbsol );
}
