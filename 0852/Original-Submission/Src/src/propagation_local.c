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
 * propagation_local.c                                                      *
 ****************************************************************************/

#include "propagation_local.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


extern IBPropagationList    *IBpropaglist;        /* global propagation list */
extern IBPropagationGlobal  *IBpropaglistglobal;
extern IBPropagationListCtr *IBpropaglistctr;

extern IBVariables           variables;         /* global array of variables */
extern IBConstraints         constraints;       /* array of constraints */


extern int    IBComputableIntervalNewton;       /* interval newton callable ? */
extern double IBPragmaPrecisionShrink;          /* precision of BC3revise */
extern double IBPragmaImprovement;              /* threshold for propagation */
extern long   IBPragmaMaxTime;      /* stop after IBPragmaMaxTime milliseconds */



IBDmodified *IBDMnew(int n)
/***************************************************************************
*  Allocation of new structure for managing at most n modified domains
*/
{
  IBDmodified *d;
  d = (IBDmodified *)malloc(sizeof(IBDmodified));
  d->dom = (int *)malloc(n*sizeof(int));
  d->N = 0;
  return( d );
}

void IBDMfree(IBDmodified *d)
/***************************************************************************
*  Desallocation
*/
{
  free(d->dom);
  free(d);
}


IBPropagationList *IBPLnew(int n)
/***************************************************************************
*  Allocation of a new propagation list of at most n projections
*/
{
  IBPropagationList *l;
  l = (IBPropagationList *)malloc(sizeof(IBPropagationList));
  l->proj = (IBProjection **)malloc(n*sizeof(IBProjection *));
  l->N = l->first = l->end = 0;
  l->size = n;
  return( l );
}

void IBPLfree(IBPropagationList *l)
/***************************************************************************
*  Desallocation
*/
{
  free(l->proj);
  free(l);
}


IBPropagationGlobal *IBPGlobalNew(int n)
/***************************************************************************
*  Allocation of a new propagation list given a set of at most n modified domains
*/
{
  IBPropagationGlobal *l;
  int i;
  int nb = (int)ceil((double)n / (double)sizeof(unsigned long));

  l = (IBPropagationGlobal *)malloc(sizeof(IBPropagationGlobal));
  l->a = (unsigned long *)malloc(nb*sizeof(unsigned long));
  l->N = nb;

  for( i=0; i<nb; i++ ) IBPGlobalValue(l,i) = 0;

  return( l );
}


void IBPGlobalFree(IBPropagationGlobal *l)
/***************************************************************************
*  Desallocation
*/
{
  free( l->a );
  free(l);
}


IBPropagationListCtr *IBPLCnew(int n)
/***************************************************************************
*  Allocation of a new propagation list of at most n constraints
*/
{
  IBPropagationListCtr *l;
  l = (IBPropagationListCtr *)malloc(sizeof(IBPropagationListCtr));
  l->ctr = (IBConstraint **)malloc(n*sizeof(IBConstraint *));
  l->N = l->first = l->end = 0;
  l->size = n;
  return( l );
}

void IBPLCfree(IBPropagationListCtr *l)
/***************************************************************************
*  Desallocation
*/
{
  free(l->ctr);
  free(l);
}


void IBFPropagationBC3(IBPropagationList *l, IBDmodified *d, int init, int allproj)
/***************************************************************************
*  Creation of propagation list l wrt. the modified domains in d
*  used in Algorithm BC3
*
*  init=1 if l is empty => propagation implemented by the
*                          concatenation of the array of projections
*
*  allproj=1 if all the projections are considered by propagation, 0 if only
*  the projections with multiple occurrences of the variables are considered
*/
{
  IBDependencyV *dep;
  IBProjection **proj;
  IBConstraint *ctr;
  int i, j, index, bitpos, n;
  unsigned long val;

  if( init )  /* initialization of propagation in BC3 */
  {
    /* Flags asleep for all the projections are initialized to 1 */
    for( i=0; i<IBCNbCtr(constraints); i++ )
    {
      ctr = IBCCtr(constraints,i);
      for( j=0; j<IBCNbProjMul(ctr); j++ )
      {
        IBCPasleep(IBCmulprj(ctr,j)) = 1;
      }
      if( allproj )
      {
        for( j=0; j<IBCNbProjOne(ctr); j++ )
        {
          IBCPasleep(IBConeprj(ctr,j)) = 1;
        }
      }
    }
 
    if( (IBDMnb(d)>1) || ((IBDMnb(d)==1) && (IBVnb(variables)==1)) )
         /* propagation wrt. all domains */
    {
      /* Concatenation of the arrays of all constraint projections */
      IBPLfirst(l) = IBPLnbelem(l) = 0;

      for( i=0; i<IBCNbCtr(constraints); i++ )
      {
        ctr = IBCCtr(constraints,i);

        /* Propagation over projections (ctr,x) s.t. x occurs more than once in ctr */
        n = IBCNbProjMul(ctr);
        if( n>0 )
        {
          proj = IBCProjMul(ctr);
          memcpy(&(IBPLproj(l,IBPLnbelem(l))),proj,n*sizeof(IBProjection *));
          IBPLnbelem(l) += n;
        }

        if( allproj )
        {
          /* Moreover, propagation over projections (ctr,x) s.t. x occurs once in ctr */
          n = IBCNbProjOne(ctr);
          if( n>0 )
	  {
            proj = IBCProjOne(ctr);
            memcpy(&(IBPLproj(l,IBPLnbelem(l))),proj,n*sizeof(IBProjection *));
            IBPLnbelem(l) += n;
	  }
        }
      }
      IBPLend(l) = IBPLnbelem(l) - 1;

      for( i=IBPLfirst(l); i<=IBPLend(l); i++ )
      {
        IBCPasleep(IBPLproj(l,i)) = 0;
      }

      return;  /* END OF PROPAGATION IN THIS CASE */
    }
  }


  /*--------------------------------------------------------------------------*/
  /* Otherwise, two cases : initialization with one domain or modification
     of one domain in BC3 -> propagation is the same with domain IBDMdom(d,0) */

  dep = IBDepV(variables,IBDMdom(d,0));  /* dependency list of variable */
  IBRoundUp();
  for( i=0; i<IBDVnb(dep); i++ )         /* decoding of the dependency list */
  {
    val   = IBDVval(dep,i);
    index = IBDVindex(dep,i);

    while( val!=0 )
    {
      bitpos = (int)floor(log((double)val)/0.6931471805599453);
      val &= ~(1<<bitpos);
      ctr = IBCCtr(constraints,(bitpos+sizeof(unsigned long)*index));

      /* Propagation over the projections (ctr,x) s.t. x occurs more than once in ctr */
      proj = IBCProjMul(ctr);
      n    = IBCNbProjMul(ctr);
      for( j=0; j<n; j++ )
      {
        if( IBCPasleep(proj[j]) )
        {
          IBPLend(l) = (IBPLend(l)+1) % IBPLsize(l);
          IBPLnbelem(l)++;

          IBPLproj(l,IBPLend(l)) = proj[j];
          IBCPasleep(proj[j]) = 0;
        }
      }

      if( allproj )
      {
        /* Moreover, propagation over the projections (ctr,x) s.t. x occurs once in ctr */
        proj = IBCProjOne(ctr);
        n = IBCNbProjOne(ctr);
        for( j=0; j<n; j++ )
        {
          if( IBCPasleep(proj[j]) )
          {
            IBPLend(l) = (IBPLend(l)+1) % IBPLsize(l);
            IBPLnbelem(l)++;
            IBPLproj(l,IBPLend(l)) = proj[j];
            IBCPasleep(proj[j]) = 0;
	  }
        }
      }
    }
  }
}


int IBFilteringBC3in(IBDomains d, IBDmodified *dmodified, int allproj)
/***************************************************************************
*  Propagation algorithm BC3 for box consistency: this algorithm is
*  an adaptation of AC3 using box consistency narrowing operators
*
*  IBDmodified contains the variable domains previously modified and
*  it is used to initialize propagation
*/
{
  IBItv out;
  int locvar, globvar;
  IBConstraint *ctr;

#if SOFTWARE_PROFILE
       IBClockBegin(IBClockBC3);
#endif

  IBFPropagationBC3(IBpropaglist,dmodified,1,allproj);
  IBDMnb(dmodified) = 1;  /* during propagation in BC3, only one domain
                             is modified at each step */

  while( IBPLnbelem(IBpropaglist) && (IBClockObserve(IBClockSolve)<=IBPragmaMaxTime) )
  {
    ctr     = IBCCtr(constraints,IBCPctr(IBPLproj(IBpropaglist,IBPLfirst(IBpropaglist))));
    locvar  = IBCPvar(IBPLproj(IBpropaglist,IBPLfirst(IBpropaglist)));
    globvar = IBCVglobvar(ctr,locvar);

    /* BC3revise over projection (ctr,globvar/locvar) */
    if( IBNarrowBC3revise(ctr,locvar,globvar,d,out,IBPragmaPrecisionShrink) )
    {
      if( IBIdiffI(IBDomV(d,globvar),out) )          /* the domain has changed */
      {
        /* Improvement factor, e.g. 10% (1-0.9) ; can be used only if out is finite */
        if( (IBWidthI(out) < IBPragmaImprovement*IBWidthI(IBDomV(d,globvar))) ||
            (IBInfiniteI(out)) )
        {
          IBDMdom(dmodified,0) = globvar;                 /* then propagation */
          IBFPropagationBC3(IBpropaglist,dmodified,0,allproj);
        }
        IBCopyI(IBDomV(d,globvar),out);    /* copy of new domain */
	if (IBPragmaPrecisionShrink==0.0)  /* idempotent narrowing operator */
        {
          IBCPasleep(IBPLproj(IBpropaglist,IBPLfirst(IBpropaglist))) = 1;
	}
      }
      else  /* no contraction */
      {
        IBCPasleep(IBPLproj(IBpropaglist,IBPLfirst(IBpropaglist))) = 1;
      }
      IBPLnbelem(IBpropaglist)--;
      IBPLfirst(IBpropaglist) = (IBPLfirst(IBpropaglist)+1) % IBPLsize(IBpropaglist);
    }
    else
    {
#if SOFTWARE_PROFILE
       IBClockEnd(IBClockBC3);
#endif

      return( 0 ); /* FAILURE OF FILTERING */
    }
  }

#if SOFTWARE_PROFILE
       IBClockEnd(IBClockBC3);
#endif

  return( 1 );        /* SUCCESS OF FILTERING */
}


int IBFilteringBC3(IBDomains d, IBDmodified *dmodified)
/***************************************************************************
*  Propagation algorithm BC3
*/
{
  return( IBFilteringBC3in(d,dmodified,1) );
}


void IBFPropagationHC4(IBPropagationListCtr *l, IBDmodified *d, int init, int onlyoccone)
/***************************************************************************
*  Creation of propagation list l wrt. the modified domains in d
*  used in Algorithm HC4
*
*  init=1 if l is empty => propagation is implemented by the
*                          concatenation of the arrays of constraints
*
*  onlyoccone=1 => propagation only wrt. constraints that contain at
*  least one variable occurring once
*      - Algorithm HC4 => no
*      - Algorithm BC5 => yes (since HC4 is combined with BC3 that efficiently 
*                              handles multiple occurrences of variables)
*/
{
  IBDependencyV *dep;
  IBProjection **proj;
  IBConstraint *ctr;
  int i, j, index, bitpos, n;
  unsigned long val;

  if( init )  /* initialization of propagation in HC4 */
  {
    /* Flags asleep for all constraints are initialized to 1 */
    for( i=0; i<IBCNbCtr(constraints); i++ )
    {
      IBCctrAsleep(IBCCtr(constraints,i)) = 1;
    }
 
    if( (IBDMnb(d)>1) || ((IBDMnb(d)==1) && (IBVnb(variables)==1)) )
         /* propagation wrt. all domains */
     {
      /* All the constraints are added in the list */
      IBPLCfirst(l) = IBPLCnbelem(l) = 0;

      for( i=0; i<IBCNbCtr(constraints); i++ )
      {
	/* propagation wrt. the i-th constraint ? */
        if( (!onlyoccone) ||
            (onlyoccone && IBCNbProjOne(IBCCtr(constraints,i))) )
	{
          IBPLCctr(l,IBPLCnbelem(l)) = IBCCtr(constraints,i);
          IBPLCnbelem(l)++;
          IBCctrAsleep(IBCCtr(constraints,i)) = 0;
	}
      }
      IBPLCend(l) = IBPLCnbelem(l) - 1;

     return;  /* END OF PROPAGATION IN THIS CASE */
    }
  }


  /*--------------------------------------------------------------------------*/
  /* Otherwise, two cases : initialization with one domain or modification
     of several domains in HC4 */
  IBRoundUp();

  if( IBDMnb(d)==1 )  /* only one modified domain */
  {
    dep = IBDepV(variables,IBDMdom(d,0));  /* dependency list of variable */
    for( i=0; i<IBDVnb(dep); i++ )         /* decoding of the dependency list */
    {
      val   = IBDVval(dep,i);
      index = IBDVindex(dep,i);

      while( val!=0 )
      {
        bitpos = (int)floor(log((double)val)/0.6931471805599453);
        val &= ~(1<<bitpos);
        ctr = IBCCtr(constraints,(bitpos+sizeof(unsigned long)*index));

	/* propagation wrt. ctr ? */
        if( (IBCctrAsleep(ctr)) &&
            ((!onlyoccone) ||
             (onlyoccone && IBCNbProjOne(IBCCtr(constraints,i)))) )
	{
          IBCctrAsleep(ctr) = 0;
          IBPLCend(l) = (IBPLCend(l)+1) % IBPLCsize(l);
          IBPLCnbelem(l)++;
          IBPLCctr(l,IBPLCend(l)) = ctr;
	}
      }
    }
  }
  /*----------------------------------------------------------------*/
  else                                      /* several domains in d */
  {
    /* Filtering of the constraints depending on the modified domains
       Consider the dependencies of all the variables in d
       and set propaglistglobal at the right place */

    for( i=0; i<IBDMnb(d); i++ )       /* for each modified domain */
    {
      dep = IBDepV(variables,IBDMdom(d,i));     /* dependency list */

      /* propaglistglobal represents all the constraints depending
         on a modified domain */
      for( j=0; j<IBDVnb(dep); j++ ) 
         IBPGlobalValue(IBpropaglistglobal,IBDVindex(dep,j)) |= IBDVval(dep,j);
    }

    /* Decoding of propaglistglobal and then propagation */
    for( i=0; i<IBPGlobalNb(IBpropaglistglobal); i++ ) 
    {
      val = IBPGlobalValue(IBpropaglistglobal,i);
      while( val!=0 )
      {
        bitpos = (int)floor(log((double)val)/0.6931471805599453);
        val &= ~(1<<bitpos);
        ctr = IBCCtr(constraints,(bitpos+sizeof(unsigned long)*i));

	/* Propagation wrt. ctr ? */
        if( (IBCctrAsleep(ctr)) &&
            ((!onlyoccone) ||
             (onlyoccone && IBCNbProjOne(IBCCtr(constraints,i)))) )
	{
          IBCctrAsleep(ctr) = 0;
          IBPLCend(l) = (IBPLCend(l)+1) % IBPLCsize(l);
          IBPLCnbelem(l)++;
          IBPLCctr(l,IBPLCend(l)) = ctr;
	}
      }
      IBPGlobalValue(IBpropaglistglobal,i) = 0;
    }
  }
}


int IBFilteringHC4in(IBDomains d, IBDmodified *dmodified, int onlyoccone)
/***************************************************************************
*  Propagation algorithm HC4 that enforces hull consistency
*
*  IBDmodified contains the variable domains previously modified
*  and it is used to initialize the propagation
*
*  onlyoccone=1 if only the constraints with at least one variable occurring
*  once are considered => in this case, propagation is stopped when no domain
*  of such a variable is contracted

*  onlyoccone=0 if all constraints are considered
*/
{
  IBConstraint *ctr;
  int i, j, globvar;
  IBDomains dsave = IBNewD(IBVnb(variables));

#if SOFTWARE_PROFILE
       IBClockBegin(IBClockHC4);
#endif

  IBFPropagationHC4(IBpropaglistctr,dmodified,1,onlyoccone);

  while( IBPLCnbelem(IBpropaglistctr) && (IBClockObserve(IBClockSolve)<=IBPragmaMaxTime) )
  {
    ctr = IBPLCctr(IBpropaglistctr,IBPLCfirst(IBpropaglistctr));

    IBCopyD(dsave,d,IBVnb(variables));   /* domains are saved */

    /* HC4revise narrowing operator */
    if( IBNarrowHC4revise(ctr,d) )
    {
      /* which domains have been modified ? */
      i = j = 0;
      for( i=0; i<IBCNbVar(ctr); i++ )
      {
        globvar = IBCVglobvar(ctr,i);
        if( IBIdiffI(IBDomV(dsave,globvar),IBDomV(d,globvar)) )
	{
          /* propagation only if the width of domain has been reduced enough, i.e.,
             improvement factor = e.g. 10%  => domain 10% tighter
             the improvement factor can be used only if the domain is finite, otherwise
             the width is infinite and the test below does not succeed */
          if( ((IBWidthI(IBDomV(d,globvar)) < IBPragmaImprovement*IBWidthI(IBDomV(dsave,globvar)))
               || (IBInfiniteI(IBDomV(d,globvar))))
              && ( (!onlyoccone) || (onlyoccone && (IBCVnbocc(ctr,i)==1))) )
	  {
            IBDMdom(dmodified,j) = globvar;
            j++;
	  }
	}
      }
      IBDMnb(dmodified) = j;   /* number of modified domains */

      /* Propagation */
      if( IBDMnb(dmodified) )
      {
        /* ctr is not added in the list since it is already in the list */
        IBFPropagationHC4(IBpropaglistctr,dmodified,0,onlyoccone);

        /* ctr is removed from the list */
        IBPLCfirst(IBpropaglistctr) =
              (IBPLCfirst(IBpropaglistctr)+1) % IBPLCsize(IBpropaglistctr);

	/* whether ctr is kept in the list ? */
        if( IBCNbProjMul(ctr)>0 )  /* HC4revise is not idempotent */
	{
          IBCctrAsleep(ctr) = 0;
          IBPLCend(IBpropaglistctr) = (IBPLCend(IBpropaglistctr)+1) % IBPLCsize(IBpropaglistctr);
          IBPLCctr(IBpropaglistctr,IBPLCend(IBpropaglistctr)) = ctr;
	}
        else                       /* HC4revise is idempotent */
	{
          IBPLCnbelem(IBpropaglistctr)--;
          IBCctrAsleep(ctr) = 1;
	}
       /* Remark: the number of elements does not change */

      }
      else  /* no domain modified, then no propagation */
      {
        IBPLCnbelem(IBpropaglistctr)--;
        IBPLCfirst(IBpropaglistctr) =
              (IBPLCfirst(IBpropaglistctr)+1) % IBPLCsize(IBpropaglistctr);

        IBCctrAsleep(ctr) = 1;
      }
    }
    else
    {
#if SOFTWARE_PROFILE
       IBClockEnd(IBClockHC4);
#endif

       IBFreeD(dsave);
       return( 0 ); /* FAILURE OF FILTERING */
    }
  }

#if SOFTWARE_PROFILE
       IBClockEnd(IBClockHC4);
#endif

  IBFreeD(dsave);
  return( 1 );      /* SUCCESS OF FILTERING */
}


int IBFilteringHC4(IBDomains d, IBDmodified *dmodified)
/***************************************************************************
*  Propagation algorithm HC4
*/
{
  return( IBFilteringHC4in(d,dmodified,0) );
}


int IBFilteringHC3(IBDomains d, IBDmodified *dmodified)
/***************************************************************************
*  Propagation algorithm HC4 that enforces hull consistency
*
*  IBDmodified contains the variable domains previously modified
*  and it is used to intialize the propagation
*
*  Equivalent to HC4 where HC4revise is replaced with HC3revise
*  HC4revise uses unions of intervals while HC3revise uses only intervals
*/
{
  IBConstraint *ctr;
  int i, j, globvar;
  int onlyoccone = 0;
  IBDomains dsave = IBNewD(IBVnb(variables));

#if SOFTWARE_PROFILE
       IBClockBegin(IBClockHC3);
#endif

  /* Propagation wrt. all the constraints */
  IBFPropagationHC4(IBpropaglistctr,dmodified,1,onlyoccone);

  while( IBPLCnbelem(IBpropaglistctr) && (IBClockObserve(IBClockSolve)<=IBPragmaMaxTime) )
  {
    ctr = IBPLCctr(IBpropaglistctr,IBPLCfirst(IBpropaglistctr));

    IBCopyD(dsave,d,IBVnb(variables));   /* domains are saved */

    /* HC3revise narrowing operator */
    if( IBNarrowHC3revise(ctr,d) )
    {
      /* which domains have been modified ? */
      i = j = 0;
      for( i=0; i<IBCNbVar(ctr); i++ )
      {
        globvar = IBCVglobvar(ctr,i);
        if( IBIdiffI(IBDomV(dsave,globvar),IBDomV(d,globvar)) )
	{
          if( ((IBWidthI(IBDomV(d,globvar)) < IBPragmaImprovement*IBWidthI(IBDomV(dsave,globvar)))
               || (IBInfiniteI(IBDomV(d,globvar))))
              && ( (!onlyoccone) || (onlyoccone && (IBCVnbocc(ctr,i)==1))) )
	  {
            IBDMdom(dmodified,j) = globvar;
            j++;
	  }
	}
      }
      IBDMnb(dmodified) = j;

      /* PROPAGATION */
      if( IBDMnb(dmodified) )
      {
        /* ctr is not added in the list since is is already in the list */
        IBFPropagationHC4(IBpropaglistctr,dmodified,0,onlyoccone);

        /* ctr is removed from the list */
        IBPLCfirst(IBpropaglistctr) =
              (IBPLCfirst(IBpropaglistctr)+1) % IBPLCsize(IBpropaglistctr);

        /* ctr is added at the end of the list since HC4revise is not idempotent */

        IBCctrAsleep(ctr) = 0;
        IBPLCend(IBpropaglistctr) = (IBPLCend(IBpropaglistctr)+1) % IBPLCsize(IBpropaglistctr);
        IBPLCctr(IBpropaglistctr,IBPLCend(IBpropaglistctr)) = ctr;

       /* Remark: the number of elements does not change */

      }
      else  /* no domain modified, then no propagation */
      {
        IBPLCnbelem(IBpropaglistctr)--;
        IBPLCfirst(IBpropaglistctr) =
              (IBPLCfirst(IBpropaglistctr)+1) % IBPLCsize(IBpropaglistctr);

        IBCctrAsleep(ctr) = 1;
      }
    }
    else
    {
#if SOFTWARE_PROFILE
       IBClockEnd(IBClockHC3);
#endif

       IBFreeD(dsave);
       return( 0 ); /* FAILURE OF FILTERING */
    }
  }
#if SOFTWARE_PROFILE
       IBClockEnd(IBClockHC3);
#endif

  IBFreeD(dsave);
  return( 1 );      /* SUCCESS OF FILTERING */
}


int IBFilteringHC3decomp(IBDomains d, IBDmodified *dmodified)
/***************************************************************************
*  Propagation algorithm HC3 than enforces hull consistency
*  over a decomposition of the user's constraints ; the decomposition
*  is already computed
*/
{
  return( IBFilteringHC3(d,dmodified) );
}


int IBFilteringBC4(IBDomains d, IBDmodified *dmodified)
/***************************************************************************
*  Propagation algorithm BC4 that enforces box consistency
*    => combination of BC3revise and HC4revise narrowing operators
*
*  IBDmodified contains the variable domains previously modified
*  and it is used to intialize the propagation
*/
{
  int notfixedpoint;
  int dom, i;
  IBDomains dsave = IBNewD(IBVnb(variables));
  
  IBCopyD(dsave,d,IBVnb(variables));    /* domains are saved */
  if( IBDMnb(dmodified)==1 )
       dom = IBDMdom(dmodified,0);      /* one domain modified */
  else dom = -1;                        /* all domains modified */

  /* first application of HC4 (3rd argument=1 since only the constraints
     containing at least one variable occurring once are considered) )*/
  if( IBFilteringHC4in(d,dmodified,1) )
  {
    if( dom==-1 )
    {
      IBDMnb(dmodified) = IBVnb(variables);  /* initialization used for BC3 */
    }
    else
    {
      IBDMnb(dmodified) = 1;
      IBDMdom(dmodified,0) = 0;
    }

  /* first application of BC3 (3rd argument=0 since only the constraint
     projections (c,x) s.t. x occurs one in c are considered */
    if( !IBFilteringBC3in(d,dmodified,0) )
    {
     IBFreeD(dsave);
     return( 0 );  /* FAILURE OF FILTERING from BC3 */
    }
  }
  else
  {
    IBFreeD(dsave);
    return( 0 );  /* FAILURE OF FILTERING from HC4 */
  }

  /* Creation of the list of modified domains */
  IBDMnb(dmodified) = 0;
  for( i=0; i<IBVnb(variables); i++ )
  {
    if( IBIdiffI(IBDomV(d,i),IBDomV(dsave,i)) )
    {
      IBDMdom(dmodified,IBDMnb(dmodified)) = i;
      IBDMnb(dmodified)++;
    }
  }
  
  while( IBDMnb(dmodified) && (IBClockObserve(IBClockSolve)<=IBPragmaMaxTime) )
  {
    IBCopyD(dsave,d,IBVnb(variables));     /* domains are saved */

    /* Call of HC4 */
    if( IBFilteringHC4in(d,dmodified,1) )
    {
      /* Creation of dmodified */
      IBDMnb(dmodified) = 0;
      for( i=0; i<IBVnb(variables); i++ )
      {
        if( IBIdiffI(IBDomV(d,i),IBDomV(dsave,i)) )
        {
          IBDMdom(dmodified,IBDMnb(dmodified)) = i;
          IBDMnb(dmodified)++;
        }
      }

      if( IBDMnb(dmodified) )
      {
        IBCopyD(dsave,d,IBVnb(variables));

	/* Call of BC3 */
        if( IBFilteringBC3in(d,dmodified,0) )
        {
          /* Creation of dmodified */
          IBDMnb(dmodified) = 0;
          for( i=0; i<IBVnb(variables); i++ )
          {
            if( IBIdiffI(IBDomV(d,i),IBDomV(dsave,i)) )
            {
              IBDMdom(dmodified,IBDMnb(dmodified)) = i;
              IBDMnb(dmodified)++;
            }
          }
        }
        else
        {
          IBFreeD(dsave);
          return( 0 );  /* FAILURE OF FILTERING from BC3 */
        }
      }
    }
    else
    {
     IBFreeD(dsave);
     return( 0 );      /* FAILURE OF FILTERING from HC4 */
    }
  }
  IBFreeD(dsave);
  return( 1 );          /* SUCCESS OF FILTERING */
}


int IBFilteringBC5(IBDomains d, IBDmodified *dmodified)
/***************************************************************************
*  Propagation algorithm BC5 than enforces box consistency
*  and the interval Newton method
*
*  It does not compute a fixed-point
*
*  IBDmodified contains the variable domains previously modified
*  and it is used to intialize the propagation
*/
{
  int dom, i, isin;
  IBDomains dsave = IBNewD(IBVnb(variables));

  IBCopyD(dsave,d,IBVnb(variables));
  if( IBDMnb(dmodified)==1 )
       dom = IBDMdom(dmodified,0);      /* one domain modified */
  else dom = -1;                        /* all domains modified */


  /* HEURISTICS: (HC4,BC3,NEWTON) MORE EFFICIENT THAN (BC4,NEWTON),
     i.e., the computation of a fixed-point of (HC4,BC3) in BC4
     may be too long, and a bisection step may be more efficient */
  /* Call of BC4 */
  /*
  if( IBFilteringBC4(d,dmodified) )
  {
    if( IBComputableIntervalNewton )
    {
      if( !IBNarrowIntervalNewton(d) )
      {
         IBFreeD(dsave);
         return( 0 );
      }
    }
  }
  else
  {
     IBFreeD(dsave);
     return( 0 );
  }
  IBFreeD(dsave);
  return( 1 );
  */


  /* Call of HC4 */
  if( IBFilteringHC4in(d,dmodified,1) )
  {
    if( dom==-1 ) IBDMnb(dmodified) = IBVnb(variables);
    else
    {
       /* Creation of dmodified before BC3 */
       IBDMnb(dmodified) = 0; isin = 0;
       for( i=0; i<IBVnb(variables); i++ )
       {
         if( IBIdiffI(IBDomV(d,i),IBDomV(dsave,i)) )
	 {
           IBDMdom(dmodified,IBDMnb(dmodified)) = i;
           IBDMnb(dmodified)++;
           if( i==dom ) isin = 1;
         }
       }
       if( !isin )  /* domain dom has not been modified by HC4
                       but dom belongs to the argument dmodified */
       {
         IBDMdom(dmodified,IBDMnb(dmodified)) = dom;
         IBDMnb(dmodified)++;
       }
    }

    /* Call of BC3 */
    if( IBFilteringBC3in(d,dmodified,0) )
    {
      if( IBComputableIntervalNewton )
      {
        /* Call of Interval Newton */
        if( !IBNarrowIntervalNewton(d) )
        {
           IBFreeD(dsave);
           return( 0 );
        }
      }
    }
    else
    {
       IBFreeD(dsave);
       return( 0 );    /* failure of filtering in BC3 */
    }
  }
  else
  {
     IBFreeD(dsave);
     return( 0 );      /* failure of filtering in HC4 */
  }
  IBFreeD(dsave);
  return( 1 );         /* success of filtering */
}


int IBFilteringHC4Newton(IBDomains d, IBDmodified *dmodified)
/***************************************************************************
*  Propagation algorithm enforcing HC4 and Interval Newton
*/
{
  /* Call of HC4 and propagation wrt. all constraints (3rd argument=0) */
  if( IBFilteringHC4in(d,dmodified,0) )
  {
    if( IBComputableIntervalNewton )
    {
      if( !IBNarrowIntervalNewton(d) )
      {
	return( 0 );    /* failure of filtering in Interval Newton */
      }
    }
  }
  else
  {
     return( 0 );        /* failure of filtering in HC4 */
  }
  return( 1 );           /* success of filtering */
} 



int IBFilteringBC3Newton(IBDomains d, IBDmodified *dmodified)
/***************************************************************************
*  Propagation algorithm enforcing BC3 and Interval Newton
*/
{
  /* Call of BC3 wrt. all constraint projections */
  if( IBFilteringBC3in(d,dmodified,1) )
  {
    if( IBComputableIntervalNewton )
    {
      if( !IBNarrowIntervalNewton(d) )
      {
	return( 0 );    /* failure of filtering in Interval Newton */
      }
    }
  }
  else
  {
     return( 0 );        /* failure of filtering in HC4 */
  }
  return( 1 );           /* success of filtering */
}
