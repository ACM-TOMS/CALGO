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
 * narrowing_hullbox.c                                                      *
 ****************************************************************************/

#include "narrowing_hullbox.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern IBVariables variables;      /* global array of constrained variables */
extern IBOperations operations;    /* global array of operations */
extern IBConstraints constraints;  /* global array of constraints */


int IBNarrowShrinkLeft(IBConstraint *c, int locvar, int globvar, IBDomains d,
                       IBInterval *in, IBInterval *out, double precision)
/***************************************************************************
*  To search for the leftmost zero of the functional part of c
*  precision is used to compute box\phi consistency
*/
{
  IBDomains stack = IBNewD(IBAllocDichotomicSearch);
  IBItv bound;
  IBTree *f;
  struct IBListDepNodes *lnodes;  /* list of nodes depending on this variable */
  int first,                      /* index of first domain in the stack */
      nballoc,                    /* size of the stack: number of allocated domains */
      found,
      resultnarrow;
  double center,
         third,
         twothirds;

  /* Copy of the input domain in the stack */
  first = 0;
  nballoc = IBAllocDichotomicSearch;
  IBCopyI(IBDomV(stack,0),in);
  f = IBCVfunc(c,locvar);
  lnodes = IBCVdepnodes(c,locvar);

  found = 0;
  while( first>=0 && !found )
  {
    IBTevalOnevar(f,IBDomV(stack,first),globvar,lnodes,d);
    if( IBEmptyI(IBTfwd(f)) )
    {
      first--;
    }
    else if( !IBDoubleInI(IBTfwd(f),0.0) )          /* 0 not in f(X) ? */
    {
      first--;
    }
    else
    {
      if( IBCanonicalI(IBDomV(stack,first)) )
      {
        found = 1;
        IBCopyI(out,IBDomV(stack,first));
      }
      else
      {
        /* Left bounded: 0 in f([a,a+]) ? */
        IBMinI(bound) = IBMinI(IBDomV(stack,first));
        IBRoundUp();
        IBMaxI(bound) = IBMin(IBMaxI(IBDomV(stack,first)),IBNextDouble(IBMinI(IBDomV(stack,first))+precision));
        IBTevalOnevar(f,bound,globvar,lnodes,d);

        if( (!IBEmptyI(IBTfwd(f))) && IBDoubleInI(IBTfwd(f),0.0) )
        {
          found = 1;                           /* left bounded */
          IBCopyI(out,IBDomV(stack,first));
        }
        else
        {
          IBMinI(IBDomV(stack,first)) = IBMaxI(bound);       /* [a,a+[ eliminated */

	  /* after the elimination of [a,a+[, the new domain can be canonical */
          if( IBCanonicalI(IBDomV(stack,first)) )
          {
            found = 1;
            IBCopyI(out,IBDomV(stack,first));
          }
	  else {
            third     = IBThirdI(IBDomV(stack,first));
            twothirds = IBTwoThirdsI(IBDomV(stack,first));
 
            if( (third>IBMinI(IBDomV(stack,first))) &&
                (twothirds>third) &&
                (IBMaxI(IBDomV(stack,first))>twothirds) )
	    {
              /* split in 3 parts */
              if( first>=nballoc-2 )
              {
                nballoc += IBAllocDichotomicSearch;
                stack = (IBDomains)realloc(stack,nballoc*sizeof(IBItv));
              }
              IBMinI(IBDomV(stack,first+2)) = IBMinI(IBDomV(stack,first));
              IBMaxI(IBDomV(stack,first+2)) = IBMinI(IBDomV(stack,first+1)) = third;
              IBMinI(IBDomV(stack,first)) = IBMaxI(IBDomV(stack,first+1)) = twothirds;
              first += 2;
	    }
            else
	    {
              /* split in 2 parts */
              if( first==nballoc-1 )
              {
                nballoc += IBAllocDichotomicSearch;
                stack = (IBDomains)realloc(stack,nballoc*sizeof(IBItv));
              }
              center = IBMidI(IBDomV(stack,first));
              IBCopyI(IBDomV(stack,first+1),IBDomV(stack,first));
              IBMaxI(IBDomV(stack,first+1)) = IBMinI(IBDomV(stack,first)) = center;
              first ++;
	    }
	  }
	}
      }
    }
  }
  IBFreeD(stack);
  if( !found )
  {
    IBSetEmptyI(out);
    return( IBNarrowFailure );
  }
  else return( IBNarrowSuccess );
}


int IBNarrowShrinkRight(IBConstraint *c, int locvar, int globvar, IBDomains d,
                        IBInterval *in, IBInterval *out, double precision)
/***************************************************************************
*  To search for the rightmost zero of the functional part of c
*  precision is used to compute box\phi consistency
*/
{
  IBDomains stack = IBNewD(IBAllocDichotomicSearch);
  IBItv bound;
  struct IBListDepNodes *lnodes;  /* list of nodes depending on this variable */
  IBTree *f;
  int first,                      /* index of first domain in the stack */
      nballoc,                    /* size of the stack: number of allocated domains */
      found,
      resultnarrow;
  double center,
         third,
         twothirds;

  /* Copy of the input domain in the stack */
  first = 0;
  nballoc = IBAllocDichotomicSearch;
  IBCopyI(IBDomV(stack,0),in);
  f = IBCVfunc(c,locvar);
  lnodes = IBCVdepnodes(c,locvar);

  found = 0;
  while( first>=0 && !found )
  {
    IBTevalOnevar(f,IBDomV(stack,first),globvar,lnodes,d);
    if( IBEmptyI(IBTfwd(f)) )
    {
      first--;
    }
    else if( !IBDoubleInI(IBTfwd(f),0.0) )          /* 0 not in f(X) ? */
    {
      first--;
    }
    else
    {
      if( IBCanonicalI(IBDomV(stack,first)) )
      {
        found = 1;
        IBCopyI(out,IBDomV(stack,first));
      }
      else
      {
        /* Right bounded: 0 in f([b-,b]) ? */
        IBRoundDown();
        IBMinI(bound) = IBMax(IBMinI(IBDomV(stack,first)),IBPrevDouble(IBMaxI(IBDomV(stack,first))-precision));
        IBMaxI(bound) = IBMaxI(IBDomV(stack,first));
        IBTevalOnevar(f,bound,globvar,lnodes,d);

        if( (!IBEmptyI(IBTfwd(f))) && IBDoubleInI(IBTfwd(f),0.0) )
        {
          found = 1;
          IBCopyI(out,IBDomV(stack,first));
        }
        else
        {
          IBMaxI(IBDomV(stack,first)) = IBMinI(bound);     /* ]b-,b] eliminated */

	  /* after the elimination of ]b-,b], the new domain can be canonical */
          if( IBCanonicalI(IBDomV(stack,first)) )
          {
            found = 1;
            IBCopyI(out,IBDomV(stack,first));
          }
	  else
          {
            third     = IBThirdI(IBDomV(stack,first));
            twothirds = IBTwoThirdsI(IBDomV(stack,first));

            if( (third>IBMinI(IBDomV(stack,first))) &&
                (twothirds>third) &&
                (IBMaxI(IBDomV(stack,first))>twothirds) )
	    {
              /* split in 3 parts */
              if( first>=nballoc-2 )
              {
                nballoc += IBAllocDichotomicSearch;
                stack = (IBDomains)realloc(stack,nballoc*sizeof(IBItv));
              }
              IBMaxI(IBDomV(stack,first+2)) = IBMaxI(IBDomV(stack,first));
              IBMaxI(IBDomV(stack,first+1)) = IBMinI(IBDomV(stack,first+2)) = twothirds;
              IBMaxI(IBDomV(stack,first)) = IBMinI(IBDomV(stack,first+1)) = third;
              first += 2;
	    }
            else
	    {
              /* split in two parts */
              if( first==nballoc-1 )
              {
                nballoc += IBAllocDichotomicSearch;
                stack = (IBDomains)realloc(stack,nballoc*sizeof(IBItv));
              }
              center = IBMidI(IBDomV(stack,first));
              IBCopyI(IBDomV(stack,first+1),IBDomV(stack,first));
              IBMinI(IBDomV(stack,first+1)) = IBMaxI(IBDomV(stack,first)) = center;
              first++;
	    }
	  }
	}
      }
    }
  }
  IBFreeD(stack);
  if( !found )
  {
    IBSetEmptyI(out);
    return( IBNarrowFailure );
  }
  else return( IBNarrowSuccess );
}


int IBNarrowBC3revise(IBConstraint *c, int locvar, int globvar, IBDomains d,
                      IBInterval *out, double precision)
/***************************************************************************
*  Narrowing operator searching for the outermost quasi-zeros
*  of IBCVfunc(c,locvar)
*/
{
  IBInterval *in = IBDomV(d,globvar);
  struct IBListDepNodes *lnodes;  /* list of nodes depending on this variable */
  IBTree *f;
  IBItv out2, in3, bound;
  int n;

  if( IBCrelfunc(c)==IBRelationINT )       /* c is equivalent to integer(x) */
  {
    IBCopyI(out,IBDomV(d,globvar));
    IBToIntegerI(out);                     /* rounding of bounds to integers */

    if( IBEmptyI(out) )
    {
       return( IBNarrowFailure );
    }
    return( IBNarrowSuccess );
  }
  else if( (IBCrelfunc(c)==IBRelationEQU) || (IBCrelfunc(c)==IBRelationSET) )       /* c is equivalent to "func = 0" */
  {
     /* first evaluation: all nodes are considered */
     IBTeval(IBCVfunc(c,locvar),in,globvar,d);
     if( IBEmptyI(IBTfwd(IBCVfunc(c,locvar))) )
     {
       IBSetEmptyI(out);
       return( IBNarrowFailure );
     }

     n = IBNarrowShrinkLeft(c,locvar,globvar,d,in,out,precision);

     if( n==IBNarrowSuccess ) /* IBMinI(out) is the leftmost zero */
     {
       IBSetI(in3,IBMinI(out),IBMaxI(in));

       if( IBNarrowShrinkRight(c,locvar,globvar,d,in3,out2,precision)
           == IBNarrowSuccess )
	 /* Necessary test since IBNarrowShrinkLeft does not compute
            a fixed-point and then may stop before the detection of
            inconsistency => hence, inconsistency can be detected
            by IBNarrowShrinkRight which takes a tighter domain as input */
       {
         IBMaxI(out) = IBMaxI(out2);
         return( IBNarrowSuccess );
       }
       else
       {
         IBSetEmptyI(out);
         return( IBNarrowFailure );
       }
     }
     else
     {
       IBSetEmptyI(out);
       return( IBNarrowFailure );
     }
  }
  else if( IBCrelfunc(c)==IBRelationSUP )  /* c is equivalent to "func >= 0" */
  {
     f = IBCVfunc(c,locvar);
     lnodes = IBCVdepnodes(c,locvar);

     /* ALL nodes of f are considered using IBTeval */
     /* right(f([left,left+])) >= 0 ? */
     IBSetI(bound,IBMinI(in),IBNextDouble(IBMinI(in)));
     IBTeval(f,bound,globvar,d);

     if( IBEmptyI(IBTfwd(f)) )
     {
       IBSetEmptyI(out);
       return( IBNarrowFailure );
     }

     if( IBMaxI(IBTfwd(f)) >= 0.0 )
     {
       IBCopyI(out,in);
     }
     else
     {
       n = IBNarrowShrinkLeft(c,locvar,globvar,d,in,out,precision);
       if( n!=IBNarrowSuccess )
       {
         IBSetEmptyI(out);
         return( IBNarrowFailure );
       }
       else
       {
         IBMaxI(out) = IBMaxI(in);
       }
     }

     /* right(f([right,right+])) >= 0 ? */
     IBSetI(bound,IBPrevDouble(IBMaxI(out)),IBMaxI(out));
     IBTevalOnevar(f,bound,globvar,lnodes,d);
     if( IBEmptyI(IBTfwd(f)) )
     {
       IBSetEmptyI(out);
       return( IBNarrowFailure );
     }
     else if( IBMaxI(IBTfwd(f)) >= 0.0 )
     {
       return( IBNarrowSuccess );
     }
     else
     {
       n = IBNarrowShrinkRight(c,locvar,globvar,d,out,out2,precision);
       if( n!=IBNarrowSuccess )
       {
         IBSetEmptyI(out);
         return( IBNarrowFailure );
       }
       else
       {
         IBMaxI(out) = IBMaxI(out2);
         return( IBNarrowSuccess );
       }
     }
  }
  else                                     /* c is equivalent to "func <= 0" */
  {
     f = IBCVfunc(c,locvar);
     lnodes = IBCVdepnodes(c,locvar);

     /* ALL nodes of f are considered using IBTeval */
     /* left(f([left,left+])) <= 0 ? */
     IBSetI(bound,IBMinI(in),IBNextDouble(IBMinI(in)));
     IBTeval(f,bound,globvar,d);
     if( IBEmptyI(IBTfwd(f)) )
     {
       IBSetEmptyI(out);
       return( IBNarrowFailure );
     }

     if( IBMinI(IBTfwd(f)) <= 0.0 )
     {
       IBCopyI(out,in);
     }
     else
     {
       if( (n=IBNarrowShrinkLeft(c,locvar,globvar,d,in,out,precision))!=IBNarrowSuccess )
       {
         IBSetEmptyI(out);
         return( IBNarrowFailure );
       }
       else
       {
         IBMaxI(out) = IBMaxI(in);
       }
     }

     /* left(f([right,right+])) <= 0 ? */
     IBSetI(bound,IBPrevDouble(IBMaxI(out)),IBMaxI(out));
     IBTevalOnevar(f,bound,globvar,lnodes,d);
     if( IBEmptyI(IBTfwd(f)) )
     {
       IBSetEmptyI(out);
       return( IBNarrowFailure );
     }
     else if( IBMinI(IBTfwd(f)) <= 0.0 )
     {
       return( IBNarrowSuccess );
     }
     else
     {
       n = IBNarrowShrinkRight(c,locvar,globvar,d,out,out2,precision);

       if( n!=IBNarrowSuccess )
       {
         IBSetEmptyI(out);
         return( IBNarrowFailure );
       }
       else
       {
         IBMaxI(out) = IBMaxI(out2);
         return( IBNarrowSuccess );
       }
     }
   }
}


int IBNarrowHC4reviseRec(IBTree *f, IBDomains d)
/***************************************************************************
*  HC4revise on a tree
*/
{
  if( f==NULL ) return( 1 );

  switch( IBTtype(f) )
  {
    case IBTNodeVar:
         if( IBInterIU(IBThc4U(f),IBDomV(d,IBTglobvar(f))) )
	 {
           IBHullU(IBThc4U(f),IBDomV(d,IBTglobvar(f)));
           IBResetU(IBThc4U(f));
           return( 1 );
	 }
         else
         {
           return( 0 );
	 }
         break;
    case IBTNodeOp:
         if( !(* IBTevalHC4(operations,IBTop(f)))(f) )
         {
	   return( 0 );
         }

         if( IBNarrowHC4reviseRec(IBTleft(f),d) )
            return( IBNarrowHC4reviseRec(IBTright(f),d) );
         else return( 0 );
         break;
    default:
         IBResetU(IBThc4U(f));
         return( 1 );
         break;
  }
}


int IBNarrowHC4revise(IBConstraint *c, IBDomains d)
/***************************************************************************
*  HC4revise narrowing operator
*/
{
  IBItv itv;
  int result, globvar;

  if( IBCrelfunc(c)==IBRelationINT )       /* c is equivalent to integer(x) */
  {
    globvar = IBCVglobvar(c,0);            /* only one variable at index 0 */
    IBToIntegerI(IBDomV(d,globvar));       /* rounding of bounds to integers */
    if( IBEmptyI(IBDomV(d,globvar)) )
    {
       return( IBNarrowFailure );
    }
    return( IBNarrowSuccess );
  }

  /* Interval evaluations of left and right terms */
  IBTevalAll(IBCleft(c),d);
  if( IBEmptyI(IBTfwd(IBCleft(c))) )
  {
    return IBNarrowFailure;
  }

  IBTevalAll(IBCright(c),d);
  if( IBEmptyI(IBTfwd(IBCright(c))) )
  {
    return IBNarrowFailure;
  }

  /* Interpretation of relation symbol */
  if( (IBCrel(c)==IBRelationEQU)  || (IBCrelfunc(c)==IBRelationSET) )
  {
    IBInterII(itv,IBTfwd(IBCleft(c)),IBTfwd(IBCright(c)));
    if( IBEmptyI(itv) )
    {
       return( IBNarrowFailure );
    }

    IBUnionIU(IBThc4U(IBCleft(c)),itv);
    IBUnionIU(IBThc4U(IBCright(c)),itv);
  }
  else if( IBCrel(c)==IBRelationINF )
  {
    /* Left term */
    IBMinI(itv) = IBMinI(IBTfwd(IBCleft(c)));
    IBMaxI(itv) = IBMin(IBMaxI(IBTfwd(IBCleft(c))),IBMaxI(IBTfwd(IBCright(c))));
    if( IBEmptyI(itv) )
    {
       return( IBNarrowFailure );
    }

    IBUnionIU(IBThc4U(IBCleft(c)),itv);

    /* right term */
    IBMaxI(itv) = IBMaxI(IBTfwd(IBCright(c)));
    IBMinI(itv) = IBMax(IBMinI(IBTfwd(IBCleft(c))),IBMinI(IBTfwd(IBCright(c))));
    if( IBEmptyI(itv) )
    {
       return( IBNarrowFailure );
    }

    IBUnionIU(IBThc4U(IBCright(c)),itv);
  }
  else        /* IBRelationSUP */
  {
    /* left term */
    IBMaxI(itv) = IBMaxI(IBTfwd(IBCleft(c)));
    IBMinI(itv) = IBMax(IBMinI(IBTfwd(IBCleft(c))),IBMinI(IBTfwd(IBCright(c))));
    if( IBEmptyI(itv) )
    {
       return( IBNarrowFailure );
    }

    IBUnionIU(IBThc4U(IBCleft(c)),itv); 

    /* right term */
    IBMinI(itv) = IBMinI(IBTfwd(IBCright(c)));
    IBMaxI(itv) = IBMin(IBMaxI(IBTfwd(IBCleft(c))),IBMaxI(IBTfwd(IBCright(c))));
    if( IBEmptyI(itv) )
    {
       return( IBNarrowFailure );
    }

    IBUnionIU(IBThc4U(IBCright(c)),itv);
  }

  /* Backward propagation */
  if( IBNarrowHC4reviseRec(IBCleft(c),d) )
  {
    result = IBNarrowHC4reviseRec(IBCright(c),d);
    return( result );
  }
  else
  {
     return( IBNarrowFailure );
  }
}


int IBNarrowHC3reviseRec(IBTree *f, IBDomains d)
/***************************************************************************
*  HC3revise on a tree
*/
{
  if( f==NULL ) return( 1 );

  switch( IBTtype(f) )
  {
    case IBTNodeVar:
         IBInterII(IBDomV(d,IBTglobvar(f)),IBDomV(d,IBTglobvar(f)),IBTbwd(f));

         if( IBEmptyI(IBDomV(d,IBTglobvar(f))) )
	 {
           return( 0 );
	 }
         else
         {
          return( 1 );
	 }
         break;
    case IBTNodeOp:
         if( (* IBTevalHC3(operations,IBTop(f)))(f) )
	 {
           if( IBNarrowHC3reviseRec(IBTleft(f),d) )
              return( IBNarrowHC3reviseRec(IBTright(f),d) );
           else return( 0 );
	 }
         else
	 {
           return( 0 );
         }
         break;
    default:
         return( 1 );
         break;
  }
}


int IBNarrowHC3revise(IBConstraint *c, IBDomains d)
/***************************************************************************
*  HC3revise narrowing operator
*/
{
  IBItv itv;
  int globvar;

  if( IBCrelfunc(c)==IBRelationINT )       /* c is equivalent to integer(x) */
  {
    globvar = IBCVglobvar(c,0);            /* only one variable at index 0 */
    IBToIntegerI(IBDomV(d,globvar));       /* rounding of bounds to integers */
    if( IBEmptyI(IBDomV(d,globvar)) )
    {
       return( IBNarrowFailure );
    }
    return( IBNarrowSuccess );
  }

  /* Interval evaluations of left and right terms */
  IBTevalAll(IBCleft(c),d);
  if( IBEmptyI(IBTfwd(IBCleft(c))) )
  {
    return IBNarrowFailure;
  }

  IBTevalAll(IBCright(c),d);
  if( IBEmptyI(IBTfwd(IBCright(c))) )
  {
    return IBNarrowFailure;
  }

  /* Interpretation of relation symbol */
  if( (IBCrel(c)==IBRelationEQU) || (IBCrelfunc(c)==IBRelationSET) )
  {
    IBInterII(itv,IBTfwd(IBCleft(c)),IBTfwd(IBCright(c)));
    if( IBEmptyI(itv) ) return( IBNarrowFailure );

    IBCopyI(IBTbwd(IBCleft(c)),itv);
    IBCopyI(IBTbwd(IBCright(c)),itv);
  }
  else if( IBCrel(c)==IBRelationINF )
  {
    /* Left term */
    IBMinI(itv) = IBMinI(IBTfwd(IBCleft(c)));
    IBMaxI(itv) = IBMin(IBMaxI(IBTfwd(IBCleft(c))),IBMaxI(IBTfwd(IBCright(c))));
    if( IBEmptyI(itv) ) return( IBNarrowFailure );

    IBCopyI(IBTbwd(IBCleft(c)),itv);

    /* right term */
    IBMaxI(itv) = IBMaxI(IBTfwd(IBCright(c)));
    IBMinI(itv) = IBMax(IBMinI(IBTfwd(IBCleft(c))),IBMinI(IBTfwd(IBCright(c))));
    if( IBEmptyI(itv) ) return( IBNarrowFailure );

    IBCopyI(IBTbwd(IBCright(c)),itv);
  }
  else        /* IBRelationSUP */
  {
    /* left term */
    IBMaxI(itv) = IBMaxI(IBTfwd(IBCleft(c)));
    IBMinI(itv) = IBMax(IBMinI(IBTfwd(IBCleft(c))),IBMinI(IBTfwd(IBCright(c))));
    if( IBEmptyI(itv) ) return( IBNarrowFailure );

    IBCopyI(IBTbwd(IBCleft(c)),itv);


    /* right term */
    IBMinI(itv) = IBMinI(IBTfwd(IBCright(c)));
    IBMaxI(itv) = IBMin(IBMaxI(IBTfwd(IBCleft(c))),IBMaxI(IBTfwd(IBCright(c))));
    if( IBEmptyI(itv) ) return( IBNarrowFailure );

    IBCopyI(IBTbwd(IBCright(c)),itv);
  }

  /* Backward propagation */
  if( IBNarrowHC3reviseRec(IBCleft(c),d) )
      return( IBNarrowHC3reviseRec(IBCright(c),d) );
  else return( IBNarrowFailure );
}
