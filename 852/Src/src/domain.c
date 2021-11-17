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
 * domain.c                                                                 *
 ****************************************************************************/

#include <memory.h>
#include <stdio.h>
#include "domain.h"


inline IBDomains IBNewD(int ndom)
/***************************************************************************
*  Allocation of an array of ndom variables' domains
*/
{
  return( (IBDomains)malloc(ndom*sizeof(IBDom)) );
}


inline void IBFreeD(IBDomains d)
/***************************************************************************
*  Desallocation of d
*/
{
  free(d);
}


inline void IBCopyD(IBDomains dcopy, IBDomains d, int ndom)
/***************************************************************************
*  dcopy <- d
*/
{
  memcpy(dcopy,d,ndom*sizeof(IBDom));
}


void IBSetDtoR(IBDomains d, double x, int ndom)
/***************************************************************************
*  all the domains in d become [x,x]
*/
{
  int i;
  for( i=0; i<ndom; i++ )
  {
    IBSetI(d[i],x,x);
  }
}


int IBEmptyD(IBDomains d, int ndom)
/***************************************************************************
*  returns 1 if d is empty
*/
{
  int i;
  for( i=0; i<ndom; i++ )
  {
    if( IBEmptyI(d[i]) ) return( 1 );
  }
  return( 0 );
}


inline IBDomains IBNewCopyD(IBDomains d, int ndom)
/***************************************************************************
*  Return a copy of d
*/
{
  return( memcpy((IBDomains)malloc(ndom*sizeof(IBDom)),
                 d,
                 ndom*sizeof(IBDom)) );
}


inline void IBSetD(IBDomains d, int var, double x1, double x2)
/***************************************************************************
*  d[var] <- i
*/
{
  IBSetI(d[var],x1,x2);
}

void IBPrintD(IBDomains d, int ndom, int digits)
/***************************************************************************
*  To print d on stdout (profiling function)
*/
{
  int i;

  for( i=0; i<ndom; i++ )
  {
    IBWriteI(stdout,IBDomV(d,i),digits,IBPrintIntervalBounds);
    printf(" ");
  }
}
