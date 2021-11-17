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
 * constant.c                                                               *
 ****************************************************************************/

#include "constant.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


IBConstants IBNewConstants()
/***************************************************************************
*  Allocation of an array of constants
*/
{
  IBConstants a;
  IBTree *f;
  IBItv itv, itv2;

  a    = (struct IBConst *)malloc(sizeof(struct IBConst));
  a->a = (struct IBCo *)malloc(IBConstAllocUnit*sizeof(struct IBCo));
  a->N = 0;
  a->Nfree = IBConstAllocUnit;

  /* Creation of mathematical constants */
  IBSetToPiI(itv);                /* itv=pi */
  f = IBTNewItv(itv);
  IBAddConstant(a,"@pi",f);       /* pi */

  IBMulRIinternal(itv2,2.0,itv);
  f = IBTNewItv(itv2);
  IBAddConstant(a,"@2_pi",f);     /* 2*pi */

  IBMulRIinternal(itv2,3.0,itv);
  f = IBTNewItv(itv2);
  IBAddConstant(a,"@3_pi",f);     /* 3*pi */

  IBMulRIinternal(itv2,4.0,itv);
  f = IBTNewItv(itv2);
  IBAddConstant(a,"@4_pi",f);     /* 4*pi */

  IBSetToHalfPiI(itv2);
  f = IBTNewItv(itv2);
  IBAddConstant(a,"@pi_2",f);     /* pi/2 */

  IBMulRIinternal(itv2,1.5,itv);
  f = IBTNewItv(itv2);
  IBAddConstant(a,"@3_pi_2",f);   /* 3*pi/2 */

  IBMulRIinternal(itv2,2.5,itv);
  f = IBTNewItv(itv2);
  IBAddConstant(a,"@5_pi_2",f);   /* 5*pi/2 */

  IBMulRIinternal(itv2,3.5,itv);
  f = IBTNewItv(itv2);
  IBAddConstant(a,"@7_pi_2",f);   /* 7*pi/2 */

  IBDivRposIinternal(itv2,1.0,itv);
  f = IBTNewItv(itv2);
  IBAddConstant(a,"@inv_pi",f);   /* 1/pi */

  IBSetToEI(itv);
  f = IBTNewItv(itv);
  IBAddConstant(a,"@e",f);        /* exp(1) */

  IBSetToLn2I(itv);
  f = IBTNewItv(itv);
  IBAddConstant(a,"@log2",f);      /* ln(2) */

  IBSetI(itv2,2.0,2.0);
  IBSqrtI(itv,itv2,itv2);
  f = IBTNewItv(itv);
  IBAddConstant(a,"@sqrt2",f);    /* sqrt(2) */

  return( a );
}


void IBFreeConstants(IBConstants a)
/***************************************************************************
*  Desallocation of an array of constants
*/
{
  int i;
  for( i=0; i<a->N; i++ )
  {
    free( a->a[i].name );   /* name of constant i */
  }
  free(a->a);
  free(a);
}


int IBAddConstant(IBConstants a, char *name, IBTree *f)
/***************************************************************************
*  To add constant name in a with value eval(f)
*  Desallocation of f
*  Returns 0 if f evaluates in the empty interval, 1 otherwise
*/
{
  if( a->Nfree==0 )
  {
    a->a = realloc(a->a,(a->N+IBConstAllocUnit)*sizeof(struct IBCo));
    a->Nfree = IBConstAllocUnit;
  }

  a->a[a->N].name = (char *)malloc((1+strlen(name))*sizeof(char));
  strcpy(a->a[a->N].name,name);
  IBTevalConstant(f);                     /* evaluation of constant */

  if( IBEmptyI(IBTfwd(f)) )
  {
    return( 0 );
  }

  IBCopyI(a->a[a->N].value,IBTfwd(f));
  IBTFree(f);                             /* desallocation of expression f */
  a->N++;
  a->Nfree--;
  return( 1 );
}


IBInterval *IBGetConstant(IBConstants a, char *name)
/***************************************************************************
*  To get the value of constant name
*  Returns NULL if name is not in a
*/
{
  int i;
  for( i=0; i<a->N; i++ )
  {
    if( strcmp(name,a->a[i].name)==0 ) return( a->a[i].value );
  }
  return( NULL );
}
