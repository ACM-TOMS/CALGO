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
 * union_interval.c                                                         *
 ****************************************************************************/

#include "union_interval.h"
#include <stdio.h>


extern IBItv IBConstantPi;
extern IBItv IBConstant_2_Pi;
extern IBItv IBConstant_Pi_2;

IBUnion *IBNewEmptyU()
/***************************************************************************
*  Allocation of a union of intervals of size 0
*/
{
  IBUnion *u;
  u = (IBUnion *)malloc(sizeof(IBUnion));
  IBUnionNf(u) = 0;
  IBUnionN(u) = 0;
  IBUnionI(u) = NULL;
  return( u );
}


IBUnion *IBNewCopyU(IBUnion *u)
/***************************************************************************
*  Allocation of a copy of u
*/
{
  IBUnion *copy;
  copy = (IBUnion *)malloc(sizeof(IBUnion));
  IBUnionNf(copy) = 0;
  IBUnionN(copy)  = IBUnionN(u);
  IBUnionI(copy) = (IBItv *)malloc(IBUnionN(u)*sizeof(IBItv));
  memcpy(IBUnionI(copy),IBUnionI(u),IBUnionN(u)*sizeof(IBItv));
  return( copy );
}


IBUnion *IBNewU(int N)
/***************************************************************************
*  Allocation of a union of intervals of size N
*/
{
  IBUnion *u;
  u = (IBUnion *)malloc(sizeof(IBUnion));
  IBUnionNf(u) = N;
  IBUnionN(u) = 0;
  IBUnionI(u) = (IBItv *)malloc(N*sizeof(IBItv));
  return( u );
}


void IBFreeU(IBUnion *u)
/***************************************************************************
*  Desallocation of a union of intervals
*/
{
  if( u!=NULL ) {
    if( IBUnionI(u)!=NULL ) free(IBUnionI(u));
    free(u);
  }
}


void IBResetU(IBUnion *u)
/***************************************************************************
*  The size of union u becomes 0
*/
{
  if( IBUnionN(u)>2 )
  {
     free(IBUnionI(u));
     IBUnionI(u) = (IBItv *)malloc(2*sizeof(IBItv));
  }
  IBUnionN(u)  = 0;
  IBUnionNf(u) = 2;
}


void IBReallocU(IBUnion *u)
/***************************************************************************
*  Reallocation of a union of intervals with IBUnionAllocUnit intervals more
*  u is supposed to be non empty
*/
{
  if( IBUnionI(u)==NULL )
  {
    IBUnionI(u) = (IBItv *)malloc(IBUnionAllocUnit*sizeof(IBItv));
    IBUnionNf(u) = IBUnionAllocUnit;
    IBUnionN(u) = 0;
  }
  else
  {
    IBUnionNf(u) += IBUnionAllocUnit;
    IBUnionI(u) = realloc(IBUnionI(u),(IBUnionN(u)+IBUnionNf(u))*sizeof(IBItv));
  }
}


void IBHullU(IBUnion *u, IBItv i)
/***************************************************************************
*  i := hull(u), u non empty
*/
{
  IBMinI(i) = IBMinI(IBUnionItv(u,0));
  IBMaxI(i) = IBMaxI(IBUnionItv(u,IBUnionN(u)-1));
}


int IBDoubleInU(IBUnion *u, double x)
/***************************************************************************
*  Returns 1 if x is in u
*/
{
  int k;
  for( k=0; k<IBUnionN(u); k++ )
  {
    if( IBDoubleInI(IBUnionItv(u,k),x) ) return 1;
  }
  return 0;
}


double IBPrevDoubleInU(IBUnion *u, double x)
/***************************************************************************
*  Returns the greatest number in u that is <= x
*  hypothesis 1 : x not in u
*  hypothesis 2 : such a number is supposed to exist, not (x < u)
*/
{
  int k;
  for( k=IBUnionN(u)-1; k>=0; k-- )
  {
    if( x>=IBMaxI(IBUnionItv(u,k)) ) {
      return( IBMaxI(IBUnionItv(u,k)) );
    }
  }
}


double IBNextDoubleInU(IBUnion *u, double x)
/***************************************************************************
*  Returns the smallest number in u that is >= x
*  hypothesis 1 : x not in u
*  hypothesis 2 : such a number is supposed to exist, not (x > u)
*/
{
  int k;
  for( k=0; k<IBUnionN(u); k++ )
  {
    if( x<=IBMinI(IBUnionItv(u,k)) ) {
      return( IBMinI(IBUnionItv(u,k)) );
    }
  }
}


void IBShiftRightU(IBUnion *u, int j, IBItv i)
/***************************************************************************
*  To insert interval i in u[j]
*    - shift to the right all intervals u[j],...,u[N-1]
*    - copy of interval i
*/
{
  int k;

  if( IBUnionNf(u)==0 ) IBReallocU(u);

  for( k=IBUnionN(u); k>j; k-- )
  {
    IBCopyI(IBUnionItv(u,k),IBUnionItv(u,k-1));
  }
  IBUnionN(u)++;
  IBUnionNf(u)--;
  IBCopyI(IBUnionItv(u,j),i);
}


void IBShiftLeftU(IBUnion *u, int j)
/***************************************************************************
*  to delete interval u[j]
*/
{
  if( j < IBUnionN(u)-1 )
    memcpy(&(IBUnionItv(u,j)),
           &(IBUnionItv(u,j+1)),
           (IBUnionN(u)-j-1)*sizeof(IBItv));
  IBUnionN(u)--;
  IBUnionNf(u)++;
}


void IBUnionIU(IBUnion *u, IBItv i)
/***************************************************************************
*  u := u union i
*  i is modified by this operation
*/
{
  int j;

  if( IBEmptyI(i) ) return;

  for( j=0 ;; )
  {
    if( j==IBUnionN(u) )                                 /* insertion in last position */
    {
      if( IBUnionNf(u)==0 )
      {
        IBReallocU(u);
      }
      IBCopyI(IBUnionItv(u,j),i);
      IBUnionN(u) ++;
      IBUnionNf(u)--;
      return;
    }
    else if( IBMaxI(i) < IBMinI(IBUnionItv(u,j)) )       /* insertion in position j */
    {
      IBShiftRightU(u,j,i);
      return;
    }
    else if( IBIncludedII(i,IBUnionItv(u,j)) )           /* i is included in u */
    {
      return;
    }
    else if( IBIncludedII(IBUnionItv(u,j),i) )           /* i contains u[j] => u[j] is deleted */
    {
      IBShiftLeftU(u,j);
    }
    else if( IBMinI(i) > IBMaxI(IBUnionItv(u,j)) )       /* i is after u[j] */
    {
     j++;
    }
    else if( (IBMinI(IBUnionItv(u,j)) >= IBMinI(i)) &&   /*    i : |--------|    */
             (IBMaxI(i) <= IBMaxI(IBUnionItv(u,j))) )    /* u[j] :    |--------| */
    {
      IBMinI(IBUnionItv(u,j)) = IBMinI(i);               /* u[j] : |-----------| */
      return;
    }
    else if( (IBMinI(IBUnionItv(u,j)) <= IBMinI(i)) &&   /*    i :    |--------| */
             (IBMaxI(i) >= IBMaxI(IBUnionItv(u,j))) )    /* u[j] : |--------|    */
    {
      IBMinI(i) = IBMinI(IBUnionItv(u,j));               /*    i : |-----------| */
      IBShiftLeftU(u,j);                                 /* u[j] is removed */
    }
  }
}


int IBInterIU(IBUnion *u, IBItv i)
/***************************************************************************
*  u := u inter i
*  i is modified by this operation
*  Returns 1 if the intersection is not empty
*/
{
  int j;

  for( j=0 ;; )
  {
    if( j==IBUnionN(u) )
    {
      return( IBUnionN(u)>0 );
    }
    else if( IBMaxI(i)<IBMinI(IBUnionItv(u,j)) )   /* i does not intersect u[j], u[j+1], ... */
    {
      IBUnionNf(u) += IBUnionN(u)-j;
      IBUnionN(u) = j;
      return( IBUnionN(u)>0 );
    }

    else if( IBMinI(i)>IBMaxI(IBUnionItv(u,j)) )   /* i and u[j] are disjoint */
    {
      IBShiftLeftU(u,j);
    }

    else if( IBIncludedII(i,IBUnionItv(u,j)) )     /* the intersection is i */
    {
      IBUnionNf(u) += IBUnionN(u)-1;
      IBUnionN(u) = 1;
      IBCopyI(IBUnionItv(u,0),i);
      return( 1 );
    }
    else if( IBIncludedII(IBUnionItv(u,j),i) )     /* u[j] is in the intersection */
    {
      j++;
    }
    else if( IBMaxI(i)<=IBMaxI(IBUnionItv(u,j)) )  /*    i : |--------|    */
    {                                              /* u[j] :    |--------| */
      IBMaxI(IBUnionItv(u,j)) = IBMaxI(i);         /* u[j] :    |-----|    */

      IBUnionNf(u) += IBUnionN(u)-j-1;
      IBUnionN(u) = j+1;
      return( 1 );
    }

    else if( IBMaxI(i)>=IBMaxI(IBUnionItv(u,j)) )  /*    i :    |--------| */
    {                                              /* u[j] : |--------|    */
      IBMinI(IBUnionItv(u,j)) = IBMinI(i);         /* u[j] :    |-----|    */
      j++;
    }
  }
}


int IBInterUU(IBUnion *u, IBUnion *v)
/***************************************************************************
*  u := u inter v
*  v is modified by this operation
*  Returns 1 if the intersection is not empty
*/
{
  int i, j;
  IBUnion *new, *save, *copy;

  save = IBNewCopyU(u);
  IBResetU(u);

  for( i=0; i<IBUnionN(v); i++ )
  {
     copy = IBNewCopyU(save);
     IBInterIU(copy,IBUnionItv(v,i));    /* copy := u inter v[i] */

     for( j=0; j<IBUnionN(copy); j++ )
     {
       IBUnionIU(u,IBUnionItv(copy,j));  /* save u inter v[i] */
     }
     IBFreeU(copy);
  }

  IBFreeU(save);
  return( IBUnionN(u)>0 );
}


void IBWriteU(FILE *out, IBUnion *u, int digits, int mode)
/***************************************************************************
*  To write u on out
*/
{
  int i;

  fprintf(out,"{");
  if( IBUnionN(u)==0 )
  {
    fprintf(out,"}");
    return;
  }

  IBWriteI(out,IBUnionItv(u,0),digits,mode);

  for( i=1; i<IBUnionN(u); i++ )
  {
    fprintf(out,",");
    IBWriteI(out,IBUnionItv(u,i),digits,mode);
  }
  fprintf(out,"}");
}


IBUnion * IBAddUU(IBUnion *u1, IBUnion *u2)
/***************************************************************************
*  Returns u1+u2
*/
{
  IBUnion *u = IBNewU(IBMin(IBUnionAllocUnit,IBUnionN(u1)*IBUnionN(u2)));
  int j, k;
  IBItv itv;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    for( k=0; k<IBUnionN(u2); k++ )
    {
      IBAddII(itv,IBUnionItv(u1,j),IBUnionItv(u2,k));
      IBUnionIU(u,itv);
    }
  }
  return( u );
}


IBUnion * IBSubUU(IBUnion *u1, IBUnion *u2)
/***************************************************************************
*  Returns u1-u2
*/
{
  IBUnion *u = IBNewU(IBMin(IBUnionAllocUnit,IBUnionN(u1)*IBUnionN(u2)));
  int j, k;
  IBItv itv;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    for( k=0; k<IBUnionN(u2); k++ )
    {
      IBSubII(itv,IBUnionItv(u1,j),IBUnionItv(u2,k));
      IBUnionIU(u,itv);
    }
  }
  return( u );
}


IBUnion * IBNegU(IBUnion *u1, IBUnion *useless)
/***************************************************************************
*  Returns -u1
*/
{
  IBUnion *u = IBNewU(IBUnionN(u1));
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    IBNegI(itv,IBUnionItv(u1,j),i1);
    IBCopyI(IBUnionItv(u,IBUnionN(u1)-j-1),itv);
  }
  IBUnionN(u) = IBUnionN(u1);
  IBUnionNf(u) = 0;
  return( u );
}


IBUnion * IBMulUU(IBUnion *u1, IBUnion *u2)
/***************************************************************************
*  Returns u1*u2
*/
{
  IBUnion *u = IBNewU(IBMin(IBUnionAllocUnit,IBUnionN(u1)*IBUnionN(u2)));
  int j, k;
  IBItv itv;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    for( k=0; k<IBUnionN(u2); k++ )
    {
      IBMulII(itv,IBUnionItv(u1,j),IBUnionItv(u2,k));
      IBUnionIU(u,itv);
    }
  }
  return( u );
}


IBUnion * IBDivUU(IBUnion *u1, IBUnion *u2)
/***************************************************************************
*  Returns u1/u2
*/
{
  IBUnion *u = IBNewU(IBMin(IBUnionAllocUnit,IBUnionN(u1)*IBUnionN(u2)));
  int j, k;
  IBItv itv;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    for( k=0; k<IBUnionN(u2); k++ )
    {
      IBDivII(itv,IBUnionItv(u1,j),IBUnionItv(u2,k));
      IBUnionIU(u,itv);
    }
  }
  return( u );
}


IBUnion * IBSqrU(IBUnion *u1, IBUnion *useless)
/***************************************************************************
*  Returns u1^2
*/
{
  IBUnion *u = IBNewU(IBUnionN(u1));
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    IBSqrI(itv,IBUnionItv(u1,j),i1);
    IBUnionIU(u,itv);
  }
  return( u );
}


IBUnion * IBSqrtU(IBUnion *u1, IBUnion *useless)
/***************************************************************************
*  Returns u1^2
*/
{
  IBUnion *u = IBNewU(IBUnionN(u1));
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    IBSqrtI(itv,IBUnionItv(u1,j),i1);
    IBUnionIU(u,itv);
  }
  return( u );
}

IBUnion *IBPowU(IBUnion *u1, IBUnion *n)
/***************************************************************************
*  Returns u1^n
*/
{
  IBUnion *u = IBNewU(IBUnionN(u1));
  int j;
  IBItv itv;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    IBPowI(itv,IBUnionItv(u1,j),IBUnionItv(n,0));
    IBUnionIU(u,itv);
  }
  return( u );
}

IBUnion * IBExpU(IBUnion *u1, IBUnion *useless)
/***************************************************************************
*  Returns exp(u1)
*/
{
  IBUnion *u = IBNewU(IBUnionN(u1));
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    IBExpI(itv,IBUnionItv(u1,j),i1);
    IBUnionIU(u,itv);
  }
  return( u );
}

IBUnion * IBLogU(IBUnion *u1, IBUnion *useless)
/***************************************************************************
*  Returns log(u1)
*/
{
  IBUnion *u = IBNewU(IBUnionN(u1));
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    IBLogI(itv,IBUnionItv(u1,j),i1);
    IBUnionIU(u,itv);
  }
  return( u );
}

IBUnion * IBDivRelOneII(IBItv i1, IBItv i2)
/***************************************************************************
*  Returns i1 divrel i2
*/
{
  IBUnion *u;
  IBItv j, k;

  if (IBExtDivII(j,k,i1,i2)==1) {
    u = IBNewU(1);
    IBCopyI(IBUnionItv(u,0),j);
    IBUnionN(u) = 1;
    IBUnionNf(u) = 0;
  }
  else {
    u = IBNewU(2);
    IBCopyI(IBUnionItv(u,0),j);
    IBCopyI(IBUnionItv(u,1),k);
    IBUnionN(u) = 2;
    IBUnionNf(u) = 0;
  }
  return( u );
}

IBUnion * IBSinRelOneI(IBItv i1, IBItv dom)
/***************************************************************************
*  Suppose that y = sin(x) <=> x = sin-1(y)
*  let dom be the domain of x and i1 be the domain of y
*  Returns the evaluation of sin-1 over i1, the result being contained in dom
*/
{
  /* Computation of sin-1 over i1 in [-pi,+pi] => u1 */
  IBUnion *u1, *u;
  IBItv j, k, l, offset;
  int nb;
  double leftB, rightB, nstep;

  /* Computation of sin-1 over i1 in [-pi,+pi] => u1 */
  if ((nb = IBSinRelI(j,k,l,i1)) == 0) {
    u1 = IBNewU(1);
    IBUnionN(u1) = 0;
    IBUnionNf(u1) = 1;   /* empty union */
    return u1;
  }
  else if (nb==1) {
    u1 = IBNewU(1);
    IBCopyI(IBUnionItv(u1,0),j);
    IBUnionN(u1) = 1;
    IBUnionNf(u1) = 0;
  }
  else if (nb==2) {
    u1 = IBNewU(2);
    IBCopyI(IBUnionItv(u1,0),j);
    IBCopyI(IBUnionItv(u1,1),k);
    IBUnionN(u1) = 2;
    IBUnionNf(u1) = 0;
  }
  else {   /* nb==3 */
    u1 = IBNewU(3);
    IBCopyI(IBUnionItv(u1,0),j);
    IBCopyI(IBUnionItv(u1,1),k);
    IBCopyI(IBUnionItv(u1,2),l);
    IBUnionN(u1) = 3;
    IBUnionNf(u1) = 0;
  }

  /* Note: only one interval is computed in the result = simplification of the problem
     of infinite unions of intervals if dom=[-oo,+oo] */
  u = IBNewU(1);
  IBUnionN(u) = 1;
  IBUnionNf(u) = 0;

  /* reduction of left bound */
  if( IBMinI(dom)<-IBMaxI(IBConstantPi) ) { /* i.inf < -pi */
    IBRoundDown();
    nstep = floor((IBMinI(dom) - IBMinI(IBConstantPi))/IBMaxI(IBConstant_2_Pi)) + 1.0;
    IBMulRIinternal(offset,nstep,IBConstant_2_Pi);
    IBRoundDown();
    leftB = IBMinI(dom) - IBMaxI(offset);

    if( leftB<(-IBMaxI(IBConstantPi)) ) leftB = -IBMaxI(IBConstantPi);
    else if( leftB>IBMaxI(IBConstantPi) ) leftB = IBMaxI(IBConstantPi);
  }
  else if( IBMinI(dom)>IBMaxI(IBConstantPi) ) {
    IBRoundUp();
    nstep = floor((IBMinI(dom) - IBMinI(IBConstantPi))/IBMinI(IBConstant_2_Pi)) + 1.0;
    IBMulRposIinternal(offset,nstep,IBConstant_2_Pi);
    IBRoundDown();
    leftB = IBMinI(dom) - IBMaxI(offset);

    if( leftB<(-IBMaxI(IBConstantPi)) ) leftB = -IBMaxI(IBConstantPi);
    else if( leftB>IBMaxI(IBConstantPi) ) leftB = IBMaxI(IBConstantPi);
  }
  else {
    nstep = 0;
    leftB = IBMinI(dom);
  }

  if( IBDoubleInU(u1,leftB) ) {  /* no reduction of left bound */
    IBMinI(IBUnionItv(u,0)) = IBMinI(dom);
  }
  else {
    double c;
    if( leftB<=IBMaxI(IBUnionItv(u1,IBUnionN(u1)-1)) ) {
      c = IBNextDoubleInU(u1,leftB);
    }
    else {
      c = IBMinI(IBUnionItv(u1,0));
      nstep += 1.0;
    }

    if( nstep==0.0 ) {
      IBMinI(IBUnionItv(u,0)) = c;
    }
    else {
      IBMulRposIinternal(offset,nstep,IBConstant_2_Pi);
      IBRoundDown();
      IBMinI(IBUnionItv(u,0)) = c + IBMaxI(offset);
    }
  }

  /* reduction of right bound */
  if( IBMaxI(dom)<-IBMaxI(IBConstantPi) ) { /* i.sup < -pi */
    IBRoundDown();
    nstep = floor((IBMaxI(dom) - IBMinI(IBConstantPi))/IBMaxI(IBConstant_2_Pi)) + 1.0;
    IBMulRIinternal(offset,nstep,IBConstant_2_Pi);
    IBRoundUp();
    rightB = IBMaxI(dom) - IBMinI(offset);

    if( rightB>IBMaxI(IBConstantPi) )    rightB = IBMaxI(IBConstantPi);
    else if( rightB<-(IBMaxI(IBConstantPi)) ) rightB = -IBMaxI(IBConstantPi);
  }
  else if( IBMaxI(dom)>IBMaxI(IBConstantPi) ) {
    IBRoundUp();
    nstep = floor((IBMaxI(dom) - IBMinI(IBConstantPi))/IBMinI(IBConstant_2_Pi)) + 1.0;
    IBMulRposIinternal(offset,nstep,IBConstant_2_Pi);
    IBRoundUp();
    rightB = IBMaxI(dom) - IBMinI(offset);

    if( rightB>IBMaxI(IBConstantPi) )    rightB = IBMaxI(IBConstantPi);
    else if( rightB<-(IBMaxI(IBConstantPi)) ) rightB = -IBMaxI(IBConstantPi);
  }
  else {
    nstep = 0;
    rightB = IBMaxI(dom);
  }

  if( IBDoubleInU(u1,rightB) ) {  /* no reduction of right bound */
    IBMaxI(IBUnionItv(u,0)) = IBMaxI(dom);
  }
  else {
    double c;
    if( rightB>=IBMinI(IBUnionItv(u1,0)) ) {
      c = IBPrevDoubleInU(u1,rightB);
    }
    else {
      c = IBMaxI(IBUnionItv(u1,IBUnionN(u1)-1));
      nstep -= 1.0;
    }

    if( nstep==0.0 ) {
      IBMaxI(IBUnionItv(u,0)) = c;
    }
    else {
      IBMulRposIinternal(offset,nstep,IBConstant_2_Pi);
      IBRoundUp();
      IBMaxI(IBUnionItv(u,0)) = c + IBMaxI(offset);
    }
  }

  if( IBEmptyI(IBUnionItv(u,0)) ) {
    IBUnionN(u) = 0;
    IBUnionNf(u) = 1;
  }

  IBFreeU(u1);
  return( u );
}

IBUnion * IBTanRelOneI(IBItv i1, IBItv dom)
/***************************************************************************
*  Suppose that y = tan(x) <=> x = tan-1(y)
*  let dom be the domain of x and i1 be the domain of y
*  Returns the evaluation of tan-1 over i1, the result being contained in dom
*/
{
  IBUnion *u;
  IBItv atan_i1, offset;
  double leftB, rightB, nstep;

  /* inverse of i1 in the interval [-pi/2,pi/2] */
  IBAtanI(atan_i1,i1,i1);

  /* Note: only one interval is computed in the result = simplification of the problem
     of infinite unions of intervals if dom=[-oo,+oo] */
  u = IBNewU(1);
  IBUnionN(u) = 1;
  IBUnionNf(u) = 0;

  /* reduction of left bound */
  if( IBMinI(dom)<-IBMaxI(IBConstant_Pi_2) ) { /* i.inf < -pi/2 */
    IBRoundDown();
    nstep = floor((IBMinI(dom) - IBMinI(IBConstant_Pi_2))/IBMaxI(IBConstantPi)) + 1.0;
    IBMulRIinternal(offset,nstep,IBConstantPi);
    IBRoundDown();
    leftB = IBMinI(dom) - IBMaxI(offset);

    if( leftB<(-IBMaxI(IBConstant_Pi_2)) ) leftB = -IBMaxI(IBConstant_Pi_2);
    else if( leftB>IBMaxI(IBConstant_Pi_2) ) leftB = IBMaxI(IBConstant_Pi_2);
  }
  else if( IBMinI(dom)>IBMaxI(IBConstant_Pi_2) ) {
    IBRoundUp();
    nstep = floor((IBMinI(dom) - IBMinI(IBConstant_Pi_2))/IBMinI(IBConstantPi)) + 1.0;
    IBMulRposIinternal(offset,nstep,IBConstantPi);
    IBRoundDown();
    leftB = IBMinI(dom) - IBMaxI(offset);

    if( leftB<(-IBMaxI(IBConstant_Pi_2)) ) leftB = -IBMaxI(IBConstant_Pi_2);
    else if( leftB>IBMaxI(IBConstant_Pi_2) ) leftB = IBMaxI(IBConstant_Pi_2);
  }
  else {
    nstep = 0;
    leftB = IBMinI(dom);
  }

  if( IBDoubleInI(atan_i1,leftB) ) {  /* no reduction of left bound */
    IBMinI(IBUnionItv(u,0)) = IBMinI(dom);
  }
  else {
    if( leftB>IBMinI(atan_i1) ) {
      nstep += 1.0;
    }
    if (nstep==0.0) {
      IBMinI(IBUnionItv(u,0)) = IBMinI(atan_i1);
    }
    else {
      IBMulRposIinternal(offset,nstep,IBConstantPi);
      IBRoundDown();
      IBMinI(IBUnionItv(u,0)) = IBMinI(atan_i1) + IBMaxI(offset);
    }
  }

  /* reduction of right bound */
  if( IBMaxI(dom)<-IBMaxI(IBConstant_Pi_2) ) { /* i.sup < -pi/2 */
    IBRoundDown();
    nstep = floor((IBMaxI(dom) - IBMinI(IBConstant_Pi_2))/IBMaxI(IBConstantPi)) + 1.0;
    IBMulRIinternal(offset,nstep,IBConstantPi);
    IBRoundUp();
    rightB = IBMaxI(dom) - IBMinI(offset);

    if( rightB<(-IBMaxI(IBConstant_Pi_2)) ) rightB = -IBMaxI(IBConstant_Pi_2);
    else if( rightB>IBMaxI(IBConstant_Pi_2) ) rightB = IBMaxI(IBConstant_Pi_2);
  }
  else if( IBMaxI(dom)>IBMaxI(IBConstant_Pi_2) ) {
    IBRoundUp();
    nstep = floor((IBMaxI(dom) - IBMinI(IBConstant_Pi_2))/IBMinI(IBConstantPi)) + 1.0;
    IBMulRposIinternal(offset,nstep,IBConstantPi);
    IBRoundUp();
    rightB = IBMaxI(dom) - IBMinI(offset);

    if( rightB<(-IBMaxI(IBConstant_Pi_2)) ) rightB = -IBMaxI(IBConstant_Pi_2);
    else if( rightB>IBMaxI(IBConstant_Pi_2) ) rightB = IBMaxI(IBConstant_Pi_2);
  }
  else {
    nstep = 0;
    rightB = IBMaxI(dom);
  }

  if( IBDoubleInI(atan_i1,rightB) ) {  /* no reduction of right bound */
    IBMaxI(IBUnionItv(u,0)) = IBMaxI(dom);
  }
  else {
    if( rightB<IBMinI(atan_i1) ) {
      nstep -= 1.0;
    }
    if (nstep==0.0) {
      IBMaxI(IBUnionItv(u,0)) = IBMaxI(atan_i1);
    }
    else {
      IBMulRposIinternal(offset,nstep,IBConstantPi);
      IBRoundUp();
      IBMaxI(IBUnionItv(u,0)) = IBMaxI(atan_i1) + IBMaxI(offset);
    }
  }

  if( IBEmptyI(IBUnionItv(u,0)) ) {
    IBUnionN(u) = 0;
    IBUnionNf(u) = 1;
  }

  return( u );
}

IBUnion * IBDivRelUU(IBUnion *u1, IBUnion *u2)
/***************************************************************************
*  Returns u1 divrel u2
*/
{
  IBUnion *u = IBNewU(IBMin(IBUnionAllocUnit,IBUnionN(u1)*IBUnionN(u2)));
  IBUnion *u3;
  int i, j, k;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    for( k=0; k<IBUnionN(u2); k++ )
    {
      u3 = IBDivRelOneII(IBUnionItv(u1,j),IBUnionItv(u2,k));
      for( i=0; i<IBUnionN(u3); i++ )
      {
        IBUnionIU(u,IBUnionItv(u3,i));
      }
      IBFreeU(u3);
    }
  }
  return( u );
}


IBUnion *IBNthRootOneI(IBItv i1, IBItv n)
/***************************************************************************
*  Returns the n-th root of i1
*/
{
  IBUnion *u;
  IBItv j, k;
  int r;

  if ((r=IBNthRootRelI(j,k,i1,n))==1) {
    u = IBNewU(1);
    IBUnionN(u) = 1;
    IBUnionNf(u) = 0;
    IBCopyI(IBUnionItv(u,0),j);
  }
  else if (r==2) {
    u = IBNewU(2);
    IBUnionN(u) = 2;
    IBUnionNf(u) = 0;
    IBCopyI(IBUnionItv(u,0),j);
    IBCopyI(IBUnionItv(u,1),k);
  }
  else {
    u = IBNewU(1);    /* u is empty, this case must not happen */
    IBUnionN(u) = 0;
    IBUnionNf(u) = 1;
  }
  return( u );
}


IBUnion *IBCoshRelOneI(IBItv i1)
/***************************************************************************
*  Returns the relational hyperbolic cosine of i1
*/
{
  IBUnion *u;
  IBItv j, k;
  int r;

  if ((r=IBCoshRelI(j,k,i1))==1) {
    u = IBNewU(1);
    IBUnionN(u) = 1;
    IBUnionNf(u) = 0;
    IBCopyI(IBUnionItv(u,0),j);
  }
  else if (r==2) {
    u = IBNewU(2);
    IBUnionN(u) = 2;
    IBUnionNf(u) = 0;
    IBCopyI(IBUnionItv(u,0),j);
    IBCopyI(IBUnionItv(u,1),k);
  }
  else {
    u = IBNewU(1);    /* u is empty, this case must not happen */
    IBUnionN(u) = 0;
    IBUnionNf(u) = 1;
  }
  return( u );
}


IBUnion *IBAcosRelOneI(IBItv i1)
/***************************************************************************
*  Returns the relational arc-cosine of i1
*/
{
  IBItv i2;
  IBUnion *u;
  u = IBNewU(1);
  IBUnionN(u) = 1;
  IBUnionNf(u) = 0;

  if( (IBMinI(i1)>=IBMaxI(IBConstantPi)) || (IBMaxI(i1)<0.0) ) {
    IBUnionN(u) = 0;
    IBUnionNf(u) = 1;   /* i1 inter [0,pi] = empty => empty result */
  }
  else {
    if (IBMinI(i1)<0.0) {
      IBMinI(i2) = 0.0;
    }
    else {
      IBMinI(i2) = IBMinI(i1);
    }
    if (IBMaxI(i1)>=IBMaxI(IBConstantPi)) {
      IBMaxI(i2) = IBMaxI(IBConstantPi);
    }
    else {
      IBMaxI(i2) = IBMaxI(i1);
    }
    IBCosI(IBUnionItv(u,0),i2,i2);
  }
  return( u );
}


IBUnion *IBAsinRelOneI(IBItv i1)
/***************************************************************************
*  Returns the relational arc-sine of i1
*/
{
  IBItv i2;
  IBUnion *u;
  u = IBNewU(1);
  IBUnionN(u) = 1;
  IBUnionNf(u) = 0;

  if( (IBMinI(i1)>=IBMaxI(IBConstant_Pi_2)) || (IBMaxI(i1)<=-IBMaxI(IBConstant_Pi_2)) ) {
    IBUnionN(u) = 0;
    IBUnionNf(u) = 1;   /* i1 inter [-pi/2,pi/2] = empty => empty result */
  }
  else {
    if (IBMinI(i1)<-IBMaxI(IBConstant_Pi_2)) {
      IBMinI(i2) = -IBMaxI(IBConstant_Pi_2);
    }
    else {
      IBMinI(i2) = IBMinI(i1);
    }
    if (IBMaxI(i1)>=IBMaxI(IBConstant_Pi_2)) {
      IBMaxI(i2) = IBMaxI(IBConstant_Pi_2);
    }
    else {
      IBMaxI(i2) = IBMaxI(i1);
    }
    IBSinI(IBUnionItv(u,0),i2,i2);
  }
  return( u );
}


IBUnion *IBAtanRelOneI(IBItv i1)
/***************************************************************************
*  Returns the relational arc-tangent of i1
*/
{
  IBItv i2;
  IBUnion *u;
  u = IBNewU(1);
  IBUnionN(u) = 1;
  IBUnionNf(u) = 0;

  if( (IBMinI(i1)>=IBMaxI(IBConstant_Pi_2)) || (IBMaxI(i1)<=-IBMaxI(IBConstant_Pi_2)) ) {
    IBUnionN(u) = 0;
    IBUnionNf(u) = 1;   /* i1 inter [-pi/2,pi/2] = empty => empty result */
  }
  else {
    if (IBMinI(i1)<-IBMaxI(IBConstant_Pi_2)) {
      IBMinI(i2) = -IBMaxI(IBConstant_Pi_2);
    }
    else {
      IBMinI(i2) = IBMinI(i1);
    }
    if (IBMaxI(i1)>=IBMaxI(IBConstant_Pi_2)) {
      IBMaxI(i2) = IBMaxI(IBConstant_Pi_2);
    }
    else {
      IBMaxI(i2) = IBMaxI(i1);
    }
    IBTanI(IBUnionItv(u,0),i2,i2);
  }
  return( u );
}


IBUnion *IBNthRootRelU(IBUnion *u1, IBUnion *n)
/***************************************************************************
*  Returns the n-th root of u1
*/
{
  IBUnion *u = IBNewU(IBUnionN(u1));
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u1); j++ )
  {
    if( (IBEven((int)(IBMinI(IBUnionItv(n,0))))) &&
        (IBMaxI(IBUnionItv(u1,j))<0.0) )
    {
    }
    else
    {
      u3 = IBNthRootOneI(IBUnionItv(u1,j),IBUnionItv(n,0));
      for( i=0; i<IBUnionN(u3); i++ ) IBUnionIU(u,IBUnionItv(u3,i));
      IBFreeU(u3);
    }
  }
  if( IBUnionN(u)==0 )  /* this case must not happen */
  {
    IBUnionN(u) = 1;
    IBUnionNf(u) = 0;
    IBMinI(IBUnionItv(u,0)) = IBNegInfinity;
    IBMaxI(IBUnionItv(u,0)) = IBPosInfinity;
  }
  return( u );
}


void IBhc4AddUI(IBUnion *result, IBUnion *u, IBItv itv)
/***************************************************************************
*  result := u + itv
*/
{
  int j;
  IBItv i2;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBAddII(i2,IBUnionItv(u,j),itv);
    IBUnionIU(result,i2);
  }
}


void IBhc4SubUI(IBUnion *result, IBUnion *u, IBItv itv)
/***************************************************************************
*  result := u - itv
*/
{
  int j;
  IBItv i2;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBSubII(i2,IBUnionItv(u,j),itv);
    IBUnionIU(result,i2);
  }
}


void IBhc4SubIU(IBUnion *result, IBItv itv, IBUnion *u)
/***************************************************************************
*  result := itv - u
*/
{
  int j;
  IBItv i2;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBSubII(i2,itv,IBUnionItv(u,j));
    IBUnionIU(result,i2);
  }
}


void IBhc4NegU(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := -u
*/
{
  int j;
  IBItv i2;
  IBInterval *i1;   /* useless */

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBNegI(i2,IBUnionItv(u,j),i1);
    IBUnionIU(result,i2);
  }
}


void IBhc4MulUI(IBUnion *result, IBUnion *u, IBItv itv)
/***************************************************************************
*  result := u * itv
*/
{
  int j;
  IBItv i2;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBMulII(i2,IBUnionItv(u,j),itv);
    IBUnionIU(result,i2);
  }
}


void IBhc4DivRelUI(IBUnion *result, IBUnion *u, IBItv itv)
/***************************************************************************
*  result := u1 divrel itv
*/
{
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u); j++ )
  {
    u3 = IBDivRelOneII(IBUnionItv(u,j),itv);
    for( i=0; i<IBUnionN(u3); i++ )
    {
      IBUnionIU(result,IBUnionItv(u3,i));
    }
    IBFreeU(u3);
  }
}


void IBhc4DivRelIU(IBUnion *result, IBItv itv, IBUnion *u)
/***************************************************************************
*  result := itv divrel u1
*/
{
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u); j++ )
  {
    u3 = IBDivRelOneII(itv,IBUnionItv(u,j));
    for( i=0; i<IBUnionN(u3); i++ )
    {
      IBUnionIU(result,IBUnionItv(u3,i));
    }
    IBFreeU(u3);
  }
}


void IBhc4NthRootRelUI(IBUnion *result, IBUnion *u, IBItv itv)
/***************************************************************************
*  result := n-th root of u, itv=[n,n]
*/
{
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u); j++ )
  {
    if( (IBEven((int)(IBMinI(itv)))) &&
        (IBMaxI(IBUnionItv(u,j))<0.0) )
    {
      /* nothing to add in result */
    }
    else
    {
      u3 = IBNthRootOneI(IBUnionItv(u,j),itv);

      for( i=0; i<IBUnionN(u3); i++ )
      {
        IBUnionIU(result,IBUnionItv(u3,i));
      }
      IBFreeU(u3);
    }
  }

  if( IBUnionN(result)==0 )  /* this case must not happen */
  {
    if( IBUnionI(result)==NULL )
    {
      IBUnionI(result) = (IBItv *)malloc(sizeof(IBItv));
    }
    IBUnionN(result) = 1;
    IBUnionNf(result) = 0;
    IBMinI(IBUnionItv(result,0)) = IBNegInfinity;
    IBMaxI(IBUnionItv(result,0)) = IBPosInfinity;
  }
}


void IBhc4CoshRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := relational hyperbolic cosine of u
*/
{
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u); j++ )
  {
    u3 = IBCoshRelOneI(IBUnionItv(u,j));

    for( i=0; i<IBUnionN(u3); i++ )
    {
      IBUnionIU(result,IBUnionItv(u3,i));
    }
    IBFreeU(u3);
  }

  if( IBUnionN(result)==0 )  /* this case must not happen */
  {
    if( IBUnionI(result)==NULL )
    {
      IBUnionI(result) = (IBItv *)malloc(sizeof(IBItv));
    }
    IBUnionN(result) = 1;
    IBUnionNf(result) = 0;
    IBMinI(IBUnionItv(result,0)) = IBNegInfinity;
    IBMaxI(IBUnionItv(result,0)) = IBPosInfinity;
  }
}


void IBhc4SqrUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := u^2
*/
{
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBSqrI(itv,IBUnionItv(u,j),i1);
    IBUnionIU(result,itv);
  }
}


void IBhc4LogUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := log(u)
*/
{
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBLogI(itv,IBUnionItv(u,j),i1);
    IBUnionIU(result,itv);
  }
}


void IBhc4ExpUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := exp(u)
*/
{
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBExpI(itv,IBUnionItv(u,j),i1);
    IBUnionIU(result,itv);
  }
}


void IBhc4MinUI(IBUnion *result, IBUnion *u, IBItv itv1, IBItv itv2)
/***************************************************************************
*  result := min-1
*/
{
  int j;
  IBItv i2;

  for( j=0; j<IBUnionN(u); j++ )
  {
    if( IBDisjointII(IBUnionItv(u,j),itv2) )
    {
      IBInterII(i2,itv1,IBUnionItv(u,j));
    }
    else
    {
      IBCopyI(i2,itv1);
    }
    if( IBMinI(IBUnionItv(u,j))>IBMinI(i2) )
    {
      IBMinI(i2) = IBMinI(IBUnionItv(u,j));
    }
    IBUnionIU(result,i2);
  }
}


void IBhc4MaxUI(IBUnion *result, IBUnion *u, IBItv itv1, IBItv itv2)
/***************************************************************************
*  result := max-1
*/
{
  int j;
  IBItv i2;

  for( j=0; j<IBUnionN(u); j++ )
  {
    if( IBDisjointII(IBUnionItv(u,j),itv2) )
    {
      IBInterII(i2,itv1,IBUnionItv(u,j));
    }
    else
    {
      IBCopyI(i2,itv1);
    }
    if( IBMaxI(IBUnionItv(u,j))<IBMaxI(i2) )
    {
      IBMaxI(i2) = IBMaxI(IBUnionItv(u,j));
    }
    IBUnionIU(result,i2);
  }
}

void IBhc4SinhRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := sinh-1 = asinh
*/
{
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBAsinhI(itv,IBUnionItv(u,j),i1);
    IBUnionIU(result,itv);
  }
}


void IBhc4AsinhRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := asinh-1 = sinh
*/
{
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBSinhI(itv,IBUnionItv(u,j),i1);
    IBUnionIU(result,itv);
  }
}

void IBhc4AcoshRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := acosh-1 = cosh
*/
{
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBCoshI(itv,IBUnionItv(u,j),i1);
    IBUnionIU(result,itv);
  }
}


void IBhc4TanhRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := tanh-1 = atanh
*/
{
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBAtanhI(itv,IBUnionItv(u,j),i1);
    IBUnionIU(result,itv);
  }
}


void IBhc4AtanhRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := atanh-1 = tanh
*/
{
  int j;
  IBItv itv;
  IBInterval *i1;

  for( j=0; j<IBUnionN(u); j++ )
  {
    IBTanhI(itv,IBUnionItv(u,j),i1);
    IBUnionIU(result,itv);
  }
}


void IBhc4CosRelUI(IBUnion *result, IBUnion *u, IBItv dom)
/***************************************************************************
*  result := cos-1(u)
*/
{
  IBItv copydom, i;
  int j;

  /* solving cos(x) = u, x in dom is equivalent to solving
        1) sin(y) = u, y in dom+pi/2
        2) x := y - pi/2 */

  IBAddII(copydom,dom,IBConstant_Pi_2);  /* copydom := dom + pi/2 */
  IBhc4SinRelUI(result,u,copydom);

  /* result := result - pi/2 */
  for( j=0; j<IBUnionN(result); j++ )
  {
    IBCopyI(i,IBUnionItv(result,j));
    IBSubII(IBUnionItv(result,j),i,IBConstant_Pi_2);
  }
}


void IBhc4SinRelUI(IBUnion *result, IBUnion *u, IBItv dom)
/***************************************************************************
*  result := sin-1(u) = cos-1(u + pi/2)
*/
{
  IBUnion *u3;
  int i, j;
  for( j=0; j<IBUnionN(u); j++ )
  {
    u3 = IBSinRelOneI(IBUnionItv(u,j),dom);
    for( i=0; i<IBUnionN(u3); i++ )
    {
      IBUnionIU(result,IBUnionItv(u3,i));
    }
    IBFreeU(u3);
  }
}


void IBhc4TanRelUI(IBUnion *result, IBUnion *u, IBItv dom)
/***************************************************************************
*  result := tan-1(u)
*/
{
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u); j++ )
  {
    u3 = IBTanRelOneI(IBUnionItv(u,j),dom);
    for( i=0; i<IBUnionN(u3); i++ )
    {
      IBUnionIU(result,IBUnionItv(u3,i));
    }
    IBFreeU(u3);
  }
}


void IBhc4AcosRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := tan-1(u)
*/
{
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u); j++ )
  {
    u3 = IBAcosRelOneI(IBUnionItv(u,j));
    for( i=0; i<IBUnionN(u3); i++ )
    {
      IBUnionIU(result,IBUnionItv(u3,i));
    }
    IBFreeU(u3);
  }
}


void IBhc4AsinRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := tan-1(u)
*/
{
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u); j++ )
  {
    u3 = IBAsinRelOneI(IBUnionItv(u,j));
    for( i=0; i<IBUnionN(u3); i++ )
    {
      IBUnionIU(result,IBUnionItv(u3,i));
    }
    IBFreeU(u3);
  }
}


void IBhc4AtanRelUI(IBUnion *result, IBUnion *u)
/***************************************************************************
*  result := tan-1(u)
*/
{
  IBUnion *u3;
  int i, j;

  for( j=0; j<IBUnionN(u); j++ )
  {
    u3 = IBAtanRelOneI(IBUnionItv(u,j));
    for( i=0; i<IBUnionN(u3); i++ )
    {
      IBUnionIU(result,IBUnionItv(u3,i));
    }
    IBFreeU(u3);
  }
}
