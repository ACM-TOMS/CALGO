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
 * evaluator.c                                                              *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "evaluator.h"


extern IBOperations operations;  /* global array of operations */


void IBBwdAddII(IBTree *f)
/***************************************************************************
*  Derivation of operation AddII
*/
{
  IBCopyI(IBTbwd(IBTleft(f)),IBTbwd(f));
  IBCopyI(IBTbwd(IBTright(f)),IBTbwd(f));
}


void IBBwdAddRI(IBTree *f)
/***************************************************************************
*  Derivation of operation AddRI                             d(u+v)/d(v) = 1
*/
{
  IBCopyI(IBTbwd(IBTright(f)),IBTbwd(f));
}


void IBBwdSubII(IBTree *f)
/***************************************************************************
*  Derivation of operation SubII
*/
{
  IBItv useless;
  IBCopyI(IBTbwd(IBTleft(f)),IBTbwd(f));
  IBNegI(IBTbwd(IBTright(f)),IBTbwd(f),useless);
}


void IBBwdSubIR(IBTree *f)
/***************************************************************************
*  Derivation of operation SubIR                             d(u-v)/d(u) = 1
*/
{
  IBCopyI(IBTbwd(IBTleft(f)),IBTbwd(f));
}


void IBBwdSubRI(IBTree *f)
/***************************************************************************
*  Derivation of operation SubRI                            d(u*v)/d(v) = -1 
*/
{
  IBItv useless;
  IBNegI(IBTbwd(IBTright(f)),IBTbwd(f),useless);
}


void IBBwdNegI(IBTree *f)
/***************************************************************************
*  Derivation of operation NegI                              d(-u)/d(u) = -1
*/
{
  IBItv useless;
  IBNegI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
}


void IBBwdMulII(IBTree *f)
/***************************************************************************
*  Derivation of operation MulII
*/
{
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),IBTfwd(IBTright(f)));
  IBMulII(IBTbwd(IBTright(f)),IBTbwd(f),IBTfwd(IBTleft(f)));
}


void IBBwdMulRI(IBTree *f)
/***************************************************************************
*  Derivation of operation MulRI                             d(u*v)/d(v) = u
*/
{
  IBMulRI(IBTbwd(IBTright(f)),IBTfwd(IBTleft(f)),IBTbwd(f));
}



void IBBwdDivII(IBTree *f)
/***************************************************************************
*  Derivation of operation DivII
*/
{
  IBItv i1, i2, i3;

  IBNegI(i1,IBTfwd(IBTleft(f)),i1);                /* i1 = -u */
  IBSqrI(i2,IBTfwd(IBTright(f)),i2);               /* i2 = v^2 */
  IBDivII(i3,i1,i2);                               /* i3 = -u/v^2 */
  IBMulII(IBTbwd(IBTright(f)),IBTbwd(f),i3);
  IBDivRposIinternal(i1,1.0,IBTfwd(IBTright(f)));  /* i1 = 1/v */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdDivIR(IBTree *f)
/***************************************************************************
*  Derivation of operation DivIR                           d(u/v)/d(u) = 1/v
*/
{
  IBItv i1;

  IBDivRposIinternal(i1,1.0,IBTfwd(IBTright(f)));  /* i1 = 1/v */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdDivRI(IBTree *f)
/***************************************************************************
*  Derivation of operation DivRI                        d(u/v)/d(v) = -u/v^2
*/
{
  IBItv i1, i2, i3;

  IBNegI(i1,IBTfwd(IBTleft(f)),i1);     /* i1 = -u */
  IBSqrI(i2,IBTfwd(IBTright(f)),i2);    /* i2 = v^2 */
  IBDivII(i3,i1,i2);                    /* i3 = -u/v^2 */
  IBMulII(IBTbwd(IBTright(f)),IBTbwd(f),i3);
}


void IBBwdSqrI(IBTree *f)
/***************************************************************************
*  Derivation of operation SqrI                                 d(u^2) = 2*u
*/
{
  IBItv i1;
  IBMulRI(i1,IBTitv(IBTright(f)),IBTfwd(IBTleft(f)));  /* i1 = 2*u */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdPowI(IBTree *f)
/***************************************************************************
*  Derivation of operation PowI                            d(u^n) = n*u^(n-1)
*/
{
  IBItv i1, i2;
  int n = IBMinI(IBTitv(IBTright(f)));

  IBPowIinternal(i1,IBTfwd(IBTleft(f)),n-1);  /* i1 = u^(n-1) */
  IBMulRposIinternal(i2,n,i1);                /* i2 = n*u^(n-1) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i2);
}


void IBBwdSqrtI(IBTree *f)
/***************************************************************************
*  Derivation of operation SqrtI                    d(sqrt(u)) = 0.5/sqrt(u)
*/
{
  IBItv i1;

  IBDivRposIinternal(i1,0.5,IBTfwd(f));  /* i1 = 0.5/sqrt(u) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdLogI(IBTree *f)
/***************************************************************************
*  Derivation of operation LogI                              d(log(u)) = 1/u
*/
{
  IBItv i1;
  IBDivII(IBTbwd(IBTleft(f)),IBTbwd(f),IBTfwd(IBTleft(f)));
}


void IBBwdExpI(IBTree *f)
/***************************************************************************
*  Derivation of operation ExpI                           d(exp(u)) = exp(u)
*/
{
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),IBTfwd(f));
}  


void IBBwdMinII(IBTree *f)
/***************************************************************************
*  Operation MinII is not derivable
*  This function must not be called
*/
{
  IBToLargestI(IBTbwd(IBTleft(f)));
  IBToLargestI(IBTbwd(IBTright(f)));
}


void IBBwdMaxII(IBTree *f)
/***************************************************************************
*  Operation MaxII is not derivable
*  This function must not be called
*/
{
  IBToLargestI(IBTbwd(IBTleft(f)));
  IBToLargestI(IBTbwd(IBTright(f)));
}


void IBBwdCosI(IBTree *f)
/***************************************************************************
*  Derivation of operation CosI                          d(cos(u)) = -sin(u)
*/
{
  IBItv i1, i2, useless;
  IBSinI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 = sin(u) */
  IBNegI(i2,i1,useless);                    /* i2 = -sin(u) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i2);
}  


void IBBwdSinI(IBTree *f)
/***************************************************************************
*  Derivation of operation SinI                           d(sin(u)) = cos(u)
*/
{
  IBItv i1, useless;
  IBCosI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := cos(u) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdTanI(IBTree *f)
/***************************************************************************
*  Derivation of operation TanI                     d(tan(u)) = 1 + tan^2(u)
*/
{
  IBItv i1, i2, one, useless;
  IBSqrI(i1,IBTfwd(f),useless);    /* i1 = tan^2(u) */
  IBSetI(one,1.0,1.0);
  IBAddII(i2,one,i1);              /* i2 = 1 + tan^2(u) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i2);
}


void IBBwdCoshI(IBTree *f)
/***************************************************************************
*  Derivation of operation CoshI                        d(cosh(u)) = sinh(u)
*/
{
  IBItv i1, useless;
  IBSinhI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := sinh(u) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdSinhI(IBTree *f)
/***************************************************************************
*  Derivation of operation SinhI                        d(sinh(u)) = cosh(u)
*/
{
  IBItv i1, useless;
  IBCoshI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := cosh(u) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdTanhI(IBTree *f)
/***************************************************************************
*  Derivation of operation TanhI                    d(tanh(u)) = 1-tanh^2(u)
*/
{
  IBItv i1, i2, one, useless;
  IBSqrI(i1,IBTfwd(f),useless);    /* i1 = tanh^2(u) */
  IBSetI(one,1.0,1.0);
  IBSubII(i2,one,i1);              /* i2 = 1 - tanh^2(u) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i2);
}


void IBBwdAcosI(IBTree *f)
/***************************************************************************
*  Derivation of operation AcosI            d(acos(u)) = - 1 / sqrt(1 - u^2)
*/
{
  IBItv i1, i2, one, useless;
  IBSqrI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := u^2 */
  IBSetI(one,1.0,1.0);
  IBSubII(i2,one,i1);                       /* i2 := 1 -  u^2 */
  IBSqrtI(i1,i2,useless);                   /* i1 := sqrt(1 - u^2) */
  IBDivII(i2,one,i1);                       /* i2 := 1 / sqrt(1 - u^2) */
  IBNegI(i1,i2,useless);                    /* i1 := -1 / sqrt(1 - u^2) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdAsinI(IBTree *f)
/***************************************************************************
*  Derivation of operation AsinI              d(asin(u)) = 1 / sqrt(1 - u^2)
*/
{
  IBItv i1, i2, one, useless;
  IBSqrI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := u^2 */
  IBSetI(one,1.0,1.0);
  IBSubII(i2,one,i1);                       /* i2 := 1 - u^2 */
  IBSqrtI(i1,i2,useless);                   /* i1 := sqrt(1 - u^2) */
  IBDivII(i2,one,i1);                       /* i2 := 1 / sqrt(1 - u^2) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i2);
}


void IBBwdAtanI(IBTree *f)
/***************************************************************************
*  Derivation of operation AtanhI                   d(atan(u)) = 1 / (1+u^2)
*/
{
  IBItv i1, i2, one, useless;
  IBSqrI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := u^2 */
  IBSetI(one,1.0,1.0);
  IBAddII(i2,one,i1);                       /* i2 := 1 + u^2 */
  IBDivII(i1,one,i2);                       /* i1 := 1 / (1 + u^2) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


void IBBwdAcoshI(IBTree *f)
/***************************************************************************
*  Derivation of operation AcoshI            d(acosh(u)) = 1 / sqrt(u^2 - 1)
*/
{
  IBItv i1, i2, one, useless;
  IBSqrI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := u^2 */
  IBSetI(one,1.0,1.0);
  IBSubII(i2,i1,one);                       /* i2 := u^2 - 1 */
  IBSqrtI(i1,i2,useless);                   /* i1 := sqrt(u^2 - 1) */
  IBDivII(i2,one,i1);                       /* i2 := 1 / sqrt(u^2 - 1) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i2);
}


void IBBwdAsinhI(IBTree *f)
/***************************************************************************
*  Derivation of operation AsinhI            d(asinh(u)) = 1 / sqrt(u^2 + 1)
*/
{
  IBItv i1, i2, one, useless;
  IBSqrI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := u^2 */
  IBSetI(one,1.0,1.0);
  IBAddII(i2,i1,one);                       /* i2 := u^2 + 1 */
  IBSqrtI(i1,i2,useless);                   /* i1 := sqrt(u^2 + 1) */
  IBDivII(i2,one,i1);                       /* i2 := 1 / sqrt(u^2 + 1) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i2);
}


void IBBwdAtanhI(IBTree *f)
/***************************************************************************
*  Derivation of operation AtanhI                d(atanh(u)) = 1 / (1 - u^2)
*/
{
  IBItv i1, i2, one, useless;
  IBSqrI(i1,IBTfwd(IBTleft(f)),useless);    /* i1 := u^2 */
  IBSetI(one,1.0,1.0);
  IBSubII(i2,one,i1);                       /* i2 := 1 - u^2 */
  IBDivII(i1,one,i2);                       /* i1 := 1 / (1 - u^2) */
  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),i1);
}


int IBHC4AddII(IBTree *f)
/***************************************************************************
*  HC4revise on operation AddII
*/
{
   int result = 1;
   IBhc4SubUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   else
   {
     IBhc4SubUI(IBThc4U(IBTright(f)),IBThc4U(f),IBTfwd(IBTleft(f)));
     if( IBIsEmptyU(IBThc4U(IBTright(f))) )
     {
       result = 0;
     }
   }
   IBResetU(IBThc4U(f));
   return( result );
}  


int IBHC4AddRI(IBTree *f)
/***************************************************************************
*  HC4revise on operation AddRI
*/
{
   int result = 1;
   IBhc4SubUI(IBThc4U(IBTright(f)),IBThc4U(f),IBTfwd(IBTleft(f)));
   if( IBIsEmptyU(IBThc4U(IBTright(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}  


int IBHC4SubII(IBTree *f)
/***************************************************************************
*  HC4revise on operation SubII
*/
{
   int result = 1;
   IBhc4AddUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   else
   {
     IBhc4SubIU(IBThc4U(IBTright(f)),IBTfwd(IBTleft(f)),IBThc4U(f));
     if( IBIsEmptyU(IBThc4U(IBTright(f))) )
     {
       result = 0;
     }
   }
   IBResetU(IBThc4U(f));
   return( result );
}  


int IBHC4SubIR(IBTree *f)
/***************************************************************************
*  HC4revise on operation SubIR
*/
{
   int result = 1;
   IBhc4AddUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}  


int IBHC4SubRI(IBTree *f)
/***************************************************************************
*  HC4revise on operation SubRI
*/
{
   int result = 1;
   IBhc4SubIU(IBThc4U(IBTright(f)),IBTfwd(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTright(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}  


int IBHC4NegI(IBTree *f)
/***************************************************************************
*  HC4revise on operation NegI
*/
{
   int result = 1;
   IBhc4NegU(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}  


int IBHC4MulII(IBTree *f)
/***************************************************************************
*  HC4revise on operation MulII
*/
{
   int result = 1;
   IBhc4DivRelUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   else
   {
     IBhc4DivRelUI(IBThc4U(IBTright(f)),IBThc4U(f),IBTfwd(IBTleft(f)));
     if( IBIsEmptyU(IBThc4U(IBTright(f))) )
     {
       result = 0;
     }
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4MulRI(IBTree *f)
/***************************************************************************
*  HC4revise on operation MulRI
*/
{
   int result = 1;
   IBhc4DivRelUI(IBThc4U(IBTright(f)),IBThc4U(f),IBTfwd(IBTleft(f)));
   if( IBIsEmptyU(IBThc4U(IBTright(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4PowI(IBTree *f)
/***************************************************************************
*  HC4revise on operation PowI
*/
{
   int result = 1;
   IBhc4NthRootRelUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4SqrtI(IBTree *f)
/***************************************************************************
*  HC4revise on operation SqrtI
*/
{
   int result = 1;
   IBhc4SqrUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4ExpI(IBTree *f)
/***************************************************************************
*  HC4revise on operation ExpI
*/
{
   int result = 1;
   IBhc4LogUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4LogI(IBTree *f)
/***************************************************************************
*  HC4revise on operation LogI
*/
{
   int result = 1;
   IBhc4ExpUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4DivII(IBTree *f)
/***************************************************************************
*  HC4revise on operation DivII
*/
{
   int result = 1;
   IBhc4MulUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   else
   {
     IBhc4DivRelIU(IBThc4U(IBTright(f)),IBTfwd(IBTleft(f)),IBThc4U(f));
     if( IBIsEmptyU(IBThc4U(IBTright(f))) )
     {
       result = 0;
     }
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4DivRI(IBTree *f)
/***************************************************************************
*  HC4revise on operation DivRI
*/
{
   int result = 1;
   IBhc4DivRelIU(IBThc4U(IBTright(f)),IBTfwd(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTright(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4DivIR(IBTree *f)
/***************************************************************************
*  HC4revise on operation DivIR
*/
{
   int result = 1;
   IBhc4MulUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4MinII(IBTree *f)
/***************************************************************************
*  HC4revise on operation MinII
*/
{
   int result = 1;
   IBhc4MinUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTleft(f)),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   else
   {
     IBhc4MinUI(IBThc4U(IBTright(f)),IBThc4U(f),IBTfwd(IBTright(f)),IBTfwd(IBTleft(f)));
     if( IBIsEmptyU(IBThc4U(IBTright(f))) )
     {
       result = 0;
     }
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4MaxII(IBTree *f)
/***************************************************************************
*  HC4revise on operation MaxII
*/
{
   int result = 1;
   IBhc4MaxUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTleft(f)),IBTfwd(IBTright(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   else
   {
     IBhc4MaxUI(IBThc4U(IBTright(f)),IBThc4U(f),IBTfwd(IBTright(f)),IBTfwd(IBTleft(f)));
     if( IBIsEmptyU(IBThc4U(IBTright(f))) )
     {
       result = 0;
     }
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4CosI(IBTree *f)
/***************************************************************************
*  HC4revise on operation CosI
*/
{
   int result = 1;
   IBhc4CosRelUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTleft(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4SinI(IBTree *f)
/***************************************************************************
*  HC4revise on operation SinI
*/
{
   int result = 1;
   IBhc4SinRelUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTleft(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4TanI(IBTree *f)
/***************************************************************************
*  HC4revise on operation TanI
*/
{
   int result = 1;
   IBhc4TanRelUI(IBThc4U(IBTleft(f)),IBThc4U(f),IBTfwd(IBTleft(f)));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC4CoshI(IBTree *f)
/***************************************************************************
*  HC4revise on operation CoshI
*/
{
   int result = 1;
   IBhc4CoshRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}

int IBHC4SinhI(IBTree *f)
/***************************************************************************
*  HC4revise on operation SinhI
*/
{
   int result = 1;
   IBhc4SinhRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}

int IBHC4TanhI(IBTree *f)
/***************************************************************************
*  HC4revise on operation TanhI
*/
{
   int result = 1;
   IBhc4TanhRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}

int IBHC4AcosI(IBTree *f)
/***************************************************************************
*  HC4revise on operation AcosI
*/
{
   int result = 1;
   IBhc4AcosRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}

int IBHC4AsinI(IBTree *f)
/***************************************************************************
*  HC4revise on operation AsinI
*/
{
   int result = 1;
   IBhc4AsinRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}

int IBHC4AtanI(IBTree *f)
/***************************************************************************
*  HC4revise on operation AtanI
*/
{
   int result = 1;
   IBhc4AtanRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}

int IBHC4AcoshI(IBTree *f)
/***************************************************************************
*  HC4revise on operation AcoshI
*/
{
   int result = 1;
   IBhc4AcoshRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}

int IBHC4AsinhI(IBTree *f)
/***************************************************************************
*  HC4revise on operation AsinhI
*/
{
   int result = 1;
   IBhc4AsinhRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}

int IBHC4AtanhI(IBTree *f)
/***************************************************************************
*  HC4revise on operation AtanhI
*/
{
   int result = 1;
   IBhc4AtanhRelUI(IBThc4U(IBTleft(f)),IBThc4U(f));
   if( IBIsEmptyU(IBThc4U(IBTleft(f))) )
   {
     result = 0;
   }
   IBResetU(IBThc4U(f));
   return( result );
}


int IBHC3AddII(IBTree *f)
/***************************************************************************
*  HC3revise on operation AddII
*/
{
  IBSubII(IBTbwd(IBTleft(f)),IBTbwd(f),IBTfwd(IBTright(f)));
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );


  IBSubII(IBTbwd(IBTright(f)),IBTbwd(f),IBTfwd(IBTleft(f)));
  IBInterII(IBTbwd(IBTright(f)),IBTbwd(IBTright(f)),IBTfwd(IBTright(f)));
  if( IBEmptyI(IBTbwd(IBTright(f))) ) return( 0 );

  return( 1 );
}  


int IBHC3AddRI(IBTree *f)
/***************************************************************************
*  HC3revise on operation AddRI
*/
{
  IBSubII(IBTbwd(IBTright(f)),IBTbwd(f),IBTfwd(IBTleft(f)));
  IBInterII(IBTbwd(IBTright(f)),IBTbwd(IBTright(f)),IBTfwd(IBTright(f)));
  if( IBEmptyI(IBTbwd(IBTright(f))) ) return( 0 );

  return( 1 );
}  


int IBHC3SubII(IBTree *f)
/***************************************************************************
*  HC3revise on operation SubII
*/
{
  IBAddII(IBTbwd(IBTleft(f)),IBTbwd(f),IBTfwd(IBTright(f)));
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  IBSubII(IBTbwd(IBTright(f)),IBTfwd(IBTleft(f)),IBTbwd(f));
  IBInterII(IBTbwd(IBTright(f)),IBTbwd(IBTright(f)),IBTfwd(IBTright(f)));
  if( IBEmptyI(IBTbwd(IBTright(f))) ) return( 0 );

  return( 1 );
}  


int IBHC3SubIR(IBTree *f)
/***************************************************************************
*  HC3revise on operation SubIR
*/
{
  IBAddII(IBTbwd(IBTleft(f)),IBTbwd(f),IBTfwd(IBTright(f)));
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}  


int IBHC3SubRI(IBTree *f)
/***************************************************************************
*  HC3revise on operation SubRI
*/
{
  IBSubII(IBTbwd(IBTright(f)),IBTfwd(IBTleft(f)),IBTbwd(f));
  IBInterII(IBTbwd(IBTright(f)),IBTbwd(IBTright(f)),IBTfwd(IBTright(f)));
  if( IBEmptyI(IBTbwd(IBTright(f))) ) return( 0 );

  return( 1 );
}  


int IBHC3NegI(IBTree *f)
/***************************************************************************
*  HC3revise on operation NegI
*/
{
  IBInterval *useless;

  IBNegI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}  


int IBHC3MulII(IBTree *f)
/***************************************************************************
*  HC3revise on operation MulII
*/
{
  IBUnion *u;

  u = IBDivRelOneII(IBTbwd(f),IBTfwd(IBTright(f)));
  if( IBInterIU(u,IBTfwd(IBTleft(f))) )
  {
    IBHullU(u,IBTbwd(IBTleft(f)));
    IBFreeU(u);
  }
  else
  {
    IBFreeU(u);
    return( 0 );
  }

  u = IBDivRelOneII(IBTbwd(f),IBTfwd(IBTleft(f)));
  if( IBInterIU(u,IBTfwd(IBTright(f))) )
  {
    IBHullU(u,IBTbwd(IBTright(f)));
    IBFreeU(u);
  }
  else
  {
    IBFreeU(u);
    return( 0 );
  }

  return( 1 );
}


int IBHC3MulRI(IBTree *f)
/***************************************************************************
*  HC3revise on operation MulRI
*/
{
  IBUnion *u;

  u = IBDivRelOneII(IBTbwd(f),IBTfwd(IBTleft(f)));
  if( IBInterIU(u,IBTfwd(IBTright(f))) )
  {
    IBHullU(u,IBTbwd(IBTright(f)));
    IBFreeU(u);
  }
  else
  {
    IBFreeU(u);
    return( 0 );
  }

  return( 1 );
}


int IBHC3PowI(IBTree *f)
/***************************************************************************
*  HC3revise on operation PowI
*/
{
  IBUnion *u;

  u = IBNthRootOneI(IBTbwd(f),IBTfwd(IBTright(f)));
  if( IBInterIU(u,IBTfwd(IBTleft(f))) )
  {
    IBHullU(u,IBTbwd(IBTleft(f)));
    IBFreeU(u);
  }
  else
  {
    IBFreeU(u);
    return( 0 );
  }

  return( 1 );
}


int IBHC3SqrtI(IBTree *f)
/***************************************************************************
*  HC3revise on operation SqrtI
*/
{
  IBInterval *useless;
  IBSqrI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


int IBHC3ExpI(IBTree *f)
/***************************************************************************
*  HC3revise on operation ExpI
*/
{
  IBInterval *useless;
  IBLogI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


int IBHC3LogI(IBTree *f)
/***************************************************************************
*  HC3revise on operation LogI
*/
{
  IBInterval *useless;
  IBExpI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


int IBHC3DivII(IBTree *f)
/***************************************************************************
*  HC3revise on operation DivII
*/
{
  IBUnion *u;

  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),IBTfwd(IBTright(f)));
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  u = IBDivRelOneII(IBTfwd(IBTleft(f)),IBTbwd(f));
  if( IBInterIU(u,IBTfwd(IBTright(f))) )
  {
    IBHullU(u,IBTbwd(IBTright(f)));
    IBFreeU(u);
  }
  else
  {
    IBFreeU(u);
    return( 0 );
  }
  return( 1 );
}


int IBHC3DivRI(IBTree *f)
/***************************************************************************
*  HC3revise on operation DivRI
*/
{
  IBUnion *u;

  u = IBDivRelOneII(IBTfwd(IBTleft(f)),IBTbwd(f));
  if( IBInterIU(u,IBTfwd(IBTright(f))) )
  {
    IBHullU(u,IBTbwd(IBTright(f)));
    IBFreeU(u);
  }
  else
  {
    IBFreeU(u);
    return( 0 );
  }
  return( 1 );
}


int IBHC3DivIR(IBTree *f)
/***************************************************************************
*  HC3revise on operation DivIR
*/
{

  IBMulII(IBTbwd(IBTleft(f)),IBTbwd(f),IBTfwd(IBTright(f)));
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


int IBHC3MinII(IBTree *f)
/***************************************************************************
*  HC3revise on operation MinII
*/
{
  if( IBDisjointII(IBTbwd(f),IBTfwd(IBTright(f))) )
  {
    IBInterII(IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)),IBTbwd(f));
  }
  else
  {
    IBCopyI(IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  }
  if( IBMinI(IBTbwd(f))>IBMinI(IBTbwd(IBTleft(f))) )
  {
    IBMinI(IBTbwd(IBTleft(f))) = IBMinI(IBTbwd(f));
  }
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );


  if( IBDisjointII(IBTbwd(f),IBTfwd(IBTleft(f))) )
  {
    IBInterII(IBTbwd(IBTright(f)),IBTfwd(IBTright(f)),IBTbwd(f));
  }
  else
  {
    IBCopyI(IBTbwd(IBTright(f)),IBTfwd(IBTright(f)));
  }
  if( IBMinI(IBTbwd(f))>IBMinI(IBTbwd(IBTright(f))) )
  {
    IBMinI(IBTbwd(IBTright(f))) = IBMinI(IBTbwd(f));
  }
  if( IBEmptyI(IBTbwd(IBTright(f))) ) return( 0 );

  return( 1 );
}  


int IBHC3MaxII(IBTree *f)
/***************************************************************************
*  HC3revise on operation MaxII
*/
{
  if( IBDisjointII(IBTbwd(f),IBTfwd(IBTright(f))) )
  {
    IBInterII(IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)),IBTbwd(f));
  }
  else
  {
    IBCopyI(IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  }
  if( IBMaxI(IBTbwd(f))<IBMaxI(IBTbwd(IBTleft(f))) )
  {
    IBMaxI(IBTbwd(IBTleft(f))) = IBMaxI(IBTbwd(f));
  }
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );


  if( IBDisjointII(IBTbwd(f),IBTfwd(IBTleft(f))) )
  {
    IBInterII(IBTbwd(IBTright(f)),IBTfwd(IBTright(f)),IBTbwd(f));
  }
  else
  {
    IBCopyI(IBTbwd(IBTright(f)),IBTfwd(IBTright(f)));
  }
  if( IBMaxI(IBTbwd(f))<IBMaxI(IBTbwd(IBTright(f))) )
  {
    IBMaxI(IBTbwd(IBTright(f))) = IBMaxI(IBTbwd(f));
  }
  if( IBEmptyI(IBTbwd(IBTright(f))) ) return( 0 );

  return( 1 );
}  


int IBHC3CosI(IBTree *f)
/***************************************************************************
*  HC3revise on operation CosI
*/
{
  /* not yet implemented */
}


int IBHC3SinI(IBTree *f)
/***************************************************************************
*  HC3revise on operation CosI
*/
{
  /* not yet implemented */
}


int IBHC3TanI(IBTree *f)
/***************************************************************************
*  HC3revise on operation TanI
*/
{
  /* not yet implemented */
}


int IBHC3CoshI(IBTree *f)
/***************************************************************************
*  HC3revise on operation CoshI
*/
{
  IBUnion *u;

  u = IBCoshRelOneI(IBTbwd(f));
  if( IBInterIU(u,IBTfwd(IBTleft(f))) )
  {
    IBHullU(u,IBTbwd(IBTleft(f)));
    IBFreeU(u);
  }
  else
  {
    IBFreeU(u);
    return( 0 );
  }

  return( 1 );
}


int IBHC3SinhI(IBTree *f)
/***************************************************************************
*  HC3revise on operation SinhI
*/
{
  IBInterval *useless;
  IBAsinhI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


int IBHC3TanhI(IBTree *f)
/***************************************************************************
*  HC3revise on operation TanhI
*/
{
  IBInterval *useless;
  IBAtanhI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


int IBHC3AcosI(IBTree *f)
/***************************************************************************
*  HC3revise on operation AcosI
*/
{
  /* not yet implemented */
}


int IBHC3AsinI(IBTree *f)
/***************************************************************************
*  HC3revise on operation AsinI
*/
{
  /* not yet implemented */
}


int IBHC3AtanI(IBTree *f)
/***************************************************************************
*  HC3revise on operation AtanI
*/
{
  /* not yet implemented */
}


int IBHC3AcoshI(IBTree *f)
/***************************************************************************
*  HC3revise on operation AcoshI
*/
{
  IBInterval *useless;
  IBCoshI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


int IBHC3AsinhI(IBTree *f)
/***************************************************************************
*  HC3revise on operation AsinhI
*/
{
  IBInterval *useless;
  IBSinhI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


int IBHC3AtanhI(IBTree *f)
/***************************************************************************
*  HC3revise on operation AtanhI
*/
{
  IBInterval *useless;
  IBTanhI(IBTbwd(IBTleft(f)),IBTbwd(f),useless);
  IBInterII(IBTbwd(IBTleft(f)),IBTbwd(IBTleft(f)),IBTfwd(IBTleft(f)));
  if( IBEmptyI(IBTbwd(IBTleft(f))) ) return( 0 );

  return( 1 );
}


IBOperations IBOperationsInit()
/***************************************************************************
*  To create the array of operations containing for each operation symbol
*  the associated evaluation functions (interval, derivative...)
*/
{
  IBOperations a;
  a = (struct IBO *)malloc(IBNbOp*sizeof(struct IBO));

  IBTevalOp(a,IBOpAddII)    = IBDefAddII;
  IBTevalOp(a,IBOpAddRI)    = IBDefAddRI;
  IBTevalOp(a,IBOpSubII)    = IBDefSubII;
  IBTevalOp(a,IBOpSubRI)    = IBDefSubRI;
  IBTevalOp(a,IBOpSubIR)    = IBDefSubIR;
  IBTevalOp(a,IBOpNegI)     = IBDefNegI;
  IBTevalOp(a,IBOpMulII)    = IBDefMulII;
  IBTevalOp(a,IBOpMulRI)    = IBDefMulRI;
  IBTevalOp(a,IBOpMulRnegI) = IBDefMulRnegI;
  IBTevalOp(a,IBOpMulRposI) = IBDefMulRposI;
  IBTevalOp(a,IBOpDivII)    = IBDefDivII;
  IBTevalOp(a,IBOpDivIR)    = IBDefDivIR;
  IBTevalOp(a,IBOpDivRI)    = IBDefDivRI;
  IBTevalOp(a,IBOpDivIRneg) = IBDefDivIRneg;
  IBTevalOp(a,IBOpDivIRpos) = IBDefDivIRpos;
  IBTevalOp(a,IBOpDivRnegI) = IBDefDivRnegI;
  IBTevalOp(a,IBOpDivRposI) = IBDefDivRposI;
  IBTevalOp(a,IBOpSqrI)     = IBDefSqrI;
  IBTevalOp(a,IBOpSqrtI)    = IBDefSqrtI;
  IBTevalOp(a,IBOpPowI)     = IBDefPowI;
  IBTevalOp(a,IBOpExpI)     = IBDefExpI;
  IBTevalOp(a,IBOpLogI)     = IBDefLogI;
  IBTevalOp(a,IBOpMinII)    = IBDefMinimumII;
  IBTevalOp(a,IBOpMaxII)    = IBDefMaximumII;
  IBTevalOp(a,IBOpCosI)     = IBDefCosI;
  IBTevalOp(a,IBOpSinI)     = IBDefSinI;
  IBTevalOp(a,IBOpTanI)     = IBDefTanI;
  IBTevalOp(a,IBOpCoshI)    = IBDefCoshI;
  IBTevalOp(a,IBOpSinhI)    = IBDefSinhI;
  IBTevalOp(a,IBOpTanhI)    = IBDefTanhI;
  IBTevalOp(a,IBOpAcosI)    = IBDefAcosI;
  IBTevalOp(a,IBOpAsinI)    = IBDefAsinI;
  IBTevalOp(a,IBOpAtanI)    = IBDefAtanI;
  IBTevalOp(a,IBOpAcoshI)   = IBDefAcoshI;
  IBTevalOp(a,IBOpAsinhI)   = IBDefAsinhI;
  IBTevalOp(a,IBOpAtanhI)   = IBDefAtanhI;

  IBTevalBwd(a,IBOpAddII)    = IBBwdAddII;
  IBTevalBwd(a,IBOpAddRI)    = IBBwdAddRI;
  IBTevalBwd(a,IBOpSubII)    = IBBwdSubII;
  IBTevalBwd(a,IBOpSubRI)    = IBBwdSubRI;
  IBTevalBwd(a,IBOpSubIR)    = IBBwdSubIR;
  IBTevalBwd(a,IBOpNegI)     = IBBwdNegI;
  IBTevalBwd(a,IBOpMulII)    = IBBwdMulII;
  IBTevalBwd(a,IBOpMulRI)    = IBBwdMulRI;
  IBTevalBwd(a,IBOpMulRnegI) = IBBwdMulRI;
  IBTevalBwd(a,IBOpMulRposI) = IBBwdMulRI;
  IBTevalBwd(a,IBOpDivII)    = IBBwdDivII;
  IBTevalBwd(a,IBOpDivIR)    = IBBwdDivIR;
  IBTevalBwd(a,IBOpDivRI)    = IBBwdDivRI;
  IBTevalBwd(a,IBOpDivIRneg) = IBBwdDivIR;
  IBTevalBwd(a,IBOpDivIRpos) = IBBwdDivIR;
  IBTevalBwd(a,IBOpDivRnegI) = IBBwdDivRI;
  IBTevalBwd(a,IBOpDivRposI) = IBBwdDivRI;
  IBTevalBwd(a,IBOpSqrI)     = IBBwdSqrI;
  IBTevalBwd(a,IBOpSqrtI)    = IBBwdSqrtI;
  IBTevalBwd(a,IBOpPowI)     = IBBwdPowI;
  IBTevalBwd(a,IBOpExpI)     = IBBwdExpI;
  IBTevalBwd(a,IBOpLogI)     = IBBwdLogI;
  IBTevalBwd(a,IBOpMinII)    = IBBwdMinII;
  IBTevalBwd(a,IBOpMaxII)    = IBBwdMaxII;
  IBTevalBwd(a,IBOpCosI)     = IBBwdCosI;
  IBTevalBwd(a,IBOpSinI)     = IBBwdSinI;
  IBTevalBwd(a,IBOpTanI)     = IBBwdTanI;
  IBTevalBwd(a,IBOpCoshI)    = IBBwdCoshI;
  IBTevalBwd(a,IBOpSinhI)    = IBBwdSinhI;
  IBTevalBwd(a,IBOpTanhI)    = IBBwdTanhI;
  IBTevalBwd(a,IBOpAcosI)    = IBBwdAcosI;
  IBTevalBwd(a,IBOpAsinI)    = IBBwdAsinI;
  IBTevalBwd(a,IBOpAtanI)    = IBBwdAtanI;
  IBTevalBwd(a,IBOpAcoshI)   = IBBwdAcoshI;
  IBTevalBwd(a,IBOpAsinhI)   = IBBwdAsinhI;
  IBTevalBwd(a,IBOpAtanhI)   = IBBwdAtanhI;

  IBTevalHC4(a,IBOpAddII)    = IBHC4AddII;
  IBTevalHC4(a,IBOpAddRI)    = IBHC4AddRI;
  IBTevalHC4(a,IBOpSubII)    = IBHC4SubII;
  IBTevalHC4(a,IBOpSubRI)    = IBHC4SubRI;
  IBTevalHC4(a,IBOpSubIR)    = IBHC4SubIR;
  IBTevalHC4(a,IBOpNegI)     = IBHC4NegI;
  IBTevalHC4(a,IBOpMulII)    = IBHC4MulII;
  IBTevalHC4(a,IBOpMulRI)    = IBHC4MulRI;
  IBTevalHC4(a,IBOpMulRnegI) = IBHC4MulRI;
  IBTevalHC4(a,IBOpMulRposI) = IBHC4MulRI;
  IBTevalHC4(a,IBOpDivII)    = IBHC4DivII;
  IBTevalHC4(a,IBOpDivIR)    = IBHC4DivIR;
  IBTevalHC4(a,IBOpDivRI)    = IBHC4DivRI;
  IBTevalHC4(a,IBOpDivIRneg) = IBHC4DivIR;
  IBTevalHC4(a,IBOpDivIRpos) = IBHC4DivIR;
  IBTevalHC4(a,IBOpDivRnegI) = IBHC4DivRI;
  IBTevalHC4(a,IBOpDivRposI) = IBHC4DivRI;
  IBTevalHC4(a,IBOpSqrI)     = IBHC4PowI;
  IBTevalHC4(a,IBOpSqrtI)    = IBHC4SqrtI;
  IBTevalHC4(a,IBOpPowI)     = IBHC4PowI;
  IBTevalHC4(a,IBOpExpI)     = IBHC4ExpI;
  IBTevalHC4(a,IBOpLogI)     = IBHC4LogI;
  IBTevalHC4(a,IBOpMinII)    = IBHC4MinII;
  IBTevalHC4(a,IBOpMaxII)    = IBHC4MaxII;
  IBTevalHC4(a,IBOpCosI)     = IBHC4CosI;
  IBTevalHC4(a,IBOpSinI)     = IBHC4SinI;
  IBTevalHC4(a,IBOpTanI)     = IBHC4TanI;
  IBTevalHC4(a,IBOpCoshI)    = IBHC4CoshI;
  IBTevalHC4(a,IBOpSinhI)    = IBHC4SinhI;
  IBTevalHC4(a,IBOpTanhI)    = IBHC4TanhI;
  IBTevalHC4(a,IBOpAcosI)    = IBHC4AcosI;
  IBTevalHC4(a,IBOpAsinI)    = IBHC4AsinI;
  IBTevalHC4(a,IBOpAtanI)    = IBHC4AtanI;
  IBTevalHC4(a,IBOpAcoshI)   = IBHC4AcoshI;
  IBTevalHC4(a,IBOpAsinhI)   = IBHC4AsinhI;
  IBTevalHC4(a,IBOpAtanhI)   = IBHC4AtanhI;

  IBTevalHC3(a,IBOpAddII)    = IBHC3AddII;
  IBTevalHC3(a,IBOpAddRI)    = IBHC3AddRI;
  IBTevalHC3(a,IBOpSubII)    = IBHC3SubII;
  IBTevalHC3(a,IBOpSubRI)    = IBHC3SubRI;
  IBTevalHC3(a,IBOpSubIR)    = IBHC3SubIR;
  IBTevalHC3(a,IBOpNegI)     = IBHC3NegI;
  IBTevalHC3(a,IBOpMulII)    = IBHC3MulII;
  IBTevalHC3(a,IBOpMulRI)    = IBHC3MulRI;
  IBTevalHC3(a,IBOpMulRnegI) = IBHC3MulRI;
  IBTevalHC3(a,IBOpMulRposI) = IBHC3MulRI;
  IBTevalHC3(a,IBOpDivII)    = IBHC3DivII;
  IBTevalHC3(a,IBOpDivIR)    = IBHC3DivIR;
  IBTevalHC3(a,IBOpDivRI)    = IBHC3DivRI;
  IBTevalHC3(a,IBOpDivIRneg) = IBHC3DivIR;
  IBTevalHC3(a,IBOpDivIRpos) = IBHC3DivIR;
  IBTevalHC3(a,IBOpDivRnegI) = IBHC3DivRI;
  IBTevalHC3(a,IBOpDivRposI) = IBHC3DivRI;
  IBTevalHC3(a,IBOpSqrI)     = IBHC3PowI;
  IBTevalHC3(a,IBOpSqrtI)    = IBHC3SqrtI;
  IBTevalHC3(a,IBOpPowI)     = IBHC3PowI;
  IBTevalHC3(a,IBOpExpI)     = IBHC3ExpI;
  IBTevalHC3(a,IBOpLogI)     = IBHC3LogI;
  IBTevalHC3(a,IBOpMinII)    = IBHC3MinII;
  IBTevalHC3(a,IBOpMaxII)    = IBHC3MaxII;
  IBTevalHC3(a,IBOpCosI)     = IBHC3CosI;
  IBTevalHC3(a,IBOpSinI)     = IBHC3SinI;
  IBTevalHC3(a,IBOpTanI)     = IBHC3TanI;
  IBTevalHC3(a,IBOpCoshI)    = IBHC3CoshI;
  IBTevalHC3(a,IBOpSinhI)    = IBHC3SinhI;
  IBTevalHC3(a,IBOpTanhI)    = IBHC3TanhI;
  IBTevalHC3(a,IBOpAcosI)    = IBHC3AcosI;
  IBTevalHC3(a,IBOpAsinI)    = IBHC3AsinI;
  IBTevalHC3(a,IBOpAtanI)    = IBHC3AtanI;
  IBTevalHC3(a,IBOpAcoshI)   = IBHC3AcoshI;
  IBTevalHC3(a,IBOpAsinhI)   = IBHC3AsinhI;
  IBTevalHC3(a,IBOpAtanhI)   = IBHC3AtanhI;

  return( a );
}


void IBOperationsFree(IBOperations a)
/***************************************************************************
*  To desallocate the array of operations
*/
{
  free(a);
}


void IBTevalOnevar(IBTree *f, IBItv domvar, int globvar,
                   struct IBListDepNodes *list, IBDomains d)
/***************************************************************************
*  To evaluate the interval expression represented by tree
*  s.t. only the nodes depending on globvar (contained in list)
*  have to be (re-)evaluated
*
*    - ivar is the domain of the variable of global index globvar
*    - d contains the domains of the other variables
*    - The results are in forward intervals of each node of tree
*
*  Remark that f is not used
*/
{
  IBTree *g;

  while( list!=NULL )
  {
    g = list->t;
    if( IBTtype(g)==IBTNodeVar )   /* necessary this variable */
    {
      IBCopyI(IBTfwd(g),domvar);
    }
    else                           /* necessary an operation */
    {
      (* IBTevalOp(operations,IBTop(g)))
                     (IBTfwd(g),
                      IBTfwd(IBTleft(g)),
                      IBTfwd(IBTright(g)));

      if( IBEmptyI(IBTfwd(g)) )
      {
	IBSetEmptyI(IBTfwd(f));
        return;
      }
    }
    list = list->next;
  }
}


void IBTeval(IBTree *f, IBItv domvar, int globvar, IBDomains d)
/***************************************************************************
*  To evaluate the interval expression represented by tree
*    - ivar is the domain of the variable of global index globvar
*    - d contains the domains of the other variables
*    - The results are in forward intervals of each node of tree
*/
{
  switch( IBTtype(f) )
  {
    case IBTNodeVar:
         if( IBTglobvar(f)==globvar )
	 {
           IBCopyI(IBTfwd(f),domvar);
	 }
         else
	 {
           IBCopyI(IBTfwd(f),IBDomV(d,IBTglobvar(f)));
	 }
         break;
    case IBTNodeOp:
         /* Left sub-expression */
         IBTeval(IBTleft(f),domvar,globvar,d);

	 if( IBEmptyI(IBTfwd(IBTleft(f))) )
	 {
	   IBSetEmptyI(IBTfwd(f));
           return;
	 }

         /* Right sub-expression */
         if( IBTtype(IBTright(f))!=IBTNodeUseless )
         {
           IBTeval(IBTright(f),domvar,globvar,d);

	   if( IBEmptyI(IBTfwd(IBTright(f))) )
	   {
  	     IBSetEmptyI(IBTfwd(f));
             return;
	   }
	 }

         (* IBTevalOp(operations,IBTop(f)))
                   (IBTfwd(f),
                    IBTfwd(IBTleft(f)),
                    IBTfwd(IBTright(f)));
         break;

    /* else if IBTNodeItv the result is already in IBTfwd */
  }

}


void IBTevalAll(IBTree *f, IBDomains d)
/***************************************************************************
*  To evaluate the interval expression represented by tree
*    - d contains the domains of the variables
*    - The results are in forward intervals of each node of tree
*/
{
  switch( IBTtype(f) )
  {
    case IBTNodeVar:
         IBCopyI(IBTfwd(f),IBDomV(d,IBTglobvar(f)));
         break;
    case IBTNodeOp:
         /* Left sub-expression */
         IBTevalAll(IBTleft(f),d);

	 if( IBEmptyI(IBTfwd(IBTleft(f))) )
	 {
	   IBSetEmptyI(IBTfwd(f));
           return;
	 }

         /* Right sub-expression */
         if( IBTtype(IBTright(f))!=IBTNodeUseless )
         {
           IBTevalAll(IBTright(f),d);

	   if( IBEmptyI(IBTfwd(IBTright(f))) )
	   {
	     IBSetEmptyI(IBTfwd(f));
             return;
	   }
	 }

         /* Expression f */
         (* IBTevalOp(operations,IBTop(f)))
                   (IBTfwd(f),
                    IBTfwd(IBTleft(f)),
                    IBTfwd(IBTright(f)));
         break;

    /* else if IBTNodeItv the result is already in IBTfwd */
  }
}


void IBTevalConstant(IBTree *f)
/***************************************************************************
*  To evaluate the constant interval expression represented by f
*/
{
  if( IBTtype(f)==IBTNodeOp )
  {
    /* Left sub-expression */
    IBTevalConstant(IBTleft(f));

    if( IBEmptyI(IBTfwd(IBTleft(f))) )
    {
      IBSetEmptyI(IBTfwd(f));
      return;
    }

    /* Right sub-expression */

    if( IBTtype(IBTright(f))!=IBTNodeUseless )
    {
      IBTevalConstant(IBTright(f));

      if( IBEmptyI(IBTfwd(IBTright(f))) )
      {
        IBSetEmptyI(IBTfwd(f));
        return;
      }
    }

    (* IBTevalOp(operations,IBTop(f)))
               (IBTfwd(f),
                IBTfwd(IBTleft(f)),
                IBTfwd(IBTright(f)));
  }
  else if( IBTtype(f)==IBTNodeVar )
  {
    IBSetEmptyI(IBTfwd(f));
    return;
  }

  /* else if IBTNodeItv the result is already in IBTfwd */
}


int IBTIsConstant(IBTree *f)
/***************************************************************************
*  Returns 1 if f is free variable, 0 otherwise
*/
{
  if( f==NULL ) return( 1 );
  else if( IBTtype(f)==IBTNodeVar ) return( 0 );
  else if( IBTtype(f)==IBTNodeOp )
  {
    if( IBTIsConstant(IBTleft(f)) )
         return( IBTIsConstant(IBTright(f)) );
    else return( 0 );
  }
  else return( 1 );
}


IBTree *IBRemoveConstantSubtrees(IBTree *f)
/***************************************************************************
*  Evaluates with IBTevalConstant and replaces all free variables
*  subtrees in f 
*  Returns NULL if an empty interval is computed, e.g., sqrt(-1)
*/
{
  IBTree *t;

  if( IBTtype(f)==IBTNodeOp )
  {
    if( IBTIsConstant(f) )
    {
      IBTevalConstant(f);
      if( IBEmptyI(IBTfwd(f)) )
      {
	return NULL;
      }
      else
      {
        t = IBTNewItv(IBTfwd(f));
        IBTFree(f);
        return( t );
      }
    }
    else
    {
      if( (IBTleft(f)=IBRemoveConstantSubtrees(IBTleft(f))) != NULL )
      {
        if( (IBTright(f)=IBRemoveConstantSubtrees(IBTright(f))) !=NULL )
	{
          return( f );
	}
      }
      return( NULL );
    }
  }
  else
  {
    return( f );
  }
}


void IBCderivRec(IBConstraint *c, IBTree *f)
/***************************************************************************
*  To evaluate all the backward intervals in f
*/
{
  if( IBTtype(f)==IBTNodeOp )
  {
    (* IBTevalBwd(operations,IBTop(f)))(f);
    IBCderivRec(c,IBTleft(f));
    IBCderivRec(c,IBTright(f));
  }
  else if( IBTtype(f)==IBTNodeVar )
  {
    /* Sum for all occurrences of variables */
    IBAddII(IBCVderiv(c,IBTlocvar(f)),
            IBCVderiv(c,IBTlocvar(f)),IBTbwd(f));
  }
}


void IBCderiv(IBConstraint *c)
/***************************************************************************
*  To evaluate all the partial derivatives of function in c
*  Forward intervals are supposed to be computed
*/
{
  int i;

  /* Initialization to 0 of intervals for the partial derivatives w.r.t.
     every variable in c since they are computed as the sum of all backward
     intervals associated to the occurrences of this variable */
  for( i=0; i<IBCNbVar(c); i++ )
  {
    IBSetI(IBCVderiv(c,i),0.0,0.0);
  }
  IBCderivRec(c,IBCfunc(c));
}



int IBCInnerBoxOfConstraint(IBConstraint *c, IBDomains d)
/***************************************************************************
*  Returns 1 if d is an inner approximation of c
*/
{
  if( IBCrelfunc(c)==IBRelationINT )           /* c is equivalent to integer(x) */
  {
    return( IBIsIntegerI(IBDomV(d,IBCVglobvar(c,0))) );
    /* true if the domain of x is reduced to an integer */
  }
  else if( IBCrelfunc(c)==IBRelationEQU )      /* c is equivalent to "func = 0" */
  {
    IBTevalAll(IBCfunc(c),d);
    if( (IBMinI(IBTfwd(IBCfunc(c))) == 0.0) &&
        (IBMaxI(IBTfwd(IBCfunc(c))) == 0.0) )
    {
      return 1;
    }
    else
    {
      return 0;
    }

  }
  else if( IBCrelfunc(c)==IBRelationSET )      /* c is equivalent to "l in r" */
  {
    IBTevalAll(IBCleft(c),d);
    IBTevalAll(IBCright(c),d);

    /*
    printf("\n");
    IBWriteI(stdout,IBTfwd(IBCleft(c)),16,1);
    printf("is in ? ");
    IBWriteI(stdout,IBTfwd(IBCright(c)),8,1);
    */


    if( IBIncludedII(IBTfwd(IBCleft(c)),IBTfwd(IBCright(c))) )
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
  else if( IBCrelfunc(c)==IBRelationSUP )      /* c is equivalent to "func >= 0" */
  {
    IBTevalAll(IBCfunc(c),d);
    if( IBMinI(IBTfwd(IBCfunc(c))) >= 0.0 )
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
  else                                         /* c is equivalent to "func <= 0" */
  {
    IBTevalAll(IBCfunc(c),d);
    if( IBMaxI(IBTfwd(IBCfunc(c))) <= 0.0 )
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
}
