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
 * default_interval.c                                                       *
 ****************************************************************************/

#include "config.h"
#include "default_interval.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


double IBBasicEpsilon;   /* machine's epsilon */


void IBBasicIntervalInit() {
/***************************************************************************
*  Initialization of the interval arithmetic module
*/
  double x, test;

  /* Computation of machine's epsilon */
  IBBasicRoundNear();
  x = 1.0;
  do {
      IBBasicEpsilon = x;
      x /= 2.0;
      test = 1.0 + x;
  }
  while (test>1.0);


  IBBasicSetToPi(IBBasicItvConstPi);
  IBBasicSetToHalfPi(IBBasicItvConst_1_Pi_2);
  IBBasicMulRposIinternal(IBBasicItvConst_2_Pi,2.0,IBBasicItvConstPi);
  IBBasicMulRposIinternal(IBBasicItvConst_4_Pi,4.0,IBBasicItvConstPi);
  IBBasicMulRposIinternal(IBBasicItvConst_3_Pi_2,3.0,IBBasicItvConst_1_Pi_2);
  IBBasicMulRposIinternal(IBBasicItvConst_5_Pi_2,5.0,IBBasicItvConst_1_Pi_2);
  IBBasicMulRposIinternal(IBBasicItvConst_7_Pi_2,7.0,IBBasicItvConst_1_Pi_2);
}

inline IBBasicBounds *IBBasicNewI()
/***************************************************************************
*  Allocation of an interval
*/
{
  return( (IBBasicBounds *)malloc(sizeof(IBBasicBounds)) );
}


IBBasicBounds *IBBasicNewLargestI()
/***************************************************************************
*  Allocation of an interval which value is [-oo,+oo]
*/
{
  IBBasicBounds *i;
  i = (IBBasicBounds *)malloc(sizeof(IBBasicBounds));
  IBBasicMinI(i) = IBBasicNegInfinity;
  IBBasicMaxI(i) = IBBasicPosInfinity;
  return( i );
}


inline
void IBBasicToLargestI(IBBasicItv i)
/***************************************************************************
*  i := [-oo,+oo]
*/
{
  IBBasicMinI(i) = IBBasicNegInfinity;
  IBBasicMaxI(i) = IBBasicPosInfinity;
}


IBBasicBounds *IBBasicNewCopyI(IBBasicItv i)
/***************************************************************************
*  Allocation of an interval initialized to i (i is not empty)
*/
{
  IBBasicBounds *i1;
  i1 = (IBBasicBounds *)malloc(sizeof(IBBasicBounds));
  IBBasicCopyI(i1,i);
  return( i1 );
}


IBBasicBounds *IBBasicSetNewI(double x1, double x2)
/***************************************************************************
*  Allocation of an interval initialized to [x1,x2], x1<=x2
*/
{
  IBBasicBounds *i;
  i = (IBBasicBounds *)malloc(sizeof(IBBasicBounds));
  IBBasicSetI(i,x1,x2);
  return( i );
}


void IBBasicStringToI(char* s, IBBasicItv i)
/***************************************************************************
*  Conversion of a string representing an 'unsigned' float to an interval
*  Replace atof(s) which is not able to round correctly
*
*  Algorithm : ex. 1.157e-5  => 1157 / 10**2
*
*  - The integer part (1157) is converted as follows:
*         7 + 10*( 5 + 10( 1 + 10*( 1 ) ) )  using interval operations
*
*  - The exponent is computed as the exponent (if specified) minus
*    the number of decimals, here -5 - 3 => -8
*
*  - The quantity 10^|exponent| is evaluated using interval operations
*    here, 10^8 => 100000000
*
*  - If exponent<0, the result is (integer_part / 10^|exponent|)
*    Otherwise (integer_part * 10^|exponent|)
*
*    here 1157 / 100000000
*/
{
  int i_s=0, i_conv=0, j, ndecimals=0, expo, expoabs, dig;
  char conv[100];
  IBBasicItv ten, exponent, digit, intpart;

  /* extraction of integer part */
  while( isdigit(s[i_s]) )
  {
    conv[i_conv++]=s[i_s++];
  }
  /* extraction of decimal part */
  if( s[i_s]=='.' )
  {
    ++i_s;
    while( isdigit(s[i_s]) )
    {
      conv[i_conv++]=s[i_s++];
      ++ndecimals;
    }
  }
  conv[i_conv] = '\0';

  if( i_s<strlen(s) )
  {
    /* extraction of exponent part */
    ++i_s;  /* e/E */
    expo = atoi(&(s[i_s]));
  }
  else expo=0;
  expo -= ndecimals;

  IBBasicSetI(ten,10.0,10.0);

  /* computation of the integer part */
  IBBasicSetI(intpart,0.0,0.0);
  for( j=0; j<strlen(conv); ++j )
  {
    dig = (int)(conv[j] - '0');
    IBBasicSetI(digit,dig,dig);
    IBBasicMulRposI(intpart,ten,intpart);
    IBBasicAddRI(intpart,digit,intpart);
  }

  /* computation of 10^|expo| */
  if( expo!=0 )
  {
    IBBasicSetI(exponent,1.0,1.0);
    expoabs = ((expo<0) ? (-expo) : expo);
    for( j=0; j<expoabs; ++j )
    {
      IBBasicMulRposI(exponent,ten,exponent);
    }
  }

  if( expo<0 )
  {
    IBBasicDivII(i,intpart,exponent);
  }
  else if( expo>0 )
  {
    IBBasicMulII(i,intpart,exponent);
  }
  else
  {
    IBBasicCopyI(i,intpart);
  }
}


inline void IBBasicToIntegerI(IBBasicItv i)
/***************************************************************************
*  The bounds of i are rounded to integers
*  let i=[a,d],  i := [ceil(a),floor(b)]
*/
{
  if (IBBasicMinI(i)!=IBBasicNegInfinity)
    IBBasicMinI(i) = ceil(IBBasicMinI(i));

  if (IBBasicMaxI(i)!=IBBasicPosInfinity)
    IBBasicMaxI(i) = floor(IBBasicMaxI(i)); 
}


void IBBasicAbsI(IBBasicItv Result, IBBasicItv i)
/***************************************************************************
*  Result := absolute value of i
*/
{
  if( IBBasicMinI(i)>=0.0 )
  {
    IBBasicMinI(Result) = IBBasicMinI(i);
    IBBasicMaxI(Result) = IBBasicMaxI(i);
  }
  else if( IBBasicMaxI(i)<=0.0 )
  {
    IBBasicMinI(Result) = -(IBBasicMaxI(i));
    IBBasicMaxI(Result) = -(IBBasicMinI(i));
  }
  else
  {
    IBBasicMinI(Result) = 0.0;
    IBBasicMaxI(Result) = IBBasicMax(-(IBBasicMinI(i)),IBBasicMaxI(i));
  }
}


void IBBasicAddII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2)
/***************************************************************************
*  Result := i1 + i2
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = IBBasicMinI(i1) + IBBasicMinI(i2);
  IBBasicRoundUp();
  IBBasicMaxI(Result) = IBBasicMaxI(i1) + IBBasicMaxI(i2);
}


void IBBasicAddRI(IBBasicItv Result, IBBasicItv x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r] + i
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = IBBasicMinI(x) + IBBasicMinI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result) = IBBasicMaxI(x) + IBBasicMaxI(i);
}


void IBBasicAddRIinternal(IBBasicItv Result, double x, IBBasicItv i)
/***************************************************************************
*  Result := x + i
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = IBBasicMinI(i) + x;
  IBBasicRoundUp();
  IBBasicMaxI(Result) = IBBasicMaxI(i) + x;
}


void IBBasicSubII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2)
/***************************************************************************
*  Result := i1 - i2
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = IBBasicMinI(i1) - IBBasicMaxI(i2);
  IBBasicRoundUp();
  IBBasicMaxI(Result) = IBBasicMaxI(i1) - IBBasicMinI(i2);
}


void IBBasicSubRI(IBBasicItv Result, IBBasicItv x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r] - i
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = IBBasicMinI(x) - IBBasicMaxI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result) = IBBasicMaxI(x) - IBBasicMinI(i);
}


void IBBasicSubRIinternal(IBBasicItv Result, double x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r] - i
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = x - IBBasicMaxI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result) = x - IBBasicMinI(i);
}


void IBBasicSubIR(IBBasicItv Result, IBBasicItv i, IBBasicItv x)
/***************************************************************************
*  let x=[r,r], Result := i - [r,r]
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = IBBasicMinI(i) - IBBasicMaxI(x);
  IBBasicRoundUp();
  IBBasicMaxI(Result) = IBBasicMaxI(i) - IBBasicMinI(x);
}


void IBBasicSubIRinternal(IBBasicItv Result, IBBasicItv i, double x)
/***************************************************************************
*  let x=[r,r], Result := i - [r,r]
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = IBBasicMinI(i) - x;
  IBBasicRoundUp();
  IBBasicMaxI(Result) = IBBasicMaxI(i) - x;
}

void IBBasicNegI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := -i
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result) = -(IBBasicMaxI(i));
  IBBasicRoundUp();
  IBBasicMaxI(Result) = -(IBBasicMinI(i));
}


void IBBasicMulII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2)
/***************************************************************************
*  Result := i1*i2
*/
{
  unsigned int sig = ((IBBasicMaxI(i1)<0.0) << 3) | ((IBBasicMinI(i1)>0.0) << 2) 
                   | ((IBBasicMaxI(i2)<0.0) << 1) | (IBBasicMinI(i2)>0.0);

  double l1, l2, u1, u2;

  switch( sig )
  {
    case 0: /* 0000: 0 in i1, 0 in i2 */
      if ( IBBasicInfinite(i1) || IBBasicInfinite(i2) ) {
	if ( ((IBBasicMinI(i1)>=0.0) && (IBBasicMinI(i2)>=0.0)) ||
             ((IBBasicMaxI(i1)<=0.0) && (IBBasicMaxI(i2)<=0.0)) ) {
          IBBasicSetI(Result,0.0,IBBasicPosInfinity);
	}
	else
	if ( ((IBBasicMinI(i1)>=0.0) && (IBBasicMaxI(i2)<=0.0)) ||
             ((IBBasicMaxI(i1)<=0.0) && (IBBasicMinI(i2)>=0.0)) ) {
          IBBasicSetI(Result,IBBasicNegInfinity,0.0);
	}
        else {
	  IBBasicToLargestI(Result);
	}
      }
      else {
        IBBasicRoundDown();
        l1=IBBasicMinI(i1)*IBBasicMaxI(i2);
        l2=IBBasicMaxI(i1)*IBBasicMinI(i2);

        IBBasicRoundUp();
        u1=IBBasicMinI(i1)*IBBasicMinI(i2);
        u2=IBBasicMaxI(i1)*IBBasicMaxI(i2);

        IBBasicMinI(Result)=IBBasicMin(l1,l2);
        IBBasicMaxI(Result)=IBBasicMax(u1,u2);
      }
    break;

    case 1: /* 0001: 0 in i1, i2>0 */
      if (IBBasicMaxI(i2)==IBBasicPosInfinity) {
        if( IBBasicMinI(i1)==0.0 ) {
          IBBasicSetI(Result,0.0,IBBasicPosInfinity);
	}
        else if ( IBBasicMaxI(i1)==0.0 ) {
          IBBasicSetI(Result,IBBasicNegInfinity,0.0);
	}
        else {
	  IBBasicToLargestI(Result);
	}
      }
      else {
        IBBasicRoundDown();
        IBBasicMinI(Result)=IBBasicMinI(i1)*IBBasicMaxI(i2);
        IBBasicRoundUp();
        IBBasicMaxI(Result)=IBBasicMaxI(i1)*IBBasicMaxI(i2);
      }
    break;

    case 2: /* 0010: 0 in i1, i2<0 */
      if (IBBasicMinI(i2)==IBBasicNegInfinity) {
        if( IBBasicMinI(i1)==0.0 ) {
          IBBasicSetI(Result,IBBasicNegInfinity,0.0);
	}
        else if ( IBBasicMaxI(i1)==0.0 ) {
          IBBasicSetI(Result,0.0,IBBasicPosInfinity);
	}
        else {
  	  IBBasicToLargestI(Result);
	}
      }
      else {
        IBBasicRoundDown();
        IBBasicMinI(Result)=IBBasicMaxI(i1)*IBBasicMinI(i2);
        IBBasicRoundUp();
        IBBasicMaxI(Result)=IBBasicMinI(i1)*IBBasicMinI(i2);
      }
    break;

    case 4: /* 0100: i1>0, 0 in i2 */
      if (IBBasicMaxI(i1)==IBBasicPosInfinity) {
        if( IBBasicMinI(i2)==0.0 ) {
          IBBasicSetI(Result,0.0,IBBasicPosInfinity);
	}
        else if ( IBBasicMaxI(i2)==0.0 ) {
          IBBasicSetI(Result,IBBasicNegInfinity,0.0);
	}
        else {
	  IBBasicToLargestI(Result);
	}
      }
      else {
        IBBasicRoundDown();
        IBBasicMinI(Result)=IBBasicMaxI(i1)*IBBasicMinI(i2);
        IBBasicRoundUp();
        IBBasicMaxI(Result)=IBBasicMaxI(i1)*IBBasicMaxI(i2);
      }
    break;

    case 5: /* 0101: i1>0, i2>0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMinI(i1)*IBBasicMinI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMaxI(i1)*IBBasicMaxI(i2);
    break;

    case 6: /* 0110: i1>0, i2<0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMaxI(i1)*IBBasicMinI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMinI(i1)*IBBasicMaxI(i2);
    break;

    case 8: /* 1000: i1<0, 0 in i2 */
      if (IBBasicMinI(i1)==IBBasicNegInfinity) {
        if( IBBasicMinI(i2)==0.0 ) {
          IBBasicSetI(Result,IBBasicNegInfinity,0.0);
	}
        else if ( IBBasicMaxI(i2)==0.0 ) {
          IBBasicSetI(Result,0.0,IBBasicPosInfinity);
	}
        else {
	  IBBasicToLargestI(Result);
	}
      }
      else {
        IBBasicRoundDown();
        IBBasicMinI(Result)=IBBasicMinI(i1)*IBBasicMaxI(i2);
        IBBasicRoundUp();
        IBBasicMaxI(Result)=IBBasicMinI(i1)*IBBasicMinI(i2);
      }
    break;

    case 9: /* 1001: i1<0, i2>0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMinI(i1)*IBBasicMaxI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMaxI(i1)*IBBasicMinI(i2);
    break;

    case 10: /* 1010: i1<0, i2<0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMaxI(i1)*IBBasicMaxI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMinI(i1)*IBBasicMinI(i2);
    break;
  }
}


void IBBasicMulRI(IBBasicItv Result, IBBasicItv x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r]*i
*/
{
  if( IBBasicMinI(x)==0.0 )
  {
    if (IBBasicInfinite(i)) {    /* 0*[a,b] where a=-oo or b=+oo */
      IBBasicToLargestI(Result);
    }
    else {
      IBBasicSetI(Result,0.0,0.0);
    }
    return;
  }

  if( IBBasicMinI(x)<0.0 )
  {
    IBBasicRoundDown();
    IBBasicMinI(Result)=IBBasicMinI(x)*IBBasicMaxI(i);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=IBBasicMinI(x)*IBBasicMinI(i);
  }
  else
  {
    IBBasicRoundDown();
    IBBasicMinI(Result)=IBBasicMinI(x)*IBBasicMinI(i);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=IBBasicMinI(x)*IBBasicMaxI(i);
  }
}


void IBBasicMulRnegI(IBBasicItv Result, IBBasicItv x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r]*i with r<0
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result)=IBBasicMinI(x)*IBBasicMaxI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result)=IBBasicMinI(x)*IBBasicMinI(i);
}


void IBBasicMulRposI(IBBasicItv Result, IBBasicItv x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r]*i with r>=0
*/
{
  if( IBBasicMinI(x)==0.0 )
  {
    if (IBBasicInfinite(i)) {    /* 0*[a,b] where a=-oo or b=+oo */
      IBBasicToLargestI(Result);
    }
    else {
      IBBasicSetI(Result,0.0,0.0);
    }
    return;
  }

  IBBasicRoundDown();
  IBBasicMinI(Result)=IBBasicMinI(x)*IBBasicMinI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result)=IBBasicMinI(x)*IBBasicMaxI(i);
}


void IBBasicMulRposIinternal(IBBasicItv Result, double x, IBBasicItv i)
/***************************************************************************
*  Result := x*i with x>=0
*/
{
  if( x==0.0 )
  {
    if (IBBasicInfinite(i)) {    /* 0*[a,b] where a=-oo or b=+oo */
      IBBasicToLargestI(Result);
    }
    else {
      IBBasicSetI(Result,0.0,0.0);
    }
    return;
  }

  IBBasicRoundDown();
  IBBasicMinI(Result)=x*IBBasicMinI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result)=x*IBBasicMaxI(i);
}


void IBBasicMulRIinternal(IBBasicItv Result, double x, IBBasicItv i)
/***************************************************************************
*  Result := x*i
*/
{
  if( x==0.0 )
  {
    if (IBBasicInfinite(i)) {    /* 0*[a,b] where a=-oo or b=+oo */
      IBBasicToLargestI(Result);
    }
    else {
      IBBasicSetI(Result,0.0,0.0);
    }
    return;
  }

  if( x<0.0 )
  {
    IBBasicRoundDown();
    IBBasicMinI(Result)=x*IBBasicMaxI(i);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=x*IBBasicMinI(i);
  }
  else
  {
    IBBasicRoundDown();
    IBBasicMinI(Result)=x*IBBasicMinI(i);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=x*IBBasicMaxI(i);
  }
}


void IBBasicDivII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2)
/***************************************************************************
*  Result := i1 / i2
*/
{
  unsigned int sig;

  if (IBBasicDoubleInI(i2,0.0)) {
    if (IBBasicMinI(i2)==0.0) {
      if (IBBasicMinI(i1)<0.0) {      /* ex. [-3,-1] / [0,4] => [-oo,-1/4] */
	IBBasicMinI(Result)=IBBasicNegInfinity;
        IBBasicRoundUp();
        IBBasicMaxI(Result)=IBBasicMaxI(i1)/IBBasicMaxI(i2);
      }
      else if (IBBasicMaxI(i1)>0.0) { /* ex. [1,3] / [0,4] => [1/4,+oo] */
        IBBasicRoundDown();
        IBBasicMinI(Result)=IBBasicMinI(i1)/IBBasicMaxI(i2);
        IBBasicMaxI(Result)=IBBasicPosInfinity;
      }
      else {
        IBBasicToLargestI(Result);
      }
    }
    else if (IBBasicMaxI(i2)==0.0) {
      if (IBBasicMinI(i1)<0.0) {      /* ex. [-3,-1] / [-4,0] => [1/4,+oo] */
        IBBasicRoundDown();
        IBBasicMinI(Result)=IBBasicMaxI(i1)/IBBasicMinI(i2);
        IBBasicMaxI(Result)=IBBasicPosInfinity;
      }
      else if (IBBasicMaxI(i1)>0.0) { /* ex. [1,3] / [-4,0] => [-oo,-1/4] */
	IBBasicMinI(Result)=IBBasicNegInfinity;
        IBBasicRoundUp();
        IBBasicMaxI(Result)=IBBasicMinI(i1)/IBBasicMinI(i2);
      }
      else {
        IBBasicToLargestI(Result);
      }
    }
    else {
      IBBasicToLargestI(Result);
    }
    return;
  }

  sig = ((IBBasicMaxI(i1)<0.0) << 3) | ((IBBasicMinI(i1)>0.0) << 2) 
        | ((IBBasicMaxI(i2)<0.0) << 1) | (IBBasicMinI(i2)>0.0);

  switch( sig )
  {
    case 1: /* 0001: 0 in i1, i2>0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMinI(i1)/IBBasicMinI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMaxI(i1)/IBBasicMinI(i2);
    break;

 case 2: /* 0010: 0 in i1, i2<0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMaxI(i1)/IBBasicMaxI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMinI(i1)/IBBasicMaxI(i2);
    break;

    case 5: /* 0101: i1>0, i2>0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMinI(i1)/IBBasicMaxI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMaxI(i1)/IBBasicMinI(i2);
    break;

    case 6: /* 0110: i1>0, i2<0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMaxI(i1)/IBBasicMaxI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMinI(i1)/IBBasicMinI(i2);
    break;

    case 9: /* 1001: i1<0, i2>0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMinI(i1)/IBBasicMinI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMaxI(i1)/IBBasicMaxI(i2);
    break;

    case 10: /* 1010: i1<0, i2<0 */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMaxI(i1)/IBBasicMinI(i2);
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMinI(i1)/IBBasicMaxI(i2);
    break;
  }
}


void IBBasicDivIR(IBBasicItv Result, IBBasicItv i, IBBasicItv x)
/***************************************************************************
*  let x=[r,r], Result := i / [r,r]
*/
{
  if( IBBasicMinI(x)==0.0 )
  {
    IBBasicToLargestI(Result);
    return;
  }

  if( IBBasicMinI(x)<0.0 )
  {
    IBBasicRoundDown();
    IBBasicMinI(Result)=IBBasicMaxI(i)/IBBasicMinI(x);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=IBBasicMinI(i)/IBBasicMinI(x);
  }
  else
  {
    IBBasicRoundDown();
    IBBasicMinI(Result)=IBBasicMinI(i)/IBBasicMinI(x);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=IBBasicMaxI(i)/IBBasicMinI(x);
  }
}


void IBBasicDivIRneg(IBBasicItv Result, IBBasicItv i, IBBasicItv x)
/***************************************************************************
*  let x=[r,r], Result := i / [r,r] with r<0
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result)=IBBasicMaxI(i)/IBBasicMinI(x);
  IBBasicRoundUp();
  IBBasicMaxI(Result)=IBBasicMinI(i)/IBBasicMinI(x);
}


void IBBasicDivIRpos(IBBasicItv Result, IBBasicItv i, IBBasicItv x)
/***************************************************************************
*  let x=[r,r], Result := i / [r,r] with r>0
*/
{
  IBBasicRoundDown();
  IBBasicMinI(Result)=IBBasicMinI(i)/IBBasicMinI(x);
  IBBasicRoundUp();
  IBBasicMaxI(Result)=IBBasicMaxI(i)/IBBasicMinI(x);
}


void IBBasicDivRI(IBBasicItv Result, IBBasicItv x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r] / i
*/
{
  if (IBBasicDoubleInI(i,0.0)) {
    if (IBBasicMinI(x)!=0.0) {
      if (IBBasicMinI(i)==0.0) {
        if (IBBasicMinI(x)>0.0) {  /* 2 / [0,4] => [2/4,+oo] */
          IBBasicRoundDown();
          IBBasicMinI(Result)=IBBasicMinI(x)/IBBasicMaxI(i);
	  IBBasicMaxI(Result)=IBBasicPosInfinity;
	}
        else {                     /* -2 / [0,4] => [-oo,-2/4] */
	  IBBasicMinI(Result)=IBBasicNegInfinity;
          IBBasicRoundUp();
          IBBasicMaxI(Result)=IBBasicMinI(x)/IBBasicMaxI(i);
	}
      }
      else if (IBBasicMaxI(i)==0.0) {
        if (IBBasicMinI(x)>0.0) {  /* 2 / [-4,0] => [-oo,-2/4] */
	  IBBasicMinI(Result)=IBBasicNegInfinity;
          IBBasicRoundUp();
          IBBasicMaxI(Result)=IBBasicMinI(x)/IBBasicMinI(i);
	}
        else {                     /* -2 / [-4,0] => [2/4,+oo] */
          IBBasicRoundDown();
          IBBasicMinI(Result)=IBBasicMinI(x)/IBBasicMinI(i);
	  IBBasicMaxI(Result)=IBBasicPosInfinity;
	}
      }
      else {
        IBBasicToLargestI(Result);
      }
    }
    else {
      IBBasicToLargestI(Result);
    }
    return;
  }

  if( IBBasicMinI(x)==0.0 ) {
    IBBasicMinI(Result)=0.0;
    IBBasicMaxI(Result)=0.0;
  }
  else if (IBBasicMinI(x)>=0.0) {
    IBBasicRoundDown();
    IBBasicMinI(Result)=IBBasicMinI(x)/IBBasicMaxI(i);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=IBBasicMinI(x)/IBBasicMinI(i);
  }
  else {
    IBBasicRoundDown();
    IBBasicMinI(Result)=IBBasicMinI(x)/IBBasicMinI(i);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=IBBasicMinI(x)/IBBasicMaxI(i);
  }
}


void IBBasicDivRnegI(IBBasicItv Result, IBBasicItv x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r] / i with r<0
*/
{
  if (IBBasicDoubleInI(i,0.0)) {
    if (IBBasicMinI(i)==0.0) {       /* -2 / [0,4] => [-oo,-2/4] */
      IBBasicMinI(Result)=IBBasicNegInfinity;
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMinI(x)/IBBasicMaxI(i);
    }
    else if (IBBasicMaxI(i)==0.0) {  /* -2 / [-4,0] => [2/4,+oo] */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMinI(x)/IBBasicMinI(i);
      IBBasicMaxI(Result)=IBBasicPosInfinity;
    }
    else {
      IBBasicToLargestI(Result);
    }
    return;
  }

  IBBasicRoundDown();
  IBBasicMinI(Result)=IBBasicMinI(x)/IBBasicMinI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result)=IBBasicMinI(x)/IBBasicMaxI(i);
}


void IBBasicDivRposI(IBBasicItv Result, IBBasicItv x, IBBasicItv i)
/***************************************************************************
*  let x=[r,r], Result := [r,r] / i with r >=0
*/
{
  if (IBBasicDoubleInI(i,0.0)) {
    if (IBBasicMinI(i)==0.0) {       /* 2 / [0,4] => [2/4,+oo] */
      IBBasicRoundDown();
      IBBasicMinI(Result)=IBBasicMinI(x)/IBBasicMaxI(i);
      IBBasicMaxI(Result)=IBBasicPosInfinity;
    }
    else if (IBBasicMaxI(i)==0.0) {  /* 2 / [-4,0] => [-oo,-2/4] */
      IBBasicMinI(Result)=IBBasicNegInfinity;
      IBBasicRoundUp();
      IBBasicMaxI(Result)=IBBasicMinI(x)/IBBasicMinI(i);
    }
    else {
      IBBasicToLargestI(Result);
    }
    return;
  }

  IBBasicRoundDown();
  IBBasicMinI(Result)=IBBasicMinI(x)/IBBasicMaxI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result)=IBBasicMinI(x)/IBBasicMinI(i);
}


void IBBasicDivRposIinternal(IBBasicItv Result, double x, IBBasicItv i)
/***************************************************************************
*  Result = [x,x] / i with x >=0
*/
{
  if (IBBasicDoubleInI(i,0.0)) {
    if (IBBasicMinI(i)==0.0) {       /* 2 / [0,4] => [2/4,+oo] */
      IBBasicRoundDown();
      IBBasicMinI(Result)=x/IBBasicMaxI(i);
      IBBasicMaxI(Result)=IBBasicPosInfinity;
    }
    else if (IBBasicMaxI(i)==0.0) {  /* 2 / [-4,0] => [-oo,-2/4] */
      IBBasicMinI(Result)=IBBasicNegInfinity;
      IBBasicRoundUp();
      IBBasicMaxI(Result)=x/IBBasicMinI(i);
    }
    else {
      IBBasicToLargestI(Result);
    }
    return;
  }

  IBBasicRoundDown();
  IBBasicMinI(Result)=x/IBBasicMaxI(i);
  IBBasicRoundUp();
  IBBasicMaxI(Result)=x/IBBasicMinI(i);
}


void IBBasicSqrI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := i^2
*/
{
  if( IBBasicMinI(i)>=0.0 ) {
    IBBasicRoundDown();
    IBBasicMinI(Result)=IBBasicMinI(i)*IBBasicMinI(i);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=IBBasicMaxI(i)*IBBasicMaxI(i);
  }
  else if( IBBasicMaxI(i)<=0.0 ) {
    IBBasicRoundDown();
    IBBasicMinI(Result)=IBBasicMaxI(i)*IBBasicMaxI(i);
    IBBasicRoundUp();
    IBBasicMaxI(Result)=IBBasicMinI(i)*IBBasicMinI(i);
  }
  else if ((-IBBasicMinI(i))<IBBasicMaxI(i)) {
    IBBasicMinI(Result) = 0.0;
    IBBasicRoundUp();
    IBBasicMaxI(Result) = IBBasicMaxI(i)*IBBasicMaxI(i);
  }
  else {
    IBBasicMinI(Result) = 0.0;
    IBBasicRoundUp();
    IBBasicMaxI(Result) = IBBasicMinI(i)*IBBasicMinI(i);
  }
}


void IBBasicSqrtI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := sqrt(i)
*/
{
  double x;

  if( IBBasicMaxI(i)<0.0 ) {
    IBBasicSetEmptyI(Result);
  }
  else if( IBBasicMinI(i)>=0.0 ) {
    IBBasicRoundDown();
    IBBasicMinI(Result) = sqrt(IBBasicMinI(i));
    IBBasicRoundUp();
    IBBasicMaxI(Result) = sqrt(IBBasicMaxI(i));
  }
  else {
    IBBasicMinI(Result) = 0.0;
    IBBasicRoundUp();
    IBBasicMaxI(Result) = sqrt(IBBasicMaxI(i));
  }
}


void IBBasicPowI(IBBasicItv Result, IBBasicItv i, IBBasicItv n)
/***************************************************************************
*  let n=[j,j], Result := i^j
*/
{
  IBBasicItv z;

  if( IBBasicIsEven((int)IBBasicMinI(n)) )   /* n even */
  {
    IBBasicAbsI(z,i);
    if (IBBasicMinI(z)==0.0) {
      IBBasicMinI(Result) = 0.0;
    }
    else {
      IBBasicRoundNear();
      IBBasicMinI(Result) = IBBasicPrevDouble(IBBasicPrevDouble(pow(IBBasicMinI(z),IBBasicMinI(n))));
    }
    if (IBBasicMaxI(z)==IBBasicPosInfinity) {
      IBBasicMaxI(Result) = IBBasicPosInfinity;
    }
    else {
      IBBasicRoundNear();
      IBBasicMaxI(Result) = IBBasicNextDouble(IBBasicNextDouble(pow(IBBasicMaxI(z),IBBasicMinI(n))));
    }
  }
  else {
    if (IBBasicMinI(i)==IBBasicNegInfinity) {
      IBBasicMinI(Result) = IBBasicNegInfinity;
    }
    else {
      IBBasicRoundNear();
      IBBasicMinI(Result) = IBBasicPrevDouble(IBBasicPrevDouble(pow(IBBasicMinI(i),IBBasicMinI(n))));
    }
    if (IBBasicMaxI(i)==IBBasicPosInfinity) {
      IBBasicMaxI(Result) = IBBasicPosInfinity;
    }
    else {
      IBBasicRoundNear();
      IBBasicMaxI(Result) = IBBasicNextDouble(IBBasicNextDouble(pow(IBBasicMaxI(i),IBBasicMinI(n))));
    }
  } 
}


void IBBasicPowIinternal(IBBasicItv Result, IBBasicItv i, int n)
/***************************************************************************
*  Result := i^n
*/
{
  IBBasicItv z;

  if( IBBasicIsEven(n) )   /* n even */
  {
    IBBasicAbsI(z,i);
    if (IBBasicMinI(z)==0.0) {
      IBBasicMinI(Result) = 0.0;
    }
    else {
      IBBasicRoundNear();
      IBBasicMinI(Result) = IBBasicPrevDouble(IBBasicPrevDouble(pow(IBBasicMinI(z),n)));
    }
    if (IBBasicMaxI(z)==IBBasicPosInfinity) {
      IBBasicMaxI(Result) = IBBasicPosInfinity;
    }
    else {
      IBBasicRoundNear();
      IBBasicMaxI(Result) = IBBasicNextDouble(IBBasicNextDouble(pow(IBBasicMaxI(z),n)));
    }
  }
  else {
    if (IBBasicMinI(i)==IBBasicNegInfinity) {
      IBBasicMinI(Result) = IBBasicNegInfinity;
    }
    else {
      IBBasicRoundNear();
      IBBasicMinI(Result) = IBBasicPrevDouble(IBBasicPrevDouble(pow(IBBasicMinI(i),n)));
    }
    if (IBBasicMaxI(i)==IBBasicPosInfinity) {
      IBBasicMaxI(Result) = IBBasicPosInfinity;
    }
    else {
      IBBasicRoundNear();
      IBBasicMaxI(Result) =  IBBasicNextDouble(IBBasicNextDouble(pow(IBBasicMaxI(i),n)));
    }
  }
}


void IBBasicExpI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := exp(i)
*/
{
  double l, u;
  IBBasicRoundNear();
  l = exp(IBBasicMinI(i));
  u = exp(IBBasicMaxI(i));

  IBBasicRoundDown();
  IBBasicMinI(Result) = IBBasicMax(0.0,l-IBBasicEpsilon);

  IBBasicRoundUp();
  IBBasicMaxI(Result) = IBBasicMax(0.0,u+IBBasicEpsilon);

  if( IBBasicMinI(Result)==IBBasicPosInfinity ) {
    IBBasicMinI(Result) = IBBasicMaxDouble;
  }

  if( IBBasicMaxI(Result)==IBBasicNegInfinity ) {
    IBBasicMaxI(Result) = IBBasicMinDouble;
  }
}


void IBBasicLogI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := log(i)
*/
{
  if( IBBasicMaxI(i)<=0.0 ) {
    IBBasicSetEmptyI(Result);
  }
  else if( IBBasicMinI(i)<=0.0 )
  {
    IBBasicMinI(Result) = IBBasicNegInfinity;
    IBBasicRoundNear();
    IBBasicMaxI(Result) = log(IBBasicMaxI(i));
    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
  else
  {
    IBBasicRoundNear();
    IBBasicMinI(Result) = log(IBBasicMinI(i));
    IBBasicMaxI(Result) = log(IBBasicMaxI(i));
    IBBasicRoundDown();
    IBBasicMinI(Result) -= IBBasicEpsilon;
    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
}


void IBBasicMinimumII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2)
/***************************************************************************
*  Result := min(i1,i2)
*/
{
  IBBasicMinI(Result) = IBBasicMin(IBBasicMinI(i1),IBBasicMinI(i2));
  IBBasicMaxI(Result) = IBBasicMin(IBBasicMaxI(i1),IBBasicMaxI(i2));
}


void IBBasicMaximumII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2)
/***************************************************************************
*  Result := max(i1,i2)
*/
{
  IBBasicMinI(Result) = IBBasicMax(IBBasicMinI(i1),IBBasicMinI(i2));
  IBBasicMaxI(Result) = IBBasicMax(IBBasicMaxI(i1),IBBasicMaxI(i2));
}


void IBBasicSinI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := sin(i)
*/
{
  double k, leftB, rightB;
  IBBasicItv offset;

  /*
  printf("sin( ");
  IBBasicWriteI(stdout,i,8,1,0);
  printf(" ) = ");
  */

  IBBasicRoundUp();
  if( IBBasicWidthI(i)>=IBBasicMinI(IBBasicItvConst_2_Pi) )
  {  /* use of minI(2pi): security */
    IBBasicSetI(Result,-1.0,1.0);
  }
  else
  {
    if( IBBasicMinI(i)<0.0 )
    {
      IBBasicRoundDown();
      k = floor(IBBasicMinI(i)/IBBasicMaxI(IBBasicItvConst_2_Pi));
      IBBasicMulRIinternal(offset,k,IBBasicItvConst_2_Pi);
      IBBasicRoundDown();
      leftB = IBBasicMinI(i) - IBBasicMaxI(offset);
      if( leftB<0.0 ) leftB = 0.0;
      IBBasicRoundUp();
      rightB = IBBasicMaxI(i) - IBBasicMinI(offset);
    }
    else if( IBBasicMinI(i)>IBBasicMinI(IBBasicItvConst_2_Pi) )
    {
      IBBasicRoundUp();
      k = floor(IBBasicMinI(i)/IBBasicMinI(IBBasicItvConst_2_Pi));
      IBBasicMulRposIinternal(offset,k,IBBasicItvConst_2_Pi);
      IBBasicRoundDown();
      leftB = IBBasicMinI(i) - IBBasicMaxI(offset);
      if( leftB<0.0 ) leftB = 0.0;
      IBBasicRoundUp();
      rightB = IBBasicMaxI(i) - IBBasicMinI(offset);
    }
    else
    {
      leftB = IBBasicMinI(i);
      rightB = IBBasicMaxI(i);
    }

    if( ((leftB<=IBBasicMinI(IBBasicItvConst_3_Pi_2)) &&
         (rightB>=IBBasicMaxI(IBBasicItvConst_3_Pi_2))) ||
         ((leftB<=IBBasicMinI(IBBasicItvConst_7_Pi_2)) &&
	  (rightB>=IBBasicMaxI(IBBasicItvConst_7_Pi_2))) )
    {
      IBBasicMinI(Result) = -1.0;
    }
    else
    {
      IBBasicRoundNear();
      IBBasicMinI(Result) = IBBasicMin(sin(leftB),sin(rightB));
      IBBasicRoundDown();
      IBBasicMinI(Result) = IBBasicMax(-1.0,IBBasicMinI(Result)-IBBasicEpsilon);
    }

    if( ((leftB<=IBBasicMinI(IBBasicItvConst_1_Pi_2)) &&
         (rightB>=IBBasicMaxI(IBBasicItvConst_1_Pi_2))) ||
         ((leftB<=IBBasicMinI(IBBasicItvConst_5_Pi_2)) &&
	  (rightB>=IBBasicMaxI(IBBasicItvConst_5_Pi_2))) )
    {
      IBBasicMaxI(Result) = 1.0;
    }
    else
    {
      IBBasicRoundNear();
      IBBasicMaxI(Result) = IBBasicMax(sin(leftB),sin(rightB));
      IBBasicRoundUp();
      IBBasicMaxI(Result) = IBBasicMin(1.0,IBBasicMaxI(Result)+IBBasicEpsilon);
    }
  }

  /*
  IBBasicWriteI(stdout,Result,8,1,0);
  printf("\n");
 */
}


void IBBasicCosI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := cos(i)
*/
{
  IBBasicItv j;
  IBBasicAddII(j,i,IBBasicItvConst_1_Pi_2);  /* j := i + pi/2 */
  IBBasicSinI(Result,j,useless);             /* cos(i) = sin(i + pi/2) */
}


void IBBasicTanI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := tan(i)
*/
{
  /*
  IBBasicItv j, k;
  IBBasicSinI(j,i,useless);
  IBBasicCosI(k,i,useless);
  IBBasicDivII(Result,j,k);
  */

  double k, leftB, rightB;
  IBBasicItv offset;

  
  //printf("tan( ");
  //IBBasicWriteI(stdout,i,8,1,0);
  //printf(" ) => ");
  


  IBBasicRoundUp();
  if( IBBasicWidthI(i)>=IBBasicMinI(IBBasicItvConstPi) )
    {  /* use of minI(pi): security */
    IBBasicSetI(Result,IBBasicNegInfinity,IBBasicPosInfinity);
  }
  else
  {
    if( IBBasicMinI(i)<=(-IBBasicMaxI(IBBasicItvConst_1_Pi_2)) )
    {
      IBBasicRoundDown();
      k = floor(IBBasicMinI(i)/IBBasicMaxI(IBBasicItvConstPi));
      IBBasicMulRIinternal(offset,k,IBBasicItvConstPi);
      IBBasicRoundDown();
      leftB = IBBasicMinI(i) - IBBasicMaxI(offset);
      if( leftB<0.0 ) leftB = 0.0;
      IBBasicRoundUp();
      rightB = IBBasicMaxI(i) - IBBasicMinI(offset);
    }
    else if( IBBasicMinI(i)>=IBBasicMaxI(IBBasicItvConst_1_Pi_2) )
    {
      IBBasicRoundUp();
      k = floor(IBBasicMinI(i)/IBBasicMinI(IBBasicItvConstPi));
      IBBasicMulRposIinternal(offset,k,IBBasicItvConstPi);
      IBBasicRoundDown();
      leftB = IBBasicMinI(i) - IBBasicMaxI(offset);
      if( leftB<0.0 ) leftB = 0.0;
      IBBasicRoundUp();
      rightB = IBBasicMaxI(i) - IBBasicMinI(offset);
    }
    else
    {
      leftB = IBBasicMinI(i);
      rightB = IBBasicMaxI(i);
    }

    if( ((leftB<=IBBasicMinI(IBBasicItvConst_1_Pi_2)) &&
         (rightB>=IBBasicMaxI(IBBasicItvConst_1_Pi_2))) ||
         ((leftB<=IBBasicMinI(IBBasicItvConst_3_Pi_2)) &&
	  (rightB>=IBBasicMaxI(IBBasicItvConst_3_Pi_2))) )
    {
      IBBasicSetI(Result,IBBasicNegInfinity,IBBasicPosInfinity);
    }
    else
    {
      IBBasicRoundNear();
      IBBasicMinI(Result) = tan(leftB);
      IBBasicMaxI(Result) = tan(rightB);
      IBBasicRoundDown();
      IBBasicMinI(Result) -= IBBasicEpsilon;
      IBBasicRoundUp();
      IBBasicMaxI(Result) += IBBasicEpsilon;
    }
  }


  //  IBBasicWriteI(stdout,Result,8,1,0);
  //printf("\n");


}


void IBBasicCoshI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := cosh(i) = 0.5 * ( exp(i) + exp(-i) )
*  decreasing function in (-oo,0]
*  increasing function in [0,+oo)
*/
{
  double l, u;

  if( IBBasicMinI(i)>=0.0 ) {
    IBBasicRoundNear();
    IBBasicMinI(Result) = cosh(IBBasicMinI(i));
    IBBasicMaxI(Result) = cosh(IBBasicMaxI(i));
    IBBasicRoundDown();
    IBBasicMinI(Result) -= IBBasicEpsilon;
    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
  else if( IBBasicMaxI(i)<=0.0 ) {
    IBBasicRoundNear();
    IBBasicMinI(Result) = cosh(IBBasicMaxI(i));
    IBBasicMaxI(Result) = cosh(IBBasicMinI(i));
    IBBasicRoundDown();
    IBBasicMinI(Result) -= IBBasicEpsilon;
    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
  else if ((-IBBasicMinI(i))<IBBasicMaxI(i)) {
    IBBasicMinI(Result) = 1.0;
    IBBasicRoundNear();
    IBBasicMaxI(Result) = cosh(IBBasicMaxI(i));
    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
  else {
    IBBasicMinI(Result) = 1.0;
    IBBasicRoundNear();
    IBBasicMaxI(Result) = cosh(IBBasicMinI(i));
    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
}


void IBBasicSinhI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := sinh(i) = 0.5 * ( exp(i) - exp(-i) )
*  increasing function in (-oo,+oo)
*/
{
  IBBasicRoundNear();
  IBBasicMinI(Result) = sinh(IBBasicMinI(i));
  IBBasicMaxI(Result) = sinh(IBBasicMaxI(i));

  IBBasicRoundDown();
  IBBasicMinI(Result) -= IBBasicEpsilon;

  IBBasicRoundUp();
  IBBasicMaxI(Result) += IBBasicEpsilon;
}


void IBBasicTanhI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := tanh(i) = sinh(i) / cosh(i)
*  increasing function in (-oo,+oo)
*/
{
  IBBasicRoundNear();
  IBBasicMinI(Result) = tanh(IBBasicMinI(i));
  IBBasicMaxI(Result) = tanh(IBBasicMaxI(i));

  IBBasicRoundDown();
  IBBasicMinI(Result) -= IBBasicEpsilon;

  IBBasicRoundUp();
  IBBasicMaxI(Result) += IBBasicEpsilon;
}


void IBBasicAsinI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := asin(i)      (increasing function in [-1,1])
*/
{
  IBBasicItv j;
  IBBasicMinI(j) = IBBasicMax(-1.0,IBBasicMinI(i));
  IBBasicMaxI(j) = IBBasicMin(1.0,IBBasicMaxI(i));

  if( IBBasicEmptyI(j) ) {
    IBBasicSetEmptyI(Result);
  }
  else {
    IBBasicRoundNear();
    IBBasicMinI(Result) = asin(IBBasicMinI(j));
    IBBasicMaxI(Result) = asin(IBBasicMaxI(j));
    IBBasicRoundDown();
    IBBasicMinI(Result) -= IBBasicEpsilon;
    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
}


void IBBasicAcosI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := acos(i)   (decreasing function in [-1,1])
*/
{
  IBBasicItv j;
  IBBasicMinI(j) = IBBasicMax(-1.0,IBBasicMinI(i));
  IBBasicMaxI(j) = IBBasicMin(1.0,IBBasicMaxI(i));

  if( IBBasicEmptyI(j) ) {
    IBBasicSetEmptyI(Result);
  }
  else {
    IBBasicRoundNear();
    IBBasicMinI(Result) = acos(IBBasicMaxI(j));
    IBBasicMaxI(Result) = acos(IBBasicMinI(j));
    IBBasicRoundDown();
    IBBasicMinI(Result) -= IBBasicEpsilon;
    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
}


void IBBasicAtanI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := atan(i)       (increasing function in (-oo,+oo))
*/
{
  IBBasicRoundNear();
  IBBasicMinI(Result) = atan(IBBasicMinI(i));
  IBBasicMaxI(Result) = atan(IBBasicMaxI(i));
  IBBasicRoundDown();
  IBBasicMinI(Result) -= IBBasicEpsilon;
  IBBasicRoundUp();
  IBBasicMaxI(Result) += IBBasicEpsilon;
}


void IBBasicAsinhI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := asinh(i) = log(i + sqrt(1+i^2))
*  increasing function in (-oo,+oo)
*/
{
  IBBasicRoundNear();
  IBBasicMinI(Result) = asinh(IBBasicMinI(i));
  IBBasicMaxI(Result) = asinh(IBBasicMaxI(i));

  IBBasicRoundDown();
  IBBasicMinI(Result) -= IBBasicEpsilon;

  IBBasicRoundUp();
  IBBasicMaxI(Result) += IBBasicEpsilon;
}


void IBBasicAcoshI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := acosh(i) = log(i + sqrt(i^2-1))
*  increasing function in [1,+oo)
*/
{
  if( IBBasicMaxI(i)<1.0 ) {
    IBBasicSetEmptyI(Result);
  }
  else if( IBBasicMinI(i)>1.0 ) {
    IBBasicRoundNear();
    IBBasicMinI(Result) = acosh(IBBasicMinI(i));
    IBBasicMaxI(Result) = acosh(IBBasicMaxI(i));

    IBBasicRoundDown();
    IBBasicMinI(Result) -= IBBasicEpsilon;

    IBBasicRoundUp();
    IBBasicMaxI(Result) += IBBasicEpsilon;
  }
  else {   /* i contains 1.0 */
    IBBasicMinI(Result) = 0.0;
    if (IBBasicMaxI(i)==1.0) {
      IBBasicMaxI(Result) = 0.0;
    }
    else {
      IBBasicRoundNear();
      IBBasicMaxI(Result) = acosh(IBBasicMaxI(i));
      IBBasicRoundUp();
      IBBasicMaxI(Result) += IBBasicEpsilon;      
    }
  }
}


void IBBasicAtanhI(IBBasicItv Result, IBBasicItv i, IBBasicItv useless)
/***************************************************************************
*  Result := atanh(i) = 0.5 * log( (1+i) / (1-i) )
*  increasing function in (-1,1)
*/
{
  if ( (IBBasicMaxI(i)<=-1.0) || (IBBasicMinI(i)>=1.0) ) {
    IBBasicSetEmptyI(Result);
  }
  else {
    if (IBBasicMinI(i)<=-1.0) {
      IBBasicMinI(Result) = IBBasicNegInfinity;
    }
    else {
      IBBasicRoundNear();
      IBBasicMinI(Result) = atanh(IBBasicMinI(i));
      IBBasicRoundDown();
      IBBasicMinI(Result) -= IBBasicEpsilon;
    }

    if (IBBasicMaxI(i)>=1.0) {
      IBBasicMaxI(Result) = IBBasicPosInfinity;
    }
    else {
      IBBasicRoundNear();
      IBBasicMaxI(Result) = atanh(IBBasicMaxI(i));
      IBBasicRoundUp();
      IBBasicMaxI(Result) += IBBasicEpsilon;      
    }
  }
}


int IBBasicExtendedDivisionInterII(IBBasicItv Result, IBBasicItv num, IBBasicItv den)
/***************************************************************************
*  Result := Result inter (num/den)
*  Returns 1 if Result is not modified, 0 otherwise
*
*  let num=[a,b] and  den = [c,d]
*/
{
  double min, max;
  int notmodified = 1;
  IBBasicItv itv;

  if( !IBBasicDoubleInI(den,0.0) )  /* standard division */
  {
    IBBasicDivII(itv,num,den);
    if( IBBasicMinI(itv)>IBBasicMinI(Result) )
    {
      notmodified = 0;
      IBBasicMinI(Result) = IBBasicMinI(itv);
    }
    if( IBBasicMaxI(itv)<IBBasicMaxI(Result) )
    {
      notmodified = 0;
      IBBasicMaxI(Result) = IBBasicMaxI(itv);
    }
  }
  else if( IBBasicIsDoubleI(den) )        /* it is necessary 0.0 */
  {
      /* num/den = [-oo,+oo] => not modified */
  }
  else if( IBBasicMinI(num)>=0.0 )
  {
    if( IBBasicMinI(den)==0.0 )      /* b>=a>=0.0, c=0 */
    {
      /* num/den = [(a/d),+oo] */
      /* only inf bound can be modified */
      IBBasicRoundDown();
      min = IBBasicMinI(num)/IBBasicMaxI(den);    /* a/d */

      if( min > IBBasicMinI(Result) )
      {
        IBBasicMinI(Result) = min;
        notmodified = 0;
      }
    }
    else if( IBBasicMaxI(den)==0.0 ) /* b>=a>=0.0, d=0 */
    {
      /* num/den = [-oo,(a/c)] */
      /* only sup bound can be modified */
      IBBasicRoundUp();
      max = IBBasicMinI(num)/IBBasicMinI(den);  /* a/c */
 
      if( max < IBBasicMaxI(Result) )
      {
        IBBasicMaxI(Result) = max;
        notmodified = 0;
      }
    }
    else                          /* b>=a>=0.0, d>0>c */
    {
      /* num/den = [-oo,(a/c)] union [(b/c),+oo] */

      if ((IBBasicMaxI(num)==IBBasicPosInfinity) && (IBBasicMinI(den)==IBBasicNegInfinity))
      {
        /* num/den = [-oo,+oo] => not modified */
      }
      else
      {
        IBBasicRoundDown();
        min = IBBasicMaxI(num)/IBBasicMinI(den);  /* b/c */
        IBBasicRoundUp();
        max = IBBasicMinI(num)/IBBasicMinI(den);  /* a/c */

        /* num/den = [-oo,max] union [min,+oo] */
        /* Result is modified if at least one of its bounds is in ]max..min[ */

        if( (IBBasicMinI(Result)>max) && (IBBasicMinI(Result)<min) )
        {
          if( IBBasicMaxI(Result)<min )  /* Result in ]max..min[ */
  	  {
            IBBasicSetEmptyI(Result);
  	  }
          else
	  {
            IBBasicMinI(Result) = min;
	  }
          notmodified = 0;
        }
        else if( (IBBasicMaxI(Result)>max) && (IBBasicMaxI(Result)<min) )
        {
          if( IBBasicMinI(Result)>max )  /* Result in ]max..min[ */	
          {
            IBBasicSetEmptyI(Result);
	  }
          else
	  {
            IBBasicMaxI(Result) = max;
	  }
          notmodified = 0;
        }
      }
    }
  }
  else if( IBBasicMaxI(num)<=0.0 )
  {
    if( IBBasicMinI(den)==0.0 )      /* 0.0>=b>=a, c=0 */
    {
      /* num/den = [-oo,(b/d)] */
      /* only sup bound is modified */
      IBBasicRoundUp();
      max = IBBasicMaxI(num)/IBBasicMaxI(den);  /* b/d */

      if( max < IBBasicMaxI(Result) )
      {
        IBBasicMaxI(Result) = max;
        notmodified = 0;
      }
    }
    else if( IBBasicMaxI(den)==0.0 ) /* 0.0>=b>=a, d=0 */
    {
      /* num/den = [(b/c),+oo] */
      /* only inf bound can be modified */
      IBBasicRoundDown();
      min = IBBasicMaxI(num)/IBBasicMinI(den);    /* b/c */

      if( min > IBBasicMinI(Result) )
      {
        IBBasicMinI(Result) = min;
        notmodified = 0;
      }
    }
    else                          /* 0.0>=b>=a, d>0>c */
    {
      /* num/den = [-oo,(b/d)] union [(b/c),+oo] */
      IBBasicRoundDown();
      min = IBBasicMaxI(num)/IBBasicMinI(den);  /* b/c */
      IBBasicRoundUp();
      max = IBBasicMaxI(num)/IBBasicMaxI(den);  /* b/d */

      /* num/den = [-oo,max] union [min,+oo] */
      /* Result is modified if at least one of its bounds is in ]max..min[ */

      if( (IBBasicMinI(Result)>max) && (IBBasicMinI(Result)<min) )
      {
        if( IBBasicMaxI(Result)<min )  /* Result in ]max..min[ */
	{
          IBBasicSetEmptyI(Result);
	}
        else
	{
          IBBasicMinI(Result) = min;
	}
        notmodified = 0;
      }
      else if( (IBBasicMaxI(Result)>max) && (IBBasicMaxI(Result)<min) )
      {
        if( IBBasicMinI(Result)>max )  /* Result in ]max..min[ */	
        {
          IBBasicSetEmptyI(Result);
	}
        else
	{
          IBBasicMaxI(Result) = max;
	}
        notmodified = 0;
      }
    }
  }
  else            /* a < 0 < b */
  {
    /* num/den = [-oo,+oo] => not modified */
  }

  return( notmodified );
}


int IBBasicExtendedDivisionII(IBBasicItv Result1, IBBasicItv Result2,
                              IBBasicItv num, IBBasicItv den)
/***************************************************************************
*  Computes (num/den) using the extended division
*  Returns:
*        1 if num/den = Result1
*        2 if num/den = Result1 union Result2
*
*  let num=[a,b] and  den = [c,d]
*/
{
  double min, max;
  IBBasicItv itv;

  if( !IBBasicDoubleInI(den,0.0) )  /* standard division */
  {
    IBBasicDivII(Result1,num,den);
    return 1;
  }
  else if( IBBasicIsDoubleI(den) )        /* den=0.0 */
  {
    /* num/den = [-oo,+oo] */
    IBBasicToLargestI(Result1);
    return 1;
  }
  else if( IBBasicMinI(num)>=0.0 )
  {
    if( IBBasicMinI(den)==0.0 )      /* b>=a>=0.0, c=0 */
    {
      /* num/den = [(a/d),+oo] */
      IBBasicRoundDown();
      IBBasicMinI(Result1) = IBBasicMinI(num)/IBBasicMaxI(den);
      IBBasicMaxI(Result1) = IBBasicPosInfinity;
      return 1;
    }
    else if( IBBasicMaxI(den)==0.0 ) /* b>=a>=0.0, d=0 */
    {
      /* num/den = [-oo,(a/c)] */
      IBBasicMinI(Result1) = IBBasicNegInfinity;
      IBBasicRoundUp();
      IBBasicMaxI(Result1) = IBBasicMinI(num)/IBBasicMinI(den);
      return 1;
    }
    else                          /* b>=a>=0.0, d>0>c */
    {
      /* num/den = [-oo,(a/c)] union [(b/c),+oo] */

      if ((IBBasicMaxI(num)==IBBasicPosInfinity) && (IBBasicMinI(den)==IBBasicNegInfinity))
      {
        /* num/den = [-oo,+oo] */
        IBBasicToLargestI(Result1);
        return 1;
      }
      else
      {
        /* num/den = [-oo,a/c] union [b/c,+oo] */
	IBBasicMinI(Result1) = IBBasicNegInfinity;
        IBBasicRoundUp();
        IBBasicMaxI(Result1) = IBBasicMinI(num)/IBBasicMinI(den);
        IBBasicRoundDown();
        IBBasicMinI(Result2) = IBBasicMaxI(num)/IBBasicMinI(den);
        IBBasicMaxI(Result2) = IBBasicPosInfinity;
	return 2;
      }
    }
  }
  else if( IBBasicMaxI(num)<=0.0 )
  {
    if( IBBasicMinI(den)==0.0 )      /* 0.0>=b>=a, c=0 */
    {
      /* num/den = [-oo,(b/d)] */
      IBBasicMinI(Result1) = IBBasicNegInfinity;
      IBBasicRoundUp();
      IBBasicMaxI(Result1) = IBBasicMaxI(num)/IBBasicMaxI(den);
      return 1;
    }
    else if( IBBasicMaxI(den)==0.0 ) /* 0.0>=b>=a, d=0 */
    {
      /* num/den = [(b/c),+oo] */
      IBBasicRoundDown();
      IBBasicMinI(Result1) = IBBasicMaxI(num)/IBBasicMinI(den);
      IBBasicMaxI(Result1) = IBBasicPosInfinity;
      return 1;
    }
    else                          /* 0.0>=b>=a, d>0>c */
    {
      /* num/den = [-oo,(b/d)] union [(b/c),+oo] */
      IBBasicMinI(Result1) = IBBasicNegInfinity;
      IBBasicRoundUp();
      IBBasicMaxI(Result1) = IBBasicMaxI(num)/IBBasicMaxI(den);
      IBBasicRoundDown();
      IBBasicMinI(Result2) = IBBasicMaxI(num)/IBBasicMinI(den);
      IBBasicMaxI(Result2) = IBBasicPosInfinity;
      return 2;
    }
  }
  else            /* a < 0 < b */
  {
    /* num/den = [-oo,+oo] */
    IBBasicToLargestI(Result1);
    return 1;
  }
}


int IBBasicNthRootRelI(IBBasicItv Result1, IBBasicItv Result2,
                       IBBasicItv i, IBBasicItv n)
/***************************************************************************
* Computes the n-root of i
* Returns:
*           0 if the result is empty
*           1 if the n-root of i = Result1
*           2 if the n-root of i = Result1 union Result2
*/
{
  if( IBBasicIsOdd(IBBasicMinI((int)n)) )
  {
    IBBasicMinI(Result1) = IBBasicNthRoot(IBBasicMinI(i),IBBasicMinI(n),1,IBBasicEpsilon);
    IBBasicMaxI(Result1) = IBBasicNthRoot(IBBasicMaxI(i),IBBasicMinI(n),2,IBBasicEpsilon);
    return 1;
  }
  else
  {
    if( IBBasicDoubleInI(i,0.0) )
    {
      /* ex. y = x^2, y in [-2,5], then x in [-sqrt(5),sqrt(5)] */

      IBBasicMaxI(Result1) = IBBasicNthRoot(IBBasicMaxI(i),IBBasicMinI(n),2,IBBasicEpsilon);
      IBBasicMinI(Result1) = -IBBasicMaxI(Result1);
      return 1;
    }
    else if( IBBasicMinI(i)>0.0 )
    {
      /* ex. y = x^2, y in [1,2], then x in [-sqrt(2),-1] union [1,sqrt(2)] */

      IBBasicMinI(Result2) = IBBasicNthRoot(IBBasicMinI(i),IBBasicMinI(n),1,IBBasicEpsilon);
      IBBasicMaxI(Result2) = IBBasicNthRoot(IBBasicMaxI(i),IBBasicMinI(n),2,IBBasicEpsilon);
      IBBasicMinI(Result1) = -IBBasicMaxI(Result2);
      IBBasicMaxI(Result1) = -IBBasicMinI(Result2);
      return 2;
    }
    else
    {
      IBBasicSetEmptyI(Result1);
      return 0;
    }
  }
}


int IBBasicCoshRelI(IBBasicItv Result1, IBBasicItv Result2, IBBasicItv i)
/***************************************************************************
* Computes the relational hyperbolic cosine of i
* Returns:
*           0 if the result is empty
*           1 if the result = Result1
*           2 if the result = Result1 union Result2
*/
{
  if( IBBasicDoubleInI(i,1.0) )
  {
    /* ex. y = cosh(x), y in [0,5], then x in [-acosh(5),acosh(5)] */
    IBBasicRoundUp();
    IBBasicMaxI(Result1) = acosh(IBBasicMaxI(i));
    IBBasicMaxI(Result1) += IBBasicEpsilon;
    IBBasicMinI(Result1) = -IBBasicMaxI(Result1);
    return 1;
  }
  else if( IBBasicMinI(i)>1.0 )
  {
    /* ex. y = cosh(x), y in [2,3], then x in [-acosh(3),-acosh(2)] union [acosh(2),acosh(3)] */

    IBBasicRoundDown();
    IBBasicMinI(Result2) = acosh(IBBasicMinI(i));
    IBBasicMinI(Result2) -= IBBasicEpsilon;

    IBBasicRoundUp();
    IBBasicMaxI(Result2) = acosh(IBBasicMaxI(i));
    IBBasicMaxI(Result2) += IBBasicEpsilon;

    IBBasicMinI(Result1) = -IBBasicMaxI(Result2);
    IBBasicMaxI(Result1) = -IBBasicMinI(Result2);
    return 2;
  }
  else
  {
    IBBasicSetEmptyI(Result1);
    return 0;
  }
}


int IBBasicSinRelI(IBBasicItv Result1, IBBasicItv Result2, IBBasicItv Result3, IBBasicItv i)
/***************************************************************************
* Computes the relational sine of i => result in [-pi, +pi]
* Returns:
*           0 if the result is empty
*           1 if the result = Result1
*           2 if the result = Result1 union Result2
*           3 if the result = Result1 union Result2 union Result3
*/
{
                                             /*           -1        0        +1               */
  if ( (IBBasicMinI(i)>1.0) ||               /*            |        |         |    |---|      */
       (IBBasicMaxI(i)<-1.0) ) {             /*  |---|     |        |         |               */
    return 0;   /* empty result */
  }
  else if (IBBasicMinI(i)<=-1.0) {           /*  |---                                         */
    if (IBBasicMaxI(i)<=0.0) {               /*  |-------------|                              */
      IBBasicRoundNear();
      IBBasicMaxI(Result1) = asin(IBBasicMaxI(i));
      IBBasicRoundUp();
      IBBasicMaxI(Result1) += IBBasicEpsilon;                       /* asin(i.sup) rounded up */

      IBBasicMinI(Result1) = -IBBasicMaxI(IBBasicItvConstPi);             /* -pi rounded down */
      IBBasicRoundDown();
      IBBasicMinI(Result1) -= IBBasicMaxI(Result1);                      /* -pi - asin(i.sup) */
      return 1;
    }
    else if (IBBasicMaxI(i)<1.0) {           /*  |-----------------------|                    */
      IBBasicMinI(Result1) = -IBBasicMaxI(IBBasicItvConstPi);
      IBBasicRoundNear();
      IBBasicMaxI(Result1) = asin(IBBasicMaxI(i));
      IBBasicRoundUp();
      IBBasicMaxI(Result1) += IBBasicEpsilon;                       /* asin(i.sup) rounded up */

      IBBasicRoundDown();
      IBBasicMinI(Result2) = IBBasicMinI(IBBasicItvConstPi);
      IBBasicMinI(Result2) -= IBBasicMaxI(Result1);                       /* pi - asin(i.sup) */
      IBBasicMaxI(Result2) = IBBasicMaxI(IBBasicItvConstPi);                 /* pi rounded up */
      return 2;
    }
    else {                                   /*  |-------------------------------------|      */
      /* whole domain valid */
      IBBasicMinI(Result1) = -IBBasicMaxI(IBBasicItvConstPi);
      IBBasicMaxI(Result1) =  IBBasicMaxI(IBBasicItvConstPi);
      return 1;   /* Result1 := [-pi, +pi] */
    }
  }
  else if (IBBasicMinI(i)<=0.0) {           /*              |---                             */
    if (IBBasicMaxI(i)<0.0) {              /*              |-----|                          */
      IBBasicRoundNear();
      IBBasicMinI(Result2) = asin(IBBasicMinI(i));
      IBBasicMaxI(Result2) = asin(IBBasicMaxI(i));
      IBBasicRoundDown();
      IBBasicMinI(Result2) -= IBBasicEpsilon;                    /* asin(i.inf) rounded down */
      IBBasicRoundUp();
      IBBasicMaxI(Result2) += IBBasicEpsilon;                      /* asin(i.sup) rounded up */

      IBBasicRoundDown();
      IBBasicMinI(Result1) = -IBBasicMaxI(IBBasicItvConstPi);            /* -pi rounded down */
      IBBasicMinI(Result1) -= IBBasicMaxI(Result2);                     /* -pi - asin(i.sup) */
      IBBasicRoundUp();
      IBBasicMaxI(Result1) = -IBBasicMinI(IBBasicItvConstPi);              /* -pi rounded up */
      IBBasicMaxI(Result1) -= IBBasicMinI(Result2);                     /* -pi - asin(i.inf) */
      return 2;
    }
    else if (IBBasicMaxI(i)<1.0) {          /*              |-----------|                    */
      IBBasicRoundNear();
      IBBasicMinI(Result2) = asin(IBBasicMinI(i));
      IBBasicMaxI(Result2) = asin(IBBasicMaxI(i));
      IBBasicRoundDown();
      IBBasicMinI(Result2) -= IBBasicEpsilon;                    /* asin(i.inf) rounded down */
      IBBasicRoundUp();
      IBBasicMaxI(Result2) += IBBasicEpsilon;                      /* asin(i.sup) rounded up */

      IBBasicMinI(Result1) = -IBBasicMaxI(IBBasicItvConstPi);
      IBBasicMaxI(Result1) = -IBBasicMinI(IBBasicItvConstPi);              /* -pi rounded up */
      IBBasicRoundUp();
      IBBasicMaxI(Result1) -= IBBasicMinI(Result2);                     /* -pi - asin(i.inf) */

      IBBasicRoundDown();
      IBBasicMinI(Result3) = IBBasicMinI(IBBasicItvConstPi);
      IBBasicMinI(Result3) -= IBBasicMaxI(Result2);                      /* pi - asin(i.sup) */
      IBBasicMaxI(Result3) = IBBasicMaxI(IBBasicItvConstPi);                /* pi rounded up */
      return 3;
    }
    else {                                  /*              |-------------------------|      */
      IBBasicRoundNear();
      IBBasicMinI(Result2) = asin(IBBasicMinI(i));
      IBBasicRoundDown();
      IBBasicMinI(Result2) -= IBBasicEpsilon;                    /* asin(i.inf) rounded down */
      IBBasicMaxI(Result2) = IBBasicMaxI(IBBasicItvConstPi);                /* pi rounded up */

      IBBasicMinI(Result1) = -IBBasicMaxI(IBBasicItvConstPi);
      IBBasicMaxI(Result1) = -IBBasicMinI(IBBasicItvConstPi);              /* -pi rounded up */
      IBBasicRoundUp();
      IBBasicMaxI(Result1) -= IBBasicMinI(Result2);                     /* -pi - asin(i.inf) */
      return 2;
    }
  }
  else {                                    /*                       |---                    */
    if (IBBasicMaxI(i)<1.0) {               /*                       |------|                */
      IBBasicRoundNear();
      IBBasicMinI(Result1) = asin(IBBasicMinI(i));
      IBBasicMaxI(Result1) = asin(IBBasicMaxI(i));
      IBBasicRoundDown();
      IBBasicMinI(Result1) -= IBBasicEpsilon;                    /* asin(i.inf) rounded down */
      IBBasicRoundUp();
      IBBasicMaxI(Result1) += IBBasicEpsilon;                      /* asin(i.sup) rounded up */

      IBBasicMinI(Result2) = IBBasicMinI(IBBasicItvConstPi);              /* pi rounded down */
      IBBasicRoundDown();
      IBBasicMinI(Result2) -= IBBasicMaxI(Result1);                      /* pi - asin(i.sup) */
      IBBasicMaxI(Result2) = IBBasicMaxI(IBBasicItvConstPi);                /* pi rounded up */
      IBBasicRoundUp();
      IBBasicMaxI(Result2) -= IBBasicMinI(Result1);                      /* pi - asin(i.inf) */
      return 2;
    }
    else {                                  /*                       |----------------|      */
      IBBasicRoundNear();
      IBBasicMinI(Result1) = asin(IBBasicMinI(i));
      IBBasicRoundDown();
      IBBasicMinI(Result1) -= IBBasicEpsilon;                    /* asin(i.inf) rounded down */

      IBBasicMaxI(Result1) = IBBasicMaxI(IBBasicItvConstPi);                /* pi rounded up */
      IBBasicRoundUp();
      IBBasicMaxI(Result1) -= IBBasicMinI(Result1);                      /* pi - asin(i.inf) */
      return 1;
    }
  }
}


void IBBasicInterII(IBBasicItv Result, IBBasicItv i1, IBBasicItv i2)
/***************************************************************************
*  Result := i1 intersection i2
*/
{
  IBBasicMinI(Result) = IBBasicMax(IBBasicMinI(i1),IBBasicMinI(i2));
  IBBasicMaxI(Result) = IBBasicMin(IBBasicMaxI(i1),IBBasicMaxI(i2));
}


void IBBasicWriteI(FILE *out, IBBasicItv i, int digits, int mode, int verbose)
/***************************************************************************
*  Writes i on out
*  mode = IBBasicPrintIntervalBounds      -> [a,b]
*  mode = IBBasicPrintIntervalMidError    -> mid + [-error,error]
*/
{
  double mid, minerror, maxerror;
  if( IBBasicEmptyI(i) )
  {
    fprintf(out,"empty");
    return;
  }
  if( IBBasicIsDoubleI(i) )
  {
    if( IBBasicMinI(i)>=0 ) fprintf(out,"%.*g",digits,IBBasicMinI(i));
    else fprintf(out,"%+.*g",digits,IBBasicMinI(i));
  }
  else
  {
    if( mode==IBBasicPrintIntervalBounds )
    {
      if( IBBasicMinI(i)>=0 )
      {
        IBBasicRoundDown();
        fprintf(out,"[%.*g , ",digits,IBBasicMinI(i));
        IBBasicRoundUp();
        if( IBBasicMaxI(i)==IBBasicPosInfinity ) fprintf(out,"+oo[");
        else fprintf(out,"%.*g]",digits,IBBasicMaxI(i));
      }
      else
      {
        IBBasicRoundDown();
        if( IBBasicMinI(i)==IBBasicNegInfinity ) fprintf(out,"]-oo , ");
        else fprintf(out,"[%+.*g , ",digits,IBBasicMinI(i));
        IBBasicRoundUp();
        if( IBBasicMaxI(i)==IBBasicPosInfinity ) fprintf(out,"+oo[");
        else fprintf(out,"%+.*g]",digits,IBBasicMaxI(i));
      }
    }
    else
    {
      if( (IBBasicMinI(i)==IBBasicNegInfinity) && (IBBasicMaxI(i)==IBBasicPosInfinity) )
      {
        fprintf(out,"0.0 + ]-oo,+oo[");
        return;
      }
      if( IBBasicMinI(i)==IBBasicNegInfinity )
      {
        IBBasicRoundDown();
        mid = IBBasicCenter(IBBasicMinDouble,IBBasicMaxI(i));
        minerror = IBBasicNegInfinity;
        IBBasicRoundUp();
        maxerror = IBBasicMaxI(i) - mid;
      }
      else if( IBBasicMaxI(i)==IBBasicPosInfinity )
      {
        IBBasicRoundDown();
        mid = IBBasicCenter(IBBasicMinI(i),IBBasicMaxDouble);
        minerror = IBBasicMinI(i) - mid;
        IBBasicRoundUp();
        maxerror = IBBasicPosInfinity;
      }
      else
      {
        IBBasicRoundDown();
        mid = IBBasicMidI(i);
        minerror = IBBasicMinI(i) - mid;
        IBBasicRoundUp();
        maxerror = IBBasicMaxI(i) - mid;
      }

      if( mid>=0 )
      {
        fprintf(out,"%.*g + ",digits,mid);
      }
      else
      {
        fprintf(out,"%+.*g + ",digits,mid);
      }

      if( minerror==IBBasicNegInfinity )
      {
        fprintf(out,"]-oo,%+.4g]",maxerror);
      }
      else if( maxerror==IBBasicPosInfinity )
      {
        IBBasicRoundDown();
        fprintf(out,"[%+.4g,+oo[",minerror);
      }
      else
      {
        IBBasicRoundDown();
        fprintf(out,"[%+.4g,",minerror);
        IBBasicRoundUp();
        fprintf(out,"%+.4g]",maxerror);
      }
    }
  }

#if SOFTWARE_PROFILE
  if (verbose) {
    if (IBBasicIsDoubleI(i)) {
      fprintf(out," **point**");
    }
    else if (IBBasicCanonicalI(i)) {
      fprintf(out," **canonical**");
    }
  }
#endif
}


void IBBasicSetToPi(IBBasicItv i) {
/***************************************************************************
*  i := smallest interval containing Pi
*/
  IBBasicMinI(i) = IBBasicPrevDouble(IBBasicConstPi);
  IBBasicMaxI(i) = IBBasicNextDouble(IBBasicConstPi);
}


void IBBasicSetToLn2(IBBasicItv i) {
/***************************************************************************
*  i := smallest interval containing Ln(2)
*/
  IBBasicMinI(i) = IBBasicPrevDouble(IBBasicConstLn2);
  IBBasicMaxI(i) = IBBasicNextDouble(IBBasicConstLn2);
}


void IBBasicSetToHalfPi(IBBasicItv i) {
/***************************************************************************
*  i := smallest interval containing Pi/2
*/
  IBBasicMinI(i) = IBBasicPrevDouble(IBBasicConstHalfPi);
  IBBasicMaxI(i) = IBBasicNextDouble(IBBasicConstHalfPi);
}

void IBBasicSetToE(IBBasicItv i) {
/***************************************************************************
*  i := smallest interval containing e
*/
  IBBasicMinI(i) = IBBasicPrevDouble(IBBasicConstE);
  IBBasicMaxI(i) = IBBasicNextDouble(IBBasicConstE);
}



int IBBasicNewtonNonzeroII(IBBasicItv Result, IBBasicItv mid, IBBasicItv eval, IBBasicItv deriv)
/***************************************************************************
*  Result := Result inter (mid - eval/deriv)   s.t. 0 not in deriv
*  Returns 1 if Result is not modified, 0 otherwise
*/
{
  IBBasicItv i, j;
  int notmodified = 1;

  IBBasicDivII(i,eval,deriv);  /* i := eval/deriv */
  IBBasicSubII(j,mid,i);      /* j := mid - eval/deriv */

  /* Result <- Result inter j */
  if( IBBasicMinI(i) > IBBasicMinI(Result) )
  {
    IBBasicMinI(Result) = IBBasicMinI(i);
    notmodified = 0;
  }

  if( IBBasicMaxI(i) < IBBasicMaxI(Result) )
  {
    IBBasicMaxI(Result) = IBBasicMaxI(i);
    notmodified = 0;
  }

  return( notmodified );
}


int IBBasicNewtonZeroII(IBBasicItv Result, IBBasicItv mid, IBBasicItv eval, IBBasicItv deriv)
/***************************************************************************
*  Result := Result inter (mid - eval/deriv)   s.t. 0 in deriv
*  Returns 1 if Result is not modified, 0 otherwise
*
*  let eval=[a,b] and  deriv = [c,d]
*/
{
  double min, max;
  int notmodified = 1;
  IBBasicItv j, k, r1, r2;


  if (IBBasicExtendedDivisionII(j,k,eval,deriv)==1) {
    IBBasicSubII(r1,mid,j);   /* r1 := mid - eval/deriv */

    if( IBBasicMinI(r1) > IBBasicMinI(Result) )
    {
      IBBasicMinI(Result) = IBBasicMinI(r1);
      notmodified = 0;
    }

    if( IBBasicMaxI(r1) < IBBasicMaxI(Result) )
    {
      IBBasicMaxI(Result) = IBBasicMaxI(r1);
      notmodified = 0;
    }
  }
  else {
    IBBasicSubII(r1,mid,j);   /* (r1 union r2) := mid - eval/deriv */
    IBBasicSubII(r2,mid,k);

    /* let (r1 union r2) = [-oo,max] union [min,+oo] */
    /* Result is modified if at least one of its bounds is in ]max..min[ */

    if( (IBBasicMinI(Result)>IBBasicMaxI(r1)) && (IBBasicMinI(Result)<IBBasicMinI(r2)) )
    {
      if( IBBasicMaxI(Result)<IBBasicMinI(r2) )  /* Result in ]max..min[ */
      {
        IBBasicSetEmptyI(Result);
      }
      else
      {
        IBBasicMinI(Result) = IBBasicMinI(r2);
      }
      notmodified = 0;
    }
    else if( (IBBasicMaxI(Result)>IBBasicMaxI(r1)) && (IBBasicMaxI(Result)<IBBasicMinI(r2)) )
    {
      if( IBBasicMinI(Result)>IBBasicMaxI(r1) )  /* Result in ]max..min[ */	
      {
        IBBasicSetEmptyI(Result);
      }
      else
      {
        IBBasicMaxI(Result) = IBBasicMaxI(r1);
      }
      notmodified = 0;
    }
  }
  return( notmodified );
}
