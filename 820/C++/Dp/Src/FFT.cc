#include "common.h"
#include "Interval.h"
#include "FFT.h"
#include <math.h>
#include <assert.h>

//****************************************************************************
// FFT.cc                                
//****************************************************************************




///////////////////////////////////////////////////////////////////////////////

integer rev(integer In, integer J)
{
  integer Out = In&1, j;
  
  for(j=1; j<J; j++){
    In  >>= 1;
    Out <<= 1;
    Out += In&1;
  }
  return Out;
}

///////////////////////////////////////////////////////////////////////////////

void bitReverse(Interval &In)
{
  real temp;
  real *dataptr = In.origin;
  integer N = In.length;
  integer J = Log2(N), j, k;

  for(j=0; j<N; j++){
    k=rev(j,J);
    if(k > j){
      temp = *dataptr;
      *dataptr = In.origin[k];
      In.origin[k] = temp;
    }
    dataptr++;
  }
}

//////////////////////////////////////////////////////////////////////////

void Weights(Interval& Rew, Interval& Imw, integer flag)
{
  integer j, k, N = Rew.length<<1, N4 = N>>2; // N is the number of data points in the FFT algrithm
  real    ReW, ImW, K = flag * (2 * M_PI)/N, ReTemp, ImTemp;
  real    *Rewptr = Rew.origin, *Imwptr = Imw.origin, *Rew1ptr = Rew.origin + 1, *Imw1ptr = Imw.origin + 1,
          *ReEndwptr = Rew.origin+( ( N>>1 ) - 1 ), *ImEndwptr = Imw.origin+((N>>1)-1), *RewTemp, *ImwTemp;

  *(Rewptr) = 1;
  *(Imwptr) = 0;
  *(ReEndwptr) = -(*(++Rewptr) = cos(K));
  *(ImEndwptr) = *(++Imwptr) = sin(K);
  for(j=2; j<N4; j<<=1){
    Rewptr++;
    Imwptr++;
    ReEndwptr--;
    ImEndwptr--;
    RewTemp = Rew1ptr;
    ImwTemp = Imw1ptr;
    *(ReEndwptr) = -( ReW = *(Rewptr) = cos(K*j) );
    *(ImEndwptr) = ImW = *(Imwptr) = sin(K*j);
    for(k=1; k<j; k++){
      ReTemp = *(RewTemp++);
      ImTemp = *(ImwTemp++);
      *(--ReEndwptr) = -( *(++Rewptr) = ReW*ReTemp-ImW*ImTemp );
      *(--ImEndwptr) = *(++Imwptr) = ReW*ImTemp + ImW*ReTemp;
    }
  }
  if( N > 4 ){
    *(++Rewptr)=0;
    *(++Imwptr)=flag;
  }
  return;
}

//////////////////////////////////////////////////////////////////////////

void CoreFFT(Interval& ReA, Interval& ImA, Interval& C, Interval& S)
{
  
  real ReTemp0, ImTemp0, ReTemp1, ImTemp1;
  real *RePtr0, *ImPtr0, *RePtr1, *ImPtr1, *Cptr, *Sptr;
  integer s, n, l, b, j;
  integer J = Log2( ReA.length );

  for(l = 1; l <= J; l++){

    s = 1 << ( l - 1 );
    n = 1 << ( J - l );
    RePtr0 = ReA.origin;
    ImPtr0 = ImA.origin;
    RePtr1 = (ReA.origin + s);
    ImPtr1 = (ImA.origin + s);
    for(b = 0; b < n; b++){
      Cptr = C.origin;
      Sptr = S.origin;
      for(j = 0; j < s; j++){
	ReTemp0 = *(RePtr0);
	ImTemp0 = *(ImPtr0);
	ReTemp1 = *(Cptr) * ( *(RePtr1)) - *(Sptr) * (*(ImPtr1));
	ImTemp1 = *(Sptr) * ( *(RePtr1)) + *(Cptr) * (*(ImPtr1));
	*(RePtr0++) = ReTemp0 + ReTemp1;
	*(ImPtr0++) = ImTemp0 + ImTemp1;
	*(RePtr1++) = ReTemp0 - ReTemp1;
	*(ImPtr1++) = ImTemp0 - ImTemp1;
	Cptr += n;
	Sptr += n;

      } // end j-loop
      RePtr0 += s;
      ImPtr0 += s;
      RePtr1 += s;
      ImPtr1 += s;
      
    } // end b-loop
  } // end l-loop
  return;
}

////////////////////////////////////////////////////////////////////////////

void FFT(Interval &ReIn, Interval &ImIn){

  assert(ReIn.length == ImIn.length);
  integer flag = -1, N = ReIn.length, J = Log2(N);
  assert(N == (1<<J));
  Interval COS(J-1), SIN(J-1);
  bitReverse(ReIn);
  bitReverse(ImIn);
  Weights(COS,SIN,flag);
  CoreFFT(ReIn, ImIn, COS, SIN);
}
