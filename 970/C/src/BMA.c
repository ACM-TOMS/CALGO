/* --------------------------------------------------------------------------

The following code is distributed under the following BSD-style license:

Copyright © 2013-2014 Marek Sys (syso@fi.muni.cz) & Zdenek Riha (zriha@fi.muni.cz).
All Rights Reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

3. The name of the author may not be used to endorse or promote products derived
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AUTHORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------- */

#include <string.h>
#include "../include/BMA.h"

#ifdef _MSC_VER
#pragma warning(disable:4146)
#endif
int log2debruins(unsigned int c)
{
  

	 static const int MultiplyDeBruijnBitPosition[32] = 
    {
       32, 2, 29, 3, 30, 15, 25, 4, 31, 23, 21, 16, 26, 18, 5, 9, 
       32, 28, 14, 24, 22, 20, 17, 8, 27, 13, 19, 7, 12, 6, 11, 10
    };

    return MultiplyDeBruijnBitPosition[((unsigned)((c & -c) * 0x077CB531U)) >> 27];
}
#ifdef _MSC_VER
#pragma warning(default:4146)
#endif

void XORT(BMAint *a, BMAint *b, int Tsize)
{
	int i;
	for(i = 0; i < Tsize; i++)
	{
		a[i] ^= b[i];
	}
}

void CPYT(BMAint *a, BMAint *b, int Tsize)
{
	int i;
	for(i = 0; i < Tsize; i++)
	{
		a[i] = b[i];
	}
}


void SETT(BMAint *a, int Tsize)
{
	int i;
	BMAint zero = 0;
	for(i = 0; i < Tsize; i++)
	{
		a[i] ^= zero;
	}
}


void LSHIFTT(BMAint *a, int shift, int Tsize)
{
	int i;
	int com = sizeof(BMAint)* 8 - shift;
	for(i = 0; i < Tsize; i++)
	{
		a[i] = (a[i] >> shift) ^ (a[i+1] << com);
	}
}


int BM_JOURNAL(BMAint  *d_b, BMAint  *d_c, BMAint  *d_t, BMAint *S, int
	n){

	static const int numbitsT = sizeof(BMAint)* 8, Tbytes = sizeof(BMAint);
	int L, Tsize, shift = 2, i;

	Tsize = (n + numbitsT - 1) / numbitsT;

	memset(d_b, 0, Tsize*Tbytes);
	memcpy(d_c, S, Tsize*Tbytes);

	L = 0;

	for (i = 0; i < n; i++) {

		if (d_c[0] & 1)
		{
			memcpy(d_t, d_c, Tsize*Tbytes);
			//CPY(d_t,d_c,Tsize);

			XORT(d_c, d_b, Tsize);
			if (L <= i / 2)
			{
				memcpy(d_b, d_t, Tsize*Tbytes);
				//CPY(d_b,d_t,Tsize);
				L = i + 1 - L;
			}
		}


		shift = log2debruins(((*(unsigned int*)d_c) | 1) ^ 1);
		LSHIFTT(d_c, shift - 1, Tsize);
		i += shift - 2;
		Tsize = (n + numbitsT - i) / numbitsT;

	}
	return L;
}
