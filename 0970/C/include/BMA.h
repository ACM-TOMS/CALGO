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

#ifndef BMA__H
#define BMA__H

typedef unsigned int BMAint;
//typedef unsigned __int64 BMAint;

extern BMAint *S;
extern BMAint *b, *c, *t;
extern BMAint *d_b, *d_c, *d_t;

int log2debruins(unsigned int c);
void XORT(BMAint *dst, BMAint *a, int num_bits);
void ANDT(BMAint  *dst, BMAint *a, int num_bits);
void LSHIFTT(BMAint *source, int size, int bit_shift);
int BM_JOURNAL(BMAint *b, BMAint *c, BMAint  *t, BMAint  *S, int n);

#endif  
