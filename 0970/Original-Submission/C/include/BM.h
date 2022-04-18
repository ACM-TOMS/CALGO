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


#ifndef _BM_H_
#define _BM_H_
#include "../include/tools.h"


/*typedef unsigned int type;


typedef struct mybitset{
	type *array;
	int first, last;
	int size,array_size;
} mybitset;



void set_bit(mybitset* bitset,int index);
void clear(mybitset* bitset);
void resize(mybitset* bitset, int new_size);
void left_shift(mybitset* bitset, int bits, int from, int to);

void XOR(mybitset* a, mybitset* b, int from, int to);
void copy(mybitset* a, mybitset* b, int from, int to);

void rand_array(mybitset* bitset);
int BM(mybitset bitstream, int N);
*/
typedef unsigned int type;
void print_bits(unsigned char c);
void array_as_bits(unsigned char* c,int byte_size);

void left_shift(type* array, int array_size, int from, int to);
void copy(type* array, type* b, int array_size, int from, int to);
void XOR(type* array, type* b, int array_size, int from, int to);
int BM_c(type* bitstream,int N,type* c, type *b,type* t);

#endif