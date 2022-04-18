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

#ifdef _MSC_VER
#if _MSC_VER <= 1700
typedef __int64 int64_t;
typedef int int32_t;
typedef unsigned __int64 uint64_t;
typedef unsigned int uint32_t;
#else
#include <stdint.h>
#endif
#else
#include <stdint.h>
#endif

// New functions
void bits(unsigned char*array, int byte_size);
unsigned int get_nth_block4(unsigned char* array, int offset);
unsigned int get_nth_block_effect(unsigned char* array, int offset);
void test_blocks(unsigned char* array, int size);
int Mirrored_int(unsigned int val, int m);

unsigned int get_block_fast(unsigned char* array, int offset);
unsigned int get_2bytes(unsigned char* array, int byte_offset);
unsigned int get_mask(int size);


// Allin ideas
unsigned int popCountLUT16_64(uint64_t* addr, uint64_t* endAddr);
unsigned int popCountLUT16_32(uint32_t* addr, uint32_t* endAddr);



//BITHACKS
int bitcount(unsigned int n);

unsigned int popCountBITHACK_32(uint32_t* addr, uint32_t* endAddr);
unsigned int runsLUT16_32(uint32_t* addr, uint32_t* endAddr);


//for whole array
void Histogram(int bitstart, int* P, int m, int bits_covered);

void LSHIFT32(int *a, int shift, int Tsize);
void LSHIFT32_p(int* a, int shift, int Tsize);
void LSHIFT64(int64_t *a, int shift, int Tsize);
