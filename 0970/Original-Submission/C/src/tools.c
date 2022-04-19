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

#include <stdio.h>
#include "../include/tools.h"
#include "../include/externs.h"

void bits(unsigned char*arr, int byte_size)
{
	int i,j;
	unsigned char val;
	for(i = 0; i < byte_size; i++)
	{
		val = arr[i];
		for(j = 0; j < 8; j++)
		{
			printf("%d",val&1);
			val >>= 1;
		}
		printf(" ");
	}
	printf("\n");
}

unsigned int get_nth_block4(unsigned char* arr, int offset)
{
	return (*(unsigned int*)(arr+(offset>>3))) >> (offset & 7);//(array2[(offset >> 3)&3][(offset >> 3)] >> (offset & 7));
}
unsigned int get_nth_block_effect(unsigned char* arr, int offset)
{	
	int shift = (offset & 7);
	int byte = (offset >> 3);
	if (shift == 0) return (*(unsigned int*)(arr + byte) >> shift);
	else return (*(unsigned int*)(arr + byte) >> shift)^(*(unsigned int*)(arr + byte+4) << (32 - shift));
}

/*
void test_blocks(unsigned char* array, int size){
	int i,sum = 0,value = 0,counter = 7;
	unsigned char *copy_array = array - 1;
	

	timings();
	for(i = 0; i < size; i++)
	{
		
		sum  += get_nth_block4(array,i);
		
	}
	timings();
	printf("%d",sum);
}
*/

int Mirrored_int(unsigned int val, int m){
	int res = 0,i;
	for(i=0; i < m; i++)
	{
		if(val & (1 << i)) res += (1 << (m - 1 - i));
	}
	return res;
}
/*
int test(unsigned char*array, int n)
{
	int sum, i;

	for(i = 0; i < n; i++)
	{
		sum += get_nth_block4(array,i);
	}
	return sum;
}
*/

unsigned int get_block_fast(unsigned char* arr, int byte_offset)
{
	unsigned int res = (*((unsigned int*)(arr + byte_offset)));

	return res;
}

unsigned int get_2bytes(unsigned char* arr, int byte_offset)
{

	return arr[byte_offset] ^ (arr[byte_offset + 1] << 8);
}

unsigned int get_mask(int size){
	return (1 << size) - 1;
}



unsigned int popCountLUT16_64(uint64_t *addr, uint64_t *endAddr)
{
	unsigned int NumberOne1 = 0;
	uint64_t val1;

	for (; addr<endAddr;)
	{
		val1 = *addr;
		addr++;
		NumberOne1 += LUT_HW_16[val1 >> 48] + LUT_HW_16[(val1 >> 32) & 0xFFFF] + LUT_HW_16[(val1 >> 16) & 0xFFFF] + LUT_HW_16[val1 & 0xFFFF];
	}
	return NumberOne1;
}

unsigned int popCountLUT16_32(uint32_t* addr, uint32_t* endAddr)
{
	unsigned int NumberOne1 = 0;
	uint32_t val1;

	for (; addr<endAddr;)
	{
		val1 = *addr;
		addr++;
		NumberOne1 += LUT_HW_16[val1 >> 16] + LUT_HW_16[val1 & 0xFFFF];
	}
	return NumberOne1;
}

//BITHACKS
int bitcount(unsigned int n) {
	/* works for 32-bit numbers only    */
	/* fix last line for 64-bit numbers */

	register unsigned int tmp;

	tmp = n - ((n >> 1) & 033333333333)
		- ((n >> 2) & 011111111111);
	return ((tmp + (tmp >> 3)) & 030707070707) % 63;
}

unsigned int popCountBITHACK_32(uint32_t* addr, uint32_t* endAddr){
	unsigned int NumberOne1 = 0;

	for (; addr<endAddr; addr++)
	{
		NumberOne1 += bitcount(*addr);
	}
	return NumberOne1;
}

unsigned int runsLUT16_32(uint32_t* addr, uint32_t* endAddr)
{
	unsigned int NumberOne1 = 0;
	unsigned int val1, val2;

	for (; addr<endAddr;)
	{
		val1 = *addr;
		addr++;

		val2 = val1 ^ ((val1 >> 1)) ^ (*addr & 1)*(1 << 31);

		NumberOne1 += LUT_HW_16[val2 >> 16] + LUT_HW_16[val2 & 0xFFFF];
	}
	return NumberOne1;
}




void Histogram(int bitstart, int* P, int m, int bitend){
	int help, mask, i;
	unsigned char* pbyte;

	mask = (1 << m) - 1;

	i = bitstart;
	while (i % 8 != 0 && i <  bitend - m + 1){
		help = get_nth_block4(array, i);
		++P[help & mask];
		i++;
	}

	pbyte = array + i/8;
	help = get_nth_block4(array, i);

	for (; i < bitend - m + 1 - 8; i += 8) {
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;
		++P[help & mask]; help >>= 1;

		++P[help & mask];
		help = *(unsigned int*)(++pbyte);
	}

	for (; i < bitend - m + 1; i++) {
		help = get_nth_block4(array, i);
		++P[help & mask];
	}
}


void LSHIFT32(int *a, int shift, int Tsize)
{
	int i;
	int com = sizeof(int)* 8 - shift;

	for (i = 0; i < Tsize; i++)
	{
		a[i] = (a[i] >> shift) ^ (a[i + 1] << com);
	}
}

/*
void LSHIFT32_p(int* a, int shift, int Tsize)
{
	unsigned int* p[4];
	int i;
	p[0] = ((unsigned char*)a + 0);
	p[1] = ((unsigned char*)a + 1);
	p[2] = ((unsigned char*)a + 2);
	p[3] = ((unsigned char*)a + 3);

	for (i = 0; i < Tsize; i++)
	{

	}
}
*/

void LSHIFT64(int64_t *a, int shift, int Tsize)
{
	int i;
	int com = sizeof(int64_t)* 8 - shift;

	for (i = 0; i < Tsize; i++)
	{

		a[i] = (a[i] >> shift) ^ (a[i + 1] << com);
	}
}


