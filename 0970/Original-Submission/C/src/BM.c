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
#include <string.h>
#include <stdlib.h>
#include "../include/BM.h"

void print_bits(unsigned char c)
{
	int i;
	for(i = 0; i < 8; i++)
	{
		if ( (c >> i) & 1) printf("1");
		else printf("0");
	}
}
void array_as_bits(unsigned char* c,int byte_size)
{
	int i;
	for(i = 0; i < byte_size; i++)
	{
		print_bits(c[i]);
		printf(" ");
	}
	printf("\n");
}


const int type_size_bytes = 4,type_size_bits = 32;

void left_shift(type* array, int array_size, int from, int to){
	
	int i;
	//if(to == -1 || to > array_size) to = array_size;
	
	for(i = to-1; i > from; i--){
		array[i] = (array[i] << 1) ^ (array[i-1] >> 31);
	}
	array[from] = (array[from]  << 1);
	
	
}
void copy(type* array, type* b, int array_size, int from, int to){
	//memcpy(array,b,to*type_size_bytes);
	int i;
	
	//if(to == -1 || to > array_size) to = array_size;
	for( i = from; i < to; i++)
	{
		array[i] = b[i];
	}
}
void XOR(type* array, type* b, int array_size, int from, int to){
	int i;
	//if(to == -1 || to > array_size) to = array_size;
	for(i = from; i < to; i++)
	{
		array[i] ^= b[i];
	}
}

/*type& bit(type* array, int index){
	return array[index / type_size / 8]
}*/


int LU_byte_weight2[256] = 
	{0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,
	1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1
	,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,0,
	1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,1,0,
	0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,
	0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,
	1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1,0,
	1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0};


int BM_c(type* bitstream,int N,type* c, type *b,type* t){
	
	int low_bound = 0,up_bound, l = 0 /*, m = -1 */,n,i;
	type d;
	int array_size = N / type_size_bits + (N % type_size_bits != 0);
	unsigned char *pw,w;
	//unsigned int help;

	memset(c,0,array_size*type_size_bytes);
	memset(b,0,array_size*type_size_bytes); 
	//memset(t,0,array_size*type_size_bytes);
	
	b[0] = c[0] = 1;
	pw = (unsigned char*)&d;
	
	for ( n = 0; n < N; n++) {
		d = 0;
		
		if( (n % type_size_bits) != 0) up_bound = n / type_size_bits + 2;
		else up_bound = (n / type_size_bits) + 1;
		
		//up_bound = (n / type_size_bits) + 2;
		if(up_bound > array_size)up_bound = array_size;
		
		//up_bound = (n/type_size_bits + 2);
		//help = up_bound - array_size;
		//up_bound = (help >> 31 )*help+array_size;


		for (i = 0; i < up_bound; i++) {
			d ^= c[i] & bitstream[i];
		}

		
		w = (pw[0]^pw[1])^(pw[2]^pw[3]);
		
		left_shift(c,array_size,low_bound,up_bound+1);
		

		//if ( (LU_byte_weight[w] & 1) == 1) {
		if ( LU_byte_weight2[w] == 1) {
			//copy(t,c,array_size,low_bound,up_bound);
			memcpy(t,c,up_bound*type_size_bytes);
			//memcpy(&t[low_bound],&c[low_bound],(up_bound-low_bound)*type_size_bytes);
			if( (b[low_bound] | c[low_bound]) == 0 && low_bound < n/type_size_bits )low_bound++;
			XOR(c,b,array_size,low_bound,up_bound+1); // c was shifted to next block
			
			if (l <= n / 2) {
				l = n + 1 - l;
				// m = n; // nepouzito
				//copy(b,t,array_size,low_bound,up_bound);
				//memcpy(b,t,up_bound*type_size_bytes);
				memcpy(&b[low_bound],&t[low_bound],(up_bound-low_bound)*type_size_bytes);
			}
		}
		
	}
	return l;
}
