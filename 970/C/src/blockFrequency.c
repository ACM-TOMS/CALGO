#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"
#include "../include/tools.h"
#include "../include/stat_fncs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                    B L O C K  F R E Q U E N C Y  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
BlockFrequency(int M, int n)
{
	
	int		i, j, N, blockSum;
	double	p_value, sum, pi, v, chi_squared;
	
	

	N = n/M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	
	for ( i=0; i<N; i++ ) {
		blockSum = 0;
		for ( j=0; j<M; j++ )
			blockSum += epsilon[j+i*M];
		
		pi = (double)blockSum/(double)M;
		v = pi - 0.5;
		sum += v*v;
		//printf("v*v: %lf\n",v*v);
		// added for printage printf("%d ",blockSum);
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N/2.0, chi_squared/2.0);
#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.blockfrequency.chi_squared=chi_squared;
	R_.blockfrequency.p_value=p_value;	
	if(BlockFrequency_v1 == BlockFrequency) R1 = R_;
	else R2 = R_;
#endif
	//printf("%lf",sum);
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
		fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
	}
#endif
#ifdef KS
	pvals.blockfrequency_pvals[pvals.seq_counter] = p_value;
#endif
}

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

// works with global variables: array, R2
// supported rangle of n: 1 - limited by size of int
// supported range of M: 8 - limited by size of int
// not dependent on little/big endian
// bits in array filled from the lowest bits towards the highest bits

void
BlockFrequency2(int M, int n)
{
	/*
	short int LU_byte_weight[256] = 
	{0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5
	,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
	,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
	,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
	,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
	,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
	,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
	,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
	*/
	int		i, N, blockSum,NM,counter;
	double	p_value, sum, pi, v, chi_squared;
	

	N = n/M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	blockSum = 0;
	NM=(N*M / 8);
	if(NM*8<N*M)NM++; // new
	counter = 0;

	for ( i=0; i < NM /*+ 1*/; i++ ) { /// removed +1
		
		blockSum = blockSum + LU_byte_weight[array[i]];
		counter = counter + 8;
		if(counter >= M)
		{
			counter = counter - M ;
			blockSum -= LU_byte_weight[array[i] >> (8 - counter)];
			pi = (double)blockSum/(double)M;
			v = pi - 0.5;
			sum = sum + v*v;
			//printf("v*v: %lf\n",v*v);
			//printf("%d ",blockSum);
			blockSum = LU_byte_weight[array[i] >> (8 - counter)];
		}	
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N/2.0, chi_squared/2.0);
#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.blockfrequency.chi_squared=chi_squared;
	R_.blockfrequency.p_value=p_value;	
	if(BlockFrequency_v1 == BlockFrequency2) R1 = R_;
	else R2 = R_;
#endif
	//printf("%lf ",sum);

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
		fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
	}
#endif
#ifdef KS
	pvals.blockfrequency_pvals[pvals.seq_counter] = p_value;
#endif
}

void
BlockFrequency3(int M, int n)
{
	int		mask, N,blockSum,block, /* n_ , */ processed_bits, padd_sum;
	double  p_value, sum, pi, v, chi_squared;

	int LUT_HW_size = 16;
	int LUT_HW_Bsize = 2;
	signed char *LUT_HW = LUT_HW_16;
	unsigned char *p_tmp, *p_end;

	
	if(0)
	{
		LUT_HW_size = 8;
		LUT_HW_Bsize = 1;
		LUT_HW = LUT_HW_8;
	}

	N = n / M; 		/* # OF SUBSTRING BLOCKS      */
	// n_ = N*M;
	sum = 0.0;
	mask = get_mask(LUT_HW_size);
	blockSum = 0;
	processed_bits = 0;
	
	//bits + LUT_HW_size < n_

	p_end = array + (n- (n%LUT_HW_size))/8;
	for (p_tmp = array; p_tmp < p_end; p_tmp += LUT_HW_Bsize){
		block = *((unsigned int*)p_tmp) & mask;
		blockSum += LUT_HW[block];
		processed_bits += LUT_HW_size;

		if(processed_bits >= M)
		{
			processed_bits -= M;
			padd_sum = LUT_HW[ block & (get_mask(processed_bits) << ( LUT_HW_size - processed_bits))];
			blockSum -= padd_sum;
			pi = (double)blockSum / (double)M;

			v = pi - 0.5;
			sum = sum + v*v;
			blockSum = padd_sum;
		}
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N / 2.0, chi_squared / 2.0);
#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.blockfrequency.chi_squared=chi_squared;
	R_.blockfrequency.p_value=p_value;	
	if(BlockFrequency_v1 == BlockFrequency3) R1 = R_;
	else R2 = R_;
#endif
	//printf("%lf ",sum);

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
		fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
	}
#endif
#ifdef KS
	pvals.blockfrequency_pvals[pvals.seq_counter] = p_value;
#endif
}

//based on Alin Suciu suggestions
void
BlockFrequency4(int M, int n)
{
	int		N, blockSum, block, /* n_, */ i,j, bitend;
	double  p_value, sum, pi, v, chi_squared;

	int LUT_HW_size = 16;
	signed char *LUT_HW = LUT_HW_16;
	const int Tsize = 64;
	uint64_t* pblock;
	uint64_t help;
	const unsigned int mask = (1 << LUT_HW_size) - 1;




	N = n / M; 		/* # OF SUBSTRING BLOCKS      */
	//n_ = N*M;
	sum = 0.0;
	


	for (block = 0; block < N; block++){
		i = block*M;
		bitend = (block+1)*M;
		blockSum = 0;

		while (i % 8 != 0 && i < bitend){
			help = get_nth_block4(array, i);
			blockSum += help & 1;
			i++;
		}

		pblock = (uint64_t*)(array+i/8);
		help = *pblock;
		for (; i < bitend + 1 - Tsize; i += Tsize) {
			for (j = 0; j < Tsize / LUT_HW_size; j++){
				blockSum += LUT_HW[help & mask];
				help >>= LUT_HW_size;
			}
			help = *(++pblock);
		}

		for (; i < bitend; i++) {
			help = get_nth_block4(array, i);
			blockSum += help & 1;
		}
		pi = (double)blockSum / (double)M;
		v = pi - 0.5;
		sum = sum + v*v;
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N / 2.0, chi_squared / 2.0);
#ifdef SPEED
	dummy_result = p_value;
#endif

#ifdef VERIFY_RESULTS
	R_.blockfrequency.chi_squared = chi_squared;
	R_.blockfrequency.p_value = p_value;
	if (BlockFrequency_v1 == BlockFrequency4) R1 = R_;
	else R2 = R_;
#endif
	//printf("%lf ",sum);

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
		fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
		fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
	}
#endif
#ifdef KS
	pvals.blockfrequency_pvals[pvals.seq_counter] = p_value;
#endif
}

