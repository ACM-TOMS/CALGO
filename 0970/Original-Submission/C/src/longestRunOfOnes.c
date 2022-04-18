/* got rid of unused 'k' */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"  
#include "../include/tools.h"
#include "../include/stat_fncs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                      L O N G E S T  R U N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
LongestRunOfOnes(int n)
{
	double			p_value, chi2, pi[7];
	int				run, v_n_obs, N, i, j, K, M, V[7];
	unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };

	if (n < 128) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
			fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_LONGEST_RUN], "\t\t   n=%d is too short\n", n);
	}
#endif
		return;
	}
	if ( n < 6272 ) {
		K = 3;
		M = 8;
		V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
		pi[0] = 0.21484375;
		pi[1] = 0.3671875;
		pi[2] = 0.23046875;
		pi[3] = 0.1875;
	}
	else if ( n < 750000 ) {
		K = 5;
		M = 128;
		V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
		pi[0] = 0.1174035788;
		pi[1] = 0.242955959;
		pi[2] = 0.249363483;
		pi[3] = 0.17517706;
		pi[4] = 0.102701071;
		pi[5] = 0.112398847;
	}
	else {
		K = 6;
		M = 10000;
			V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
		pi[0] = 0.0882;
		pi[1] = 0.2092;
		pi[2] = 0.2483;
		pi[3] = 0.1933;
		pi[4] = 0.1208;
		pi[5] = 0.0675;
		pi[6] = 0.0727;
	}
	
	N = n/M;
	for ( i=0; i<N; i++ ) {
		v_n_obs = 0;
		run = 0;
		for ( j=0; j<M; j++ ) {
			if ( epsilon[i*M+j] == 1 ) {
				run++;
				if ( run > v_n_obs )
					v_n_obs = run;
			}
			else
				run = 0;
		}
		if ( v_n_obs < V[0] )
			nu[0]++;
		for ( j=0; j<=K; j++ ) {
			if ( v_n_obs == V[j] )
				nu[j]++;
		}
		if ( v_n_obs > V[K] )
			nu[K]++;
		//printf("%d ", v_n_obs);
	}

	chi2 = 0.0;
	for ( i=0; i<=K; i++ )
		chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

	p_value = cephes_igamc((double)(K / 2.0), chi2 / 2.0);

	//for(i = 0; i < K; i++)
	//	printf("%d ",nu[i]);

#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.longestrunofones.N=N;
	R_.longestrunofones.M=M;
	R_.longestrunofones.chi2=chi2;
	for(i = 0; i < 7; i++) R_.longestrunofones.nu[i]=nu[i];
	R_.longestrunofones.p_value=p_value;
	if(LongestRunOfOnes_v1 == LongestRunOfOnes) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(a) N (# of substrings)  = %d\n", N);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(b) M (Substring Length) = %d\n", M);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(c) Chi^2                = %f\n", chi2);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t      F R E Q U E N C Y\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");

		if (K == 3) {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t  <=1     2     3    >=4   P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
		}
		else if (K == 5) {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t<=4  5  6  7  8  >=9 P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5]);
		}
		else {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t<=10  11  12  13  14  15 >=16 P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5], nu[6]);
		}
		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_LONGEST_RUN], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

		fprintf(stats[TEST_LONGEST_RUN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_LONGEST_RUN]);
		fprintf(results[TEST_LONGEST_RUN], "%f\n", p_value); fflush(results[TEST_LONGEST_RUN]);
	}
#endif
#ifdef KS
	pvals.longestrunofones_pvals[pvals.seq_counter] = p_value;
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

void
LongestRunOfOnes2(int n)
{
	unsigned char  start_run[256]={0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,8};
	unsigned char  end_run[256]  ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,6,6,7,8};
	unsigned char  max_run[256]  ={0,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,1,1,1,2,1,1,2,3,2,2,2,2,3,3,4,5,1,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,2,2,2,2,2,2,2,3,3,3,3,3,4,4,5,6,1,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,1,1,1,2,1,1,2,3,2,2,2,2,3,3,4,5,2,2,2,2,2,2,2,3,2,2,2,2,2,2,3,4,3,3,3,3,3,3,3,3,4,4,4,4,5,5,6,7,1,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,1,1,1,2,1,1,2,3,2,2,2,2,3,3,4,5,1,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,2,2,2,2,2,2,2,3,3,3,3,3,4,4,5,6,2,2,2,2,2,2,2,3,2,2,2,2,2,2,3,4,2,2,2,2,2,2,2,3,2,2,2,2,3,3,4,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,5,5,5,5,6,6,7,8};


	double			p_value, chi2, pi[7];
	int				run, v_n_obs, N, i, j, K, M, /*Mbytes, */ V[7],rest_bits;
	unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };
//	unsigned int end_mask[8] = {255,254,252,248,240,224,192,128};
	if (n < 128) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
			fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_LONGEST_RUN], "\t\t   n=%d is too short\n", n);
	}
#endif
		return;
	}
	if ( n < 6272 ) {
		K = 3;
		M = 8;
		V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
		pi[0] = 0.21484375;
		pi[1] = 0.3671875;
		pi[2] = 0.23046875;
		pi[3] = 0.1875;
	}
	else if ( n < 750000 ) {
		K = 5;
		M = 128;
		V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
		pi[0] = 0.1174035788;
		pi[1] = 0.242955959;
		pi[2] = 0.249363483;
		pi[3] = 0.17517706;
		pi[4] = 0.102701071;
		pi[5] = 0.112398847;
	}
	else {
		K = 6;
		M = 10000;
			V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
		pi[0] = 0.0882;
		pi[1] = 0.2092;
		pi[2] = 0.2483;
		pi[3] = 0.1933;
		pi[4] = 0.1208;
		pi[5] = 0.0675;
		pi[6] = 0.0727;
	}
	
	N = n/M;
	//Mbytes = M/8;
	//bits(array, 8);
	for ( i=0; i<N; i++ ) {

		//start byte of block
		rest_bits = i*M % 8;
		if(rest_bits != 0 )
		{
			j = i*M/8 + 1;
			v_n_obs = max_run[array[j] >> (8 - rest_bits)];
			run = end_run[((array[j] >> (8 - rest_bits)) << (8 - rest_bits))];
			//run = end_run[ array[j] & end_mask[rest_bits] ];
		}
		else
		{
			j = i*M/8;
			v_n_obs = max_run[array[j]];
			run = end_run[array[j]];
		}
		//j = i*M/8 + (rest_bits != 0);
		//v_n_obs = max_run[array[j] >> (8 - rest_bits)];
		//run = end_run[ array[j] & end_mask[rest_bits] ];
		
		

		j++;
		//whole bytes
		for (; j < M*(i+1)/8; j++ ) {
			run += start_run[array[j]];
			
			if ( run > v_n_obs ){
				v_n_obs = run;
			}
			if (  max_run[array[j]] > v_n_obs ){
				v_n_obs =  max_run[array[j]];
			}
			
			//v_n_obs = v_n_obs+( (v_n_obs-run) >> 31)*(v_n_obs-run);
		
			
			if(array[j] != 255) run = end_run[array[j]];
			//if(v_n_obs > V[K])break;		
		

		}
		//enb byte of block
		rest_bits = (i+1)*M % 8;
		if(rest_bits != 0 )
		{
			run += start_run[array[j] & (255 >> (8 - rest_bits))];;
			if ( run > v_n_obs ){
				v_n_obs = run;
			}
			
		}
		
		
		//printf("%d ", v_n_obs);
		if(v_n_obs < V[0])nu[0]++;
		else if(v_n_obs > V[K])nu[K]++;
		else nu[v_n_obs-V[0]]++;
		


	}


	chi2 = 0.0;
	for ( i=0; i<=K; i++ )
		chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

	p_value = cephes_igamc((double)(K / 2.0), chi2 / 2.0);
	//for(i = 0; i < K; i++)
	//	printf("%d ",nu[i]);

#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.longestrunofones.N=N;
	R_.longestrunofones.M=M;
	R_.longestrunofones.chi2=chi2;
	for(i = 0; i < 7; i++) R_.longestrunofones.nu[i]=nu[i];
	R_.longestrunofones.p_value=p_value;
	if(LongestRunOfOnes_v1 == LongestRunOfOnes2) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(a) N (# of substrings)  = %d\n", N);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(b) M (Substring Length) = %d\n", M);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(c) Chi^2                = %f\n", chi2);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t      F R E Q U E N C Y\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");

		if (K == 3) {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t  <=1     2     3    >=4   P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
		}
		else if (K == 5) {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t<=4  5  6  7  8  >=9 P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5]);
		}
		else {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t<=10  11  12  13  14  15 >=16 P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5], nu[6]);
		}
		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_LONGEST_RUN], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

		fprintf(stats[TEST_LONGEST_RUN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_LONGEST_RUN]);
		fprintf(results[TEST_LONGEST_RUN], "%f\n", p_value); fflush(results[TEST_LONGEST_RUN]);
	}
#endif
#ifdef KS
	pvals.longestrunofones_pvals[pvals.seq_counter] = p_value;
#endif
}

void
LongestRunOfOnes3(int n)
{
	unsigned int    mask, tmp, processed_bits, block;

	double			p_value, chi2, pi[7];
	int				run, v_n_obs, N, i, K, M, V[7];
	unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };

	/*
	int LUT_Lrun_size = 8;
	int LUT_Lrun_Bsize = 1;

	signed char *LUT_Lrun_start = LUT_Lrun_start_8;
	signed char *LUT_Lrun_end = LUT_Lrun_end_8;
	signed char *LUT_Lrun_max = LUT_Lrun_max_8;
	*/
	
	int LUT_Lrun_size = 16;
	int LUT_Lrun_Bsize = 2;

	signed char *LUT_Lrun_start = LUT_Lrun_start_16;
	signed char *LUT_Lrun_end = LUT_Lrun_end_16;
	signed char *LUT_Lrun_max = LUT_Lrun_max_16;
	
	unsigned char *p_tmp, *p_end;

	if (n < 128) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
			fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_LONGEST_RUN], "\t\t   n=%d is too short\n", n);
	}
#endif
		return;
	}
	if (n < 6272) {
		K = 3;
		M = 8;
		V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
		pi[0] = 0.21484375;
		pi[1] = 0.3671875;
		pi[2] = 0.23046875;
		pi[3] = 0.1875;
	}
	else if (n < 750000) {
		K = 5;
		M = 128;
		V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
		pi[0] = 0.1174035788;
		pi[1] = 0.242955959;
		pi[2] = 0.249363483;
		pi[3] = 0.17517706;
		pi[4] = 0.102701071;
		pi[5] = 0.112398847;
	}
	else {
		K = 6;
		M = 10000;
		V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
		pi[0] = 0.0882;
		pi[1] = 0.2092;
		pi[2] = 0.2483;
		pi[3] = 0.1933;
		pi[4] = 0.1208;
		pi[5] = 0.0675;
		pi[6] = 0.0727;
	}




	if (M == 8)
	{
		LUT_Lrun_size = 8;
		LUT_Lrun_Bsize = 1;

		LUT_Lrun_start = LUT_Lrun_start_8;
		LUT_Lrun_end = LUT_Lrun_end_8;
		LUT_Lrun_max = LUT_Lrun_max_8;
	}

	
	N = n / M;
	run = 0;
	block = 1;
	processed_bits = 0;
	mask = get_mask(LUT_Lrun_size);
	v_n_obs = 0;
	
	p_end = array + M*N/8;
	for (p_tmp = array; p_tmp < p_end; p_tmp += LUT_Lrun_Bsize)
	{
		tmp =  *((unsigned int*)p_tmp) & mask;
		run += LUT_Lrun_start[tmp];

		if (run > v_n_obs){
			v_n_obs = run;
		}
		if (LUT_Lrun_max[tmp] > v_n_obs){
			v_n_obs = LUT_Lrun_max[tmp];
		}

		if (tmp != mask) run = LUT_Lrun_end[tmp];
		processed_bits += LUT_Lrun_size;

		if (processed_bits == M*block){
			if (v_n_obs < V[0])nu[0]++;
			else if (v_n_obs > V[K])nu[K]++;
			else nu[v_n_obs - V[0]]++;
			block++;
			v_n_obs = run = 0;
		}
	}

	
	chi2 = 0.0;
	for (i = 0; i <= K; i++)
		chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

	p_value = cephes_igamc((double)(K / 2.0), chi2 / 2.0);

#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.longestrunofones.N=N;
	R_.longestrunofones.M=M;
	R_.longestrunofones.chi2=chi2;
	for(i = 0; i < 7; i++) R_.longestrunofones.nu[i]=nu[i];
	R_.longestrunofones.p_value=p_value;
	if(LongestRunOfOnes_v1 == LongestRunOfOnes3) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(a) N (# of substrings)  = %d\n", N);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(b) M (Substring Length) = %d\n", M);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t(c) Chi^2                = %f\n", chi2);
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t      F R E Q U E N C Y\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");

		if (K == 3) {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t  <=1     2     3    >=4   P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
		}
		else if (K == 5) {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t<=4  5  6  7  8  >=9 P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5]);
		}
		else {
			fprintf(stats[TEST_LONGEST_RUN], "\t\t<=10  11  12  13  14  15 >=16 P-value  Assignment");
			fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5], nu[6]);
		}
		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_LONGEST_RUN], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

		fprintf(stats[TEST_LONGEST_RUN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_LONGEST_RUN]);
		fprintf(results[TEST_LONGEST_RUN], "%f\n", p_value); fflush(results[TEST_LONGEST_RUN]);
	}
#endif
#ifdef KS
	pvals.longestrunofones_pvals[pvals.seq_counter] = p_value;
#endif
}
