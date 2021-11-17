#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"
#include "../include/tools.h"  
#include "../include/stat_fncs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		    C U M U L A T I V E  S U M S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
CumulativeSums(int n)
{
	int		S, sup, inf, z, zrev, k;
	double	sum1, sum2, p_value;
	

	S = 0;
	sup = 0;
	inf = 0;
	for ( k=0; k<n; k++ ) {
		epsilon[k] ? S++ : S--;
		if ( S > sup )
			sup++;
		if ( S < inf )
			inf--;
		z = (sup > -inf) ? sup : -inf;
		zrev = (sup-S > S-inf) ? sup-S : S-inf;
	}
	
	//printf("%d %d", z, zrev);
#ifdef VERIFY_RESULTS
	R_.cusum.z=z;
	R_.cusum.zrev=zrev;
#endif


	// forward
	sum1 = 0.0;
	for ( k=(-n/z+1)/4; k<=(n/z-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*z)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*z)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/z-3)/4; k<=(n/z-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*z)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*z)/sqrt(n));
	}

	p_value = 1.0 - sum1 + sum2;
#ifdef SPEED
	dummy_result = p_value;
	//printf("%lf %lf %lf", sum1, sum2, p_value);
#endif
#ifdef VERIFY_RESULTS
	R_.cusum.p_valueA=p_value;
	R_.cusum.sum1A=sum1;
	R_.cusum.sum2A=sum2;
#endif
	
#ifdef KS
	pvals.cusum_pvals[0][pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
 if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (FORWARD) TEST\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", z);
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

		fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
		fprintf(results[TEST_CUSUM], "%f\n", p_value);
	}
#endif	
	// backwards
	sum1 = 0.0;
	for ( k=(-n/zrev+1)/4; k<=(n/zrev-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*zrev)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/zrev-3)/4; k<=(n/zrev-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*zrev)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt(n));
	}
	p_value = 1.0 - sum1 + sum2;
#ifdef SPEED
	dummy_result += p_value;
	//printf("%lf %lf %lf", sum1, sum2, p_value);
#endif
#ifdef VERIFY_RESULTS
	R_.cusum.p_valueB=p_value;
	R_.cusum.sum1B=sum1;
	R_.cusum.sum2B=sum2;
	if(CumulativeSums_v1 == CumulativeSums) R1 = R_;
	else R2 = R_;
#endif

#ifdef KS
	pvals.cusum_pvals[1][pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
 if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (REVERSE) TEST\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", zrev);
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

		fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_CUSUM]);
		fprintf(results[TEST_CUSUM], "%f\n", p_value); fflush(results[TEST_CUSUM]);
	}
#endif
	//printf("sum1 %lf   sum2 %lf",sum1,sum2);
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
CumulativeSums2(int n)
{
	int		S, sup, inf, z, zrev, k, byte_ind;
	double	sum1, sum2, p_value;
	int max_minus[256] = {-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-7,-5,-5,-3,-5,-3,-3,-1,-5,-3,-3,-1,-3,-1,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0};
	int max_plus[256] = {0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,1,3,1,3,3,5,1,3,3,5,3,5,5,7,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8};
	int byte_sum[256] = {-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8};
unsigned int window;

	S = 0;
	sup = 0;
	inf = 0;
	for ( byte_ind = 0; byte_ind < n/8; byte_ind += 1 ) {

		if(S + max_plus[array[byte_ind]] > sup) sup = S + max_plus[array[byte_ind]];
		if(S + max_minus[array[byte_ind]] < inf) inf = S + max_minus[array[byte_ind]];
		S += byte_sum[array[byte_ind]];
	}
	//LAST BYTE
	if( n % 8 != 0)
	{
		for ( k=byte_ind*8; k<n; k++ ) {
			
			window = get_nth_block4(array,k);
			(window & 1) ? S++ : S--;

			if ( S > sup )
				sup++;
			if ( S < inf )
				inf--;
		}
	}


	z = (sup > -inf) ? sup : -inf;

	zrev = (sup-S > S-inf) ? sup-S : S-inf;

#ifdef VERIFY_RESULTS
	R_.cusum.z=z;
	R_.cusum.zrev=zrev;
#endif

	//printf("z = %d zrev=%d S= %d", z, zrev,S);
	
	// forward
	sum1 = 0.0;
	for ( k=(-n/z+1)/4; k<=(n/z-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*z)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*z)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/z-3)/4; k<=(n/z-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*z)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*z)/sqrt(n));
	}

	p_value = 1.0 - sum1 + sum2;
#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS	
	R_.cusum.p_valueA=p_value;
	R_.cusum.sum1A=sum1;
	R_.cusum.sum2A=sum2;
#endif


#ifdef KS
	pvals.cusum_pvals[0][pvals.seq_counter] = p_value;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
 if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (FORWARD) TEST\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", z);
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

		fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
		fprintf(results[TEST_CUSUM], "%f\n", p_value);
	}
#endif
	// backwards
	sum1 = 0.0;
	for ( k=(-n/zrev+1)/4; k<=(n/zrev-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*zrev)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/zrev-3)/4; k<=(n/zrev-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*zrev)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt(n));
	}
	p_value = 1.0 - sum1 + sum2;
#ifdef SPEED
	dummy_result += p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.cusum.p_valueB=p_value;
	R_.cusum.sum1B=sum1;
	R_.cusum.sum2B=sum2;
	if(CumulativeSums_v1 == CumulativeSums2) R1 = R_;
	else R2 = R_;
#endif

#ifdef KS
	pvals.cusum_pvals[1][pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
 if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (REVERSE) TEST\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", zrev);
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

		fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_CUSUM]);
		fprintf(results[TEST_CUSUM], "%f\n", p_value); fflush(results[TEST_CUSUM]);
	}
#endif
	//printf("sum1 %lf   sum2 %lf\n\n",sum1,sum2);
}

void
CumulativeSums3(int n)
{
	int		S, sup, inf, z, zrev, k, plus, minus, tmp, i, mask;
	double	sum1, sum2, p_value;

	int LUT_Cusum_size = 16;
	int LUT_Cusum_Bsize = 2;
	signed char *LUT_Cusum = LUT_Cusum_16;
	signed char *LUT_Cusum_max_positiv = LUT_Cusum_max_positiv_16;
	signed char *LUT_Cusum_max_negativ = LUT_Cusum_max_negativ_16;
	unsigned char *p_tmp, *p_end;

	if(0)
	{
		LUT_Cusum_size = 8;
		LUT_Cusum_Bsize = 1;
		LUT_Cusum = LUT_Cusum_8;
		LUT_Cusum_max_positiv = LUT_Cusum_max_positiv_8;
		LUT_Cusum_max_negativ = LUT_Cusum_max_negativ_8;
	}
	S = 0;
	sup = 0;
	inf = 0;

	mask = get_mask(LUT_Cusum_size);

	
	p_end = array + (n- (n%LUT_Cusum_size))/8;
	for (p_tmp = array; p_tmp < p_end; p_tmp += LUT_Cusum_Bsize){
		tmp = *((unsigned int*)p_tmp) & mask & mask;

		plus = LUT_Cusum_max_positiv[tmp];
		minus = LUT_Cusum_max_negativ[tmp];
		if (S + plus > sup) sup = S + plus;
		if (S + minus < inf) inf = S + minus;

		S += LUT_Cusum[tmp];
	}




	//LAST BLOCK
	if ( n % LUT_Cusum_size != 0)
	{
		for (i = n - (n % LUT_Cusum_size)  ; i < n; i++) {

			if (get_nth_block4(array, i) & 1) S++;
			else  S--;

			if (S > sup)
				sup++;
			if (S < inf)
				inf--;
		}
	}


	z = (sup > -inf) ? sup : -inf;

	zrev = (sup - S > S - inf) ? sup - S : S - inf;

#ifdef VERIFY_RESULTS
	R_.cusum.z = z;
	R_.cusum.zrev = zrev;
#endif

	//printf("z = %d zrev=%d S= %d", z, zrev,S);

	// forward
	sum1 = 0.0;
	for (k = (-n / z + 1) / 4; k <= (n / z - 1) / 4; k++) {
		sum1 += cephes_normal(((4 * k + 1)*z) / sqrt(n));
		sum1 -= cephes_normal(((4 * k - 1)*z) / sqrt(n));
	}
	sum2 = 0.0;
	for (k = (-n / z - 3) / 4; k <= (n / z - 1) / 4; k++) {
		sum2 += cephes_normal(((4 * k + 3)*z) / sqrt(n));
		sum2 -= cephes_normal(((4 * k + 1)*z) / sqrt(n));
	}

	p_value = 1.0 - sum1 + sum2;
#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS	
	R_.cusum.p_valueA = p_value;
	R_.cusum.sum1A = sum1;
	R_.cusum.sum2A = sum2;
#endif

#ifdef KS
	pvals.cusum_pvals[0][pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
 if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (FORWARD) TEST\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", z);
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

		fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
		fprintf(results[TEST_CUSUM], "%f\n", p_value);
	}
#endif
	// backwards
	sum1 = 0.0;
	for (k = (-n / zrev + 1) / 4; k <= (n / zrev - 1) / 4; k++) {
		sum1 += cephes_normal(((4 * k + 1)*zrev) / sqrt(n));
		sum1 -= cephes_normal(((4 * k - 1)*zrev) / sqrt(n));
	}
	sum2 = 0.0;
	for (k = (-n / zrev - 3) / 4; k <= (n / zrev - 1) / 4; k++) {
		sum2 += cephes_normal(((4 * k + 3)*zrev) / sqrt(n));
		sum2 -= cephes_normal(((4 * k + 1)*zrev) / sqrt(n));
	}
	p_value = 1.0 - sum1 + sum2;
#ifdef SPEED
	dummy_result += p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.cusum.p_valueB=p_value;
	R_.cusum.sum1B=sum1;
	R_.cusum.sum2B=sum2;
	if(CumulativeSums_v1 == CumulativeSums3) R1 = R_;
	else R2 = R_;
#endif

#ifdef KS
	pvals.cusum_pvals[1][pvals.seq_counter] = p_value;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
 if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (REVERSE) TEST\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", zrev);
		fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

		fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_CUSUM]);
		fprintf(results[TEST_CUSUM], "%f\n", p_value); fflush(results[TEST_CUSUM]);
	}
#endif
	//printf("sum1 %lf   sum2 %lf\n\n",sum1,sum2);
}
