#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/cephes.h"
#include "../include/matrix.h"
#include "../include/stat_fncs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                              R A N K  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	
void
Rank(int n)
{
	int			N, i, k, r;
	double		p_value, product, chi_squared, arg1, p_32, p_31, p_30, R, F_32, F_31, F_30;
	BitSequence	**matrix = create_matrix(32, 32);
	
	N = n/(32*32);
	if (isZero(N)) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RANK], "\t\t\t\tRANK TEST\n");
			fprintf(stats[TEST_RANK], "\t\tError: Insuffucient # Of Bits To Define An 32x32 (%dx%d) Matrix\n", 32, 32);
	}
#endif
		p_value = 0.00;
	}
	else {
		r = 32;					/* COMPUTE PROBABILITIES */
		product = 1;
		for ( i=0; i<=r-1; i++ )
			product *= ((1.e0-pow(2, i-32))*(1.e0-pow(2, i-32)))/(1.e0-pow(2, i-r));
		p_32 = pow(2, r*(32+32-r)-32*32) * product;
		
		r = 31;
		product = 1;
		for ( i=0; i<=r-1; i++ )
			product *= ((1.e0-pow(2, i-32))*(1.e0-pow(2, i-32)))/(1.e0-pow(2, i-r));
		p_31 = pow(2, r*(32+32-r)-32*32) * product;
		
		p_30 = 1 - (p_32+p_31);
		
		F_32 = 0;
		F_31 = 0;
		for ( k=0; k<N; k++ ) {			/* FOR EACH 32x32 MATRIX   */
			def_matrix(32, 32, matrix, k);
#if (DISPLAY_MATRICES == 1)
			display_matrix(32, 32, matrix);
#endif
			R = computeRank(32, 32, matrix);
			if ( R == 32 )
				F_32++;			/* DETERMINE FREQUENCIES */
			if ( R == 31 )
				F_31++;
		}
		F_30 = (double)N - (F_32+F_31);
		
		chi_squared =(pow(F_32 - N*p_32, 2)/(double)(N*p_32) +
					  pow(F_31 - N*p_31, 2)/(double)(N*p_31) +
					  pow(F_30 - N*p_30, 2)/(double)(N*p_30));
		
		arg1 = -chi_squared/2.e0;

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RANK], "\t\t\t\tRANK TEST\n");
			fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_RANK], "\t\tCOMPUTATIONAL INFORMATION:\n");
			fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_RANK], "\t\t(a) Probability P_%d = %f\n", 32, p_32);
			fprintf(stats[TEST_RANK], "\t\t(b)             P_%d = %f\n", 31, p_31);
			fprintf(stats[TEST_RANK], "\t\t(c)             P_%d = %f\n", 30, p_30);
			fprintf(stats[TEST_RANK], "\t\t(d) Frequency   F_%d = %d\n", 32, (int)F_32);
			fprintf(stats[TEST_RANK], "\t\t(e)             F_%d = %d\n", 31, (int)F_31);
			fprintf(stats[TEST_RANK], "\t\t(f)             F_%d = %d\n", 30, (int)F_30);
			fprintf(stats[TEST_RANK], "\t\t(g) # of matrices    = %d\n", N);
			fprintf(stats[TEST_RANK], "\t\t(h) Chi^2            = %f\n", chi_squared);
			fprintf(stats[TEST_RANK], "\t\t(i) NOTE: %d BITS WERE DISCARDED.\n", n % (32 * 32));
			fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
		}
#endif

		p_value = exp(arg1);
#ifdef SPEED
		dummy_result = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			if (isNegative(p_value) || isGreaterThanOne(p_value))
				fprintf(stats[TEST_RANK], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
		}
#endif
		for ( i=0; i<32; i++ )				/* DEALLOCATE MATRIX  */
			free(matrix[i]);
		free(matrix);

#ifdef VERIFY_RESULTS
		R_.rank.p_30=p_30;
		R_.rank.p_31=p_31;
		R_.rank.p_32=p_32;
		R_.rank.F_30=F_30;
		R_.rank.F_31=F_31;
		R_.rank.F_32=F_32;
		R_.rank.chi_squared=chi_squared;
		R_.rank.p_value=p_value;
		R_.rank.N=N;
		if(Rank_v1 == Rank) R1 = R_;
		else R2 = R_;

#endif
	}
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RANK], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RANK]);
		fprintf(results[TEST_RANK], "%f\n", p_value); fflush(results[TEST_RANK]);
	}
#endif

#ifdef KS
	pvals.rank_pvals[pvals.seq_counter] = p_value;
#endif
	//printf("%lf %lf %lf \n",F_30,F_31,F_32);
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
Rank2(int n){

	unsigned int matrix[32];
	//unsigned int * p_array = (unsigned int*)array;
	unsigned char *p_array = array;
	int N, i, r;
	double p_value, product, chi_squared, arg1, p_32, p_31, p_30, R, F_32, F_31, F_30;

	N = n / 1024;

	if (isZero(N)) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RANK], "\t\t\t\tRANK TEST\n");
			fprintf(stats[TEST_RANK], "\t\tError: Insuffucient # Of Bits To Define An 32x32 (%dx%d) Matrix\n", 32, 32);
		}
#endif
		p_value = 0.00;
		return;
	}


	// COMPUTE PROBABILITIES
	r = 32;
	product = 1.0;
	for (i = 0; i <= r - 1; ++i)
		product *= ((1.e0 - pow((double)2.0, (double)i - 32))*(1.e0 - pow((double)2, (double)i - 32))) / (1.e0 - pow((double)2, (double)i - r));
	p_32 = pow((double)2, (double)r*(32 + 32 - r) - 32 * 32) * product;

	r = 31;
	product = 1;
	for (i = 0; i <= r - 1; ++i)
		product *= ((1.e0 - pow((double)2, (double)i - 32))*(1.e0 - pow((double)2, (double)i - 32))) / (1.e0 - pow((double)2, (double)i - r));
	p_31 = pow((double)2, (double)r*(32 + 32 - r) - 32 * 32) * product;

	p_30 = 1 - (p_32 + p_31);

	F_32 = 0;
	F_31 = 0;
	// FOR EACH 32x32 MATRIX
	for (i = 0; i < N; i++)
	{
		//if(i % 10000 == 0)cout << i << endl;
		memcpy(matrix, p_array, 128);
		/*for( i = 0; i < 32; i++)
		{
		matrix[i] = *(p_array++);
		}*/
		p_array = p_array + 128;

		R = Mrank(matrix);
		if (R == 32)F_32++;
		if (R == 31)F_31++;
	}

	F_30 = (double)N - (F_32 + F_31);

	chi_squared = (pow((double)F_32 - (double)(N)*p_32, 2.0) / ((double)(N)*p_32) +
		pow((double)F_31 - (double)(N)*p_31, 2.0) / ((double)(N)*p_31) +
		pow((double)F_30 - (double)(N)*p_30, 2.0) / ((double)(N)*p_30));

	arg1 = -chi_squared / 2.e0;
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RANK], "\t\t\t\tRANK TEST\n");
		fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_RANK], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_RANK], "\t\t(a) Probability P_%d = %f\n", 32, p_32);
		fprintf(stats[TEST_RANK], "\t\t(b)             P_%d = %f\n", 31, p_31);
		fprintf(stats[TEST_RANK], "\t\t(c)             P_%d = %f\n", 30, p_30);
		fprintf(stats[TEST_RANK], "\t\t(d) Frequency   F_%d = %d\n", 32, (int)F_32);
		fprintf(stats[TEST_RANK], "\t\t(e)             F_%d = %d\n", 31, (int)F_31);
		fprintf(stats[TEST_RANK], "\t\t(f)             F_%d = %d\n", 30, (int)F_30);
		fprintf(stats[TEST_RANK], "\t\t(g) # of matrices    = %d\n", N);
		fprintf(stats[TEST_RANK], "\t\t(h) Chi^2            = %f\n", chi_squared);
		fprintf(stats[TEST_RANK], "\t\t(i) NOTE: %d BITS WERE DISCARDED.\n", n % (32 * 32));
		fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
	}
#endif

	p_value = exp(arg1);
#ifdef SPEED
	dummy_result = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_RANK], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
	}
#endif

#ifdef VERIFY_RESULTS
	R_.rank.p_30=p_30;
	R_.rank.p_31=p_31;
	R_.rank.p_32=p_32;
	R_.rank.F_30=F_30;
	R_.rank.F_31=F_31;
	R_.rank.F_32=F_32;
	R_.rank.chi_squared=chi_squared;
	R_.rank.p_value=p_value;
	R_.rank.N=N;
	if(Rank_v1 == Rank2) R1 = R_;
	else R2 = R_;

#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RANK], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RANK]);
		fprintf(results[TEST_RANK], "%f\n", p_value); fflush(results[TEST_RANK]);
}
#endif

#ifdef KS
	pvals.rank_pvals[pvals.seq_counter] = p_value;
#endif
    //printf("%lf %lf %lf \n",F_30,F_31,F_32);
}
