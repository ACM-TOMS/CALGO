#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"  
#include "../include/BM.h"
#include "../include/BMA.h"
#include "../include/stat_fncs.h"

//100MB - M=1000 927sec
//		  M=5000 4500sec
void
LinearComplexity(int M, int n)
{
	int       i, ii, j, d, N, L, m, N_, parity, sign, K = 6, nu[7];
	double    p_value, T_, mean, chi2;
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	BitSequence  *T = NULL, *P = NULL, *B_ = NULL, *C = NULL;

	N = (int)floor(n / M);
	if (((B_ = (BitSequence *)calloc(M, sizeof(BitSequence))) == NULL) ||
		((C = (BitSequence *)calloc(M, sizeof(BitSequence))) == NULL) ||
		((P = (BitSequence *)calloc(M, sizeof(BitSequence))) == NULL) ||
		((T = (BitSequence *)calloc(M, sizeof(BitSequence))) == NULL)) {
		printf("Insufficient Memory for Work Space:: Linear Complexity Test\n");
		if (B_ != NULL)
			free(B_);
		if (C != NULL)
			free(C);
		if (P != NULL)
			free(P);
		if (T != NULL)
			free(T);
		return;
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);
	}
#endif

	for (i = 0; i < K + 1; i++)
		nu[i] = 0;
	for (ii = 0; ii < N; ii++) {
		for (i = 0; i < M; i++) {
			B_[i] = 0;
			C[i] = 0;
			T[i] = 0;
			P[i] = 0;
		}
		L = 0;
		m = -1;
		d = 0;
		C[0] = 1;
		B_[0] = 1;

		/* DETERMINE LINEAR COMPLEXITY */
		N_ = 0;
		while (N_ < M) {
			d = (int)epsilon[ii*M + N_];
			for (i = 1; i <= L; i++)
				d += C[i] * epsilon[ii*M + N_ - i];


			d = d % 2;
			if (d == 1) {
				for (i = 0; i < M; i++) {
					T[i] = C[i];
					P[i] = 0;
				}
				for (j = 0; j < M; j++)
					if (B_[j] == 1)
						P[j + N_ - m] = 1;
				for (i = 0; i < M; i++)
					C[i] = (C[i] + P[i]) % 2;
				if (L <= N_ / 2) {
					L = N_ + 1 - L;
					m = N_;
					for (i = 0; i < M; i++)
						B_[i] = T[i];
				}
			}
			N_++;
		}
		//printf("L = %d \n",L);

		if ((parity = (M + 1) % 2) == 0)
			sign = -1;
		else
			sign = 1;
		mean = M / 2.0 + (9.0 + sign) / 36.0 - 1.0 / pow(2, M) * (M / 3.0 + 2.0 / 9.0);
		if ((parity = M % 2) == 0)
			sign = 1;
		else
			sign = -1;
		T_ = sign * (L - mean) + 2.0 / 9.0;

		if (T_ <= -2.5)
			nu[0]++;
		else if (T_ > -2.5 && T_ <= -1.5)
			nu[1]++;
		else if (T_ > -1.5 && T_ <= -0.5)
			nu[2]++;
		else if (T_ > -0.5 && T_ <= 0.5)
			nu[3]++;
		else if (T_ > 0.5 && T_ <= 1.5)
			nu[4]++;
		else if (T_ > 1.5 && T_ <= 2.5)
			nu[5]++;
		else
			nu[6]++;
	}
	chi2 = 0.00;
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		for ( i=0; i<K+1; i++ ) 
			fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
}
#endif
	for ( i=0; i<K+1; i++ )
	{
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
		//printf("%d ",nu[i]);
#ifdef VERIFY_RESULTS
		R_.linear_complexity.nu[i]=nu[i];
#endif
	}

	p_value = cephes_igamc(K/2.0, chi2/2.0);
#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.linear_complexity.chi2=chi2;
	R_.linear_complexity.p_value=p_value;
	if(LinearComplexity_v1 == LinearComplexity) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_LINEARCOMPLEXITY]);
		fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value); fflush(results[TEST_LINEARCOMPLEXITY]);
	}
#endif

#ifdef KS
	pvals.linear_complexity_pvals[pvals.seq_counter] = p_value;
#endif
	free(B_);
	free(P);
	free(C);
	free(T);
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
LinearComplexity2(int M, int n)
{
	int       i, ii, j, N, L, parity, sign, K = 6, nu[7];
	double    p_value, T_, mean, chi2;
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	
	int array_size, type_size_bits,type_size_bytes;
	type *type_array,*c,*t,*b;

	type_size_bytes = sizeof(type);
	type_size_bits = sizeof(type)*8;

	array_size = M / type_size_bits + (M % type_size_bits != 0);

	type_array = (type*)malloc(array_size*type_size_bytes);
	c = (type*)malloc(array_size*type_size_bytes);
	t = (type*)malloc(array_size*type_size_bytes);
	b = (type*)malloc(array_size*type_size_bytes);

	if (type_array==NULL||c==NULL||t==NULL||b==NULL) {
		printf("Insufficient Memory for Work Space: Linear Complexity Test\n");
		if ( type_array!= NULL )
			free(type_array);
		if ( c != NULL )
			free(c);
		if ( t != NULL )
			free(t);
		if ( b != NULL )
			free(b);
		return;
	}
	
	N = (int)floor(n/M);

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);
	}
#endif
	for ( i=0; i<K+1; i++ )
		nu[i] = 0;
	for ( ii=0; ii<N; ii++ ) {
		
		j = 0;
		for(i = M*ii; i < M*(ii+1); i+= type_size_bits)
		{
			type_array[j] = get_nth_block_effect(array,i);
			j++;			
		}
		if((M % type_size_bits) != 0) type_array[j] &= ((M % type_size_bits) << 1);


		L = BM_c(type_array,M,c,b,t);
		if ( (parity = (M+1)%2) == 0 ) 
			sign = -1;
		else 
			sign = 1;
		mean = M/2.0 + (9.0+sign)/36.0 - 1.0/pow(2, M) * (M/3.0 + 2.0/9.0);
		if ( (parity = M%2) == 0 )
			sign = 1;
		else 
			sign = -1;
		T_ = sign * (L - mean) + 2.0/9.0;
		
		if ( T_ <= -2.5 )
			nu[0]++;
		else if ( T_ > -2.5 && T_ <= -1.5 )
			nu[1]++;
		else if ( T_ > -1.5 && T_ <= -0.5 )
			nu[2]++;
		else if ( T_ > -0.5 && T_ <= 0.5 )
			nu[3]++;
		else if ( T_ > 0.5 && T_ <= 1.5 )
			nu[4]++;
		else if ( T_ > 1.5 && T_ <= 2.5 )
			nu[5]++;
		else
			nu[6]++;
	}
	chi2 = 0.00;
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		for (i = 0; i < K + 1; i++)
			fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
	}
#endif
	for ( i=0; i<K+1; i++ )
	{
		chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
		//printf("%d ",nu[i]);
#ifdef VERIFY_RESULTS
		R_.linear_complexity.nu[i]=nu[i];
#endif
	}
	
	p_value = cephes_igamc(K/2.0, chi2/2.0);
#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.linear_complexity.chi2=chi2;
	R_.linear_complexity.p_value=p_value;
	if(LinearComplexity_v1 == LinearComplexity2) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_LINEARCOMPLEXITY]);
		fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value); fflush(results[TEST_LINEARCOMPLEXITY]);
	}
#endif

#ifdef KS
	pvals.linear_complexity_pvals[pvals.seq_counter] = p_value;
#endif
}

void
LinearComplexity3(int M, int n)
{
	int       i, ii, j, N, L, parity, sign, K = 6, nu[7];
	double    p_value, T_, mean, chi2;
	double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };
	int size;
	//double low, up;
	BMAint *d_b, *d_c, *d_t, *S;

	size = (M + sizeof(BMAint)* 8) / (sizeof(BMAint)* 8) + 4 /* pro jistotu: */ + 100;

	d_b = (BMAint*)malloc(sizeof(BMAint)*size);
	d_c = (BMAint*)malloc(sizeof(BMAint)*size);
	d_t = (BMAint*)malloc(sizeof(BMAint)*size);
	S = (BMAint*)malloc(sizeof(BMAint)*size);

	N = (int)floor(n / M);

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
		fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);
	}
#endif

	if ((parity = (M + 1) % 2) == 0)
		sign = -1;
	else
		sign = 1;
	mean = M / 2.0 + (9.0 + sign) / 36.0 - 1.0 / pow(2, M) * (M / 3.0 + 2.0 / 9.0);
	if ((parity = M % 2) == 0)
		sign = 1;
	else
		sign = -1;

	//low = (2.5 - 2.0 / 9.0) / sign + mean;
	//up = (-2.5 - 2.0 / 9.0) / sign + mean;


	for (i = 0; i<K + 1; i++)
		nu[i] = 0;

	for (ii = 0; ii<N; ii++) {

		j = 0;
		for (i = M*ii; i < M*(ii + 1); i += 32)
		{
			S[j] = get_nth_block_effect(array, i);
			j++;
		}
		//if((M % 32) != 0) S[j] &= ((M % 32) << 1);

		L = BM_JOURNAL(d_b, d_c, d_t, S, M);
		T_ = sign * (L - mean) + 2.0 / 9.0;


		if (T_ <= -2.5)
			nu[0]++;
		else if (T_ > -2.5 && T_ <= -1.5)
			nu[1]++;
		else if (T_ > -1.5 && T_ <= -0.5)
			nu[2]++;
		else if (T_ > -0.5 && T_ <= 0.5)
			nu[3]++;
		else if (T_ > 0.5 && T_ <= 1.5)
			nu[4]++;
		else if (T_ > 1.5 && T_ <= 2.5)
			nu[5]++;
		else
			nu[6]++;
	}
	chi2 = 0.00;
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		for (i = 0; i < K + 1; i++)
			fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
	}
#endif
	for (i = 0; i<K + 1; i++)
	{
		chi2 += pow(nu[i] - N*pi[i], 2) / (N*pi[i]);
		//printf("%d ",nu[i]);
#ifdef VERIFY_RESULTS
		R_.linear_complexity.nu[i] = nu[i];
#endif
	}

	p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
#ifdef SPEED
	dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
	R_.linear_complexity.chi2=chi2;
	R_.linear_complexity.p_value=p_value;
	if(LinearComplexity_v1 == LinearComplexity3) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value); fflush(stats[TEST_LINEARCOMPLEXITY]);
		fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value); fflush(results[TEST_LINEARCOMPLEXITY]);
	}
#endif

#ifdef KS
	pvals.linear_complexity_pvals[pvals.seq_counter] = p_value;
#endif
}