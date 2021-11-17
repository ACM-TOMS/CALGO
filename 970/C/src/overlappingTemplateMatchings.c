#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  
#include "../include/tools.h"  
#include "../include/stat_fncs.h"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
               O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
double	Pr(int u, double eta);

void
OverlappingTemplateMatchings(int m, int n)
{
	int				i, k, match;
	double			W_obs, eta, sum, chi2, p_value, lambda;
	int				M, N, j, K = 5;
	unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 };
	//double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 }; //old incorrect
	double			pi[6] = { 0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865 };

	
	BitSequence		*sequence;

	M = 1032;
	N = n/M;
	
	if ((sequence = (BitSequence *)calloc(m, sizeof(BitSequence))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
			fprintf(stats[TEST_OVERLAPPING], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_OVERLAPPING], "\t\tTEMPLATE DEFINITION:  Insufficient memory, Overlapping Template Matchings test aborted!\n");
	}
#endif
		return;
	}
	else
		for ( i=0; i<m; i++ )
			sequence[i] = 1;
	
	lambda = (double)(M-m+1)/pow(2,m);
	eta = lambda/2.0;
	sum = 0.0;
	for ( i=0; i<K; i++ ) {			/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;

	for ( i=0; i<N; i++ ) {
		W_obs = 0;
		for ( j=0; j<M-m+1; j++ ) {
			match = 1;
			for ( k=0; k<m; k++ ) {
				if ( sequence[k] != epsilon[i*M+j+k] )
					match = 0;
			}
			if ( match == 1 )
				W_obs++;
		}
		if ( W_obs <= 4 )
			nu[(int)W_obs]++;
		else
			nu[K]++;
	}
	sum = 0;
	chi2 = 0.0;                                   /* Compute Chi Square */
	for ( i=0; i<K+1; i++ ) {
		chi2 += pow((double)nu[i] - (double)N*pi[i], 2)/((double)N*pi[i]);
		sum += nu[i];
		//printf("%d ",nu[i]); // TODO: comment
	}
	p_value = cephes_igamc(K/2.0, chi2/2.0);
	//printf("chi2:%lf value: %lf",chi2,p_value); // TODO: comment
#ifdef SPEED
	dummy_result = p_value + chi2 + nu[0] + nu[1] + nu[2] + nu[3] + nu[4] + nu[5];
#endif
#ifdef VERIFY_RESULTS
	R_.overlapping.chi2=chi2;
	R_.overlapping.p_value=p_value;
	for(i=0;i<6;i++) R_.overlapping.nu[i]=nu[i];
	if(OverlappingTemplateMatchings_v1 == OverlappingTemplateMatchings) R1 = R_;
	else R2 = R_;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
			nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
	}
#endif

	free(sequence);
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_OVERLAPPING], "%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS"); fflush(stats[TEST_OVERLAPPING]);
		fprintf(results[TEST_OVERLAPPING], "%f\n", p_value); fflush(results[TEST_OVERLAPPING]);
	}
#endif
#ifdef KS
	pvals.overlapping_pvals[pvals.seq_counter] = p_value;
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
OverlappingTemplateMatchings2(int m, int n)
{
	int				i;
	double			W_obs, eta, sum, chi2, p_value, lambda;
	int				M, K = 5;
	unsigned int	N, nu[6] = { 0, 0, 0, 0, 0, 0 },sequence;
	//double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 }; //old incorrect
	double			pi[6] = { 0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865 };

	unsigned int    window, /* **Wj = NULL, */ mask;
	unsigned int	bit_ind,byte_ind,block, byte_size;

	//byte_size=n/8 + (n%8!=0);
	byte_size=(unsigned int)ceil(n/8)+4;
	M = 1032;
	N = n/M;
	sequence = (1 << m)-1;
	mask = (1 << m)-1;

	lambda = (double)(M-m+1)/pow(2,m);
	eta = lambda/2.0;
	sum = 0.0;
	for ( i=0; i<K; i++ ) {			/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;

	window = 0;
	for(byte_ind = 0; (byte_ind < 4) && (byte_ind < byte_size); byte_ind++ )
	{
		window ^= (array[byte_ind] << byte_ind*8) ;
	}
	//bits(&window,4);
	for(block = 0; block < N; block++)
	{
		W_obs = 0;
		for(bit_ind = block*M; bit_ind < (block+1)*M; bit_ind++)
		{
			//bits(&window,4);
			if(((mask & window)== sequence) && ((block+1)*M-m+1 > bit_ind)){
					 ++W_obs;
					//printf("%d",window & mask);
			}
			window = window >>  1;
			if( ((bit_ind & 7) == 7) && (byte_ind < byte_size) ){
				window ^= (array[byte_ind++] << 24) ;
				//printf(" ");
			}			
		}
		if ( W_obs <= 4 )
			nu[(int)W_obs]++;
		else
			nu[K]++;
	}


	sum = 0;
	chi2 = 0.0;                                   /* Compute Chi Square */
	for ( i=0; i<K+1; i++ ) {
		chi2 += pow((double)nu[i] - (double)N*pi[i], 2)/((double)N*pi[i]);
		sum += nu[i];
		//printf("%d ",nu[i]);
	}
	p_value = cephes_igamc(K/2.0, chi2/2.0);
	//printf("chi2:%lf value: %lf",chi2,p_value);
#ifdef SPEED
	dummy_result = p_value + chi2 + nu[0] + nu[1] + nu[2] + nu[3] + nu[4] + nu[5];
#endif
#ifdef VERIFY_RESULTS
	R_.overlapping.chi2=chi2;
	R_.overlapping.p_value=p_value;
	for(i=0;i<6;i++) R_.overlapping.nu[i]=nu[i];
	if(OverlappingTemplateMatchings_v1 == OverlappingTemplateMatchings2) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
			nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

		fprintf(stats[TEST_OVERLAPPING], "%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS"); fflush(stats[TEST_OVERLAPPING]);
		fprintf(results[TEST_OVERLAPPING], "%f\n", p_value); fflush(results[TEST_OVERLAPPING]);
	}
#endif
#ifdef KS
	pvals.overlapping_pvals[pvals.seq_counter] = p_value;
#endif
}



void
OverlappingTemplateMatchings3(int m, int n) // formerly _effective
{
	int				i;
	double			W_obs, eta, sum, chi2, p_value, lambda;
	int				M, N, K = 5;
	unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 }, sequence;
	//double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 }; //old incorrect
	double			pi[6] = { 0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865 };
	unsigned int    window, /* **Wj = NULL, */ mask;
	int				bit_ind, block;

	M = 1032;
	N = n / M;
	sequence = (1 << m) - 1;
	mask = (1 << m) - 1;

	lambda = (double)(M - m + 1) / pow(2, m);
	eta = lambda / 2.0;
	sum = 0.0;
	for (i = 0; i<K; i++) {			/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;


	//bits(&window,4);
	for (block = 0; block < N; block++)
	{
		W_obs = 0;
		for (bit_ind = block*M; bit_ind < (block + 1)*M - m + 1; bit_ind++)
		{
			window = get_nth_block4(array, bit_ind);
			//bits(&window,4);
			if ((mask & window) == sequence){
				++W_obs;
				//printf("%d",window & mask);
			}

		}
		if (W_obs <= 4)
			nu[(int)W_obs]++;
		else
			nu[K]++;
	}


	sum = 0;
	chi2 = 0.0;                                   /* Compute Chi Square */
	for (i = 0; i<K + 1; i++) {
		chi2 += pow((double)nu[i] - (double)N*pi[i], 2) / ((double)N*pi[i]);
		sum += nu[i];
		//printf("%d ",nu[i]);
	}
	p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
	//printf("chi2:%lf value: %lf",chi2,p_value);
#ifdef SPEED
	dummy_result = p_value + chi2 + nu[0] + nu[1] + nu[2] + nu[3] + nu[4] + nu[5];
#endif
#ifdef VERIFY_RESULTS
	R_.overlapping.chi2 = chi2;
	R_.overlapping.p_value = p_value;
	for (i = 0; i<6; i++) R_.overlapping.nu[i] = nu[i];
	if (OverlappingTemplateMatchings_v1 == OverlappingTemplateMatchings3) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
			nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

		fprintf(stats[TEST_OVERLAPPING], "%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS"); fflush(stats[TEST_OVERLAPPING]);
		fprintf(results[TEST_OVERLAPPING], "%f\n", p_value); fflush(results[TEST_OVERLAPPING]);
	}
#endif
#ifdef KS
	pvals.overlapping_pvals[pvals.seq_counter] = p_value;
#endif
}


void
OverlappingTemplateMatchings4(int m, int n)
{
	int				i;
	double			W_obs, eta, sum, chi2, p_value, lambda;
	int				M, N, K = 5;
	unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 }, sequence;
	//double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 }; //old incorrect
	double			pi[6] = { 0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865 };

	unsigned int    window, /* **Wj = NULL, */ mask;
	int				block;
	unsigned char* pbyte;

	M = 1032;
	N = n / M;
	sequence = (1 << m) - 1;
	mask = (1 << m) - 1;

	lambda = (double)(M - m + 1) / pow(2, m);
	eta = lambda / 2.0;
	sum = 0.0;
	for (i = 0; i<K; i++) {			/* Compute Probabilities */
		pi[i] = Pr(i, eta);
		sum += pi[i];
	}
	pi[K] = 1 - sum;


	//bits(&window,4);
	for (block = 0; block < N; block++)
	{
		W_obs = 0;
		/*for(bit_ind = block*M; bit_ind < (block+1)*M-m+1; bit_ind++)
		{
		window = get_nth_block4(array,bit_ind);
		//bits(&window,4);
		if((mask & window)== sequence){
		++W_obs;
		//printf("%d",window & mask);
		}

		}*/
		i = block*M;
		while (i % 8 != 0 && i < (block+1)*M - m + 1){
			window = get_nth_block4(array, i);
			if ((window & mask) == sequence) ++W_obs;
			i++;
		}

		pbyte = array + i / 8;
		window = get_nth_block4(array, i);

		for (; i < (block + 1)*M - m + 1 - 8; i += 8) {
			if ((window & mask) == sequence) ++W_obs; window >>= 1;
			if ((window & mask) == sequence) ++W_obs; window >>= 1;
			if ((window & mask) == sequence) ++W_obs; window >>= 1;
			if ((window & mask) == sequence) ++W_obs; window >>= 1;
			if ((window & mask) == sequence) ++W_obs; window >>= 1;
			if ((window & mask) == sequence) ++W_obs; window >>= 1;
			if ((window & mask) == sequence) ++W_obs; window >>= 1;
			if ((window & mask) == sequence) ++W_obs;
			window = *(unsigned int*)(++pbyte);
		}

		for (; i < (block + 1)*M - m + 1; i++) {
			window = get_nth_block4(array, i);
			if ((window & mask) == sequence) ++W_obs;
		}

		if (W_obs <= 4)
			nu[(int)W_obs]++;
		else
			nu[K]++;
	}


	sum = 0;
	chi2 = 0.0;                                   /* Compute Chi Square */
	for (i = 0; i<K + 1; i++) {
		chi2 += pow((double)nu[i] - (double)N*pi[i], 2) / ((double)N*pi[i]);
		sum += nu[i];
		//printf("%d ",nu[i]);
	}
	p_value = cephes_igamc(K / 2.0, chi2 / 2.0);
	//printf("chi2:%lf value: %lf",chi2,p_value);
#ifdef SPEED
	dummy_result = p_value + chi2 + nu[0] + nu[1] + nu[2] + nu[3] + nu[4] + nu[5];
#endif
#ifdef VERIFY_RESULTS
	R_.overlapping.chi2 = chi2;
	R_.overlapping.p_value = p_value;
	for (i = 0; i<6; i++) R_.overlapping.nu[i] = nu[i];
	if (OverlappingTemplateMatchings_v1 == OverlappingTemplateMatchings4) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
		fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
		fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
			nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

		if (isNegative(p_value) || isGreaterThanOne(p_value))
			fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

		fprintf(stats[TEST_OVERLAPPING], "%f %s\n\n", p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS"); fflush(stats[TEST_OVERLAPPING]);
		fprintf(results[TEST_OVERLAPPING], "%f\n", p_value); fflush(results[TEST_OVERLAPPING]);
	}
#endif
#ifdef KS
	pvals.overlapping_pvals[pvals.seq_counter] = p_value;
#endif
}

double
Pr(int u, double eta)
{
	int		l;
	double	sum, p;
	
	if ( u == 0 )
		p = exp(-eta);
	else {
		sum = 0.0;
		for ( l=1; l<=u; l++ )
			sum += exp(-eta-u*log(2)+l*log(eta)-cephes_lgam(l+1)+cephes_lgam(u)-cephes_lgam(l)-cephes_lgam(u-l+1));
		p = sum;
	}
	return p;
}
