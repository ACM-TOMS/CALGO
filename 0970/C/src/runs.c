#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/cephes.h"
#include "../include/erf.h"
#include "../include/tools.h"
#include "../include/stat_fncs.h"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                              R U N S  T E S T 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
Runs(int n)
{

	int		S, k;
	double	pi, V, erfc_arg, p_value;

	S = 0;
	for (k = 0; k < n; k++)
		if (epsilon[k])
			S++;

	pi = (double)S / (double)n;

	if (fabs(pi - 0.5) > (2.0 / sqrt(n))) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
		}
#endif
		p_value = 0.0;

#ifdef VERIFY_RESULTS
		R_.runs.p_value=p_value;
		R_.runs.pi=pi;
		R_.runs.V=0.0;
		R_.runs.erfc_arg=0.0;
#endif
	}
	else {

		V = 1;
		for (k = 1; k < n; k++)
			if (epsilon[k] != epsilon[k - 1])
				V++;

		erfc_arg = fabs(V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * sqrt(2 * n));
		p_value = erfc(erfc_arg);
#ifdef SPEED
		dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
		R_.runs.p_value=p_value;
		R_.runs.pi=pi;
		R_.runs.V=V;
		R_.runs.erfc_arg=erfc_arg;
		if(Runs_v1 == Runs) R1 = R_;
		else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tCOMPUTATIONAL INFORMATION:\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\t(a) Pi                        = %f\n", pi);
			fprintf(stats[TEST_RUNS], "\t\t(b) V_n_obs (Total # of runs) = %d\n", (int)V);
			fprintf(stats[TEST_RUNS], "\t\t(c) V_n_obs - 2 n pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t    -----------------------   = %f\n", erfc_arg);
			fprintf(stats[TEST_RUNS], "\t\t      2 sqrt(2n) pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			if (isNegative(p_value) || isGreaterThanOne(p_value))
				fprintf(stats[TEST_RUNS], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_RUNS], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RUNS]);
		}
#endif
	}

#ifdef KS
	pvals.runs_pvals[pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(results[TEST_RUNS], "%f\n", p_value); fflush(results[TEST_RUNS]);
}
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
Runs2(int n)
{
	int		S, k, i, diff;
	double	pi, erfc_arg, p_value;
	unsigned long int V;
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
		short int LU_switches[256] =
		{0,1,2,1,2,3,2,1,2,3,4,3,2,3,2,1,2,3,4,3,4,5,4,3,2,3,4,3,2,3,2,1
		,2,3,4,3,4,5,4,3,4,5,6,5,4,5,4,3,2,3,4,3,4,5,4,3,2,3,4,3,2,3,2,1
		,2,3,4,3,4,5,4,3,4,5,6,5,4,5,4,3,4,5,6,5,6,7,6,5,4,5,6,5,4,5,4,3
		,2,3,4,3,4,5,4,3,4,5,6,5,4,5,4,3,2,3,4,3,4,5,4,3,2,3,4,3,2,3,2,1
		,1,2,3,2,3,4,3,2,3,4,5,4,3,4,3,2,3,4,5,4,5,6,5,4,3,4,5,4,3,4,3,2
		,3,4,5,4,5,6,5,4,5,6,7,6,5,6,5,4,3,4,5,4,5,6,5,4,3,4,5,4,3,4,3,2
		,1,2,3,2,3,4,3,2,3,4,5,4,3,4,3,2,3,4,5,4,5,6,5,4,3,4,5,4,3,4,3,2
		,1,2,3,2,3,4,3,2,3,4,5,4,3,4,3,2,1,2,3,2,3,4,3,2,1,2,3,2,1,2,1,0} ;
		*/
	S = 0;
	for (k = 0; k * 8 < n; k++)
	{
		S += LU_byte_weight[array[k]];
		//sum += LU_byte_weight[(get_nth_block4(array,i*8) & 255)];
	}
	if (k * 8 > n)S -= LU_byte_weight[array[k - 1] >> (8 - (k * 8 - n))];


	pi = (double)S / (double)n;

	if (fabs(pi - 0.5) > (2.0 / sqrt(n))) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
		}
#endif
		p_value = 0.0;

#ifdef VERIFY_RESULTS
		R_.runs.p_value=p_value;
		R_.runs.pi=pi;
		R_.runs.V=0.0;
		R_.runs.erfc_arg=0.0;
#endif


	}
	else {

		V = 1;
		for (k = 0; (k + 1) * 8 < n; k++)
		{
			V += LU_byte_switches[array[k]];
			diff = ((array[k] >> 7) - (array[k + 1] & 1));
			V = V + diff*diff;

			//if((array[k]>>7) != (array[k] & 1) ) V++;
			//sum += LU_byte_weight[(get_nth_block4(array,i*8) & 255)];
		}
		//if((get_nth_block4(array,k*8)&1) == (get_nth_block4(array,k*8+1)&1)) V--;

		for (i = k * 8; i + 1 < n; i++)
		{
			if ((get_nth_block4(array, i) & 1) != (get_nth_block4(array, i + 1) & 1)) V++;
		}

		erfc_arg = fabs((double)V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * sqrt(2 * n));
		p_value = erfc(erfc_arg);

#ifdef SPEED
		dummy_result = p_value;
#endif
#ifdef VERIFY_RESULTS
		R_.runs.p_value=p_value;
		R_.runs.pi=pi;
		R_.runs.V=V;
		R_.runs.erfc_arg=erfc_arg;
		if(Runs_v1 == Runs2) R1 = R_;
		else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tCOMPUTATIONAL INFORMATION:\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\t(a) Pi                        = %f\n", pi);
			fprintf(stats[TEST_RUNS], "\t\t(b) V_n_obs (Total # of runs) = %d\n", (int)V);
			fprintf(stats[TEST_RUNS], "\t\t(c) V_n_obs - 2 n pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t    -----------------------   = %f\n", erfc_arg);
			fprintf(stats[TEST_RUNS], "\t\t      2 sqrt(2n) pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			if (isNegative(p_value) || isGreaterThanOne(p_value))
				fprintf(stats[TEST_RUNS], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_RUNS], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RUNS]);
		}
#endif
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(results[TEST_RUNS], "%f\n", p_value); fflush(results[TEST_RUNS]);
}
#endif
	
#ifdef KS
	pvals.runs_pvals[pvals.seq_counter] = p_value;
#endif
	//printf("R_: %d %lf \n",S,(double)V);
	
}


void
Runs3(int n)
{
	int		S, i;
	double	pi, erfc_arg, p_value;
	unsigned long int V;
	unsigned int mask;
	unsigned char *p_tmp, *p_end;

	int LUT_HW_size = 16;
	int LUT_HW_Bsize = 2;
	signed char *LUT_HW = LUT_HW_16;

	int LUT_Switches_size = 16;
	int LUT_Switches_Bsize = 2;
	signed char *LUT_Switches = LUT_Switches_16;

	if (0)
	{
		LUT_HW_size = 8;
		LUT_HW_Bsize = 1;
		LUT_HW = LUT_HW_8;

		LUT_Switches_size = 8;
		LUT_Switches_Bsize = 1;
		LUT_Switches = LUT_Switches_8;
	}

	mask = get_mask(LUT_HW_size);

	//count bits
	S = 0;
	mask = get_mask(LUT_HW_size);

	p_end = array + (n - (n%LUT_HW_size)) / 8;
	for (p_tmp = array; p_tmp < p_end; p_tmp += LUT_HW_Bsize){
		S += LUT_HW[*((unsigned int*)p_tmp) & mask];
	}
	if (n % LUT_HW_size){
		S += LUT_HW[*((unsigned int*)p_tmp) & get_mask(n % LUT_HW_size)];
	}

	pi = (double)S / (double)n;
	if (fabs(pi - 0.5) > (2.0 / sqrt(n))) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
		}
#endif
		p_value = 0.0;

#ifdef VERIFY_RESULTS
		R_.runs.p_value = p_value;
		R_.runs.pi = pi;
		R_.runs.V = 0.0;
		R_.runs.erfc_arg = 0.0;
#endif


	}
	else {


		V = 1;

		mask = get_mask(LUT_Switches_size + 1);
		p_end = array + ((n - 1) - (n - 1) % LUT_Switches_size) / 8 - LUT_Switches_Bsize + 1;

		for (p_tmp = array; p_tmp < p_end; p_tmp += LUT_Switches_Bsize){
			V += LUT_Switches[*((unsigned int*)p_tmp) & mask];
		}

		//last incomplete block
		if ((n - 1) % LUT_Switches_size != 0){
			for (i = n - ((n - 1) % LUT_Switches_size); i < n; i++)
			{
				if ((get_nth_block_effect(array, i - 1) & 1) != (get_nth_block_effect(array, i) & 1))
					V++;
			}

		}



		erfc_arg = fabs((double)V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * sqrt(2 * n));
		p_value = erfc(erfc_arg);
#ifdef SPEED
		dummy_result = p_value;
#endif

#ifdef VERIFY_RESULTS
		R_.runs.p_value=p_value;
		R_.runs.pi=pi;
		R_.runs.V=V;
		R_.runs.erfc_arg=erfc_arg;
		if(Runs_v1 == Runs3) R1 = R_;
		else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tCOMPUTATIONAL INFORMATION:\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\t(a) Pi                        = %f\n", pi);
			fprintf(stats[TEST_RUNS], "\t\t(b) V_n_obs (Total # of runs) = %d\n", (int)V);
			fprintf(stats[TEST_RUNS], "\t\t(c) V_n_obs - 2 n pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t    -----------------------   = %f\n", erfc_arg);
			fprintf(stats[TEST_RUNS], "\t\t      2 sqrt(2n) pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			if (isNegative(p_value) || isGreaterThanOne(p_value))
				fprintf(stats[TEST_RUNS], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_RUNS], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RUNS]);
		}
#endif
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(results[TEST_RUNS], "%f\n", p_value); fflush(results[TEST_RUNS]);
}
#endif
#ifdef KS
	pvals.runs_pvals[pvals.seq_counter] = p_value;
#endif
}


//based on Alin Suciu suggestions and bit tricks from http://graphics.stanford.edu/~seander/bithacks.html
void
Runs4(int n)
{
	int		S=0, i, j;
	double	pi, erfc_arg, p_value;
	double V;

	//unsigned int NumberOfRuns = 0;
	
	int LUT_HW_size = 16;
	signed char *LUT_HW = LUT_HW_16;


	const int Tsize = 64;
	uint64_t *pblock;
	uint64_t help;
	const unsigned int mask = (1 << LUT_HW_size) - 1;

	S = 0;

	pblock = (uint64_t*)array;
	help = *pblock;
	
	for (i = 0; i < n + 1 - Tsize; i += Tsize) {
		for (j = 0; j < Tsize / LUT_HW_size; j++){
			S += LUT_HW[help & mask];
			help >>= LUT_HW_size;
		}
		help = *(++pblock);
	}

	for (; i < n; i++) {
		help = get_nth_block4(array, i);
		S += help & 1;
	}

	pi = (double)S / (double)n;

	if (fabs(pi - 0.5) > (2.0 / sqrt(n))) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
	}
#endif
		p_value = 0.0;

#ifdef VERIFY_RESULTS
		R_.runs.p_value = p_value;
		R_.runs.pi = pi;
		R_.runs.V = 0.0;
		R_.runs.erfc_arg = 0.0;
#endif
	}
	else {

		V = 1;
		pblock = (uint64_t*)array;
		
		for (i = 0; i < n - Tsize - 1; i += Tsize)
		{
			help = *pblock;
			++pblock;

			help = help ^ (help >> 1) ^ (*pblock << (Tsize - 1));

			V += LUT_HW_16[help & 0xFFFF]; help >>= LUT_HW_size;
			V += LUT_HW_16[help & 0xFFFF]; help >>= LUT_HW_size;
			V += LUT_HW_16[help & 0xFFFF]; help >>= LUT_HW_size;
			V += LUT_HW_16[help & 0xFFFF]; 	
		}

		for (; i + 1 < n; i++)
		{
			if ((get_nth_block_effect(array, i ) & 1) != (get_nth_block_effect(array, i + 1) & 1))
				V++;
		}

		erfc_arg = fabs((double)V - 2.0 * n * pi * (1 - pi)) / (2.0 * pi * (1 - pi) * sqrt(2 * n));
		p_value = erfc(erfc_arg);
		//printf("Pval: %lf sum %lf", p_value, V);
#ifdef SPEED
		dummy_result = p_value;
#endif

#ifdef VERIFY_RESULTS
		R_.runs.p_value = p_value;
		R_.runs.pi = pi;
		R_.runs.V = V;
		R_.runs.erfc_arg = erfc_arg;
		if (Runs_v1 == Runs4) R1 = R_;
		else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RUNS], "\t\t\t\tRUNS TEST\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\tCOMPUTATIONAL INFORMATION:\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			fprintf(stats[TEST_RUNS], "\t\t(a) Pi                        = %f\n", pi);
			fprintf(stats[TEST_RUNS], "\t\t(b) V_n_obs (Total # of runs) = %d\n", (int)V);
			fprintf(stats[TEST_RUNS], "\t\t(c) V_n_obs - 2 n pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t    -----------------------   = %f\n", erfc_arg);
			fprintf(stats[TEST_RUNS], "\t\t      2 sqrt(2n) pi (1-pi)\n");
			fprintf(stats[TEST_RUNS], "\t\t------------------------------------------\n");
			if (isNegative(p_value) || isGreaterThanOne(p_value))
				fprintf(stats[TEST_RUNS], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_RUNS], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_RUNS]);
		}
#endif
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(results[TEST_RUNS], "%f\n", p_value); fflush(results[TEST_RUNS]);
	}
#endif

#ifdef KS
	pvals.runs_pvals[pvals.seq_counter] = p_value;
#endif
}