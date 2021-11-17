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
            R A N D O M  E X C U R S I O N S  V A R I A N T  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
RandomExcursionsVariant(int n)
{
	int		i, p, J, x, constraint, count, *S_k;
	int		stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double	p_value;

#ifdef VERIFY_RESULTS
	R_.random_excursion_variant.valid=0;
#endif


	if ((S_k = (int *)calloc(n, sizeof(int))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\tRANDOM EXCURSIONS VARIANT: Insufficent memory allocated.\n");
		}
#endif
		return;
	}
	J = 0;
	S_k[0] = 2 * (int)epsilon[0] - 1;
	for (i = 1; i < n; i++) {
		S_k[i] = S_k[i - 1] + 2 * epsilon[i] - 1;
		if (S_k[i] == 0)
			J++;
	}
	if (S_k[n - 1] != 0)
		J++;

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\tRANDOM EXCURSIONS VARIANT TEST\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(a) Number Of Cycles (J) = %d\n", J);
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) Sequence Length (n)  = %d\n", n);
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
	}
#endif

	constraint = (int)MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RND_EXCURSION_VAR], "\n\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
			fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
			fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t---------------------------------------------\n");
			for (i = 0; i < 18; i++)
				fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", 0.0);
		}
#endif
	}
	else {
		for (p = 0; p <= 17; p++) {
			x = stateX[p];
			count = 0;
			for (i = 0; i < n; i++)
				if (S_k[i] == x)
					count++;
			//PRINT
			//printf("%d [%d]",count,J);
			p_value = erfc(fabs(count - J) / (sqrt(2.0*J*(4.0*fabs(x) - 2))));
#ifdef SPEED
			dummy_result += p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				if (isNegative(p_value) || isGreaterThanOne(p_value))
					fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) WARNING: P_VALUE IS OUT OF RANGE.\n");
				fprintf(stats[TEST_RND_EXCURSION_VAR], "%s\t\t", p_value < ALPHA ? "FAILURE" : "SUCCESS");
				fprintf(stats[TEST_RND_EXCURSION_VAR], "(x = %2d) Total visits = %4d; p-value = %f\n", x, count, p_value);
				fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", p_value); fflush(results[TEST_RND_EXCURSION_VAR]);
			}
#endif

#ifdef KS
			pvals.random_excursion_variant_pvals[p][pvals.seq_counter] = p_value;
#endif
#ifdef VERIFY_RESULTS
			R_.random_excursion_variant.valid=1;
			R_.random_excursion_variant.x[p]=x;
			R_.random_excursion_variant.count[p]=count;
			R_.random_excursion_variant.p_value[p]=p_value;
			if(RandomExcursionsVariant_v1 == RandomExcursionsVariant) R1 = R_;
			else R2 = R_;	
#endif

		}
	}
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\n"); fflush(stats[TEST_RND_EXCURSION_VAR]);
}
#endif
	//printf("\n\n");
	free(S_k);
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
RandomExcursionsVariant2(int n)
{
	int		i, p, J = 0, x, constraint, count, S_k = 0, bit_ind, window;
	int		stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double	p_value;
	int counter[19] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

#ifdef VERIFY_RESULTS
	R_.random_excursion_variant.valid=0;
#endif

	for (bit_ind = 0; bit_ind < n; bit_ind++)
	{
		window = get_nth_block4(array, bit_ind);
		S_k += (window & 1) * 2 - 1;
		if (S_k == 0)
		{
			++J;
		}
		if ((S_k >= -9) && (S_k <= 9)) counter[S_k + 9]++;
	}
	////Last bit
	if (S_k)
	{
		J++;
	}
	for (i = 9; i < 18; i++)
	{
		counter[i] = counter[i + 1];
	}
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\tRANDOM EXCURSIONS VARIANT TEST\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(a) Number Of Cycles (J) = %d\n", J);
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) Sequence Length (n)  = %d\n", n);
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
	}
#endif

	constraint = (int)MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RND_EXCURSION_VAR], "\n\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
			fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
			fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t---------------------------------------------\n");
			for (i = 0; i < 18; i++)
				fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", 0.0);
		}
#endif
#ifdef KS
		for (i = 0; i < 18; i++)
			pvals.random_excursion_variant_pvals[i][pvals.seq_counter] = 0.0;
#endif
	}
	else {
		for (p = 0; p <= 17; p++) {
			x = stateX[p];
			count = counter[p];
			//PRINT
			//printf("%d [%d]",count,J);
			p_value = erfc(fabs(count - J) / (sqrt(2.0*J*(4.0*fabs(x) - 2))));
#ifdef SPEED
			dummy_result += p_value;
#endif
#ifdef VERIFY_RESULTS
			R_.random_excursion_variant.valid=1;
			R_.random_excursion_variant.x[p]=x;
			R_.random_excursion_variant.count[p]=count;
			R_.random_excursion_variant.p_value[p]=p_value;
			if(RandomExcursionsVariant_v1 == RandomExcursionsVariant2) R1 = R_;
			else R2 = R_;	
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				if (isNegative(p_value) || isGreaterThanOne(p_value))
					fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) WARNING: P_VALUE IS OUT OF RANGE.\n");
				fprintf(stats[TEST_RND_EXCURSION_VAR], "%s\t\t", p_value < ALPHA ? "FAILURE" : "SUCCESS");
				fprintf(stats[TEST_RND_EXCURSION_VAR], "(x = %2d) Total visits = %4d; p-value = %f\n", x, count, p_value);
				fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", p_value); fflush(results[TEST_RND_EXCURSION_VAR]);
			}
#endif
#ifdef KS
			pvals.random_excursion_variant_pvals[p][pvals.seq_counter] = p_value;
#endif
		}
	}
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\n"); fflush(stats[TEST_RND_EXCURSION_VAR]);
}
#endif
}

