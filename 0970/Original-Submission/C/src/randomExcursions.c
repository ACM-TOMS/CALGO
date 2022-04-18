#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/cephes.h"  
#include "../include/tools.h" 
#include "../include/stat_fncs.h"

#define MAXCYCLES 1000

// random excursion in original version does not work with i=1999 and j=4

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                     R A N D O M  E X C U R S I O N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
RandomExcursions(int n)
{
	int		b, i, j, k, J, x;
	int		cycleStart, cycleStop, *cycle = NULL, *S_k = NULL;
	int		stateX[8] = { -4, -3, -2, -1, 1, 2, 3, 4 };
	int		counter[8] = { 0, 0, 0, 0, 0, 0, 0, 0 },nu[6][8];
	double	p_value, sum, constraint;
	double	pi[5][6] = { {0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000}, 
						 {0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000},
						 {0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625},
						 {0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143},
						 {0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051} };
	//int maxcyc;
	
#ifdef VERIFY_RESULTS
	R_.random_excursion.valid=0;
#endif
	//maxcyc=MAX(MAXCYCLES, n/100);

	if ( ((S_k = (int *)calloc(n, sizeof(int))) == NULL) ||
		 ((cycle = (int *)calloc(1+MAX(MAXCYCLES, n/100), sizeof(int))) == NULL) ) {  // fixed bug accessing later one more element :-)
		printf("Random Excursions Test:  Insufficent Work Space Allocated.\n");
		if ( S_k != NULL )
			free(S_k);
		if ( cycle != NULL )
			free(cycle);
#ifdef VERIFY_RESULTS
		if (RandomExcursions_v1 == RandomExcursions) R1 = R_;
		else R2 = R_;
#endif
		return;
	}
	
	J = 0; 					/* DETERMINE CYCLES */
	S_k[0] = 2*(int)epsilon[0] - 1;
	for( i=1; i<n; i++ ) {
		S_k[i] = S_k[i-1] + 2*epsilon[i] - 1;
		if ( S_k[i] == 0 ) {
			J++;
			if ( J >= MAX(MAXCYCLES, n/100) ) {
				printf("ERROR IN FUNCTION randomExcursions:  EXCEEDING THE MAX NUMBER OF CYCLES EXPECTED\n.");
				free(S_k);
				free(cycle);
#ifdef VERIFY_RESULTS
				if (RandomExcursions_v1 == RandomExcursions) R1 = R_;
				else R2 = R_;
#endif
				return;
			}
			//if(J>=maxcyc)printf("STOP1: accesing wrong element (n=%i).\n",n);
			cycle[J] = i;
		}
	}
	if ( S_k[n-1] != 0 )
		J++;
	//if(J>=maxcyc)printf("STOP2: accesing wrong element (n=%i) [%i vs %i].\n",n,J,maxcyc);
	cycle[J] = n;

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  RANDOM EXCURSIONS TEST\n");
		fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_RND_EXCURSION], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
		fprintf(stats[TEST_RND_EXCURSION], "\t\t(a) Number Of Cycles (J) = %04d\n", J);
		fprintf(stats[TEST_RND_EXCURSION], "\t\t(b) Sequence Length (n)  = %d\n", n);
	}
#endif

	constraint = MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
			for(i = 0; i < 8; i++)
				fprintf(results[TEST_RND_EXCURSION], "%f\n", 0.0);
	}
#endif
		//printf("WARNING: INSUFFICIENT NUMBER OF CYCLES (n=%i, J=%i, constraint=%i).\n",n,J,(int)constraint);
	}
	else {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_RND_EXCURSION], "\t\t(c) Rejection Constraint = %f\n", constraint);
			fprintf(stats[TEST_RND_EXCURSION], "\t\t-------------------------------------------\n");
	}
#endif
		cycleStart = 0;
		//if(1>=maxcyc)printf("STOP3: accesing wrong element (n=%i).\n",n);
		cycleStop  = cycle[1];
		for ( k=0; k<6; k++ )
			for ( i=0; i<8; i++ )
				nu[k][i] = 0;
		for ( j=1; j<=J; j++ ) {                           /* FOR EACH CYCLE */
			for ( i=0; i<8; i++ )
				counter[i] = 0;
			for ( i=cycleStart; i<cycleStop; i++ ) {
				if ( (S_k[i] >= 1 && S_k[i] <= 4) || (S_k[i] >= -4 && S_k[i] <= -1) ) {
					if ( S_k[i] < 0 )
						b = 4;
					else
						b = 3;
					counter[S_k[i]+b]++;
				}
			}
			//if(j>=maxcyc)printf("STOP4: accesing wrong element (n=%i)  [%i vs %i].\n",n,j,maxcyc);
			cycleStart = cycle[j]+1;
			if ( j < J )
			{
				//if((j+1)>=maxcyc)printf("STOP5: accesing wrong element (n=%i)  [%i vs %i].\n",n,j+1,maxcyc);
				cycleStop = cycle[j+1];
			}
			
			for ( i=0; i<8; i++ ) {
				if ( (counter[i] >= 0) && (counter[i] <= 4) )
					nu[counter[i]][i]++;
				else if ( counter[i] >= 5 )
					nu[5][i]++;
			}
		}
		
		//PRINT
		
		/*
		printf("Cycle count = %d \n",J);
		for( i = 0; i < 8; i++)
		{
			for ( k=0; k<6; k++ )
			{
				printf("%5d",nu[k][i]);
			}
			printf("\n");
		}
		*/
		

		for (i = 0; i < 8; i++) {
			x = stateX[i];
			sum = 0.;
			for (k = 0; k < 6; k++)
				sum += pow(nu[k][i] - J*pi[(int)fabs(x)][k], 2) / (J*pi[(int)fabs(x)][k]);
			p_value = cephes_igamc(2.5, sum / 2.0);
#ifdef SPEED
			dummy_result += p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				if ( isNegative(p_value) || isGreaterThanOne(p_value) )
					fprintf(stats[TEST_RND_EXCURSION], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

				fprintf(stats[TEST_RND_EXCURSION], "%s\t\tx = %2d chi^2 = %9.6f p_value = %f\n",
					p_value < ALPHA ? "FAILURE" : "SUCCESS", x, sum, p_value);
				fprintf(results[TEST_RND_EXCURSION], "%f\n", p_value); fflush(results[TEST_RND_EXCURSION]);
		}
#endif
#ifdef KS
			pvals.random_excursion_pvals[i][pvals.seq_counter] = p_value;
#endif
#ifdef VERIFY_RESULTS
			R_.random_excursion.valid=1;
			R_.random_excursion.J[i]=J;
			R_.random_excursion.x[i]=x;
			R_.random_excursion.p_value[i]=p_value;
			R_.random_excursion.sum[i]=sum;
#endif
		}
		}
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RND_EXCURSION], "\n"); fflush(stats[TEST_RND_EXCURSION]);
	}
#endif
	
	free(S_k);
	free(cycle);

#ifdef VERIFY_RESULTS
	if (RandomExcursions_v1 == RandomExcursions) R1 = R_;
	else R2 = R_;
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
RandomExcursions2(int n)
{
	int		i, k, J = 0, x;
	int     bit_ind, window;
	int		stateX[8] = { -4, -3, -2, -1, 1, 2, 3, 4 };
	int		counter[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 }, nu[6][9];
	int		S_k = 0;
	double	p_value, sum, constraint;
	double	pi[5][6] = { { 0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000 },
	{ 0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000 },
	{ 0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625 },
	{ 0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143 },
	{ 0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051 } };

	//int max_minus[256] = {-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-7,-5,-5,-3,-5,-3,-3,-1,-5,-3,-3,-1,-3,-1,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0};
	//int max_plus[256] = {0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,1,3,1,3,3,5,1,3,3,5,3,5,5,7,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8};
	//int byte_sum[256] = {-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8};

#ifdef VERIFY_RESULTS
	R_.random_excursion.valid=0;
#endif


	//move to
	/*
	#ifdef FILE_OUTPUT
	fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  RANDOM EXCURSIONS TEST\n");
	fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_RND_EXCURSION], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_RND_EXCURSION], "\t\t(a) Number Of Cycles (J) = %04d\n", J);
	fprintf(stats[TEST_RND_EXCURSION], "\t\t(b) Sequence Length (n)  = %d\n", n);
	#endif
	*/

	for (k = 0; k < 6; k++)
		for (i = 0; i < 9; i++)
			nu[k][i] = 0;

	for (bit_ind = 0; bit_ind < n; bit_ind++)
	{
		window = get_nth_block4(array, bit_ind);
		S_k += (window & 1) * 2 - 1;
		if (S_k == 0)
		{
			++J;
			for (i = 0; i < 9; i++)
			{
				if (counter[i] > 4)nu[5][i]++;
				else nu[counter[i]][i]++;
				counter[i] = 0;
			}
		}
		else
		{
			if ((S_k >= -4) && (S_k <= 4)) counter[S_k + 4]++;

		}
	}
	////Last cycle

	if (S_k)
	{
		++J;
		for (i = 0; i < 9; i++)
		{
			if (counter[i] > 4)nu[5][i]++;
			else nu[counter[i]][i]++;
			counter[i] = 0;
		}
	}

	////
	for (i = 5; i < 9; i++)
		for (k = 0; k < 6; k++)
		{
		nu[k][i - 1] = nu[k][i];
		}

	constraint = MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){


			fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
			for (i = 0; i < 8; i++)
				fprintf(results[TEST_RND_EXCURSION], "%f\n", 0.0);
		}
#endif
#ifdef KS
		for(i=0;i<8;i++)
			pvals.random_excursion_pvals[i][pvals.seq_counter] = 0.0;
#endif
	}
	else {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			//new
			fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  RANDOM EXCURSIONS TEST\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\tCOMPUTATIONAL INFORMATION:\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
			fprintf(stats[TEST_RND_EXCURSION], "\t\t(a) Number Of Cycles (J) = %04d\n", J);
			fprintf(stats[TEST_RND_EXCURSION], "\t\t(b) Sequence Length (n)  = %d\n", n);
			///
			fprintf(stats[TEST_RND_EXCURSION], "\t\t(c) Rejection Constraint = %f\n", constraint);
			fprintf(stats[TEST_RND_EXCURSION], "\t\t-------------------------------------------\n");
		}
#endif

		//PRINT
		/*
		printf("Cycle count = %d \n",J);
		for( i = 0; i < 8; i++)
		{
		for ( k=0; k<6; k++ )
		{
		printf("%5d",nu[k][i]);
		}
		printf("\n");
		}*/


		for (i = 0; i < 8; i++) {
			x = stateX[i];
			sum = 0.;
			for (k = 0; k < 6; k++)
				sum += pow(nu[k][i] - J*pi[(int)fabs(x)][k], 2) / (J*pi[(int)fabs(x)][k]);
			p_value = cephes_igamc(2.5, sum / 2.0);
#ifdef SPEED
			dummy_result += p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				if (isNegative(p_value) || isGreaterThanOne(p_value))
					fprintf(stats[TEST_RND_EXCURSION], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

				fprintf(stats[TEST_RND_EXCURSION], "%s\t\tx = %2d chi^2 = %9.6f p_value = %f\n",
					p_value < ALPHA ? "FAILURE" : "SUCCESS", x, sum, p_value);
				fprintf(results[TEST_RND_EXCURSION], "%f\n", p_value); fflush(results[TEST_RND_EXCURSION]);
			}
#endif

#ifdef KS
			pvals.random_excursion_pvals[i][pvals.seq_counter] = p_value;
#endif
#ifdef VERIFY_RESULTS
			R_.random_excursion.valid=1;
			R_.random_excursion.J[i]=J;
			R_.random_excursion.x[i]=x;
			R_.random_excursion.p_value[i]=p_value;
			R_.random_excursion.sum[i]=sum;
#endif
		}
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_RND_EXCURSION], "\n"); fflush(stats[TEST_RND_EXCURSION]);
}
#endif

#ifdef VERIFY_RESULTS
		if (RandomExcursions_v1 == RandomExcursions2) R1 = R_;
		else R2 = R_;
#endif

}
