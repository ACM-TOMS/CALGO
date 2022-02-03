#include "../include/config.h"
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
          N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
NonOverlappingTemplateMatchings(int m, int n)
{
	int		numOfTemplates[100] = { 0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
		2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152 };
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must
	first be constructed, saved into files and then the corresponding
	number of nonperiodic templates for that file be stored in the m-th
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	unsigned int	bit, W_obs, /* nu[6],*/ *Wj = NULL;
	FILE			*fp = NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int				i, j, jj, k, match, SKIP, M, N, K = 5;
	char			directory[100];
	BitSequence		*sequence = NULL;

	N = 8;
	M = n / N;

	if ((Wj = (unsigned int*)calloc(N, sizeof(unsigned int))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		}
#endif
		printf("\tNONOVERLAPPING TEMPLATES TESTS: Insufficient memory for required work space.\n");
		return;
	}
	lambda = (M - m + 1) / pow(2, m);
	varWj = M*(1.0 / pow(2.0, m) - (2.0*m - 1.0) / pow(2.0, 2.0*m));
	sprintf(directory, "templates/template%d", m);

	if (((isNegative(lambda)) || (isZero(lambda))) ||
		((fp = fopen(directory, "r")) == NULL) ||
		((sequence = (BitSequence *)calloc(m, sizeof(BitSequence))) == NULL)) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
			fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		}
#endif
		printf("\tNONOVERLAPPING TEMPLATES TESTS ABORTED.\n");
		if (sequence != NULL)
		{
			free(sequence);
			sequence = NULL;
		}

#ifdef VERIFY_RESULTS
		R_.nonoverlapping.templates=0;
		R_.nonoverlapping.W=NULL;
		R_.nonoverlapping.chi2=NULL;
		R_.nonoverlapping.p_value=NULL;
#endif
	}
	else {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
			fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		}
#endif

		if (numOfTemplates[m] < MAXNUMOFTEMPLATES)
			SKIP = 1;
		else
			SKIP = (int)(numOfTemplates[m] / MAXNUMOFTEMPLATES);
		numOfTemplates[m] = (int)numOfTemplates[m] / SKIP;

#ifdef KS
		pvals.num_Nonoverlap_pvals = (int)numOfTemplates[m] / SKIP;
#endif

#ifdef VERIFY_RESULTS
		R_.nonoverlapping.templates=MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]);
		R_.nonoverlapping.W=(unsigned int*)malloc(sizeof(unsigned int)*R_.nonoverlapping.templates*8);
		R_.nonoverlapping.chi2=(double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		R_.nonoverlapping.p_value=(double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		if(R_.nonoverlapping.W==NULL||R_.nonoverlapping.chi2==NULL||R_.nonoverlapping.p_value==NULL)
		{ printf("NONOVERLAPPING TEMPLATES TEST: Cannot allocate memory"); return; }
#endif


		sum = 0.0;
		for (i = 0; i < 2; i++) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda + i*log(lambda) - cephes_lgam(i + 1));
			sum += pi[i];
		}
		pi[0] = sum;
		for (i = 2; i <= K; i++) {                      /* Compute Probabilities */
			pi[i - 1] = exp(-lambda + i*log(lambda) - cephes_lgam(i + 1));
			sum += pi[i - 1];
		}
		pi[K] = 1 - sum;

		for (jj = 0; jj < MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]); jj++) {
			sum = 0;

			for (k = 0; k < m; k++) {
				fscanf(fp, "%d", &bit);
				sequence[k] = (unsigned char)bit;
#if defined(FILE_OUTPUT) ||  defined(KS)
				if (cmdFlags.output == 1 || cmdFlags.output == -1){
					fprintf(stats[TEST_NONPERIODIC], "%d", sequence[k]);
				}
#endif
			}
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				fprintf(stats[TEST_NONPERIODIC], " ");
			}
#endif
			// nadbytecny kod? nu nikde nevyuzito...
			//			for ( k=0; k<=K; k++ )
			//				nu[k] = 0;
			for (i = 0; i < N; i++) {
				W_obs = 0;
				for (j = 0; j < M - m + 1; j++) {
					match = 1;
					for (k = 0; k < m; k++) {
						if ((int)sequence[k] != (int)epsilon[i*M + j + k]) {
							match = 0;
							break;
						}
					}
					if (match == 1)
					{
						W_obs++;
						j += m - 1;
					}
				}
				Wj[i] = W_obs;
			}
			sum = 0;
			chi2 = 0.0;                                   /* Compute Chi Square */
			for (i = 0; i < N; i++) {
#if defined(FILE_OUTPUT) ||  defined(KS)
				if (cmdFlags.output == 1 || cmdFlags.output == -1){
					if (m == 10)
						fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i]);
					else
						fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i]);
				}
#endif
				/*if ( m == 10 )
					printf("%3d  ", Wj[i]);
					else
					printf("%4d ", Wj[i]);*/
#ifdef VERIFY_RESULTS
				R_.nonoverlapping.W[jj*N+i]=Wj[i];
#endif
				//
				chi2 += pow(((double)Wj[i] - lambda) / pow(varWj, 0.5), 2);
			}
			p_value = cephes_igamc(N / 2.0, chi2 / 2.0);
#ifdef SPEED
			dummy_result += p_value;
#endif
#ifdef VERIFY_RESULTS
			R_.nonoverlapping.chi2[jj]=chi2;
			R_.nonoverlapping.p_value[jj]=p_value;
#endif

#ifdef KS
			pvals.nonoverlapping_pvals[jj][pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				if (isNegative(p_value) || isGreaterThanOne(p_value))
					fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

				fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
			}
#endif
			if (SKIP > 1)
				fseek(fp, (long)(SKIP - 1) * 2 * m, SEEK_CUR);
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
			}
#endif
		}
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
}
#endif
	if ( sequence != NULL )
		free(sequence);

	free(Wj);
	if(fp != NULL)
	   fclose(fp);

#ifdef VERIFY_RESULTS
	if (NonOverlappingTemplateMatchings_v1 == NonOverlappingTemplateMatchings) R1 = R_;
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


//TODO save p-value into right file
void
NonOverlappingTemplateMatchings2(int m, int n)
{
	int		numOfTemplates[100] = { 0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
		2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152 };
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must
	first be constructed, saved into files and then the corresponding
	number of nonperiodic templates for that file be stored in the m-th
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	FILE			*fp = NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int	i, j, k, M, N, K = 5, SKIP;
	char			directory[100];

	unsigned int sequence, numoftemplates, window, one_template, **Wj = NULL, *templates, mask;
	int bit, bit_ind;
	int block;

	N = 8;
	M = n / N;

	if (numOfTemplates[m] < MAXNUMOFTEMPLATES)
		SKIP = 1;
	else
		SKIP = (int)(numOfTemplates[m] / MAXNUMOFTEMPLATES);


	numoftemplates = numOfTemplates[m];


	templates = (unsigned int*)malloc(numoftemplates*sizeof(unsigned int));
	if (templates == NULL)
	{
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
		}
#endif
		printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");
		return;
	}
	//Wj = (unsigned int**)malloc((1 << m)*sizeof(unsigned int*));
	Wj = (unsigned int**)malloc(N*sizeof(unsigned int*));
	if (Wj == NULL)
	{
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
		}
#endif
		printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");
		return;
	}
	for (i = 0; i < N; i++)
	{
		Wj[i] = (unsigned int*)malloc(((size_t)1 << m)*sizeof(unsigned int));
		if (Wj[i] == NULL)
		{
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
			}
#endif
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");

			for (j = 0; j < i; j++)
			{
				free(Wj[j]);
			}
			free(Wj);
			return;
		}
		for (j = 0; j < (1 << m); j++)
		{
			Wj[i][j] = 0;
		}
	}

	mask = (1 << m) - 1;
	lambda = (M - m + 1) / pow(2, m);
	varWj = M*(1.0 / pow(2.0, m) - (2.0*m - 1.0) / pow(2.0, 2.0*m));

	sprintf(directory, "templates/template%d", m);


	if (((isNegative(lambda)) || (isZero(lambda))) ||
		((fp = fopen(directory, "r")) == NULL)
		)
	{
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
			fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		}
#endif
		//if(fp==NULL)
		//	printf("NONOVERLAPPING TEMPLATES TESTS: Cannot open templates file.\n"); 

#ifdef VERIFY_RESULTS
		R_.nonoverlapping.templates=0;
		R_.nonoverlapping.W=(unsigned int*)malloc(sizeof(unsigned int)*R_.nonoverlapping.templates*8);
		R_.nonoverlapping.chi2=(double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		R_.nonoverlapping.p_value=(double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		if(R_.nonoverlapping.W==NULL||R_.nonoverlapping.chi2==NULL||R_.nonoverlapping.p_value==NULL)
		{ printf("NONOVERLAPPING TEMPLATES TEST: Cannot allocate memory"); return; }
#endif
	}
	else
	{
#ifdef VERIFY_RESULTS
		R_.nonoverlapping.templates=MIN(numoftemplates,MAXNUMOFTEMPLATES);
		R_.nonoverlapping.W=(unsigned int*)malloc(sizeof(unsigned int)*R_.nonoverlapping.templates*8);
		if(R_.nonoverlapping.W==NULL)
		{
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 
			return; 
		}
		R_.nonoverlapping.chi2=(double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		if(R_.nonoverlapping.chi2==NULL)
		{ 
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 
			return; 
		}
		R_.nonoverlapping.p_value=(double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		if(R_.nonoverlapping.p_value==NULL) 
		{ 
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 
			return; 
		}
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
			fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		}
#endif

		for (i = 0; i < (int)numoftemplates; i++)
		{
			sequence = 0;
			for (k = 0; k < m; k++) {
				fscanf(fp, "%d", &bit);
				sequence <<= 1;
				sequence ^= bit;
			}
			templates[i] = Mirrored_int(sequence, m); //sequence;
		}

		//timings();
		for (block = 0; block < N; block++)
		{
			for (bit_ind = block*M; bit_ind < (block + 1)*M - m + 1; bit_ind++)
			{
				window = get_nth_block4(array, bit_ind);
				//bits(&window,4);
				Wj[block][window & mask]++;
				//printf("%d",window & mask);
			}
		}

		sum = 0.0;
		for (i = 0; i < 2; i++) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda + i*log(lambda) - cephes_lgam(i + 1));
			sum += pi[i];
		}

		pi[0] = sum;
		for (i = 2; i <= K; i++) {                      /* Compute Probabilities */
			pi[i - 1] = exp(-lambda + i*log(lambda) - cephes_lgam(i + 1));
			sum += pi[i - 1];
		}
		pi[K] = 1 - sum;

		///new
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fclose(fp);
			fp = fopen(directory, "r");
		}
#endif
		///
		for (j = 0; j / SKIP < (int)MIN(numoftemplates, MAXNUMOFTEMPLATES); j += SKIP)
		{
			///new
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				for (k = 0; k < m; k++) {
					fscanf(fp, "%d", &bit);
					fprintf(stats[TEST_NONPERIODIC], "%d", bit);

				}
				fprintf(stats[TEST_NONPERIODIC], " ");
			}
#endif
			///
			one_template = templates[j];
			chi2 = 0.0;                                   /* Compute Chi Square */
			for (i = 0; i < N; i++) {
				chi2 += pow(((double)Wj[i][one_template] - lambda) / pow(varWj, 0.5), 2);

#if defined(FILE_OUTPUT) ||  defined(KS)
				if (cmdFlags.output == 1 || cmdFlags.output == -1){
					if (m == 10)
						fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i][one_template]);
					else
						fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i][one_template]);
				}
#endif

#ifdef VERIFY_RESULTS
				//if((j/SKIP)>=R_.nonoverlapping.templates) printf("!!!\n");
				R_.nonoverlapping.W[(j/SKIP)*N+i]=Wj[i][one_template];
#endif

			}
			p_value = cephes_igamc(N / 2.0, chi2 / 2.0);
			//printf("(2) chi2:%lf value: %lf\n",chi2,p_value);

#ifdef SPEED
			dummy_result += p_value;
#endif

#ifdef VERIFY_RESULTS
			R_.nonoverlapping.chi2[j / SKIP] = chi2;
			R_.nonoverlapping.p_value[j / SKIP] = p_value;
#endif

#ifdef KS
			pvals.nonoverlapping_pvals[j / SKIP][pvals.seq_counter] = p_value;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				if (isNegative(p_value) || isGreaterThanOne(p_value))
					fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

				fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", j / SKIP);
				fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
			}
#endif
		}
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
}
#endif

	for(i = 0; i < N; i++)
	{
		if(Wj[i])free(Wj[i]);
	}
	free(Wj);
	if(fp)fclose(fp); //PADA PRI M == 16

#ifdef VERIFY_RESULTS
	if (NonOverlappingTemplateMatchings_v1 == NonOverlappingTemplateMatchings2) R1 = R_;
	else R2 = R_;
#endif
}

//TODO save p-value into right file

void
NonOverlappingTemplateMatchings4(int m, int n)
{
	int		numOfTemplates[100] = { 0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
		2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152 };
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must
	first be constructed, saved into files and then the corresponding
	number of nonperiodic templates for that file be stored in the m-th
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	FILE			*fp = NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int	i, j, k, M, N, K = 5, SKIP;
	char			directory[100];

	unsigned int sequence, numoftemplates, one_template, *templates /*, mask */;
	int **Wj = NULL;
	int bit;
	int block;

	N = 8;
	M = n / N;

	if (numOfTemplates[m] < MAXNUMOFTEMPLATES)
		SKIP = 1;
	else
		SKIP = (int)(numOfTemplates[m] / MAXNUMOFTEMPLATES);


	numoftemplates = numOfTemplates[m];


	templates = (unsigned int*)malloc(numoftemplates*sizeof(unsigned int));
	if (templates == NULL)
	{
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
		}
#endif
		printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");
		return;
	}
	//Wj = (unsigned int**)malloc((1 << m)*sizeof(unsigned int*));
	Wj = (int**)malloc(N*sizeof(int*));
	if (Wj == NULL)
	{
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
		}
#endif
		printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");
		return;
	}
	for (i = 0; i < N; i++)
	{
		Wj[i] = (int*)malloc(((size_t)1 << m)*sizeof(int));
		if (Wj[i] == NULL)
		{
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
			}
#endif
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");

			for (j = 0; j < i; j++)
			{
				free(Wj[j]);
			}
			free(Wj);
			return;
		}
		for (j = 0; j < (1 << m); j++)
		{
			Wj[i][j] = 0;
		}
	}

	// mask = (1 << m) - 1;
	lambda = (M - m + 1) / pow(2, m);
	varWj = M*(1.0 / pow(2.0, m) - (2.0*m - 1.0) / pow(2.0, 2.0*m));

	sprintf(directory, "templates/template%d", m);


	if (((isNegative(lambda)) || (isZero(lambda))) ||
		((fp = fopen(directory, "r")) == NULL)
		)
	{
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
			fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
			fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
			fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
		}
#endif
		//if(fp==NULL)
		//	printf("NONOVERLAPPING TEMPLATES TESTS: Cannot open templates file.\n"); 

#ifdef VERIFY_RESULTS
		R_.nonoverlapping.templates = 0;
		R_.nonoverlapping.W = (unsigned int*)malloc(sizeof(unsigned int)*R_.nonoverlapping.templates * 8);
		R_.nonoverlapping.chi2 = (double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		R_.nonoverlapping.p_value = (double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		if (R_.nonoverlapping.W == NULL || R_.nonoverlapping.chi2 == NULL || R_.nonoverlapping.p_value == NULL)
		{
			printf("NONOVERLAPPING TEMPLATES TEST: Cannot allocate memory"); return;
		}
#endif
	}
	else
	{
#ifdef VERIFY_RESULTS
		R_.nonoverlapping.templates = MIN(numoftemplates, MAXNUMOFTEMPLATES);
		R_.nonoverlapping.W = (unsigned int*)malloc(sizeof(unsigned int)*R_.nonoverlapping.templates * 8);
		if (R_.nonoverlapping.W == NULL)
		{
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");
			return;
		}
		R_.nonoverlapping.chi2 = (double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		if (R_.nonoverlapping.chi2 == NULL)
		{
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");
			return;
		}
		R_.nonoverlapping.p_value = (double *)malloc(sizeof(double)*R_.nonoverlapping.templates);
		if (R_.nonoverlapping.p_value == NULL)
		{
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n");
			return;
		}
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
			fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
			fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
			fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		}
#endif

		for (i = 0; i < (int)numoftemplates; i++)
		{
			sequence = 0;
			for (k = 0; k < m; k++) {
				fscanf(fp, "%d", &bit);
				sequence <<= 1;
				sequence ^= bit;
			}
			templates[i] = Mirrored_int(sequence, m); //sequence;
		}

		//timings();
		for (block = 0; block < N; block++)
		{
			Histogram(block*M, Wj[block], m, (block + 1)*M);
		}

		sum = 0.0;
		for (i = 0; i < 2; i++) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda + i*log(lambda) - cephes_lgam(i + 1));
			sum += pi[i];
		}

		pi[0] = sum;
		for (i = 2; i <= K; i++) {                      /* Compute Probabilities */
			pi[i - 1] = exp(-lambda + i*log(lambda) - cephes_lgam(i + 1));
			sum += pi[i - 1];
		}
		pi[K] = 1 - sum;

		///new
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fclose(fp);
			fp = fopen(directory, "r");
		}
#endif
		///
		for (j = 0; j / SKIP < (int)MIN(numoftemplates, MAXNUMOFTEMPLATES); j += SKIP)
		{
			///new
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				for (k = 0; k < m; k++) {
					fscanf(fp, "%d", &bit);
					fprintf(stats[TEST_NONPERIODIC], "%d", bit);
				}
				fprintf(stats[TEST_NONPERIODIC], " ");
			}
#endif
			///
			one_template = templates[j];
			chi2 = 0.0;                                   /* Compute Chi Square */
			for (i = 0; i < N; i++) {
				chi2 += pow(((double)Wj[i][one_template] - lambda) / pow(varWj, 0.5), 2);

#if defined(FILE_OUTPUT) ||  defined(KS)
				if (cmdFlags.output == 1 || cmdFlags.output == -1){
					if (m == 10)
						fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i][one_template]);
					else
						fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i][one_template]);
				}
#endif

#ifdef VERIFY_RESULTS
				//if((j/SKIP)>=R_.nonoverlapping.templates) printf("!!!\n");
				R_.nonoverlapping.W[(j / SKIP)*N + i] = Wj[i][one_template];
#endif

			}
			p_value = cephes_igamc(N / 2.0, chi2 / 2.0);
			//printf("(2) chi2:%lf value: %lf\n",chi2,p_value);

#ifdef SPEED
			dummy_result += p_value;
#endif

#ifdef VERIFY_RESULTS
			R_.nonoverlapping.chi2[j/SKIP] = chi2;
			R_.nonoverlapping.p_value[j/SKIP] = p_value;
#endif

#ifdef KS
			pvals.nonoverlapping_pvals[j / SKIP][pvals.seq_counter] = p_value;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
			if (cmdFlags.output == 1 || cmdFlags.output == -1){
				if (isNegative(p_value) || isGreaterThanOne(p_value))
					fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

				fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", j / SKIP);
				fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
			}
#endif
		}
	}

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
}
#endif

	for (i = 0; i < N; i++)
	{
		if (Wj[i])free(Wj[i]);
	}
	free(Wj);
	free(templates);
	if (fp)fclose(fp); //PADA PRI M == 16

#ifdef VERIFY_RESULTS
	if (NonOverlappingTemplateMatchings_v1 == NonOverlappingTemplateMatchings4) R1 = R_;
	else R2 = R_;
#endif
}
