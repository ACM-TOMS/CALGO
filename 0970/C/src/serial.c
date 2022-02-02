#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/cephes.h"  
#include "../include/tools.h" 
#include "../include/stat_fncs.h"

double psi2(int m, int n);

void
Serial(int m, int n)
{
	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;
	
	psim0 = psi2(m, n);
	psim1 = psi2(m-1, n);
	psim2 = psi2(m-2, n);
	//printf("%lf %lf\n",psim1,psim2);
#ifdef SPEED
	dummy_result =  psim1 + psim2;
#endif
	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m-1)/2, del1/2.0);
	p_value2 = cephes_igamc(pow(2, m-2)/2, del2/2.0);
	//printf("%lf %lf\n",p_value1,p_value2);
#ifdef SPEED
	dummy_result = p_value1 + p_value2; // +psim0 + psim1 + psim2 + del1 + del2;
#endif

#ifdef VERIFY_RESULTS
	R_.serial.psim0=psim0;
	R_.serial.psim1=psim1;
	R_.serial.psim2=psim2;
	R_.serial.del1=del1;
	R_.serial.del2=del2;
	R_.serial.p_value1=p_value1;
	R_.serial.p_value2=p_value2;
	if(Serial_v1 == Serial) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
		fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
		fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
		fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
		fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
		fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
		fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
		fprintf(results[TEST_SERIAL], "%f\n", p_value1);

		fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2); fflush(stats[TEST_SERIAL]);
		fprintf(results[TEST_SERIAL], "%f\n", p_value2); fflush(results[TEST_SERIAL]);
	}
#endif
#ifdef KS
	pvals.serial_pvals[0][pvals.seq_counter] = p_value1;
	pvals.serial_pvals[1][pvals.seq_counter] = p_value2;
#endif
}

double
psi2(int m, int n)
{
	int				i, j, k, powLen;
	double			sum, numOfBlocks;
	unsigned int	*P;
	
	if ( (m == 0) || (m == -1) )
		return 0.0;
	numOfBlocks = n;
	powLen = (int)pow(2, m+1)-1;
	if ((P = (unsigned int*)calloc(powLen, sizeof(unsigned int))) == NULL) {
		printf("Serial Test:  Insufficient memory available.\n");
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_SERIAL], "Serial Test:  Insufficient memory available.\n");
			fflush(stats[TEST_SERIAL]);
	}
#endif
		return 0.0;
	}
	for ( i=1; i<powLen-1; i++ )
		P[i] = 0;	  /* INITIALIZE NODES */
	for ( i=0; i<numOfBlocks; i++ ) {		 /* COMPUTE FREQUENCY */
		k = 1;
		for ( j=0; j<m; j++ ) {
			if ( epsilon[(i+j)%n] == 0 )
				k *= 2;
			else if ( epsilon[(i+j)%n] == 1 )
				k = 2*k+1;
		}
		P[k-1]++;
	}
	sum = 0.0;
	for ( i=(int)pow(2, m)-1; i<(int)pow(2, m+1)-1; i++ )
	{
		sum += pow(P[i], 2);
		//printf("%d %d ",i,P[i]);
	}

	sum = (sum * pow(2, m)/(double)n) - (double)n;
	free(P);
	
	//printf("%lf\n",sum);
	return sum;
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
Serial2(int m, int n)
{
	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;

	int				i, k, len;
	double			sum /*, numOfBlocks */;
	unsigned int	*P, mask, help;
	//deleted
	/*#ifdef FILE_OUTPUT
	fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);
	#endif
	*/

	//numOfBlocks = n;
	mask = (1 << m) - 1;
	len = (1 << m);

	if ((P = (unsigned int*)calloc(len, sizeof(unsigned int))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
		}
#endif
		printf("Serial Test:  Insufficient memory available.\n");
		return;
	}
	for (i = 0; i < len; i++)
		P[i] = 0;
	/* COMPUTE FREQUENCY */
	for (i = 0; i < n - m + 1; i++) {
		help = get_nth_block4(array, i);
		++P[/*get_nth_block4(array,i)*/ help & mask];
		//if(i < 100)printf(" %i ",(1 << m)- 1 + Mirrored_int((get_nth_block4(array,i)&mask),m));
	}

	for (i = 1; i < m; i++) {
		k = get_nth_block4(array, n - m + i)&(mask >> i);
		//printf("%d ",k);
		k ^= (((unsigned int*)array)[0] << (m - i));
		//printf("%d ",k);
		k &= mask;
		//printf("%d ",k);
		P[k]++;
	}

	if (m > 0)
	{
		sum = 0.0;
		for (i = 0; i < len; i++) {
			help = P[Mirrored_int(i, m)];
			if (help > 0)
				sum += (double)help*help;
			//printf("%d %d ",(1 << m)-1+i,help);	
		}
		sum = (sum * (1 << m) / (double)n) - (double)n;
		//printf("SUM = %lf   \n",sum);
	}
	else sum = 0.0;
	psim0 = sum;


	if (m - 1 > 0)
	{
		sum = 0.0;
		for (i = 0; i < len / 2; i++) {
			help = P[Mirrored_int(i, m - 1)] + P[Mirrored_int(i, m - 1) + len / 2];
			if (help > 0)
				sum += (double)help*help;
			//printf("%d %d ",(1 << (m-1))-1+Mirrored_int(i,m-1),help);	
		}
		sum = (sum * (1 << (m - 1)) / (double)n) - (double)n;
		//printf("SUM = %lf   \n",sum);

	}
	else sum = 0.0;
	psim1 = sum;

	if (m - 2 > 0)
	{
		sum = 0.0;
		for (i = 0; i < len / 4; i++) {
			help = P[Mirrored_int(i, m - 2)] + P[Mirrored_int(i, m - 2) + len / 4] + P[Mirrored_int(i, m - 2) + 2 * len / 4] + P[Mirrored_int(i, m - 2) + 3 * len / 4];
			if (help > 0)
				sum += (double)help*help;
			//printf("%d %d ",(1 << (m-2))-1+Mirrored_int(i,m-2),help);	
		}
		sum = (sum * (1 << (m - 2)) / (double)n) - (double)n;
		//printf("SUM = %lf   \n",sum);

	}
	else sum = 0.0;
	psim2 = sum;

	//printf("%lf %lf %lf\n",psim0,psim1,psim2);
	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m - 1) / 2, del1 / 2.0);
	p_value2 = cephes_igamc(pow(2, m - 2) / 2, del2 / 2.0);

#ifdef SPEED
	dummy_result = p_value1;
	dummy_result += p_value2;
#endif

#ifdef VERIFY_RESULTS
	R_.serial.psim0 = psim0;
	R_.serial.psim1 = psim1;
	R_.serial.psim2 = psim2;
	R_.serial.del1 = del1;
	R_.serial.del2 = del2;
	R_.serial.p_value1 = p_value1;
	R_.serial.p_value2 = p_value2;
	if (Serial_v1 == Serial2) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
		fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
		fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
		fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
		fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
		fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
		fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
		fprintf(results[TEST_SERIAL], "%f\n", p_value1);

		fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2); fflush(stats[TEST_SERIAL]);
		fprintf(results[TEST_SERIAL], "%f\n", p_value2); fflush(results[TEST_SERIAL]);
}
#endif
#ifdef KS
	pvals.serial_pvals[0][pvals.seq_counter] = p_value1;
	pvals.serial_pvals[1][pvals.seq_counter] = p_value2;
#endif
	free(P);
}

void
Serial4(int m, int n)
{
	double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;

	int				i, k, len;
	double			sum /*, numOfBlocks */;
	unsigned int	mask, help;
	int             *P;
	//deleted
	/*#ifdef FILE_OUTPUT
	fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);
	#endif
	*/

	//numOfBlocks = n;
	mask = (1 << m) - 1;
	len = (1 << m);

	if ((P = (int*)calloc(len, sizeof(int))) == NULL) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
	}
#endif
		printf("Serial Test:  Insufficient memory available.\n");
		return;
	}
	for (i = 0; i < len; i++)
		P[i] = 0;
	/* COMPUTE FREQUENCY */
	Histogram(0, P, m, n);

	for (i = 1; i<m; i++) {
		k = get_nth_block4(array, n - m + i)&(mask >> i);
		//printf("%d ",k);
		k ^= (((unsigned int*)array)[0] << (m - i));
		//printf("%d ",k);
		k &= mask;
		//printf("%d ",k);
		P[k]++;
	}

	if (m > 0)
	{
		sum = 0.0;
		for (i = 0; i < len; i++) {
			help = P[Mirrored_int(i, m)];
			if (help > 0)
				sum += (double)help*help;
			//printf("%d %d ",(1 << m)-1+i,help);	
		}
		sum = (sum * (1 << m) / (double)n) - (double)n;
		//printf("SUM = %lf   \n",sum);
	}
	else sum = 0.0;
	psim0 = sum;


	if (m - 1 > 0)
	{
		sum = 0.0;
		for (i = 0; i < len / 2; i++) {
			help = P[Mirrored_int(i, m - 1)] + P[Mirrored_int(i, m - 1) + len / 2];
			if (help > 0)
				sum += (double)help*help;
			//printf("%d %d ",(1 << (m-1))-1+Mirrored_int(i,m-1),help);	
		}
		sum = (sum * (1 << (m - 1)) / (double)n) - (double)n;
		//printf("SUM = %lf   \n",sum);

	}
	else sum = 0.0;
	psim1 = sum;

	if (m - 2 > 0)
	{
		sum = 0.0;
		for (i = 0; i < len / 4; i++) {
			help = P[Mirrored_int(i, m - 2)] + P[Mirrored_int(i, m - 2) + len / 4] + P[Mirrored_int(i, m - 2) + 2 * len / 4] + P[Mirrored_int(i, m - 2) + 3 * len / 4];
			if (help > 0)
				sum += (double)help*help;
			//printf("%d %d ",(1 << (m-2))-1+Mirrored_int(i,m-2),help);	
		}
		sum = (sum * (1 << (m - 2)) / (double)n) - (double)n;
		//printf("SUM = %lf   \n",sum);

	}
	else sum = 0.0;
	psim2 = sum;

	//printf("%lf %lf %lf\n",psim0,psim1,psim2);
	del1 = psim0 - psim1;
	del2 = psim0 - 2.0*psim1 + psim2;
	p_value1 = cephes_igamc(pow(2, m - 1) / 2, del1 / 2.0);
	p_value2 = cephes_igamc(pow(2, m - 2) / 2, del2 / 2.0);

#ifdef SPEED
	dummy_result = p_value1;
	dummy_result += p_value2;
#endif

#ifdef VERIFY_RESULTS
	R_.serial.psim0 = psim0;
	R_.serial.psim1 = psim1;
	R_.serial.psim2 = psim2;
	R_.serial.del1 = del1;
	R_.serial.del2 = del2;
	R_.serial.p_value1 = p_value1;
	R_.serial.p_value2 = p_value2;
	if (Serial_v1 == Serial4) R1 = R_;
	else R2 = R_;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
		fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
		fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
		fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
		fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
		fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
		fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
		fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

		fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
		fprintf(results[TEST_SERIAL], "%f\n", p_value1);

		fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2); fflush(stats[TEST_SERIAL]);
		fprintf(results[TEST_SERIAL], "%f\n", p_value2); fflush(results[TEST_SERIAL]);
	}
#endif
#ifdef KS
	pvals.serial_pvals[0][pvals.seq_counter] = p_value1;
	pvals.serial_pvals[1][pvals.seq_counter] = p_value2;
#endif
	free(P);
}
