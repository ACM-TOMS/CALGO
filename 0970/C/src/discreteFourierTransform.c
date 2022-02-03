#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"
#include "../include/erf.h"
#include "../include/tools.h"
#include "../include/stat_fncs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         D I S C R E T E  F O U R I E R  T R A N S F O R M  T E S T 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void  __ogg_fdrffti(int n, double *wsave, int *ifac);
void  __ogg_fdrfftf(int n, double *X, double *wsave, int *ifac);

void
DiscreteFourierTransform(int n)
{
	double  p_value, upperBound, percentile, N_l, N_o, d, *m = NULL, *X = NULL, *wsave = NULL;
	int		i, count, ifac[15];

	if ( ((X = (double*) calloc(n,sizeof(double))) == NULL) ||
		 ((wsave = (double *)calloc(2*n,sizeof(double))) == NULL) ||
		 ((m = (double*)calloc(n/2+1, sizeof(double))) == NULL) ) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[7], "\t\tUnable to allocate working arrays for the DFT.\n");
		}
#endif
			 printf("\t\tUnable to allocate working arrays for the DFT.\n");
			if( X != NULL )
				free(X);
			if( wsave != NULL )
				free(wsave);
			if( m != NULL )
				free(m);
			return;
	}
	for ( i=0; i<n; i++ ){
		X[i] = 2*(int)epsilon[i] - 1;
	}
	
	__ogg_fdrffti(n, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
	__ogg_fdrfftf(n, X, wsave, ifac);	/* APPLY FORWARD FFT */
	
	m[0] = sqrt(X[0]*X[0]);	    /* COMPUTE MAGNITUDE */
	
	for ( i=0; i<n/2; i++ )
		m[i+1] = sqrt(pow(X[2*i+1],2)+pow(X[2*i+2],2)); 
	count = 0;				       /* CONFIDENCE INTERVAL */
	upperBound = sqrt(2.995732274*n);
	for ( i=0; i<n/2; i++ )
		if ( m[i] < upperBound )
			count++;
	percentile = (double)count/(n/2)*100;
	N_l = (double) count;       /* number of peaks less than h = sqrt(3*n) */
	N_o = (double) 0.95*n/2.0;
	d = (N_l - N_o)/sqrt(n/4.0*0.95*0.05);
	p_value = erfc(fabs(d)/sqrt(2.0));
	//printf("%lf ",p_value);
#ifdef SPEED
	dummy_result = p_value + percentile;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
		fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
		//fprintf(stats[TEST_FFT], "\t\t(c) N_o        = r%f\n", N_o); replaced by
		fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o);
		fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

		fprintf(stats[TEST_FFT], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
		fprintf(results[TEST_FFT], "%f\n", p_value);
	}
#endif

#ifdef VERIFY_RESULTS
	R_.dft.p_value=p_value;
	R_.dft.percentile=percentile;
	R_.dft.N_l=N_l;
	R_.dft.N_o=N_o;
	R_.dft.d=d;
	if(DiscreteFourierTransform_v1 == DiscreteFourierTransform) R1 = R_;
	else R2 = R_;
#endif

#ifdef KS
	pvals.dft_pvals[pvals.seq_counter] = p_value;
#endif
	free(X);
	free(wsave);
	free(m);
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
DiscreteFourierTransform2(int n)
{
	double	p_value, upperBound, percentile, N_l, N_o, d, *m=NULL, *X=NULL, *wsave=NULL;
	int		i, count, ifac[15];

	if (((X = (double*)calloc(n, sizeof(double))) == NULL) ||
		((wsave = (double *)calloc(2 * n, sizeof(double))) == NULL) ||
		((m = (double*)calloc(n / 2 + 1, sizeof(double))) == NULL)) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1){
			fprintf(stats[7], "\t\tUnable to allocate working arrays for the DFT.\n");
	}
#endif
			printf("\t\tUnable to allocate working arrays for the DFT.\n");
			if( X != NULL )
				free(X);
			if( wsave != NULL )
				free(wsave);
			if( m != NULL )
				free(m);
			return;
	}
	for ( i=0; i<n; i++ )
	{
		X[i] = 2*((double)(get_nth_block4(array,i) & 1))-1;
	}
	
	__ogg_fdrffti(n, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
	__ogg_fdrfftf(n, X, wsave, ifac);	/* APPLY FORWARD FFT */
	
	m[0] = sqrt(X[0]*X[0]);	    /* COMPUTE MAGNITUDE */
	
	for ( i=0; i<n/2; i++ )
		m[i+1] = sqrt(pow(X[2*i+1],2)+pow(X[2*i+2],2)); 
	count = 0;				       /* CONFIDENCE INTERVAL */
	upperBound = sqrt(2.995732274*n);
	for ( i=0; i<n/2; i++ )
		if ( m[i] < upperBound )
			count++;
	percentile = (double)count/(n/2)*100;
	N_l = (double) count;       /* number of peaks less than h = sqrt(3*n) */
	N_o = (double) 0.95*n/2.0;
	d = (N_l - N_o)/sqrt(n/4.0*0.95*0.05);
	p_value = erfc(fabs(d)/sqrt(2.0));
	//printf("%lf\n",p_value);
#ifdef SPEED
	dummy_result = p_value + percentile;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1){
		fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
		fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
		fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o);
		fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

		fprintf(stats[TEST_FFT], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
		fprintf(results[TEST_FFT], "%f\n", p_value);
	}
#endif

#ifdef VERIFY_RESULTS
	R_.dft.p_value=p_value;
	R_.dft.percentile=percentile;
	R_.dft.N_l=N_l;
	R_.dft.N_o=N_o;
	R_.dft.d=d;
	if(DiscreteFourierTransform_v1 == DiscreteFourierTransform2) R1 = R_;
	else R2 = R_;
#endif

#ifdef KS
	pvals.dft_pvals[pvals.seq_counter] = p_value;
#endif
	free(X);
	free(wsave);
	free(m);
}

/*

#include "../include/fftw3.h"

void
DiscreteFourierTransform3(int n)
{
	double	p_value, upperBound, percentile, N_l, N_o, d, *m = NULL;
	fftw_complex *in, *out;
	fftw_plan p;
	int		i, count;

	if (((in = (fftw_complex*)fftw_malloc(n * sizeof(fftw_complex))) == NULL) ||
		((out = (fftw_complex*)fftw_malloc(n * sizeof(fftw_complex))) == NULL) ||
		((m = (double*)calloc(n / 2 + 1, sizeof(double))) == NULL)) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[7], "\t\tUnable to allocate working arrays for the DFT.\n");
		}
#endif
		printf("\t\tUnable to allocate working arrays for the DFT.\n");
		if (in != NULL)
			fftw_free(in);
		if (out != NULL)
			fftw_free(out);
		if (m != NULL)
			free(m);
		return;
	}

	p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	for (i = 0; i<n; i++)
	{
		in[i][0] = 2 * ((double)(get_nth_block4(array, i) & 1)) - 1;
		in[i][1] = 0.0;
	}

	fftw_execute(p);

	m[0] = sqrt(out[0][0] * out[0][0]);	  

	for (i = 0; i<n / 2; i++)
		m[i + 1] = sqrt(pow(out[i + 1][0], 2) + pow(out[i + 1][1], 2));

	fftw_destroy_plan(p);

	count = 0;				      
	upperBound = sqrt(2.995732274*n);
	for (i = 0; i<n / 2; i++)
	if (m[i] < upperBound)
		count++;
	percentile = (double)count / (n / 2) * 100;
	N_l = (double)count;      
	N_o = (double) 0.95*n / 2.0;
	d = (N_l - N_o) / sqrt(n / 4.0*0.95*0.05);
	p_value = erfc(fabs(d) / sqrt(2.0));
#ifdef SPEED
	dummy_result = p_value + percentile;
#endif
#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
		fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
		fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o);
		fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

		fprintf(stats[TEST_FFT], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
		fprintf(results[TEST_FFT], "%f\n", p_value);
	}
#endif

#ifdef VERIFY_RESULTS
	R_.dft.p_value=p_value;
	R_.dft.percentile=percentile;
	R_.dft.N_l=N_l;
	R_.dft.N_o=N_o;
	R_.dft.d=d;
	if(DiscreteFourierTransform_v1 == DiscreteFourierTransform3) R1 = R_;
	else R2 = R_;
#endif
#ifdef KS
	pvals.dft_pvals[pvals.seq_counter] = p_value;
#endif
	fftw_free(in);
	fftw_free(out);
	free(m);
}


void
DiscreteFourierTransform4(int n)
{
	double	p_value, upperBound, percentile, N_l, N_o, d, *m = NULL;
	fftw_complex *out;
	double *in;
	fftw_plan p;
	int		i, count;

	if (((in = (double*)fftw_malloc(n * sizeof(double))) == NULL) ||
		((out = (fftw_complex*)fftw_malloc((n / 2 + 1) * sizeof(fftw_complex))) == NULL) ||
		((m = (double*)calloc(n / 2 + 1, sizeof(double))) == NULL)) {
#if defined(FILE_OUTPUT) ||  defined(KS)
		if (cmdFlags.output == 1 || cmdFlags.output == -1) {
			fprintf(stats[7], "\t\tUnable to allocate working arrays for the DFT.\n");
		}
#endif
		printf("\t\tUnable to allocate working arrays for the DFT.\n");
		if (in != NULL)
			fftw_free(in);
		if (out != NULL)
			fftw_free(out);
		if (m != NULL)
			free(m);
		return;
	}

	p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE );

	for (i = 0; i<n; i++)
	{
		in[i] = 2 * ((double)(get_nth_block4(array, i) & 1)) - 1;
	}

	fftw_execute(p);

	m[0] = sqrt(out[0][0] * out[0][0]);

	for (i = 0; i<n / 2; i++)
		m[i + 1] = sqrt(pow(out[i + 1][0], 2) + pow(out[i + 1][1], 2));

	fftw_destroy_plan(p);

	count = 0;				     
	upperBound = sqrt(2.995732274*n);
	for (i = 0; i<n / 2; i++)
	if (m[i] < upperBound)
		count++;
	percentile = (double)count / (n / 2) * 100;
	N_l = (double)count;      
	N_o = (double) 0.95*n / 2.0;
	d = (N_l - N_o) / sqrt(n / 4.0*0.95*0.05);
	p_value = erfc(fabs(d) / sqrt(2.0));

#ifdef SPEED
	dummy_result = p_value + percentile;
#endif

#if defined(FILE_OUTPUT) ||  defined(KS)
	if (cmdFlags.output == 1 || cmdFlags.output == -1) {
		fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
		fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
		fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
		fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o);
		fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
		fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

		fprintf(stats[TEST_FFT], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
		fprintf(results[TEST_FFT], "%f\n", p_value);
	}
#endif

#ifdef VERIFY_RESULTS
	R_.dft.p_value=p_value;
	R_.dft.percentile=percentile;
	R_.dft.N_l=N_l;
	R_.dft.N_o=N_o;
	R_.dft.d=d;
	if(DiscreteFourierTransform_v1 == DiscreteFourierTransform4) R1 = R_;
	else R2 = R_;
#endif
#ifdef KS
	pvals.dft_pvals[pvals.seq_counter] = p_value;
#endif
	fftw_free(in);
	fftw_free(out);
	free(m);
}
*/