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

#include "../include/config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "../include/externs.h"
#include "../include/cephes.h"
#include "../include/utilities.h"
#include "../include/tools.h"
#include "../include/stat_fncs.h"

#ifdef _MSC_VER
#if _MSC_VER <= 1700
#define M_LOG2E 1.44269504088896340736	//log2(e)
long double
log2(const long double x)
{
	return log(x) * M_LOG2E;
}
#endif
#endif

unsigned char *
load_array(FILE * f, int bit_size)
{
	unsigned char *array;
	unsigned int byte_size;
	int n;

	if (bit_size < 1)
	{
		fseek(f, 0, SEEK_END);
		bit_size = ftell(f) * 8;
		fseek(f, 0, SEEK_SET);
	}
	n = bit_size;
	byte_size = (unsigned int)ceil(n / 8.0) + 4;
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size));
	if (array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	fread(array, 1, byte_size, f);
	//after_byte = &array[byte_size]+1;
	array[byte_size - 1] = array[byte_size - 2] = array[byte_size - 3] =
		array[byte_size - 4] = 0;
	return array;
}

void
transform(BitSequence * epsilon, unsigned char *array, int n)
{
	int i;

	for (i = 0; i < n; i++)
	{
		//printf("%d %d %d\n",array[i >> 3],(1 << (i & 7)), array[i >> 3] & (1 << (i & 7)));
		if ((array[i >> 3] & (1 << (i & 7))) != 0)
			epsilon[i] = 1;
		else
			epsilon[i] = 0;
	}
}

void
data_all_zeros(int bit_size)
{
	int i;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = 0;
	transform(epsilon, array, bit_size);
}

void
data_all_ones_padded(int bit_size)
{
	int i;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = 0xff;
	transform(epsilon, array, bit_size);
}

void
data_all_ones(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = 0xff;

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	array[i - 1] ^= padding;
	transform(epsilon, array, bit_size);
}

void
data_all_zeros_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = 0;

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	array[i - 1] |= padding;
	transform(epsilon, array, bit_size);
}

void
data_pattern_01_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = 0xAA;

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	array[i - 1] |= padding;
	transform(epsilon, array, bit_size);
}

void
data_pattern_01(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = 0xAA;

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	padding = 0xFF ^ padding;
	array[i - 1] &= padding;
	transform(epsilon, array, bit_size);
}

void
data_pattern_10_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = 0x55;

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	array[i - 1] |= padding;
	transform(epsilon, array, bit_size);
}

void
data_pattern_10(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = 0x55;

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	padding = 0xFF ^ padding;
	array[i - 1] &= padding;
	transform(epsilon, array, bit_size);
}

void
data_prandom_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	srand((unsigned)0 /*time(NULL) */);
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = (0xFF & ((unsigned)rand()));

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	array[i - 1] |= padding;
	transform(epsilon, array, bit_size);
}

void
data_prandom(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	srand((unsigned)0 /*time(NULL) */);
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = (0xFF & ((unsigned)rand()));

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	padding = 0xFF ^ padding;
	array[i - 1] &= padding;
	transform(epsilon, array, bit_size);
}

void
data_prandom_padded_random(int bit_size)
{
	int i /*, pad */;
	// unsigned char mask = 0x80;
	//unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	// pad = byte_size * 8 - bit_size;
	srand((unsigned)time(NULL));
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = (0xFF & ((unsigned)rand()));

	transform(epsilon, array, bit_size);
}

void
data_bad_prandom_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	srand((unsigned)time(NULL));
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = (0xFF & ((unsigned)rand())) | (0xFF & ((unsigned)rand()));

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	array[i - 1] |= padding;
	transform(epsilon, array, bit_size);
}

void
data_bad_prandom(int bit_size)
{
	int i, j, pad;
	unsigned char mask = 0x80;
	unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	pad = byte_size * 8 - bit_size;
	srand((unsigned)time(NULL));
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = (0xFF & ((unsigned)rand())) | (0xFF & ((unsigned)rand()));

	for (j = 0; j < pad; j++)
	{
		padding += mask;
		mask >>= 1;
	}

	padding = 0xFF ^ padding;
	array[i - 1] &= padding;
	transform(epsilon, array, bit_size);
}

void
data_bad_prandom_padded_random(int bit_size)
{
	int i /* , pad */ ;
	//unsigned char mask = 0x80;
	//unsigned char padding = 0x00;
	int byte_size = bit_size / 8;
	if (byte_size * 8 < bit_size)
		byte_size++;
	// pad = byte_size * 8 - bit_size;
	srand((unsigned)time(NULL));
	epsilon =
		(unsigned char *)malloc(sizeof(unsigned char) * (byte_size * 8));
	array = (unsigned char *)malloc(sizeof(unsigned char) * (byte_size)+4);
	if (epsilon == NULL || array == NULL)
	{
		printf("Cannot allocate memory.\n");
		exit(10);
	}
	for (i = 0; i < byte_size; i++)
		array[i] = (0xFF & ((unsigned)rand())) | (0xFF & ((unsigned)rand()));

	transform(epsilon, array, bit_size);
}

#include <float.h>

union Double_Int
{
	int64_t i;
	double d;
};

int
ReasonablyEqualDoubles(double D1, double D2, int Max)
{
	union Double_Int DD1;
	union Double_Int DD2;
	int64_t Difference;

	if (D1 == D2)
		return 1;

	assert(sizeof(double) == 8);

	DD1.d = D1;
	DD2.d = D2;
	// Compare the sign
	if (((DD1.i >> 63) != 0) != ((DD2.i >> 63) != 0))
	{
		//if(D1 == D2) return 1;
		return 0;
	}
	// Compare the integers made from doubles
	Difference = DD1.i - DD2.i;
	if (Difference < 0)
		Difference = -Difference;
	if (Difference <= Max)
		return 1;
	return 0;
}

#ifdef VERIFY_RESULTS
int
equal_frequency_results()
{
	if (R1.frequency.sum != R2.frequency.sum)
		return 0;
	if (R1.frequency.sum_n != R2.frequency.sum_n)
		return 0;
	if (R1.frequency.p_value != R2.frequency.p_value)
		return 0;
	return 1;
}

int
equal_blockfrequency_results()
{
	if (R1.blockfrequency.chi_squared != R2.blockfrequency.chi_squared)
		return 0;
	if (R1.blockfrequency.p_value != R2.blockfrequency.p_value)
		return 0;
	return 1;
}



int
equal_runs_results()
{
	if (!ReasonablyEqualDoubles(R1.runs.p_value, R2.runs.p_value, 1))
		return 0;
	if (R1.runs.V != R2.runs.V)
		return 0;
	if (R1.runs.pi != R2.runs.pi)
		return 0;
	if (!ReasonablyEqualDoubles(R1.runs.erfc_arg, R2.runs.erfc_arg, 1))
		return 0;
	return 1;
}

int
equal_longestrunofones_results()
{
	int i;
	if (R1.longestrunofones.p_value != R2.longestrunofones.p_value)
		return 0;
	if (R1.longestrunofones.chi2 != R2.longestrunofones.chi2)
		return 0;
	if (R1.longestrunofones.N != R2.longestrunofones.N)
		return 0;
	if (R1.longestrunofones.M != R2.longestrunofones.M)
		return 0;
	for (i = 0; i < 7; i++)
		if (R1.longestrunofones.nu[i] != R2.longestrunofones.nu[i])
			return 0;
	return 1;
}

int
equal_rank_results()
{
	if (R1.rank.p_30 != R2.rank.p_30)
		return 0;
	if (R1.rank.p_31 != R2.rank.p_31)
		return 0;
	if (R1.rank.p_32 != R2.rank.p_32)
		return 0;
	if (R1.rank.F_30 != R2.rank.F_30)
		return 0;
	if (R1.rank.F_31 != R2.rank.F_31)
		return 0;
	if (R1.rank.F_32 != R2.rank.F_32)
		return 0;
	if (R1.rank.chi_squared != R2.rank.chi_squared)
		return 0;
	if (R1.rank.p_value != R2.rank.p_value)
		return 0;
	if (R1.rank.N != R2.rank.N)
		return 0;
	return 1;
}

int
equal_serial_results()
{
	if (R1.serial.psim0 != R2.serial.psim0)
		return 0;
	if (R1.serial.psim1 != R2.serial.psim1)
		return 0;
	if (R1.serial.psim2 != R2.serial.psim2)
		return 0;
	if (R1.serial.del1 != R2.serial.del1)
		return 0;
	if (R1.serial.del2 != R2.serial.del2)
		return 0;
	if (R1.serial.p_value1 != R2.serial.p_value1)
		return 0;
	if (R1.serial.p_value2 != R2.serial.p_value2)
		return 0;
	return 1;
}

int
equal_nonoverlapping_results()
{
	unsigned int i, j;
	if (R1.nonoverlapping.templates != R2.nonoverlapping.templates)
		return 0;
	for (i = 0; i < R1.nonoverlapping.templates; i++)
	{
		if (R1.nonoverlapping.chi2[i] != R2.nonoverlapping.chi2[i])
			return 0;
		if (R1.nonoverlapping.p_value[i] != R2.nonoverlapping.p_value[i])
			return 0;
		for (j = 0; j < 8; j++)
			if (R1.nonoverlapping.W[i * 8 + j] != R2.nonoverlapping.W[i * 8 + j])
				return 0;
	}
	if (R1.nonoverlapping.chi2)
		free(R1.nonoverlapping.chi2);
	if (R1.nonoverlapping.p_value)
		free(R1.nonoverlapping.p_value);
	if (R1.nonoverlapping.W)
		free(R1.nonoverlapping.W);

	if (R2.nonoverlapping.chi2)
		free(R2.nonoverlapping.chi2);
	if (R2.nonoverlapping.p_value)
		free(R2.nonoverlapping.p_value);
	if (R2.nonoverlapping.W)
		free(R2.nonoverlapping.W);
	return 1;
}

int
equal_overlapping_results()
{
	int i;
	if (R1.overlapping.chi2 == R1.overlapping.chi2)
	{
		if (R1.overlapping.chi2 != R2.overlapping.chi2)
			return 0;
	}
	else
	{
		if (R2.overlapping.chi2 == R2.overlapping.chi2)
			return 0;
	}
	if (R1.overlapping.p_value == R1.overlapping.p_value)
	{
		if (R1.overlapping.p_value != R2.overlapping.p_value)
			return 0;
	}
	else
	{
		if (R2.overlapping.p_value == R2.overlapping.p_value)
			return 0;
	}
	for (i = 0; i < 6; i++)
	{
		if (R1.overlapping.nu[i] != R2.overlapping.nu[i])
			return 0;
	}
	return 1;
}

int
equal_universal_results()
{
	if (R1.universal.phi != R2.universal.phi)
		return 0;
	if (R1.universal.p_value != R2.universal.p_value)
		return 0;
	if (R1.universal.sum != R2.universal.sum)
		return 0;
	return 1;
}

int
equal_apen_results()
{
	unsigned int i;
	if (R1.approximate_entropy.pp != R2.approximate_entropy.pp)
		return 0;
	for (i = 0; i < R1.approximate_entropy.pp; i++)
		if (R1.approximate_entropy.P[i] != R2.approximate_entropy.P[i])
			return 0;
	if (R1.approximate_entropy.P)
		free(R1.approximate_entropy.P);
	if (R2.approximate_entropy.P)
		free(R2.approximate_entropy.P);

	/*
	   if(!ReasonablyEqualDoubles(R1.approximate_entropy.chi_squared, R2.approximate_entropy.chi_squared,10000000)) return 0;
	   if(!ReasonablyEqualDoubles(R1.approximate_entropy.p_value, R2.approximate_entropy.p_value, 10000000)) return 0;
	   if(!ReasonablyEqualDoubles(R1.approximate_entropy.ApEn[0], R2.approximate_entropy.ApEn[0], 10000000)) return 0;
	   if(!ReasonablyEqualDoubles(R1.approximate_entropy.ApEn[1], R2.approximate_entropy.ApEn[1], 10000000)) return 0;
	   */

	/*
	if (fabs
		(R1.approximate_entropy.chi_squared -
			R2.approximate_entropy.chi_squared) > FLT_EPSILON)
		return 0;
	if (fabs(R1.approximate_entropy.p_value - R2.approximate_entropy.p_value) >
		FLT_EPSILON)
		return 0;
	if (fabs(R1.approximate_entropy.ApEn[0] - R2.approximate_entropy.ApEn[0]) >
		FLT_EPSILON)
		return 0;
	if (fabs(R1.approximate_entropy.ApEn[1] - R2.approximate_entropy.ApEn[1]) >
		FLT_EPSILON)
		return 0;
	*/

	   if(R1.approximate_entropy.chi_squared!=R2.approximate_entropy.chi_squared)
	   return 0;
	   if(R1.approximate_entropy.p_value!=R2.approximate_entropy.p_value)
	   return 0;
	   if(R1.approximate_entropy.ApEn[0]!=R2.approximate_entropy.ApEn[0])
	   return 0;
	   if(R1.approximate_entropy.ApEn[1]!=R2.approximate_entropy.ApEn[1])
	   return 0;
	return 1;
}

int
equal_cusum_results()
{
	if (R1.cusum.z != R2.cusum.z)
		return 0;
	if (R1.cusum.zrev != R2.cusum.zrev)
		return 0;
	if (R1.cusum.sum1A != R2.cusum.sum1A)
		return 0;
	if (R1.cusum.sum1B != R2.cusum.sum1B)
		return 0;
	if (R1.cusum.sum2A != R2.cusum.sum2A)
		return 0;
	if (R1.cusum.sum2B != R2.cusum.sum2B)
		return 0;
	if (R1.cusum.p_valueA != R2.cusum.p_valueA)
		return 0;
	if (R1.cusum.p_valueB != R2.cusum.p_valueB)
		return 0;
	return 1;
}

int
equal_random_excursion_results()
{
	int i;

	// Sometimes the original code gives up...
	if (!R1.random_excursion.valid)
		return 1;

	if (R1.random_excursion.valid != R2.random_excursion.valid)
		return 0;
	if (R1.random_excursion.valid)
	{
		for (i = 0; i < 8; i++)
		{
			if (R1.random_excursion.p_value[i] !=
				R2.random_excursion.p_value[i])
				return 0;
			if (R1.random_excursion.sum[i] != R2.random_excursion.sum[i])
				return 0;
			if (R1.random_excursion.J[i] != R2.random_excursion.J[i])
				return 0;
			if (R1.random_excursion.x[i] != R2.random_excursion.x[i])
				return 0;
		}
	}
	return 1;
}

int
equal_random_excursion_var_results()
{
	int i;

	// Sometimes the original code gives up...
	if (!R1.random_excursion_variant.valid)
		return 1;
	//
	if (R1.random_excursion_variant.valid != R2.random_excursion_variant.valid)
		return 0;
	if (R1.random_excursion_variant.valid)
	{
		for (i = 0; i < 18; i++)
		{
			if (R1.random_excursion_variant.p_value[i] !=
				R2.random_excursion_variant.p_value[i])
				return 0;
			if (R1.random_excursion_variant.count[i] !=
				R2.random_excursion_variant.count[i])
				return 0;
			if (R1.random_excursion_variant.x[i] !=
				R2.random_excursion_variant.x[i])
				return 0;
		}
	}
	return 1;
}

int
equal_linear_complexity_results()
{
	int i;

	if (R1.linear_complexity.p_value != R2.linear_complexity.p_value)
		return 0;
	//printf("[1chi2]: %f [2chi2]: %f\n",R1.linear_complexity.chi2, R2.linear_complexity.chi2);
	if (R1.linear_complexity.chi2 != R2.linear_complexity.chi2)
		return 0;
	for (i = 0; i < 7; i++)
	{
		if (R1.linear_complexity.nu[i] != R2.linear_complexity.nu[i])
			return 0;
	}
	return 1;
}

int
equal_dft_results()
{
	if (R1.dft.p_value != R2.dft.p_value)
		return 0;
	// works also with NaN
	if (R1.dft.percentile == R1.dft.percentile)
	{
		if (R1.dft.percentile != R2.dft.percentile)
			return 0;
	}
	else
	{
		if (R2.dft.percentile == R2.dft.percentile)
			return 0;
	}
	if (R1.dft.N_l != R2.dft.N_l)
		return 0;
	if (R1.dft.N_o != R2.dft.N_o)
		return 0;
	if (R1.dft.d != R2.dft.d)
		return 0;
	return 1;
}

int
compare_results(int what)
{
	int result = 0;
	switch (what)
	{
	case TEST_FREQUENCY:
		result = equal_frequency_results();
		break;
	case TEST_BLOCK_FREQUENCY:
		result = equal_blockfrequency_results();
		break;
	case TEST_RUNS:
		result = equal_runs_results();
		break;
	case TEST_LONGEST_RUN:
		result = equal_longestrunofones_results();
		break;
	case TEST_RANK:
		result = equal_rank_results();
		break;
	case TEST_SERIAL:
		result = equal_serial_results();
		break;
	case TEST_NONPERIODIC:
		result = equal_nonoverlapping_results();
		break;
	case TEST_OVERLAPPING:
		result = equal_overlapping_results();
		break;
	case TEST_UNIVERSAL:
		result = equal_universal_results();
		break;
	case TEST_APEN:
		result = equal_apen_results();
		break;
	case TEST_CUSUM:
		result = equal_cusum_results();
		break;
	case TEST_RND_EXCURSION:
		result = equal_random_excursion_results();
		break;
	case TEST_RND_EXCURSION_VAR:
		result = equal_random_excursion_var_results();
		break;
	case TEST_LINEARCOMPLEXITY:
		result = equal_linear_complexity_results();
		break;
	case TEST_FFT:
		result = equal_dft_results();
		break;
	}

	if (result)
		; // printf("");
	else
		printf("!");

	return result;
}

#define DATA_ALL_ZEROS 0
#define DATA_ALL_ZEROS_PADDED 1
#define DATA_ALL_ONES 2
#define DATA_ALL_ONES_PADDED 3
#define DATA_PATTERN_01 4
#define DATA_PATTERN_01_PADDED 5
#define DATA_PATTERN_10 6
#define DATA_PATTERN_10_PADDED 7
#define DATA_PRANDOM 8
#define DATA_PRANDOM_PADDED 9
#define DATA_PRANDOM_PADDED_RANDOM 10
#define DATA_BAD_PRANDOM 11
#define DATA_BAD_PRANDOM_PADDED 12
#define DATA_BAD_PRANDOM_PADDED_RANDOM 13
#define DATA_MAX_NUMBER 13

void
prepare_data(int j, int i)
{
	switch (j)
	{
	case DATA_ALL_ZEROS:
		data_all_zeros(i);
		break;
	case DATA_ALL_ZEROS_PADDED:
		data_all_zeros_padded(i);
		break;
	case DATA_ALL_ONES:
		data_all_ones(i);
		break;
	case DATA_ALL_ONES_PADDED:
		data_all_ones_padded(i);
		break;
	case DATA_PATTERN_01:
		data_pattern_01(i);
		break;
	case DATA_PATTERN_01_PADDED:
		data_pattern_01_padded(i);
		break;
	case DATA_PATTERN_10:
		data_pattern_10(i);
		break;
	case DATA_PATTERN_10_PADDED:
		data_pattern_10_padded(i);
		break;
	case DATA_PRANDOM:
		data_prandom(i);
		break;
	case DATA_PRANDOM_PADDED:
		data_prandom_padded(i);
		break;
	case DATA_PRANDOM_PADDED_RANDOM:
		data_prandom_padded_random(i);
		break;
	case DATA_BAD_PRANDOM:
		data_bad_prandom(i);
		break;
	case DATA_BAD_PRANDOM_PADDED:
		data_bad_prandom_padded(i);
		break;
	case DATA_BAD_PRANDOM_PADDED_RANDOM:
		data_bad_prandom_padded_random(i);
		break;
	}
}


void
test(int testcase, int smallnumbers)
{
	unsigned int j, r;
	int from, to;
	int i, k;

	if (smallnumbers)
	{
		from = 1;
		to = 1000000;
	}
	else
	{
		from = 1;
		to = 1000000000;
	}

	for (i = from; i <= to; i = (smallnumbers) ? i + 1 : i * 10 + rand() % 8)
	{
		if (!smallnumbers || i % 1000 == 0)
		{
			printf("Testing %s %u bit.\n", testNames[testcase], i);
			fflush(stdout);
		}

		for (j = 0; j <= DATA_MAX_NUMBER; j++)
		{
			prepare_data(j, i);

			if (testcase == TEST_FREQUENCY)
			{
				Frequency_v1(i);
				Frequency_v2(i);
				r = compare_results(TEST_FREQUENCY);
				if (!r)
				{
					printf("i: %i, j: %i\n", i, j);
					exit(1);
				}
			}

			if (testcase == TEST_BLOCK_FREQUENCY)
			{
				for (k = 8; k <= i / 10; k++)
				{
					if (i >= k)
					{
						BlockFrequency_v1(k, i);
						BlockFrequency_v2(k, i);
						r = compare_results(TEST_BLOCK_FREQUENCY);
						if (!r)
						{
							printf("i: %i, j: %i\n", i, j);
							exit(1);

						}
					}
				}
			}

			if (testcase == TEST_CUSUM)
			{
				CumulativeSums_v1(i);
				CumulativeSums_v2(i);
				r = compare_results(TEST_CUSUM);
				if (!r)
				{
					printf("i: %i, j: %i\n", i, j);
					getchar();
					exit(1);

				}
			}

			if (testcase == TEST_RUNS)
			{
				Runs_v1(i);
				Runs_v2(i);
				r = compare_results(TEST_RUNS);
				if (!r)
				{
					printf("i: %i, j: %i\n", i, j);
					exit(1);

				}
			}

			if (testcase == TEST_LONGEST_RUN)
			{
				if (i >= 128)
				{
					LongestRunOfOnes_v1(i);
					LongestRunOfOnes_v2(i);
					r = compare_results(TEST_LONGEST_RUN);
					if (!r)
					{
						printf("i: %i, j: %i\n", i, j);
						exit(1);

					}
				}
			}

			if (testcase == TEST_RANK)
			{
				if (i > 32 * 32)
				{
					Rank_v1(i);
					Rank_v2(i);
					r = compare_results(TEST_RANK);
					if (!r)
					{
						printf("i: %i, j: %i\n", i, j);
						exit(1);

					}
				}
			}

			if (testcase == TEST_FFT)
			{
				DiscreteFourierTransform_v1(i);
				DiscreteFourierTransform_v2(i);
				r = compare_results(TEST_FFT);
				if (!r)
				{
					printf("i: %i, j: %i\n", i, j);
					exit(1);

				}
			}

			if (testcase == TEST_NONPERIODIC)
			{
				for (k = 2; k <= 16; k++)	// 2-16
				{
					if (i >= k)
					{
						NonOverlappingTemplateMatchings_v1(k, i);
						NonOverlappingTemplateMatchings_v2(k, i);
						r = compare_results(TEST_NONPERIODIC);
						if (!r)
						{
							printf("i: %i, j: %i, m: %i\n", i, j, k);
							exit(1);

						}
					}
				}
			}


			if (testcase == TEST_OVERLAPPING)
			{
				for (k = 2; k <= 10; k++)
				{
					if (i >= k)
					{
						OverlappingTemplateMatchings_v1(k, i);
						OverlappingTemplateMatchings_v2(k, i);
						r = compare_results(TEST_OVERLAPPING);
						if (!r)
						{
							printf("i: %i, j: %i, m: %i\n", i, j, k);
							exit(1);

						}
					}
				}
			}

			if (testcase == TEST_UNIVERSAL)
			{
				Universal_v1(i);
				Universal_v2(i);
				r = compare_results(TEST_UNIVERSAL);
				if (!r)
				{
					printf("i: %i, j: %i\n", i, j);
					exit(1);

				}
			}

			if (testcase == TEST_APEN)
			{
				for (k = 1; k < ((int)log2(i) - 5); k++)
				{
					if (i >= k)
					{
						ApproximateEntropy_v1(k, i);
						ApproximateEntropy_v2(k, i);
						r = compare_results(TEST_APEN);
						if (!r)
						{
							printf("i: %i, j: %i, m: %i\n", i, j, k);
							exit(1);
						}
					}
				}
			}

			if (testcase == TEST_RND_EXCURSION)
			{
				RandomExcursions_v1(i);
				RandomExcursions_v2(i);
				r = compare_results(TEST_RND_EXCURSION);
				if (!r)
				{
					printf("i: %i, j: %i\n", i, j);
					exit(1);

				}
			}


			if (testcase == TEST_RND_EXCURSION_VAR)
			{
				RandomExcursionsVariant_v1(i);
				RandomExcursionsVariant_v2(i);
				r = compare_results(TEST_RND_EXCURSION_VAR);
				if (!r)
				{
					printf("i: %i, j: %i\n", i, j);
					exit(1);

				}
			}

			if (testcase == TEST_SERIAL)
			{
				for (k = 1; k < ((int)log2(i) - 2); k++)
				{
					if (i >= k)
					{
						Serial_v1(k, i);
						Serial_v2(k, i);
						r = compare_results(TEST_SERIAL);
						if (!r)
						{
							printf("i: %i, j: %i, k: %i\n", i, j, k);
							exit(1);

						}
					}
				}
			}


			if (testcase == TEST_LINEARCOMPLEXITY)
			{
				if (i >= 1000000)
				{
					for (k = 500; k <= 5000; k++)
					{
						LinearComplexity_v1(k, i);
						LinearComplexity_v2(k, i);
						r = compare_results(TEST_LINEARCOMPLEXITY);
						if (!r)
						{
							printf("i: %i, j: %i, m: %i\n", i, j, k);
							exit(1);

						}
					}
				}
			}

			free(epsilon);
			free(array);
		}
	}
}
#endif

#ifdef SPEED

#ifdef _WIN64
// In x64 Windows systems use __rdtsc function (asm is not available)
#include <intrin.h>
int64_t GetCpuClocks()
{
	return __rdtsc();
}
#else
#ifdef _WIN32
// In win32 systems use asm
int64_t GetCpuClocks()
{
	struct
	{
		int32_t low, high;
	} clocks;
	__asm push EAX
	__asm push EDX
	__asm __emit 0x0f
	__asm __emit 0x31
	__asm mov clocks.low, EAX
	__asm mov clocks.high, EDX
	__asm pop EDX 
	__asm pop EAX 
	return *(int64_t *)(&clocks);
}
#else
// On Linux use asm
#include <stdint.h>
uint64_t
GetCpuClocks()
{
	unsigned int low, high;
	__asm__ __volatile__("rdtsc":"=a" (low), "=d" (high));
	return ((uint64_t)high << 32) | low;
}
#endif
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
clock_t
GetTickCount()
{
	return clock();
}
#endif


char *
test_name(int t)
{
	if (t <= 15)
		return testNames[t];
	switch (t)
	{
	case 16:
	case 29:
		return "BlockFrequency";
	case 17:
	case 18:
	case 30:
		return "NonOverlappingTemplate";
	case 19:
	case 20:
	case 31:
		return "OverlappingTemplate";
	case 21:
	case 22:
	case 23:
	case 34:
	case 35:
	case 36:
		return "ApproximateEntropy";
	case 24:
	case 25:
	case 32:
	case 26:
	case 33:
		return "Serial";
	case 27:
	case 28:
		return "LinearComplexity";
	}
	return "Unknown";
}

#include <locale.h>

void
flushcache()
{
	int i, j;
	const int size = 20 * 1024 * 1024;	// Allocate 20M. 
	char *c = (char *)malloc(size);
	if (c == NULL)
		exit(100);
	for (i = 0; i < 0xff; i++)
		for (j = 0; j < size; j++)
			c[j] = i * j;
	free(c);
}

void
speed(int scale, int repeat, int test_from, int test_to)
{
	int n, t, j, from, to, param;
	FILE *f;
#ifdef VERIFY_RESULTS
	int r;
#endif
	uint64_t Astart, Amiddle1, Amiddle2, Aend, Atime1, Atime2, Amintime1, Amintime2;
#ifdef _WIN32
	unsigned long Bstart, Bmiddle1, Bmiddle2, Bend, Btime1, Btime2, Bmintime1, Bmintime2;
	//      LARGE_INTEGER Cstart, Cmiddle, Cend;
	//      __int64 Ctime1, Ctime2, Dtime1, Dtime2;
	//      FILETIME Dstart,Dmiddle,Dend,Dtest1,Dtest2,Dtest3,Dvoid;
	DWORD affinity;
#else
	clock_t Bstart, Bmiddle1, Bmiddle2, Bend, Btime1, Btime2, Bmintime1, Bmintime2;
#endif

#ifdef _WIN32

	affinity = 1;
	SetThreadAffinityMask(GetCurrentThread(), affinity);
	SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
#endif

	setlocale(LC_NUMERIC, "Czech");
	f = fopen("speed.csv", "wt");
	if (!f)
	{
		printf("Cannot open file speed.csv for writing.\n");
		exit(3);
	};

	if (scale)
	{
		//from=1024*1*8+1;
		//to=1024*1024*8*100;
		//from = 1024 * 1024+1;
		//to = 1024 * 1024+1;
		//from = 1000000; to = 1000000;
		//from = 1024 * 1024; to = 1024 * 1024;
		from = 3 * 3 * 3 * 3 * 3;
		to = 1024 * 1024 * 1024;
	}
	else
	{
		from = 1024 * 1024 * 8 * 20;
		to = 1024 * 1024 * 8 * 20;
	}

	for (n = from; n <= to;
	n = (scale) ? n * 3 /*((n-1)*2)+1 *//* *10 +rand()%8 */ : n + 1)
	{
		for (t = test_from; t <= test_to; t++)
		{
			data_prandom(n);

			for (j = 1; j <= repeat; j++)
			{
				flushcache();
				Astart = GetCpuClocks();
				Bstart = GetTickCount();
#ifdef _WIN32
				/*
				   GetProcessTimes(GetCurrentProcess(),&Dvoid,&Dvoid,&Dtest1,&Dstart);
				   QueryPerformanceCounter(&Cstart);
				 */
#endif

				param = 0;

				switch (t)
				{
				case TEST_FREQUENCY:
					Frequency_v1(n);
					break;
				case TEST_BLOCK_FREQUENCY:
					if (n >= 100)
						BlockFrequency_v1(n / 100, n);
					param = n / 100;
					break;
				case TEST_CUSUM:
					CumulativeSums_v1(n);
					break;
				case TEST_RUNS:
					Runs_v1(n);
					break;
				case TEST_LONGEST_RUN:
					LongestRunOfOnes_v1(n);
					break;
				case TEST_RANK:
					if (n > 32 * 32)
						Rank_v1(n);
					break;
				case TEST_FFT:
					DiscreteFourierTransform_v1(n);
					break;
				case TEST_NONPERIODIC:
					NonOverlappingTemplateMatchings_v1(10, n);
					param = 10;
					break;
				case TEST_OVERLAPPING:
					OverlappingTemplateMatchings_v1(10, n);
					param = 10;
					break;
				case TEST_UNIVERSAL:
					Universal_v1(n);
					break;
				case TEST_APEN:
					ApproximateEntropy_v1(9, n);
					param = 9;
					break;
				case TEST_RND_EXCURSION:
					RandomExcursions_v1(n);
					break;
				case TEST_RND_EXCURSION_VAR:
					RandomExcursionsVariant_v1(n);
					break;
				case TEST_SERIAL:
					if (log2(n) >= 9)
						Serial_v1(9, n);
					param = 9;
					break;
				case TEST_LINEARCOMPLEXITY:
					LinearComplexity_v1(5000, n);
					param = 5000;
					break;
				case 16:
					BlockFrequency_v1(20, n);
					param = 20;
					break;
				case 17:
					NonOverlappingTemplateMatchings_v1(2, n);
					param = 2;
					break;
				case 18:
					NonOverlappingTemplateMatchings_v1(9, n);
					param = 9;
					break;
				case 19:
					OverlappingTemplateMatchings_v1(2, n);
					param = 2;
					break;
				case 20:
					OverlappingTemplateMatchings_v1(9, n);
					param = 9;
					break;
				case 21:
					ApproximateEntropy_v1(2, n);
					param = 2;
					break;
				case 22:
					ApproximateEntropy_v1(5, n);
					param = 5;
					break;
				case 23:
					ApproximateEntropy_v1(24, n);
					param = 24;
					break;
				case 24:
					if (log2(n) >= 2)
						Serial_v1(2, n);
					param = 2;
					break;
				case 25:
					if (log2(n) >= 24)
						Serial_v1(24, n);
					param = 24;
					break;
				case 26:
					if (log2(n) >= 5)
						Serial_v1(5, n);
					param = 5;
					break;
				case 27:
					LinearComplexity_v1(500, n);
					param = 500;
					break;
				case 28:
					LinearComplexity_v1(1000, n);
					param = 1000;
					break;
				case 29:
					BlockFrequency_v1(128, n);
					param = 128;
					break;
				case 30:
					NonOverlappingTemplateMatchings_v1(21, n);
					param = 21;
					break;
				case 31:
					OverlappingTemplateMatchings_v1(24, n);
					param = 24;
					break;
				case 32:
					if (log2(n) >= 15)
						Serial_v1(13, n);
					param = 13;
					break;
				case 33:
					if (log2(n) >= 17)
						Serial_v1(14, n);
					param = 14;
					break;
				case 34:
					ApproximateEntropy_v1(27, n);
					param = 27;
					break;
				case 35:
					ApproximateEntropy_v1(8, n);
					param = 8;
					break;
				case 36:
					ApproximateEntropy_v1(10, n);
					param = 10;
					break;
				}

				Amiddle1 = GetCpuClocks();
				Bmiddle1 = GetTickCount();

				if (dummy_result == 1.23456789)printf("\n");
#ifdef _WIN32
				/*
				   QueryPerformanceCounter(&Cmiddle);
				   GetProcessTimes(GetCurrentProcess(),&Dvoid,&Dvoid,&Dtest2,&Dmiddle);
				 */
#endif
				flushcache();
				Amiddle2 = GetCpuClocks();
				Bmiddle2 = GetTickCount();

				switch (t)
				{
				case TEST_FREQUENCY:
					Frequency_v2(n);
					break;
				case TEST_BLOCK_FREQUENCY:
					if (n >= 100)
						BlockFrequency_v2(n / 100, n);
					break;
				case TEST_CUSUM:
					CumulativeSums_v2(n);
					break;
				case TEST_RUNS:
					Runs_v2(n);
					break;
				case TEST_LONGEST_RUN:
					LongestRunOfOnes_v2(n);
					break;
				case TEST_RANK:
					if (n > 32 * 32)
						Rank_v2(n);
					break;
				case TEST_FFT:
					DiscreteFourierTransform_v2(n);
					break;
				case TEST_NONPERIODIC:
					NonOverlappingTemplateMatchings_v2(10, n);
					break;
				case TEST_OVERLAPPING:
					OverlappingTemplateMatchings_v2(10, n);
					break;
				case TEST_UNIVERSAL:
					Universal_v2(n);
					break;
				case TEST_APEN:
					ApproximateEntropy_v2(9, n);
					break;
				case TEST_RND_EXCURSION:
					RandomExcursions_v2(n);
					break;
				case TEST_RND_EXCURSION_VAR:
					RandomExcursionsVariant_v2(n);
					break;
				case TEST_SERIAL:
					if (log2(n) >= 9)
						Serial_v2(9, n);
					break;
				case TEST_LINEARCOMPLEXITY:
					LinearComplexity_v2(5000, n);
					break;
				case 16:
					BlockFrequency_v2(20, n);
					break;
				case 17:
					NonOverlappingTemplateMatchings_v2(2, n);
					break;
				case 18:
					NonOverlappingTemplateMatchings_v2(9, n);
					break;
				case 19:
					OverlappingTemplateMatchings_v2(2, n);
					break;
				case 20:
					OverlappingTemplateMatchings_v2(9, n);
					break;
				case 21:
					ApproximateEntropy_v2(2, n);
					break;
				case 22:
					ApproximateEntropy_v2(5, n);
					break;
				case 23:
					ApproximateEntropy_v2(24, n);
					break;
				case 24:
					if (log2(n) >= 2)
						Serial_v2(2, n);
					break;
				case 25:
					if (log2(n) >= 24)
						Serial_v2(24, n);
					break;
				case 26:
					if (log2(n) >= 5)
						Serial_v2(5, n);
					break;
				case 27:
					LinearComplexity_v2(500, n);
					break;
				case 28:
					LinearComplexity_v2(1000, n);
					break;
				case 29:
					BlockFrequency_v2(128, n);
					break;
				case 30:
					NonOverlappingTemplateMatchings_v2(21, n);
					break;
				case 31:
					OverlappingTemplateMatchings_v2(24, n);
					break;
				case 32:
					if (log2(n) >= 15)
						Serial_v2(13, n);
					break;
				case 33:
					if (log2(n) >= 17)
						Serial_v2(14, n);
					break;
				case 34:
					ApproximateEntropy_v2(27, n);
					break;
				case 35:
					ApproximateEntropy_v2(8, n);
					break;
				case 36:
					ApproximateEntropy_v2(10, n);
					break;
				}

				Aend = GetCpuClocks();
				Bend = GetTickCount();

				if (dummy_result == 1.23456789)printf("\n");
#ifdef _WIN32
				/*
				   QueryPerformanceCounter(&Cend);
				   GetProcessTimes(GetCurrentProcess(),&Dvoid,&Dvoid,&Dtest3,&Dend);
				 */
#endif

				Atime1 = Amiddle1 - Astart;
				Atime2 = Aend - Amiddle2;
				Btime1 = Bmiddle1 - Bstart;
				Btime2 = Bend - Bmiddle2;

				if (j == 1)
				{
					Amintime1 = Atime1;
					Amintime2 = Atime2;
					Bmintime1 = Btime1;
					Bmintime2 = Btime2;
				}

				if (Atime1 < Amintime1)
					Amintime1 = Atime1;
				if (Atime2 < Amintime2)
					Amintime2 = Atime2;
				if (Btime1 < Bmintime1)
					Bmintime1 = Btime1;
				if (Btime2 < Bmintime2)
					Bmintime2 = Btime2;


#ifdef _WIN32
				/*
				   Ctime1=Cmiddle.QuadPart-Cstart.QuadPart;
				   Ctime2=Cend.QuadPart-Cmiddle.QuadPart;
				   Dtime1=((((__int64)Dmiddle.dwHighDateTime)<<32) | (__int64)(Dmiddle.dwLowDateTime)) - ((((__int64)Dstart.dwHighDateTime)<<32) | (__int64)Dstart.dwLowDateTime);
				   Dtime2=((((__int64)Dend.dwHighDateTime)<<32) | (__int64)(Dend.dwLowDateTime)) - ((((__int64)Dmiddle.dwHighDateTime)<<32) | (__int64)Dmiddle.dwLowDateTime);
				 */
#endif

				 //printf("%s (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %ul ms vs. %ul ms], %f x faster [QueryPerformanceCounter], %f x faster [GetProcessTimes]\n",testNames[t],n,(float)Atime1/Atime2, (float)Btime1/Btime2, Btime1, Btime2, (float)Ctime1/Ctime2, (float)Dtime1/Dtime2);
#ifdef _WIN32
				printf("%s [%i] (%i bits): %f x faster [CPU clocks: %I64u vs. %I64u], %f x faster [GetTickCount - %lu ms vs. %lu ms].\n",
#else
				printf("%s [%i] (%i bits): %f x faster [CPU clocks: %" PRIu64 " vs. %" PRIu64 "], %f x faster [GetTickCount - %lu ms vs. %lu ms].\n",
#endif
					test_name(t), param, n, (float)Atime1 / Atime2, Atime1, Atime2, (float)Btime1 / Btime2, Btime1, Btime2);
				//printf("%s (%i bits): %f x faster\n",testNames[t],n,(float)Atime1/Atime2);
				//printf("%s (%i bits): %f x faster {GetProcessTimes}\n",testNames[t],n,(float)Dtime1/Dtime2);
				//fprintf(f, "%s [%i];%I64i;%I64i;%f;%lu;%lu;%f\n", test_name(t),param, Atime1, Atime2, (float)Atime1 / Atime2, Btime1, Btime2, (float)Btime1 / Btime2);
				fflush(stdout);
				//fflush(f);

#ifdef VERIFY_RESULTS
				r = compare_results(t);
				if (!r)
				{
					printf
						("%s (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %lu ms vs. %lu ms].\n",
							test_name(t), n, (float)Atime1 / Atime2,
							(float)Btime1 / Btime2, Btime1, Btime2);
					printf("Results DO NOT MATCH\n");
					//exit(2);
				}
#endif

			}
#ifdef _WIN32
			printf("MINIMUM - %s (%i bits) [%i]: %f x faster [CPU clocks: %I64u vs. %I64u], %f x faster [ms: %lu vs. %lu]\n",
#else
			printf("MINIMUM - %s (%i bits) [%i]: %f x faster [CPU clocks: %" PRIu64 " vs. %" PRIu64 "], %f x faster [ms: %lu vs. %lu]\n",
#endif
					test_name(t), n, param, (float)Amintime1 / Amintime2, Amintime1, Amintime2, 
				    (float)Bmintime1 / Bmintime2, Bmintime1, Bmintime2);
#ifdef _WIN32
			fprintf(f, "%s [%i];%I64i;%I64i;%f;%lu;%lu;%f\n",
#else
			fprintf(f, "%s [%i];%" PRIu64 ";%" PRIu64 ";%f;%lu;%lu;%f\n", 
#endif
				test_name(t), param, Amintime1, Amintime2, (float)Amintime1 / Amintime2,
				Bmintime1, Bmintime2, (float)Bmintime1 / Bmintime2);
			fflush(f);
			free(epsilon);
			free(array);
	}
		printf("Done n=%i.\n", n);
}
	fclose(f);
}
#endif

/*
void TestRank()
{

	int n, j, from, to, repeat;
	FILE *f;
	__int64 Astart, Amiddle, Aend, Atime1, Atime2, Amintime1, Amintime2;
	unsigned long Bstart, Bmiddle, Bend, Btime1, Btime2, Bmintime1, Bmintime2;
	DWORD affinity;

	affinity = 1;
	SetThreadAffinityMask(GetCurrentThread(), affinity);
	SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);

	setlocale(LC_NUMERIC, "Czech");
	f = fopen("speed_rank.csv", "wt");
	if (!f) { printf("Cannot open file speed.csv for writing.\n"); exit(3); };

	from = 1024 * 1024 * 8 * 150;
	to = 1024 * 1024 * 8 * 200;
	repeat = 10;

	for (n = from; n <= to; n += 1024*1024*8)
	{
		data_prandom(n);

		for (j = 1; j <= repeat; j++)
		{
			Astart = GetCpuClocks();
			Bstart = GetTickCount();
			if (n > 32 * 32) Rank(n);
			Amiddle = GetCpuClocks();
			Bmiddle = GetTickCount();
			if (n > 32 * 32) Rank2(n);

			Aend = GetCpuClocks();
			Bend = GetTickCount();

			Atime1 = Amiddle - Astart;
			Atime2 = Aend - Amiddle;
			Btime1 = Bmiddle - Bstart;
			Btime2 = Bend - Bmiddle;

			if (j == 1)
			{
				Amintime1 = Atime1;
				Amintime2 = Atime2;
				Bmintime1 = Btime1;
				Bmintime2 = Btime2;
			}

			if (Atime1 < Amintime1) Amintime1 = Atime1;
			if (Atime2 < Amintime2) Amintime2 = Atime2;
			if (Btime1 < Bmintime1) Bmintime1 = Btime1;
			if (Btime2 < Bmintime2) Bmintime2 = Btime2;

			printf("Rank (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %lu ms vs. %lu ms].\n", n, (float)Atime1 / Atime2, (float)Btime1 / Btime2, Btime1, Btime2);
			fflush(stdout);
		}
		printf("MINIMUM - Rank (%i bits): %f x faster [CPU clocks: %I64i vs. %I64i], %f x faster [ms: %i vs. %i]\n", n, (float)Amintime1 / Amintime2, Amintime1, Amintime2, (float)Bmintime1 / Bmintime2, Bmintime1, Bmintime2);
		fprintf(f, "Rank; %i; %i; %I64i; %I64i; %f; %lu; %lu; %f\n", n, n / 8,  Amintime1, Amintime2, (float)Amintime1 / Amintime2, Bmintime1, Bmintime2, (float)Bmintime1 / Bmintime2);
		fflush(f);
		free(epsilon); free(array);
		printf("Done n=%i.\n", n);
	}
	fclose(f);
}
*/

/*
int TestNonParam[] = { 2, 4, 6, 10, 16, 21};

void TestNonOver()
{

	int n, j, from, to, repeat, param, index;
	FILE *f;
	__int64 Astart, Amiddle, Aend, Atime1, Atime2, Amintime1, Amintime2;
	unsigned long Bstart, Bmiddle, Bend, Btime1, Btime2, Bmintime1, Bmintime2;
	DWORD affinity;

	affinity = 1;
	SetThreadAffinityMask(GetCurrentThread(), affinity);
	SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);

	setlocale(LC_NUMERIC, "Czech");
	f = fopen("speed_nonover.csv", "wt");
	if (!f) { printf("Cannot open file speed.csv for writing.\n"); exit(3); };

	from = 1024 * 1024 * 8;
	to = 1024 * 1024 * 8 * 200;
	repeat = 3;

	for (n = from; n <= to; n += 1024 * 1024 * 8)
	{
		data_prandom(n);
		//for (index = 0, param = TestNonParam[0]; index<=5; index++, param=TestNonParam[index])
		for (param=9; param<=9; param++)
		{

			for (j = 1; j <= repeat; j++)
			{
				Astart = GetCpuClocks();
				Bstart = GetTickCount();
				NonOverlappingTemplateMatchings(param, n);
				Amiddle = GetCpuClocks();
				Bmiddle = GetTickCount();
				NonOverlappingTemplateMatchings_v2(param, n);

				Aend = GetCpuClocks();
				Bend = GetTickCount();

				Atime1 = Amiddle - Astart;
				Atime2 = Aend - Amiddle;
				Btime1 = Bmiddle - Bstart;
				Btime2 = Bend - Bmiddle;

				if (j == 1)
				{
					Amintime1 = Atime1;
					Amintime2 = Atime2;
					Bmintime1 = Btime1;
					Bmintime2 = Btime2;
				}

				if (Atime1 < Amintime1) Amintime1 = Atime1;
				if (Atime2 < Amintime2) Amintime2 = Atime2;
				if (Btime1 < Bmintime1) Bmintime1 = Btime1;
				if (Btime2 < Bmintime2) Bmintime2 = Btime2;

				printf("NonOver [%i] (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %lu ms vs. %lu ms].\n", param, n, (float)Atime1 / Atime2, (float)Btime1 / Btime2, Btime1, Btime2);
				fflush(stdout);
			}
			printf("MINIMUM - NonOver [%i] (%i bits): %f x faster [CPU clocks: %I64i vs. %I64i], %f x faster [ms: %i vs. %i]\n", param, n, (float)Amintime1 / Amintime2, Amintime1, Amintime2, (float)Bmintime1 / Bmintime2, Bmintime1, Bmintime2);
			fprintf(f, "NonOver;%i; %i; %i; %I64i; %I64i; %f; %lu; %lu; %f\n", param, n, n / 8, Amintime1, Amintime2, (float)Amintime1 / Amintime2, Bmintime1, Bmintime2, (float)Bmintime1 / Bmintime2);
			fflush(f);
		}
		free(epsilon); free(array);
		printf("Done n=%i.\n", n);
	}
	fclose(f);
}

void TestLinCom()
{

	int n, j, from, to, repeat,param;
	FILE *f;
	__int64 Astart, Amiddle, Aend, Atime1, Atime2, Amintime1, Amintime2;
	unsigned long Bstart, Bmiddle, Bend, Btime1, Btime2, Bmintime1, Bmintime2;
	DWORD affinity;

	affinity = 1;
	SetThreadAffinityMask(GetCurrentThread(), affinity);
	SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);

	setlocale(LC_NUMERIC, "Czech");
	f = fopen("speed_lincom.csv", "wt");
	if (!f) { printf("Cannot open file speed.csv for writing.\n"); exit(3); };

	from = 1024 * 1024 * 8 * 10;
	to = 1024 * 1024 * 8 * 10;
	repeat = 3;

	for (n = from; n <= to; n += 1024 * 1024 * 8)
	{
		data_prandom(n);
		for (param = 500; param <= 5000; param+=100) // 500
		{

			for (j = 1; j <= repeat; j++)
			{
				Astart = GetCpuClocks();
				Bstart = GetTickCount();
				LinearComplexity(param, n);
				Amiddle = GetCpuClocks();
				Bmiddle = GetTickCount();
				LinearComplexity3(param, n);

				Aend = GetCpuClocks();
				Bend = GetTickCount();

				Atime1 = Amiddle - Astart;
				Atime2 = Aend - Amiddle;
				Btime1 = Bmiddle - Bstart;
				Btime2 = Bend - Bmiddle;

				if (j == 1)
				{
					Amintime1 = Atime1;
					Amintime2 = Atime2;
					Bmintime1 = Btime1;
					Bmintime2 = Btime2;
				}

				if (Atime1 < Amintime1) Amintime1 = Atime1;
				if (Atime2 < Amintime2) Amintime2 = Atime2;
				if (Btime1 < Bmintime1) Bmintime1 = Btime1;
				if (Btime2 < Bmintime2) Bmintime2 = Btime2;

				printf("LinCom [%i] (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %lu ms vs. %lu ms].\n", param, n, (float)Atime1 / Atime2, (float)Btime1 / Btime2, Btime1, Btime2);
				fflush(stdout);
			}
			printf("MINIMUM - LinCom [%i] (%i bits): %f x faster [CPU clocks: %I64i vs. %I64i], %f x faster [ms: %i vs. %i]\n", param, n, (float)Amintime1 / Amintime2, Amintime1, Amintime2, (float)Bmintime1 / Bmintime2, Bmintime1, Bmintime2);
			fprintf(f, "LinCom;%i; %i; %i; %I64i; %I64i; %f; %lu; %lu; %f\n", param, n, n / 8, Amintime1, Amintime2, (float)Amintime1 / Amintime2, Bmintime1, Bmintime2, (float)Bmintime1 / Bmintime2);
			fflush(f);
		}
		free(epsilon); free(array);
		printf("Done n=%i.\n", n);
	}
	fclose(f);
}
*/

#ifdef VERIFY_RESULTS
#ifndef SPEED
int
main(int argc, char **argv)
{
	int testcase;

	if (argc >= 2)
		testcase = atoi(argv[1]);
	else
		testcase = 1;

	test(testcase, 0);
	test(testcase, 1);

	return 0;
}
#endif
#endif

#ifdef SPEED
int
main(int argc, char **argv)
{

	int scale, repeat, test_from, test_to;

	if (argc >= 2)
		scale = atoi(argv[1]);
	else
		scale = 0;
	if (argc >= 3)
		repeat = atoi(argv[2]);
	else
		repeat = 10;
	if (argc >= 4)
		test_from = atoi(argv[3]);
	else
		test_from = 1;
	if (argc >= 5)
		test_to = atoi(argv[4]);
	else
		test_to = 36;

	speed(scale, repeat, test_from, test_to);

	//TestRank();
	//TestNonOver();
	//TestLinCom();
	return 0;
}
#endif

#ifdef KS
#define MAXBUFFER 2000000
#define THRESHOLD 0.000002
int fails = 0;

void check(double real, double expected)
{
	double difference = fabs(real - expected);
	if (difference < THRESHOLD)
	{
		printf("OK\n");

	}
	else
	{
		fails++;
		printf("FAIL\n");
	}
}

void selfcheck(void)
{
	FILE *data = NULL;
	unsigned char *buffer;
	size_t r, el;
	int i;

	fails = 0;
	buffer = malloc(MAXBUFFER);
	if (buffer == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.frequency_pvals = malloc(sizeof(double) * 1);
	if (pvals.frequency_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.blockfrequency_pvals = malloc(sizeof(double) * 1);
	if (pvals.blockfrequency_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.cusum_pvals[0] = malloc(sizeof(double) * 1);
	if (pvals.cusum_pvals[0] == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.cusum_pvals[1] = malloc(sizeof(double) * 1);
	if (pvals.cusum_pvals[1] == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.runs_pvals = malloc(sizeof(double) * 1);
	if (pvals.runs_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.longestrunofones_pvals = malloc(sizeof(double) * 1);
	if (pvals.longestrunofones_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.rank_pvals = malloc(sizeof(double) * 1);
	if (pvals.rank_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.dft_pvals = malloc(sizeof(double) * 1);
	if (pvals.dft_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	for (i = 0; i < MAXNUMOFTEMPLATES; i++)
	{
		pvals.nonoverlapping_pvals[i] = (double*)malloc(sizeof(double) * 1);
		if (pvals.nonoverlapping_pvals[i] == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	}
	pvals.overlapping_pvals = malloc(sizeof(double) * 1);
	if (pvals.overlapping_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.universal_pvals = malloc(sizeof(double) * 1);
	if (pvals.universal_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.approximate_entropy_pvals = malloc(sizeof(double) * 1);
	if (pvals.approximate_entropy_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	for (i = 0; i < 8; i++)
	{
		pvals.random_excursion_pvals[i] = malloc(sizeof(double) * 1);
		if (pvals.random_excursion_pvals[i] == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	}
	for (i = 0; i < 18; i++)
	{
		pvals.random_excursion_variant_pvals[i] = malloc(sizeof(double) * 1);
		if (pvals.random_excursion_variant_pvals[i] == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	}
	pvals.serial_pvals[0] = malloc(sizeof(double) * 1);
	if (pvals.serial_pvals[0] == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.serial_pvals[1] = malloc(sizeof(double) * 1);
	if (pvals.serial_pvals[1] == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.linear_complexity_pvals = malloc(sizeof(double) * 1);
	if (pvals.linear_complexity_pvals == NULL) { fprintf(stderr, "Cannot allocate memory for buffer.\n"); return; }
	pvals.seq_counter = 0;
	array = buffer;

	// data.pi
	el = 125000;
	data = fopen("data/pi.bin", "rb");
	if (data == NULL)
	{
		fprintf(stderr, "Cannot open 'data/pi.bin' file to load data for testing.\n"); 
		return;
	}
	r=fread(buffer, 1, MAXBUFFER, data);
	if (r != el)
	{
		fprintf(stderr, "Cannot read from 'data/pi.bin' (file length mismatch).\n"); 
		fclose(data);
		return;
	}
	printf("Selftest using the 'pi' binary file:\n\n");
	printf("Testing Frequency test ........................ ");
	Frequency_v2(8 * (int) el);
	check(pvals.frequency_pvals[0], 0.578211);
	printf("Testing Block Frequency test .................. ");
	BlockFrequency_v2(128, 8 * (int)el);
	check(pvals.blockfrequency_pvals[0], 0.380615);
	printf("Testing Cumulative Sum test - Forward ......... ");
	CumulativeSums_v2(8 * (int)el);
	check(*pvals.cusum_pvals[0], 0.628308);
	printf("Testing Cumulative Sum test - Reverse ......... ");	
	check(*pvals.cusum_pvals[1], 0.663369);
	printf("Testing Runs test ............................. ");
	Runs_v2(8 * (int)el);
	check(pvals.runs_pvals[0], 0.419268);
	printf("Testing Longest Run test ...................... ");
	LongestRunOfOnes_v2(8 * (int)el);
	check(pvals.longestrunofones_pvals[0], 0.024390);
	printf("Testing Rank test ............................. ");
	Rank_v2(8 * (int)el);
	check(pvals.rank_pvals[0], 0.083553);
	printf("Testing Spectral test ......................... ");
	DiscreteFourierTransform_v2(8 * (int)el);
	check(pvals.dft_pvals[0], 0.010186);
	printf("Testing Non-overlapping Templates test ........ ");
	NonOverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.nonoverlapping_pvals[0][0], 0.165757);
	printf("Testing Overlapping Templates test ............ ");
	OverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.overlapping_pvals[0], 0.296897);
	printf("Testing Universal test ........................ ");
	Universal_v2(8 * (int)el);
	check(pvals.universal_pvals[0], 0.669012);
	printf("Testing Approximate Entropy test .............. ");
	ApproximateEntropy_v2(10, 8 * (int)el);
	check(pvals.approximate_entropy_pvals[0], 0.361595);
	printf("Testing Random Excursion test ................. ");
	RandomExcursions_v2(8 * (int)el);
	check(pvals.random_excursion_pvals[4][0], 0.844143);
	printf("Testing Random Excursion Variant test ......... ");
	RandomExcursionsVariant_v2(8 * (int)el);
	check(pvals.random_excursion_variant_pvals[8][0], 0.760966);
	printf("Testing Linear Complexity test ................ ");
	LinearComplexity_v2(500, 8 * (int)el);
	check(pvals.linear_complexity_pvals[0], 0.255475);
	printf("Testing Serial test ........................... ");
	Serial_v2(16, 8 * (int)el);
	check(pvals.serial_pvals[0][0], 0.143005);
	fclose(data);

	// data.e
	el = 125000;
	data = fopen("data/e.bin", "rb");
	if (data == NULL)
	{
		fprintf(stderr, "Cannot open 'data/e.bin' file to load data for testing.\n");
		return;
	}
	r = fread(buffer, 1, MAXBUFFER, data);
	if (r != el)
	{
		fprintf(stderr, "Cannot read from 'data/e.bin' (file length mismatch).\n");
		fclose(data);
		return;
	}
	printf("\n\nSelftest using the 'e' binary file:\n\n");
	printf("Testing Frequency test ........................ ");
	Frequency_v2(8 * (int)el);
	check(pvals.frequency_pvals[0], 0.953749);
	printf("Testing Block Frequency test .................. ");
	BlockFrequency_v2(128, 8 * (int)el);
	check(pvals.blockfrequency_pvals[0], 0.211072);
	printf("Testing Cumulative Sum test - Forward ......... ");
	CumulativeSums_v2(8 * (int)el);
	check(pvals.cusum_pvals[0][0], 0.669887);
	printf("Testing Cumulative Sum test - Reverse ......... ");
	check(pvals.cusum_pvals[1][0], 0.724266);
	printf("Testing Runs test ............................. ");
	Runs_v2(8 * (int)el);
	check(pvals.runs_pvals[0], 0.561917);
	printf("Testing Longest Run test ...................... ");
	LongestRunOfOnes_v2(8 * (int)el);
	check(pvals.longestrunofones_pvals[0], 0.718945);
	printf("Testing Rank test ............................. ");
	Rank_v2(8 * (int)el);
	check(pvals.rank_pvals[0], 0.306156);
	printf("Testing Spectral test ......................... ");
	DiscreteFourierTransform_v2(8 * (int)el);
	check(pvals.dft_pvals[0], 0.847187);
	printf("Testing Non-overlapping Templates test ........ ");
	NonOverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.nonoverlapping_pvals[0][0], 0.078790);
	printf("Testing Overlapping Templates test ............ ");
	OverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.overlapping_pvals[0], 0.110434);
	printf("Testing Universal test ........................ ");
	Universal_v2(8 * (int)el);
	check(pvals.universal_pvals[0], 0.282568);
	printf("Testing Approximate Entropy test .............. ");
	ApproximateEntropy_v2(10, 8 * (int)el);
	check(pvals.approximate_entropy_pvals[0], 0.700073);
	printf("Testing Random Excursion test ................. ");
	RandomExcursions_v2(8 * (int)el);
	check(pvals.random_excursion_pvals[4][0], 0.786868);
	printf("Testing Random Excursion Variant test ......... ");
	RandomExcursionsVariant_v2(8 * (int)el);
	check(pvals.random_excursion_variant_pvals[8][0], 0.826009);
	printf("Testing Linear Complexity test ................ ");
	LinearComplexity_v2(500, 8 * (int)el);
	check(pvals.linear_complexity_pvals[0], 0.826335);
	printf("Testing Serial test ........................... ");
	Serial_v2(16, 8 * (int)el);
	check(*pvals.serial_pvals[0], 0.766182);
	fclose(data);

	// data.sha1
	el = 125000;
	data = fopen("data/sha1.bin", "rb");
	if (data == NULL)
	{
		fprintf(stderr, "Cannot open 'data/sha1.bin' file to load data for testing.\n");
		return;
	}
	r = fread(buffer, 1, MAXBUFFER, data);
	if (r != el)
	{
		fprintf(stderr, "Cannot read from 'data/sha1.bin' (file length mismatch).\n");
		fclose(data);
		return;
	}
	printf("\n\nSelftest using the 'sha1' binary file:\n\n");
	printf("Testing Frequency test ........................ ");
	Frequency_v2(8 * (int)el);
	check(pvals.frequency_pvals[0], 0.604458);
	printf("Testing Block Frequency test .................. ");
	BlockFrequency_v2(128, 8 * (int)el);
	check(pvals.blockfrequency_pvals[0], 0.091517);
	printf("Testing Cumulative Sum test - Forward ......... ");
	CumulativeSums_v2(8 * (int)el);
	check(pvals.cusum_pvals[0][0], 0.451231);
	printf("Testing Cumulative Sum test - Reverse ......... ");
	check(pvals.cusum_pvals[1][0], 0.550134);
	printf("Testing Runs test ............................. ");
	Runs_v2(8 * (int)el);
	check(pvals.runs_pvals[0], 0.309757);
	printf("Testing Longest Run test ...................... ");
	LongestRunOfOnes_v2(8 * (int)el);
	check(pvals.longestrunofones_pvals[0], 0.657812);
	printf("Testing Rank test ............................. ");
	Rank_v2(8 * (int)el);
	check(pvals.rank_pvals[0], 0.577829);
	printf("Testing Spectral test ......................... ");
	DiscreteFourierTransform_v2(8 * (int)el);
	check(pvals.dft_pvals[0], 0.163062);
	printf("Testing Non-overlapping Templates test ........ ");
	NonOverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.nonoverlapping_pvals[0][0], 0.496601);
	printf("Testing Overlapping Templates test ............ ");
	OverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.overlapping_pvals[0], 0.339426);
	printf("Testing Universal test ........................ ");
	Universal_v2(8 * (int)el);
	check(pvals.universal_pvals[0], 0.411079);
	printf("Testing Approximate Entropy test .............. ");
	ApproximateEntropy_v2(10, 8 * (int)el);
	check(pvals.approximate_entropy_pvals[0], 0.982885);
	printf("Testing Random Excursion test ................. ");
	RandomExcursions_v2(8 * (int)el);
	check(pvals.random_excursion_pvals[4][0], 0.000000);
	printf("Testing Random Excursion Variant test ......... ");
	RandomExcursionsVariant_v2(8 * (int)el);
	check(pvals.random_excursion_variant_pvals[8][0], 0.000000);
	printf("Testing Linear Complexity test ................ ");
	LinearComplexity_v2(500, 8 * (int)el);
	check(pvals.linear_complexity_pvals[0], 0.309412);
	printf("Testing Serial test ........................... ");
	Serial_v2(16, 8 * (int)el);
	check(pvals.serial_pvals[0][0], 0.760793);
	fclose(data);

	// data.sqrt2
	el = 125000;
	data = fopen("data/sqrt2.bin", "rb");
	if (data == NULL)
	{
		fprintf(stderr, "Cannot open 'data/sqrt2.bin' file to load data for testing.\n");
		return;
	}
	r = fread(buffer, 1, MAXBUFFER, data);
	if (r != el)
	{
		fprintf(stderr, "Cannot read from 'data/sqrt2.bin' (file length mismatch).\n");
		fclose(data);
		return;
	}
	printf("\n\nSelftest using the 'sqrt2' binary file:\n\n");
	printf("Testing Frequency test ........................ ");
	Frequency_v2(8 * (int)el);
	check(pvals.frequency_pvals[0], 0.811881);
	printf("Testing Block Frequency test .................. ");
	BlockFrequency_v2(128, 8 * (int)el);
	check(pvals.blockfrequency_pvals[0], 0.833222);
	printf("Testing Cumulative Sum test - Forward ......... ");
	CumulativeSums_v2(8 * (int)el);
	check(pvals.cusum_pvals[0][0], 0.879009);
	printf("Testing Cumulative Sum test - Reverse ......... ");
	check(pvals.cusum_pvals[1][0], 0.957206);
	printf("Testing Runs test ............................. ");
	Runs_v2(8 * (int)el);
	check(pvals.runs_pvals[0], 0.313427);
	printf("Testing Longest Run test ...................... ");
	LongestRunOfOnes_v2(8 * (int)el);
	check(pvals.longestrunofones_pvals[0], 0.012117);
	printf("Testing Rank test ............................. ");
	Rank_v2(8 * (int)el);
	check(pvals.rank_pvals[0], 0.823810);
	printf("Testing Spectral test ......................... ");
	DiscreteFourierTransform_v2(8 * (int)el);
	check(pvals.dft_pvals[0], 0.581909);
	printf("Testing Non-overlapping Templates test ........ ");
	NonOverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.nonoverlapping_pvals[0][0], 0.569461);
	printf("Testing Overlapping Templates test ............ ");
	OverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.overlapping_pvals[0], 0.791982);
	printf("Testing Universal test ........................ ");
	Universal_v2(8 * (int)el);
	check(pvals.universal_pvals[0], 0.130805);
	printf("Testing Approximate Entropy test .............. ");
	ApproximateEntropy_v2(10, 8 * (int)el);
	check(pvals.approximate_entropy_pvals[0], 0.884740);
	printf("Testing Random Excursion test ................. ");
	RandomExcursions_v2(8 * (int)el);
	check(pvals.random_excursion_pvals[4][0], 0.216235);
	printf("Testing Random Excursion Variant test ......... ");
	RandomExcursionsVariant_v2(8 * (int)el);
	check(pvals.random_excursion_variant_pvals[8][0], 0.566118);
	printf("Testing Linear Complexity test ................ ");
	LinearComplexity_v2(500, 8 * (int)el);
	check(pvals.linear_complexity_pvals[0], 0.317127);
	printf("Testing Serial test ........................... ");
	Serial_v2(16, 8 * (int)el);
	check(pvals.serial_pvals[0][0], 0.861925);
	fclose(data);

	// data.sqrt3
	el = 125000;
	data = fopen("data/sqrt3.bin", "rb");
	if (data == NULL)
	{
		fprintf(stderr, "Cannot open 'data/sqrt3.bin' file to load data for testing.\n");
		return;
	}
	r = fread(buffer, 1, MAXBUFFER, data);
	if (r != el)
	{
		fprintf(stderr, "Cannot read from 'data/sqrt3.bin' (file length mismatch).\n");
		fclose(data);
		return;
	}
	printf("\n\nSelftest using the 'sqrt3' binary file:\n\n");
	printf("Testing Frequency test ........................ ");
	Frequency_v2(8 * (int)el);
	check(pvals.frequency_pvals[0], 0.610051);
	printf("Testing Block Frequency test .................. ");
	BlockFrequency_v2(128, 8 * (int)el);
	check(pvals.blockfrequency_pvals[0], 0.473961);
	printf("Testing Cumulative Sum test - Forward ......... ");
	CumulativeSums_v2(8 * (int)el);
	check(pvals.cusum_pvals[0][0], 0.917121);
	printf("Testing Cumulative Sum test - Reverse ......... ");
	check(pvals.cusum_pvals[1][0], 0.689519);
	printf("Testing Runs test ............................. ");
	Runs_v2(8 * (int)el);
	check(pvals.runs_pvals[0], 0.261123);
	printf("Testing Longest Run test ...................... ");
	LongestRunOfOnes_v2(8 * (int)el);
	check(pvals.longestrunofones_pvals[0], 0.446726);
	printf("Testing Rank test ............................. ");
	Rank_v2(8 * (int)el);
	check(pvals.rank_pvals[0], 0.314498);
	printf("Testing Spectral test ......................... ");
	DiscreteFourierTransform_v2(8 * (int)el);
	check(pvals.dft_pvals[0], 0.776046);
	printf("Testing Non-overlapping Templates test ........ ");
	NonOverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.nonoverlapping_pvals[0][0], 0.532235);
	printf("Testing Overlapping Templates test ............ ");
	OverlappingTemplateMatchings_v2(9, 8 * (int)el);
	check(pvals.overlapping_pvals[0], 0.082716);
	printf("Testing Universal test ........................ ");
	Universal_v2(8 * (int)el);
	check(pvals.universal_pvals[0], 0.165981);
	printf("Testing Approximate Entropy test .............. ");
	ApproximateEntropy_v2(10, 8 * (int)el);
	check(pvals.approximate_entropy_pvals[0], 0.180481);
	printf("Testing Random Excursion test ................. ");
	RandomExcursions_v2(8 * (int)el);
	check(pvals.random_excursion_pvals[4][0], 0.783283);
	printf("Testing Random Excursion Variant test ......... ");
	RandomExcursionsVariant_v2(8 * (int)el);
	check(pvals.random_excursion_variant_pvals[8][0], 0.155066);
	printf("Testing Linear Complexity test ................ ");
	LinearComplexity_v2(500, 8 * (int)el);
	check(pvals.linear_complexity_pvals[0], 0.346469);
	printf("Testing Serial test ........................... ");
	Serial_v2(16, 8 * (int)el);
	check(pvals.serial_pvals[0][0], 0.157500);
	fclose(data);

	printf("\n\nFinal result: ");
	if (fails) printf("FAIL");
	else printf("PASS");
	printf("\n\n");

	free(buffer);
	array = NULL;
	free(pvals.frequency_pvals);
	free(pvals.blockfrequency_pvals);
	free(pvals.cusum_pvals[0]);
	free(pvals.cusum_pvals[1]);
	free(pvals.runs_pvals);
	free(pvals.longestrunofones_pvals);
	free(pvals.rank_pvals);
	free(pvals.dft_pvals);
	for (i = 0; i < MAXNUMOFTEMPLATES; i++)
		free(pvals.nonoverlapping_pvals[i]);
	free(pvals.overlapping_pvals);
	free(pvals.universal_pvals);
	free(pvals.approximate_entropy_pvals);
	for (i = 0; i < 8; i++)
		free(pvals.random_excursion_pvals[i]);
	for (i = 0; i < 18; i++)
		free(pvals.random_excursion_variant_pvals[i]);
	free(pvals.serial_pvals[0]);
	free(pvals.serial_pvals[1]);
	free(pvals.linear_complexity_pvals);
}
#endif
