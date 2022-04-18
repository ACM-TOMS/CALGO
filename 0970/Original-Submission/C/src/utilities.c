/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
U T I L I T I E S
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "../include/config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/generators.h"
#include "../include/stat_fncs.h"

int
displayGeneratorOptions()
{
	int		option = 0;

	printf("           G E N E R A T O R    S E L E C T I O N \n");
	printf("           ______________________________________\n\n");
	printf("    [0] Input File                 [1] Linear Congruential\n");
	printf("    [2] Quadratic Congruential I   [3] Quadratic Congruential II\n");
	printf("    [4] Cubic Congruential         [5] XOR\n");
	printf("    [6] Modular Exponentiation     [7] Blum-Blum-Shub\n");
	printf("    [8] Micali-Schnorr             [9] G Using SHA-1\n\n");
	printf("   Enter Choice: ");
	scanf("%d", &option);
	printf("\n\n");

	return option;
}


int
generatorOptions(char** streamFile)
{
	char	file[200];
	int		option = NUMOFGENERATORS+1;
	FILE	*fp;
	
	while ( (option < 0) || (option > NUMOFGENERATORS) ) {
		option = displayGeneratorOptions();
		switch( option ) {
			case 0:
				printf("\t\tUser Prescribed Input File: ");
				scanf("%s", file);
				*streamFile = (char*)calloc(200, sizeof(char));
				if(streamFile==NULL) { printf("Cannot allocate memory.\n"); exit(1); }
				sprintf(*streamFile, "%s", file);
				printf("\n");
				if ( (fp = fopen(*streamFile, "r")) == NULL ) {
					printf("File Error:  file %s could not be opened.\n",  *streamFile);
					exit(-1);
				}
				else
					fclose(fp);
				break;
			case 1:
				*streamFile = "Linear-Congruential";
				break;
			case 2:
				*streamFile = "Quadratic-Congruential-1";
				break;
			case 3:
				*streamFile = "Quadratic-Congruential-2";
				break;
			case 4:
				*streamFile = "Cubic-Congruential";
				break;
			case 5:
				*streamFile = "XOR";
				break;
			case 6:
				*streamFile = "Modular-Exponentiation";
				break;
			case 7:
				*streamFile = "Blum-Blum-Shub";
				break;
			case 8:
				*streamFile = "Micali-Schnorr";
				break;
			case 9:
				*streamFile = "G using SHA-1";
				break;
				
			/* INTRODUCE NEW PRNG NAMES HERE */
			/*
			case 10:  *streamFile = "myNewPRNG";
				break;
			*/
			default:
				printf("Error:  Out of range - Try again!\n");
				break;
		}
	}
	return option;
}


void
chooseTests()
{
	int		i;
	
	printf("                S T A T I S T I C A L   T E S T S\n");
	printf("                _________________________________\n\n");
	printf("    [01] Frequency                       [02] Block Frequency\n");
	printf("    [03] Cumulative Sums                 [04] Runs\n");
	printf("    [05] Longest Run of Ones             [06] Rank\n");
	printf("    [07] Discrete Fourier Transform      [08] Nonperiodic Template Matchings\n");
	printf("    [09] Overlapping Template Matchings  [10] Universal Statistical\n");
	printf("    [11] Approximate Entropy             [12] Random Excursions\n");
	printf("    [13] Random Excursions Variant       [14] Serial\n");
	printf("    [15] Linear Complexity\n\n");
	printf("         INSTRUCTIONS\n");
	printf("            Enter 0 if you DO NOT want to apply all of the\n");
	printf("            statistical tests to each sequence and 1 if you DO.\n\n");
	printf("   Enter Choice: ");
	scanf("%d", &testVector[0]);
	printf("\n");
	if ( testVector[0] == 1 )
		for( i=1; i<=NUMOFTESTS; i++ )
			testVector[i] = 1;
	else {
		printf("         INSTRUCTIONS\n");
		printf("            Enter a 0 or 1 to indicate whether or not the numbered statistical\n");
		printf("            test should be applied to each sequence.\n\n");
		printf("      123456789111111\n");
		printf("               012345\n");
		printf("      ");
		for ( i=1; i<=NUMOFTESTS; i++ ) 
			scanf("%1d", &testVector[i]);
		printf("\n\n");
	}
}


void
fixParameters()
{

	int		counter, testid;
	
	//  Check to see if any parameterized tests are selected
	if ( (testVector[TEST_BLOCK_FREQUENCY] != 1) && (testVector[TEST_NONPERIODIC] != 1) && 
		 (testVector[TEST_OVERLAPPING] != 1) && (testVector[TEST_APEN] != 1) &&
		 (testVector[TEST_SERIAL] != 1) && (testVector[TEST_LINEARCOMPLEXITY] != 1) )
			return;
		
	do {
		counter = 1;
		printf("        P a r a m e t e r   A d j u s t m e n t s\n");
		printf("        -----------------------------------------\n");
		if ( testVector[TEST_BLOCK_FREQUENCY] == 1 )
			printf("    [%d] Block Frequency Test - block length(M):         %d\n", counter++, tp.blockFrequencyBlockLength);
		if ( testVector[TEST_NONPERIODIC] == 1 )
			printf("    [%d] NonOverlapping Template Test - block length(m): %d\n", counter++, tp.nonOverlappingTemplateBlockLength);
		if ( testVector[TEST_OVERLAPPING] == 1 )
			printf("    [%d] Overlapping Template Test - block length(m):    %d\n", counter++, tp.overlappingTemplateBlockLength);
		if ( testVector[TEST_APEN] == 1 )
			printf("    [%d] Approximate Entropy Test - block length(m):     %d\n", counter++, tp.approximateEntropyBlockLength);
		if ( testVector[TEST_SERIAL] == 1 )
			printf("    [%d] Serial Test - block length(m):                  %d\n", counter++, tp.serialBlockLength);
		if ( testVector[TEST_LINEARCOMPLEXITY] == 1 )
			printf("    [%d] Linear Complexity Test - block length(M):       %d\n", counter++, tp.linearComplexitySequenceLength);
		printf("\n");
		printf("   Select Test (0 to continue): ");
		scanf("%1d", &testid);
		printf("\n");
		
		counter = 0;
		if ( testVector[TEST_BLOCK_FREQUENCY] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Block Frequency Test block length: ");
				scanf("%d", &tp.blockFrequencyBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_NONPERIODIC] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter NonOverlapping Template Test block Length: ");
				scanf("%d", &tp.nonOverlappingTemplateBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_OVERLAPPING] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Overlapping Template Test block Length: ");
				scanf("%d", &tp.overlappingTemplateBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_APEN] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Approximate Entropy Test block Length: ");
				scanf("%d", &tp.approximateEntropyBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_SERIAL] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Serial Test block Length: ");
				scanf("%d", &tp.serialBlockLength);
				printf("\n");
				continue;
			}
		}
		if ( testVector[TEST_LINEARCOMPLEXITY] == 1 ) {
			counter++;
			if ( counter == testid ) {
				printf("   Enter Linear Complexity Test block Length: ");
				scanf("%d", &tp.linearComplexitySequenceLength);
				printf("\n");
				continue;
			}
		}
	} while ( testid != 0 );
}

#if defined(FILE_OUTPUT) || defined(KS)
void
fileBasedBitStreams(char *streamFile)
{
	FILE	*fp;
	int		mode=-1;
	if (cmdFlags.fileFormat == 1 || cmdFlags.fileFormat == 0){
		mode = cmdFlags.fileFormat;
	}
	else {
		printf("   Input File Format:\n");
		printf("    [0] ASCII - A sequence of ASCII 0's and 1's\n");
		printf("    [1] Binary - Each byte in data file contains 8 bits of data\n\n");
		printf("   Select input mode:  ");
		scanf("%1d", &mode);
		printf("\n");
	}
	
	if ( mode == 0 ) {
		if ( (fp = fopen(streamFile, "r")) == NULL ) {
			printf("ERROR IN FUNCTION fileBasedBitStreams:  file %s could not be opened.\n",  streamFile);
			exit(-1);
		}
		readBinaryDigitsInASCIIFormat(fp, streamFile);
		fclose(fp);
	}
	else if ( mode == 1 ) {
		if ( (fp = fopen(streamFile, "rb")) == NULL ) {
			printf("ERROR IN FUNCTION fileBasedBitStreams:  file %s could not be opened.\n", streamFile);
			exit(-1);
		}
		readHexDigitsInBinaryFormat(fp);
		fclose(fp);
	}
}
#endif


void
readBinaryDigitsInASCIIFormat(FILE *fp, char *streamFile)
{
	int		i, j, num_0s, num_1s, bitsRead, bit;
	
	if ( (epsilon = (BitSequence *) calloc(tp.n, sizeof(BitSequence))) == NULL ) {
		printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
		printf("Statistical Testing Aborted!\n");
		return;
	}
	printf("     Statistical Testing In Progress.........\n\n");   
	for ( i=0; i<tp.numOfBitStreams; i++ ) {
		num_0s = 0;
		num_1s = 0;
		bitsRead = 0;
		for ( j=0; j<tp.n; j++ ) {
			if ( fscanf(fp, "%1d", &bit) == EOF ) {
				printf("ERROR:  Insufficient data in file %s.  %d bits were read.\n", streamFile, bitsRead);
				fclose(fp);
				free(epsilon);
				return;
			}
			else {
				bitsRead++;
				if ( bit == 0 ) 
					num_0s++;
				else 
					num_1s++;
				epsilon[j] = (unsigned char) bit;
			}
		}
		fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
		// TMP
		//convert_epsilon_to_array(tp.n);
		//FILE *f = fopen("xxx", "wb");
		//fwrite(array, 1, 125000, f);
		//fclose(f);
		//
		if(tp.fast)convert_epsilon_to_array(tp.n);
		nist_test_suite();
#ifdef KS
		pvals.seq_counter++;
#endif
		if(tp.fast)free(array);
	}
	free(epsilon);
}

void
readHexDigitsInBinaryFormat(FILE *fp)
{
	int		i, done, num_0s, num_1s, bitsRead, index, ones, bits, j;
	BYTE	buffer[4];
	
	if(tp.fast)
	{
		if ( (array = (unsigned char *) calloc((tp.n>>3)+4,sizeof(unsigned char))) == NULL ) {
			printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
			return;
		}

		printf("     Statistical Testing In Progress.........\n\n");   
		for ( i=0; i<tp.numOfBitStreams; i++ ) {
			num_0s = 0;
			num_1s = 0;
			bitsRead = 0;
			done = 0;
			index =0;
			do {
				if ( fread(buffer, sizeof(unsigned char), 1, fp) != 1 ) {
					printf("READ ERROR:  Insufficient data in file.\n");
					free(array);
					return;
				}
				array[index]=LU_byte_inverted[buffer[0]];
				index++;
				if((index<<3)<=tp.n)
				{
					bitsRead+=8; 
					ones=LU_byte_weight[array[index-1]];
					num_1s+= ones;
					num_0s+=(8-ones);
				}
				else
				{
					bits=tp.n-((index-1)<<3);
					bitsRead+=bits;
					ones=0;
					for (j=0;j<bits;j++)
						if(array[index-1] & (1<<j)) ones++;
					num_1s+=ones;
					num_0s+=(bits-ones);
				}
				if((index<<3)>=tp.n) 
					done = 1;
			} while (!done);
			fread(buffer, sizeof(unsigned char), (index%4==0)?0:4-index%4, fp);
			fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
			//printf("\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
			// TMP
			//FILE *f = fopen("zzz.bin", "wb");
			//fwrite(array, 1, 125000, f);
			//fclose(f);
			// TMP
			nist_test_suite();
#ifdef KS
			pvals.seq_counter++;
#endif
		}
		free(array);
	}
	else
	{
		if ( (epsilon = (BitSequence *) calloc(tp.n,sizeof(BitSequence))) == NULL ) {
			printf("BITSTREAM DEFINITION:  Insufficient memory available.\n");
			return;
		}

		printf("     Statistical Testing In Progress.........\n\n");   
		for ( i=0; i<tp.numOfBitStreams; i++ ) {
			num_0s = 0;
			num_1s = 0;
			bitsRead = 0;
			done = 0;
			do {
				if ( fread(buffer, sizeof(unsigned char), 4, fp) != 4 ) {
					printf("READ ERROR:  Insufficient data in file.\n");
					free(epsilon);
					return;
				}
				done = convertToBits(buffer, 32, tp.n, &num_0s, &num_1s, &bitsRead);
			} while ( !done );
			fprintf(freqfp, "\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
			//printf("\t\tBITSREAD = %d 0s = %d 1s = %d\n", bitsRead, num_0s, num_1s);
			// TMP
			//FILE *f = fopen("zzz.bin", "wb");
			//fwrite(array, 1, 125000, f);
			//fclose(f);
			// TMP
			nist_test_suite();
#ifdef KS
			pvals.seq_counter++;
#endif
		}
		free(epsilon);
	}
}


int
convertToBits(BYTE *x, int xBitLength, int bitsNeeded, int *num_0s, int *num_1s, int *bitsRead)
{
	int		i, j, count, bit;
	BYTE	mask;
	int		zeros, ones;

	count = 0;
	zeros = ones = 0;
	for ( i=0; i<(xBitLength+7)/8; i++ ) {
		mask = 0x80;
		for ( j=0; j<8; j++ ) {
			if ( *(x+i) & mask ) {
				bit = 1;
				(*num_1s)++;
				ones++;
			}
			else {
				bit = 0;
				(*num_0s)++;
				zeros++;
			}
			mask >>= 1;
			epsilon[*bitsRead] = (unsigned char) bit;
			(*bitsRead)++;
			if ( *bitsRead == bitsNeeded )
				return 1;
			if ( ++count == xBitLength )
				return 0;
		}
	}
	
	return 0;
}


#if defined(FILE_OUTPUT) || defined(KS)
void
openOutputStreams(int option)
{
	//hodit ifdef FILE OUPUT
	int		i, numOfBitStreams;

	//#ifdef FILE_OUTPUT
	int		numOfOpenFiles = 0;
	//#endif

	char	freqfn[200], summaryfn[200];

	//#ifdef FILE_OUTPUT
	char	statsDir[200], resultsDir[200];
	//#endif
	
#ifdef KS	
	int		numOfTemplates[100] = { 0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
		2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152 };
#endif
		

	sprintf(freqfn, "experiments/%s/freq.txt", generatorDir[option]);
	if ( (freqfp = fopen(freqfn, "w")) == NULL ) {
		printf("\t\tMAIN:  Could not open freq file: <%s>", freqfn);
		exit(-1);
	}
	sprintf(summaryfn, "experiments/%s/finalAnalysisReport.txt", generatorDir[option]);
	if ( (summary = fopen(summaryfn, "w")) == NULL ) {
		printf("\t\tMAIN:  Could not open stats file: <%s>", summaryfn);
		exit(-1);
	}
	if (cmdFlags.output == 1)
	{
		for (i = 1; i <= NUMOFTESTS; i++) {
			if (testVector[i] == 1) {
				sprintf(statsDir, "experiments/%s/%s/stats.txt", generatorDir[option], testNames[i]);
				sprintf(resultsDir, "experiments/%s/%s/results.txt", generatorDir[option], testNames[i]);
				if ((stats[i] = fopen(statsDir, "w")) == NULL) {	/* STATISTICS LOG */
					printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
					printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
					printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
					exit(-1);
				}
				else
					numOfOpenFiles++;
				if ((results[i] = fopen(resultsDir, "w")) == NULL) {	/* P_VALUES LOG   */
					printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
					printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
					printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
					exit(-1);
				}
				else
					numOfOpenFiles++;
			}
		}
}
	#ifdef FILE_OUTPUT
	for( i=1; i<=NUMOFTESTS; i++ ) {
		if ( testVector[i] == 1 ) {
			sprintf(statsDir, "experiments/%s/%s/stats.txt", generatorDir[option], testNames[i]);
			sprintf(resultsDir, "experiments/%s/%s/results.txt", generatorDir[option], testNames[i]);
			if ( (stats[i] = fopen(statsDir, "w")) == NULL ) {	/* STATISTICS LOG */
				printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
				printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
				printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
				exit(-1);
			}
			else
				numOfOpenFiles++;
			if ( (results[i] = fopen(resultsDir, "w")) == NULL ) {	/* P_VALUES LOG   */
				 printf("ERROR: LOG FILES COULD NOT BE OPENED.\n");
				 printf("       MAX # OF OPENED FILES HAS BEEN REACHED = %d\n", numOfOpenFiles);
				 printf("-OR-   THE OUTPUT DIRECTORY DOES NOT EXIST.\n");
				 exit(-1);
			}
			else
				numOfOpenFiles++;
		}
	}
	#endif
	if (cmdFlags.bitStreams == 0){
		printf("   How many bitstreams? ");
		scanf("%d", &numOfBitStreams);
		tp.numOfBitStreams = numOfBitStreams;
		printf("\n");
	}
	else numOfBitStreams = tp.numOfBitStreams;
	
#ifdef KS
	if(testVector[1]) pvals.frequency_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[2])pvals.blockfrequency_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[3])pvals.cusum_pvals[0] = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[3])pvals.cusum_pvals[1] = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[4])pvals.runs_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[5])pvals.longestrunofones_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[6])pvals.rank_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[7])pvals.dft_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[8])for (i = 0; i < MAXNUMOFTEMPLATES; i++) pvals.nonoverlapping_pvals[i] = (double*)malloc(sizeof(double)*numOfBitStreams);
	if (testVector[9])pvals.overlapping_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[10])pvals.universal_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[11])pvals.approximate_entropy_pvals = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[12])for (i = 0; i < 8; i++) pvals.random_excursion_pvals[i] = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[13])for (i = 0; i < 18; i++) pvals.random_excursion_variant_pvals[i] = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[14])pvals.serial_pvals[0] = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[14])pvals.serial_pvals[1] = malloc(sizeof(double)*numOfBitStreams);
	if (testVector[15])pvals.linear_complexity_pvals = malloc(sizeof(double)*numOfBitStreams);
	
	pvals.seq_counter = 0;
	//***************************
	if (testVector[8])pvals.num_Nonoverlap_pvals = numOfTemplates[tp.nonOverlappingTemplateBlockLength];
	sprintf(summaryfn, "experiments/%s/results.txt", generatorDir[option]);
	//***************************
	if ((pvals.results = fopen(summaryfn, "w")) == NULL) {
		printf("\t\tMAIN:  Could not open stats file: <%s>", summaryfn);
		exit(-1);
	}
#endif
}


void
invokeTestSuite(int option, char *streamFile)
{
	fprintf(freqfp, "________________________________________________________________________________\n\n");
	fprintf(freqfp, "\t\tFILE = %s\t\tALPHA = %6.4f\n", streamFile, ALPHA);
	fprintf(freqfp, "________________________________________________________________________________\n\n");
	if ( option != 0 )
		printf("     Statistical Testing In Progress.........\n\n");
	switch( option ) {
		case 0:
			fileBasedBitStreams(streamFile);
			break;
		case 1:
			lcg();
			break;
		case 2:
			quadRes1();
			break;
		case 3:
			quadRes2();
			break;
		case 4:
			cubicRes();
			break;
		case 5:
			exclusiveOR();
			break;
		case 6:
			modExp();
			break;
		case 7:
			bbs();
			break;
		case 8:
			micali_schnorr();
			break;
		case 9:
			SHA1();
			break;
			
		/* INTRODUCE NEW PSEUDO RANDOM NUMBER GENERATORS HERE */
			
		default:
			printf("Error in invokeTestSuite!\n");
			break;
	}
	printf("     Statistical Testing Complete!!!!!!!!!!!!\n\n");
}
#endif


void
nist_test_suite()
{
	if ((testVector[0] == 1) || (testVector[TEST_FREQUENCY] == 1))
	{
		if (tp.fast) Frequency_v2(tp.n); else Frequency_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_BLOCK_FREQUENCY] == 1))
	{
		if (tp.fast) BlockFrequency_v2(tp.blockFrequencyBlockLength, tp.n); else BlockFrequency_v1(tp.blockFrequencyBlockLength, tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_CUSUM] == 1))
	{
		if (tp.fast) CumulativeSums_v2(tp.n); else CumulativeSums_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_RUNS] == 1))
	{
		if (tp.fast) Runs_v2(tp.n); else Runs_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_LONGEST_RUN] == 1))
	{
		if (tp.fast) LongestRunOfOnes_v2(tp.n); else LongestRunOfOnes_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_RANK] == 1))
	{
		if (tp.fast) Rank_v2(tp.n); else Rank_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_FFT] == 1))
	{
		if (tp.fast) DiscreteFourierTransform_v2(tp.n); else DiscreteFourierTransform_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_NONPERIODIC] == 1))
	{
		if (tp.fast) NonOverlappingTemplateMatchings_v2(tp.nonOverlappingTemplateBlockLength, tp.n);
		else NonOverlappingTemplateMatchings_v1(tp.nonOverlappingTemplateBlockLength, tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_OVERLAPPING] == 1))
	{
		if (tp.fast) OverlappingTemplateMatchings_v2(tp.overlappingTemplateBlockLength, tp.n);
		else OverlappingTemplateMatchings_v1(tp.overlappingTemplateBlockLength, tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_UNIVERSAL] == 1))
	{
		if (tp.fast) Universal_v2(tp.n); else Universal_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_APEN] == 1))
	{
		if (tp.fast) ApproximateEntropy_v2(tp.approximateEntropyBlockLength, tp.n);
		else ApproximateEntropy_v1(tp.approximateEntropyBlockLength, tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_RND_EXCURSION] == 1))
	{
		if (tp.fast) RandomExcursions_v2(tp.n); else RandomExcursions_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_RND_EXCURSION_VAR] == 1))
	{
		if (tp.fast) RandomExcursionsVariant_v2(tp.n); else RandomExcursionsVariant_v1(tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_SERIAL] == 1))
	{
		if (tp.fast) Serial_v2(tp.serialBlockLength, tp.n); else Serial_v1(tp.serialBlockLength, tp.n);
	}
	
	if ((testVector[0] == 1) || (testVector[TEST_LINEARCOMPLEXITY] == 1))
	{
		if (tp.fast) LinearComplexity_v2(tp.linearComplexitySequenceLength, tp.n); else LinearComplexity_v1(tp.linearComplexitySequenceLength, tp.n);
	}
}
#ifdef KS
void freeMemory()
{
	int i;
	if (testVector[1]) free(pvals.frequency_pvals);
	if (testVector[2]) free(pvals.blockfrequency_pvals);
	if (testVector[3]) free(pvals.cusum_pvals[0]);
	if (testVector[3]) free(pvals.cusum_pvals[1]);
	if (testVector[4]) free(pvals.runs_pvals);
	if (testVector[5]) free(pvals.longestrunofones_pvals);
	if (testVector[6]) free(pvals.rank_pvals);
	if (testVector[7]) free(pvals.dft_pvals);
	if (testVector[8])for (i = 0; i < MAXNUMOFTEMPLATES; i++) free(pvals.nonoverlapping_pvals[i]);
	if (testVector[9]) free(pvals.overlapping_pvals);
	if (testVector[10]) free(pvals.universal_pvals);
	if (testVector[11]) free(pvals.approximate_entropy_pvals);
	if (testVector[12])for (i = 0; i < 8; i++) free(pvals.random_excursion_pvals[i]);
	if (testVector[13])for (i = 0; i < 18; i++) free(pvals.random_excursion_variant_pvals[i]);
	if (testVector[14]) free(pvals.serial_pvals[0]);
	if (testVector[14]) free(pvals.serial_pvals[1]);
	if (testVector[15]) free(pvals.linear_complexity_pvals);
}


void mMultiply(double *A, double *B, double *C, int m)
{
	int i, j, k; double s;
	for (i = 0; i<m; i++) for (j = 0; j<m; j++)
	{
		s = 0.; for (k = 0; k<m; k++) s += A[i*m + k] * B[k*m + j]; C[i*m + j] = s;
	}
}
void mPower(double *A, int eA, double *V, int *eV, int m, int n)
{
	double *B; int eB, i;
	if (n == 1) { for (i = 0; i<m*m; i++) V[i] = A[i]; *eV = eA; return; }
	mPower(A, eA, V, eV, m, n / 2);
	B = (double*)malloc((m*m) * sizeof(double));
	mMultiply(V, V, B, m); eB = 2 * (*eV);
	if (n % 2 == 0) { for (i = 0; i<m*m; i++) V[i] = B[i]; *eV = eB; }
	else { mMultiply(A, B, V, m); *eV = eA + eB; }
	if (V[(m / 2)*m + (m / 2)]>1e140) { for (i = 0; i<m*m; i++) V[i] = V[i] * 1e-140; *eV += 140; }
	free(B);
}
double K(int n, double d)
{
	int k, m, i, j, g, eH, eQ;
	double h, s, *H, *Q;

	if (d<0) return 0;

	//OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL
	s = d*d*n; if (s>7.24 || (s>3.76&&n>99)) return 1 - 2 * exp(-(2.000071 + .331 / sqrt(n) + 1.409 / n)*s);
	k = (int)(n*d) + 1; m = 2 * k - 1; h = k - n*d;
	H = (double*)malloc((m*m) * sizeof(double));
	Q = (double*)malloc((m*m) * sizeof(double));
	for (i = 0; i<m; i++) for (j = 0; j<m; j++)
		if (i - j + 1<0) H[i*m + j] = 0; else H[i*m + j] = 1;
	for (i = 0; i<m; i++) { H[i*m] -= pow(h, i + 1); H[(m - 1)*m + i] -= pow(h, (m - i)); }
	H[(m - 1)*m] += (2 * h - 1>0 ? pow(2 * h - 1, m) : 0);
	for (i = 0; i<m; i++) for (j = 0; j<m; j++)
		if (i - j + 1>0) for (g = 1; g <= i - j + 1; g++) H[i*m + j] /= g;
	eH = 0; mPower(H, eH, Q, &eQ, m, n);
	s = Q[(k - 1)*m + k - 1];
	for (i = 1; i <= n; i++) { s = s*i / n; if (s<1e-140) { s *= 1e140; eQ -= 140; } }
	s *= pow(10., eQ); free(H); free(Q); return s;
}
#endif

/*
void  prepare_optimisation(){

	op.num_of_ones = -1;
	op.LUT_HW_size = op.LUT_Cusum_size = op.LUT_Run_size = op.LUT_Lrun_size = op.LUT_Switches_size = 16;
	op.LUT_HW_Bsize = op.LUT_Cusum_Bsize = op.LUT_Run_Bsize = op.LUT_Lrun_Bsize = op.LUT_Switches_Bsize = 2;



	//LUT Hamming weigh
	set_LUT(op.LUT_HW_size, &op.LUT_HW, func_HW);

	//LUT for Runs - overlapping blocks => size +1
	set_LUT(op.LUT_Switches_size + 1, &op.LUT_Switches, func_Runs);

	//LUTs for longest run
	set_LUT(op.LUT_Lrun_size, &op.LUT_Lrun_start, func_LRun_start);
	set_LUT(op.LUT_Lrun_size, &op.LUT_Lrun_max, func_LRun_max);
	set_LUT(op.LUT_Lrun_size, &op.LUT_Lrun_end, func_LRun_end);

	//LUts for Cusum
	set_LUT(op.LUT_Cusum_size, &op.LUT_Cusum_max_positiv, func_Cusum_max_positiv);
	set_LUT(op.LUT_Cusum_size, &op.LUT_Cusum_max_negativ, func_Cusum_max_negativ);
	set_LUT(op.LUT_Cusum_size, &op.LUT_Cusum, func_Cusum);

}
*/
