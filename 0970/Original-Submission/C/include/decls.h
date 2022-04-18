
#include "../include/defs.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                   G L O B A L   D A T A  S T R U C T U R E S 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* New structures */
unsigned char *array;
#ifdef VERIFY_RESULTS
struct results R1,R2,R_;
#endif

#ifdef SPEED
double dummy_result = 0.0;
#endif

#if defined(FILE_OUTPUT) || defined(KS)
CmdFlags cmdFlags;
#endif

#ifdef KS
Pvals pvals;
#endif

/* original stuff */
BitSequence	*epsilon;				// BIT STREAM
TP			tp;						// TEST PARAMETER STRUCTURE
FILE		*stats[NUMOFTESTS+1];	// FILE OUTPUT STREAM
FILE		*results[NUMOFTESTS+1];	// FILE OUTPUT STREAM
FILE		*freqfp;				// FILE OUTPUT STREAM
FILE		*summary;				// FILE OUTPUT STREAM
int			testVector[NUMOFTESTS+1];

char	generatorDir[NUMOFGENERATORS][20] = { "AlgorithmTesting", "LCG", "QCG1", "QCG2","CCG", "XOR",
			"MODEXP", "BBS", "MS", "G-SHA1" };
				
char	testNames[NUMOFTESTS+1][32] = { " ", "Frequency", "BlockFrequency", "CumulativeSums", "Runs", "LongestRun", "Rank",
			"FFT", "NonOverlappingTemplate", "OverlappingTemplate", "Universal", "ApproximateEntropy", "RandomExcursions",
			"RandomExcursionsVariant", "Serial", "LinearComplexity" };
