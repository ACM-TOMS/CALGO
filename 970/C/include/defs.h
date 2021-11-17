/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                       D E B U G G I N G  A I D E S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "config.h"
#include <time.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                              M A C R O S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define MAX(x,y)             ((x) <  (y)  ? (y)  : (x))
#define MIN(x,y)             ((x) >  (y)  ? (y)  : (x))
#define isNonPositive(x)     ((x) <= 0.e0 ?   1  : 0)
#define isPositive(x)        ((x) >  0.e0 ?   1 : 0)
#define isNegative(x)        ((x) <  0.e0 ?   1 : 0)
#define isGreaterThanOne(x)  ((x) >  1.e0 ?   1 : 0)
#define isZero(x)            ((x) == 0.e0 ?   1 : 0)
#define isOne(x)             ((x) == 1.e0 ?   1 : 0)

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
                         G L O B A L  C O N S T A N T S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ALPHA							0.01	/* SIGNIFICANCE LEVEL */
#define MAXNUMOFTEMPLATES				148		/* APERIODIC TEMPLATES: 148=>temp_length=9 */
#define NUMOFTESTS						15		/* MAX TESTS DEFINED  */
#define NUMOFGENERATORS					10		/* MAX PRNGs */
#define MAXFILESPERMITTEDFORPARTITION	148
#define	TEST_FREQUENCY					1
#define	TEST_BLOCK_FREQUENCY			2
#define	TEST_CUSUM						3
#define	TEST_RUNS						4
#define	TEST_LONGEST_RUN				5
#define	TEST_RANK						6
#define	TEST_FFT						7
#define	TEST_NONPERIODIC				8
#define	TEST_OVERLAPPING				9
#define	TEST_UNIVERSAL					10
#define	TEST_APEN						11
#define	TEST_RND_EXCURSION				12
#define	TEST_RND_EXCURSION_VAR			13
#define	TEST_SERIAL						14
#define	TEST_LINEARCOMPLEXITY			15


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                   G L O B A L   D A T A  S T R U C T U R E S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef unsigned char	BitSequence;

typedef struct _testParameters {
	int		n;
	int		blockFrequencyBlockLength;
	int		nonOverlappingTemplateBlockLength;
	int		overlappingTemplateBlockLength;
	int		serialBlockLength;
	int		linearComplexitySequenceLength;
	int		approximateEntropyBlockLength;
	int		numOfBitStreams;
// New Flag
	int		fast;
} TP;

// New structure for results comparison
#ifdef VERIFY_RESULTS 
struct results{
	struct frequency_str{double sum, sum_n, p_value;} frequency;
	struct blockfrequency_str{double chi_squared, p_value;} blockfrequency;
	struct runs_str {double pi, V, erfc_arg, p_value;} runs;
	struct longestrunofones_str {int N,M; double chi2, p_value; unsigned int	nu[7];} longestrunofones;
	struct rank_str {double p_30,p_31,p_32, F_30, F_31, F_32, N, chi_squared, p_value;} rank;
	struct serial_str {double	p_value1, p_value2, psim0, psim1, psim2, del1, del2;} serial;
	struct nonoverlapping_str {unsigned int templates; unsigned int *W; double *chi2; double *p_value;} nonoverlapping;
	struct overlapping_str {unsigned int nu[6]; double chi2, p_value;} overlapping;
	struct universal_str { double p_value, sum, phi; } universal;
	struct approximate_entropy_str {double ApEn[2], chi_squared, p_value; unsigned int *P,pp;  } approximate_entropy;
	struct cusum_str {int z, zrev; double sum1A, sum2A, sum1B, sum2B, p_valueA, p_valueB;} cusum;
	struct random_excursion_str {double p_value[8], sum[8]; int valid,x[8],J[8]; } random_excursion;
	struct random_excursion_variant_str { double p_value[18]; int valid,x[18],count[18]; } random_excursion_variant;
	struct linear_complexity_str { double p_value, chi2; int nu[7]; } linear_complexity;
	struct dft_str {double	p_value, percentile, N_l, N_o, d;} dft;
};
#endif

// New structure for faster postprocessing
#ifdef KS 
typedef struct _Pvalues {
	double *frequency_pvals;
	double *blockfrequency_pvals;
	double *cusum_pvals[2];
	double *runs_pvals;
	double *longestrunofones_pvals;
	double *rank_pvals;
	double *dft_pvals;
	double *nonoverlapping_pvals[MAXNUMOFTEMPLATES];
	double *overlapping_pvals;
	double *universal_pvals;
	double *approximate_entropy_pvals;
	double *random_excursion_pvals[8];
	double *random_excursion_variant_pvals[18];
	double *serial_pvals[2];
	double *linear_complexity_pvals;
	

	int seq_counter;
	int num_Nonoverlap_pvals;
	FILE* results;
} Pvals;
#endif

// New structure for use of command line arguments 
#if defined(FILE_OUTPUT) ||  defined(KS)
typedef struct _CmdArgFlags {

	int output;
	int bitStreams;
	int fixArguments;
	int tests;
	int fileFormat;
	int fileGen;
	int argCounter;
} CmdFlags;
#endif
