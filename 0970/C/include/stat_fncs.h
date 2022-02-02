
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     S T A T I S T I C A L  T E S T  F U N C T I O N  P R O T O T Y P E S 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#define	Frequency_v1 Frequency
#define	Frequency_v2 Frequency4  //2,3,4
#define BlockFrequency_v1 BlockFrequency
#define BlockFrequency_v2 BlockFrequency4 //2,3,4
#define	CumulativeSums_v1 CumulativeSums
#define	CumulativeSums_v2 CumulativeSums3 // 2,3
#define	Runs_v1 Runs
#define	Runs_v2 Runs4 //2,3,4
#define	LongestRunOfOnes_v1 LongestRunOfOnes
#define	LongestRunOfOnes_v2 LongestRunOfOnes3 //2,3
#define	Rank_v1 Rank
#define	Rank_v2 Rank2 //2
#define	DiscreteFourierTransform_v1 DiscreteFourierTransform
#define	DiscreteFourierTransform_v2 DiscreteFourierTransform2 //2,3,4
#define	NonOverlappingTemplateMatchings_v1 NonOverlappingTemplateMatchings
#define	NonOverlappingTemplateMatchings_v2 NonOverlappingTemplateMatchings4 //2,4
#define	OverlappingTemplateMatchings_v1 OverlappingTemplateMatchings
#define	OverlappingTemplateMatchings_v2 OverlappingTemplateMatchings4 //2,3,4
#define	Universal_v1 Universal
#define	Universal_v2 Universal2 //2
#define	ApproximateEntropy_v1 ApproximateEntropy
#define	ApproximateEntropy_v2 ApproximateEntropy4 // 2,4
#define	RandomExcursions_v1 RandomExcursions
#define	RandomExcursions_v2 RandomExcursions2 // 2
#define	RandomExcursionsVariant_v1 RandomExcursionsVariant
#define	RandomExcursionsVariant_v2 RandomExcursionsVariant2 //2
#define	LinearComplexity_v1 LinearComplexity
#define	LinearComplexity_v2 LinearComplexity3 //2,3
#define	Serial_v1 Serial
#define	Serial_v2 Serial4 //2,4

//original functions

void	Frequency(int n);
void	BlockFrequency(int M, int n);
void	CumulativeSums(int n);
void	Runs(int n);
void	LongestRunOfOnes(int n);
void	Rank(int n);
void	DiscreteFourierTransform(int n);
void	NonOverlappingTemplateMatchings(int m, int n);
void	OverlappingTemplateMatchings(int m, int n);
void	Universal(int n);
void	ApproximateEntropy(int m, int n);
void	RandomExcursions(int n);
void	RandomExcursionsVariant(int n);
void	LinearComplexity(int M, int n);
void	Serial(int m, int n);

// New functions

//versions published in SPACE 2014 - "Faster randomness testing with the NIST STS"
void	Frequency2(int n);
void	BlockFrequency2(int M, int n);
void	CumulativeSums2(int n);
void	Runs2(int n);
void	LongestRunOfOnes2(int n);
void	Rank2(int n);
void	DiscreteFourierTransform2(int n);
void	NonOverlappingTemplateMatchings2(int m, int n);
void	OverlappingTemplateMatchings2(int m, int n);
void	Universal2(int n);
void	ApproximateEntropy2(int m, int n);
void	RandomExcursions2(int n);
void	RandomExcursionsVariant2(int n);
void	LinearComplexity2(int M, int n);
void	Serial2(int m, int n);
// Journal 
void	Frequency3(int n);
void	BlockFrequency3(int M, int n);
void	CumulativeSums3(int n);
void	Runs3(int n);
void	LongestRunOfOnes3(int n);
void	DiscreteFourierTransform3(int n);
void	OverlappingTemplateMatchings3(int m, int n);

void	LinearComplexity3(int M, int n);

//version 4 and more
void	DiscreteFourierTransform4(int n);

void	BlockFrequency4(int M, int n);  
void	Runs4();
void	Serial4(int m, int n);
void	Frequency4(int n); 
void	ApproximateEntropy4(int m, int n);

void	OverlappingTemplateMatchings4(int m, int n); 
void    NonOverlappingTemplateMatchings4(int m, int n); 





/*

//Frequency 1
void Frequency(int n);
void Frequency2(int n);
void Frequency3(int n);
void Frequency4(int n);
//BlockFrequency 2
void BlockFrequency(int M, int n);
void BlockFrequency2(int M, int n);
void BlockFrequency3(int M, int n);
void BlockFrequency4(int M, int n);
//Runs 4
void Runs(int n);
void Runs2(int n);
void Runs3(int n);
void Runs4(int n);
// LongestRunOfOnes
void LongestRunOfOnes(int n);
void LongestRunOfOnes2(int n);
void LongestRunOfOnes3(int n);
//Rank
void Rank(int n);
void Rank2(int n);
//Serial 14
void Serial(int n, int m);
void Serial2(int n, int m);
void Serial4(int n, int m); //OK correct
//NonOverlappingTemplateMatchings 8
void NonOverlappingTemplateMatchings(int m, int n);
void NonOverlappingTemplateMatchings2(int m, int n);
void NonOverlappingTemplateMatchings4(int m, int n); //OK correct
// OverlappingTemplateMatchings
void OverlappingTemplateMatchings(int m, int n);
void OverlappingTemplateMatchings2(int m, int n);
void OverlappingTemplateMatchings3(int m, int n);
void OverlappingTemplateMatchings4(int m, int n);
//Universal
void Universal(int n);
void Universal2(int n);
//ApproximateEntropy
void ApproximateEntropy4(int m, int n);
void ApproximateEntropy2(int m, int n);
void ApproximateEntropy4(int m, int n);
// CumulativeSums
void CumulativeSums(int n)
void CumulativeSums2(int n)
void CumulativeSums3(int n)
// RandomExcursions
void RandomExcursions(int n);
void RandomExcursions2(int n);
//RandomExcursionsVariant
void RandomExcursionsVariant(int n);
void RandomExcursionsVariant2(int n);
//LinearComplexity
void LinearComplexity(int M, int n);
void LinearComplexity2(int M, int n);
void LinearComplexity3(int M, int n);
// DiscreteFourierTransform
void DiscreteFourierTransform(int n);
void DiscreteFourierTransform2(int n);
void DiscreteFourierTransform3(int n);
void DiscreteFourierTransform4(int n);
















*/