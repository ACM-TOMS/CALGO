/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
              U T I L I T Y  F U N C T I O N  P R O T O T Y P E S 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int		displayGeneratorOptions();
int		generatorOptions(char** streamFile);
void	chooseTests();
void	fixParameters();
void	fileBasedBitStreams(char *streamFile);
void	readBinaryDigitsInASCIIFormat(FILE *fp, char *streamFile);
void	readHexDigitsInBinaryFormat(FILE *fp);
int		convertToBits(BYTE *x, int xBitLength, int bitsNeeded, int *num_0s, int *num_1s, int *bitsRead);
void	openOutputStreams(int option);
void	invokeTestSuite(int option, char *streamFile);
void	nist_test_suite();

// new added funtions for faster output processing
#ifdef KS
void freeMemory();
//http://www.jstatsoft.org/v08/i18/paper?ev=pub_ext_btn_xdl
//code for Kolomogorov-Smirnov p-value

void mMultiply(double *A, double *B, double *C, int m);
void mPower(double *A, int eA, double *V, int *eV, int m, int n);
double K(int n, double d);
#endif
