*************************************************************************
*        Alternative (faster) implementation of the NIST STS tests      *
*                                                                       *
*            NIST STS Version 2.1.2 - Alt Version F (v6)                *
*************************************************************************
*                                                                       *
*   Please write your bug reports, questions and comments to:           *
*                        zriha@mail.muni.cz                             *
*                        syso@mail.muni.cz                              *
*                                                                       *
*************************************************************************
*                                                                       *
*                         I N S T A L L                                 *
*                                                                       *
*************************************************************************

Just type 'make' (unix-like systems) or use NIST.sln (Visual Studio).
Result is a single executable file. For the Non-overlapping Template 
Matching Test a subdirectory "template" with the template files must be
available. For saving the results the "experiments" subtree must exist.

If you want to use the FFTW library for the FFT tests then install the
library and in stat_fncs.h define DiscreteFourierTransform_v2 as
DiscreteFourierTransform3 or DiscreteFourierTransform4.

Compilation tested on recent Linux (gcc 4.4.7), cygwin64 (gcc 4.9.3) and 
mingw64 (gcc 5.2.0). Tested on Visual Studio 2015 (also works with previous
versions).

*************************************************************************
*                                                                       *
*                           U S A G E                                   *
*                                                                       *
*************************************************************************

This is an alternative implementation of the NIST statistical randomness
tests. The source codes include both the original and the new alternative
improved variant/variants of tests. You can also compare performance of 
both implementations and correctness of the new implementation. Also 
Kolmogorov-Smirnov test is added to check uniformity (additional test 
to Chi^2 test).

Each test can be used in two versions v1, v2 (for instance Frequency_v1 
and Frequency_v2),  where v1 (Frequency_v1) defines original implementation
of the test and v2 (Frequency_v2) defines a new implementation. 
 
You can add new implementation of the arbitrary test as follows: 
	1.Add the new function (Frequency_new(int n){...}) to file frequency.h.
	2.Add its prototype (Frequency_new(int n);) to stat-fncs.h
	3.Change #define Frequency_v2 Frequency3 
	  to #define Frequency_v2 Frequency_new
 
*************************************************************************
*                                                                       *
*                           Randomness testing                          *
*                                                                       *
*************************************************************************
1. In config.h set the following:
    //#define VERIFY_RESULTS 1
    //#define SPEED 1
     #define FILE_OUTPUT 1 
    //define KS 1 
    if you want to test randomness.
    
    //#define FILE_OUTPUT 1 
    define KS 1 
    if you want also to test uniformity using Kolmogorov-Smirnov test.
    KS can compute p-value greater than 1 since it computes accurate 
    p-values only for small p-values. 
    
 
2. Usage: assess <stream length> [-fast] 
   <stream length> is the length of the individual bit stream(s) to be processed
   -fast           use the faster alternative implementation (version 2) of tests
   (set Command arguments: 10000 -fast if you want to test 10 000 bits with faster version of NIST tests)

*************************************************************************
*                                                                       *
*                           Performance testing                         *
*                                                                       *
*************************************************************************
1. In config.h set the following:
    //#define VERIFY_RESULTS 1
    #define SPEED 1
    //#define FILE_OUTPUT 1
    //define KS 1 
2. Usage: assess scale repeat test_from test_to (for instance 0 10 1 2)
	scale      - 0 (one fixed bit size = 20MB) or 1 (bitsizes are increases in steps)
	repeat     - is number of times( 10) each test is executed (minimum time is taken as the result)
	test_from test_to  - which tests are measured ( 1 - 2 = Freqency and BlockFrequency) 

You can also compare speed of your implementation with our implemenation.
If you want to compare speed of our Frequency2 and your Frequency_your function
it suffices to change defines to:
    1. #define Frequency_v1 Frequency2
    2. #define Frequency_v2 Frequency_your

*************************************************************************
*                                                                       *
*                           Correctness testing                         *
*                                                                       *
*************************************************************************
1. In config.h set the following:
    #define VERIFY_RESULTS 1
    //#define SPEED 1
    //#define FILE_OUTPUT 1
    //define KS 1 
2. Usage: assess test
	test - which test to verify (e.g. 1 = Frequency) 

*************************************************************************
*                                                                       *
*              Introduce support for m>25 (and m<=32)                   *
*              in Serial and Overlapping Template Matching tests        *
*                                                                       *
*************************************************************************

If you feel limited with the maximum value of m (m <=25) in the Serial
and Overlapping Template Matching tests then simply replace the calls
of get_nth_block4() with the calls of get_nth_block_effect() in the 
relevant parts of the source code. For the serial test search for function 
Serial4() in serial.c. For the Overlapping Template Matching test search for
the funtion OverlappingTemplateMatchings2() in overlappingTemplateMatchings.c

 