#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <stdexcept>
#include "mex.h"

using namespace std;

unsigned long long
gcd(unsigned long long x, unsigned long long y){
    while (y != 0){
        unsigned long long t = x % y;
        x = y;
        y = t;
    }
    return x;
}

unsigned long long
choose(unsigned long long n, unsigned long long k){
    if (k > n){
        return 0;//throw invalid_argument("invalid argument in choose");
    }
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d, --n){
        unsigned long long g = gcd(r, d);
        r /= g;
        unsigned long long t = n / (d / g);
        if (r > numeric_limits<unsigned long long>::max() / t)
           throw overflow_error("overflow in choose");
        r *= t;
    }
    return r;
}

bool CombWithRep( unsigned int Set, unsigned int Comb, std::vector<unsigned int> &vi ) {
    if( Set == 0 || Comb == 0 ){
        return false;
    }
    
    bool bEndReach = false;
    for( int x=Comb-1; x>=0 ; --x ){
        if( x == 0 && vi[x] == Set - 1 ){return false;}

        if( bEndReach ){
            if( vi[x] != Set - 1 ){
                unsigned int level = vi[x] + 1;
                for( unsigned int y=x; y<Comb; ++y ){
                    vi[y] = level;
                }
                    
                return true;
            }
        }
        
        // At the end of the Set
        if( vi[x] == Set - 1 ){
            bEndReach = true;
        }else if( vi[x] < Set - 1 ){
            (vi[x])++;
            return true;        
        }
    }
    return true;        
}

void mexFunction(
         int nlhs,       mxArray *plhs[],
         int nrhs, const mxArray *prhs[]
         )
{
    /* Check for proper number of input and output arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithRep:invalidNumInputs",
                "Two input arguments required.");
    }
    if (nlhs > 2){
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithRep:maxlhs",
                "Too many output arguments.");
    }

    /* Check data type of input argument  */
    if ( !(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1])) ){
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithRep:inputNotDouble",
                "Input argument must be of type double.");
    }

    /* Get the size and pointers to input data */
    mwSize m0  = mxGetM(prhs[0]);
    mwSize n0  = mxGetN(prhs[0]);
    mwSize m1  = mxGetM(prhs[1]);
    mwSize n1  = mxGetN(prhs[1]);
    if ( !(m0==1 && n0==1 && m1==1 && n1==1) ){
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithRep:inputNotDouble",
                "Input arguments must be scalar.");        
    }

    /* Get Size */
    const mwSize SET = (mwSize)*mxGetPr(prhs[0]);
    const mwSize COMB = (mwSize)*mxGetPr(prhs[1]);
    if (COMB <= 0) {
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithRep:inputNotDouble",
                "Second argument must be a positive integer.");  
    }

    /* Allocate Output Matrix */
    const mwSize numComb = choose(SET+COMB-1,COMB);
    plhs[0] = mxCreateNumericMatrix(COMB, numComb, mxDOUBLE_CLASS, mxREAL);
    double *pl = mxGetPr(plhs[0]);

    // Initialize the first combination vector to zeros
    std::vector<unsigned int> vi( COMB, 0 );

    // Display the first combination
    int head = 0;
    do{
        for( mwSize j = 0; j < COMB; ++j){
            pl[head+j] = vi[j];
        }
        head += COMB;
    } while ( CombWithRep( SET, COMB, vi ) );
}

