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
        return 0; //throw invalid_argument("invalid argument in choose");
    }
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d, --n){
        unsigned long long g = gcd(r, d);
        r /= g;
        unsigned long long t = n / (d / g);
        if (r > numeric_limits<unsigned long long>::max() / t){
           throw overflow_error("overflow in choose");
        }
        r *= t;
    }
    return r;
}

class CIdxComb{
public:
    // Constructor
    CIdxComb(){
        Init( 2, 1 );
    };

    CIdxComb( unsigned int SetSize, unsigned int CombSize ){
        Init( SetSize, CombSize ); 
    };
    
    // Destructor
    ~CIdxComb() {};

    void Init( unsigned int SetSize, unsigned int CombSize );

    bool SetSizes( unsigned int SetSize, unsigned int CombSize );

    bool GetNextComb( std::vector<unsigned int> &vi );

protected:
    unsigned int m_ArrSize;
    unsigned int m_LastIdx;
    
    unsigned int m_SetSize;
    unsigned int m_LastSetIdx;

};

void CIdxComb::Init( unsigned int SetSize, unsigned int CombSize  ){
    // Assign CombSize
    ////////////////////////
    if( CombSize == 0 ){CombSize = 1;}

    m_ArrSize = CombSize;
    m_LastIdx = CombSize - 1;

    // Assign SetSize
    ////////////////////////
    if( SetSize == 0 ){SetSize = 2;}

    if( CombSize > SetSize ){CombSize = SetSize;}
    
    m_SetSize = SetSize;
    m_LastSetIdx = SetSize - 1;
}


bool CIdxComb::SetSizes( unsigned int SetSize, unsigned int CombSize ){
    if( SetSize == 0 ){return false;}

    if( CombSize == 0 ){return false;}

    if( CombSize > SetSize ){return false;}
    
    m_ArrSize = CombSize;
    m_LastIdx = CombSize - 1;

    m_SetSize = SetSize;
    m_LastSetIdx = SetSize - 1;

    return true;
}


bool CIdxComb::GetNextComb( std::vector<unsigned int> &vi ){
    // Check if the last element is at the end
    if( vi[m_LastIdx] == m_LastSetIdx ){
        if( m_ArrSize == 1 ){return false;}

        // Check if the subsequent elements(counted from back)
        // is also at their subsequent positions
        //////////////////////////////////////////////////////
        bool Completed = true;
        // Incomplete Index, init value not used
        unsigned int IncompIdx = m_LastIdx - 1; 
        
        bool FirstIdx = false;
        unsigned int ArrIdx = m_LastIdx - 1;

        unsigned int SetIdx = m_LastSetIdx - 1;
        
        while( !FirstIdx ){
            if( vi[ArrIdx] != SetIdx ){
                Completed = false;
                IncompIdx = vi[ArrIdx] + 1;
                break;
            }

            if( SetIdx ){--SetIdx;}

            if( !ArrIdx ){
                FirstIdx = true;
            }else{
                --ArrIdx; 
            }
        }

        if( Completed ){
            return false;
        }else{
            for( unsigned int i=ArrIdx; i<=m_LastIdx; ++i, ++IncompIdx )            {
                vi[i] = IncompIdx;
            }
        }
    }else if ( vi[m_LastIdx] < m_LastSetIdx ){
        (vi[m_LastIdx])++;
    }else{
        return false;
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
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithNoRep:invalidNumInputs",
                "Two input arguments required.");
    }
    if (nlhs > 2){
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithNoRep:maxlhs",
                "Too many output arguments.");
    }

    /* Check data type of input argument  */
    if ( !(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1])) ){
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithNoRep:inputNotDouble",
                "Input argument must be of type double.");
    }

    /* Get the size and pointers to input data */
    mwSize m0  = mxGetM(prhs[0]);
    mwSize n0  = mxGetN(prhs[0]);
    mwSize m1  = mxGetM(prhs[1]);
    mwSize n1  = mxGetN(prhs[1]);
    if ( !(m0==1 && n0==1 && m1==1 && n1==1) ){
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithNoRep:inputNotDouble",
                "Input arguments must be scalar.");        
    }

    /* Get Size */
    const mwSize SET = (mwSize)*mxGetPr(prhs[0]);
    const mwSize COMB = (mwSize)*mxGetPr(prhs[1]);
    if (COMB <= 0) {
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithNoRep:inputNotDouble",
                "Second argument must be a positive integer.");  
    }
    if (SET < COMB) {
        mexErrMsgIdAndTxt( "MATLAB:enumCombWithNoRep:inputNotDouble",
                "Second argument must be smaller than first argument.");  
    }
            

    /* Allocate Output Matrix */
    const mwSize numComb = choose(SET,COMB);
    plhs[0] = mxCreateNumericMatrix(COMB, numComb, mxDOUBLE_CLASS, mxREAL);
    double *pl = mxGetPr(plhs[0]);

    // Initialize the first combination vector
    std::vector<unsigned int> vi( COMB );
    for( mwSize j = 0; j < COMB; ++j ){
        vi[j] = j;
    }
    
    CIdxComb cb;

    cb.SetSizes(SET,COMB);

    // Display the first combination
    int head = 0;
    do{
        for( mwSize j = 0; j < COMB; ++j){
            pl[head+j] = vi[j];
        }
        head += COMB;
    } while ( cb.GetNextComb( vi ) );
}