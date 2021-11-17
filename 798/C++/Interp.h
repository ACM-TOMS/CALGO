/***********************************************************************
                             CLASS Interpolator
This class contains all of the variables and functions related to the
interpolation process.  This code was directly translated from Renka's
Fortran code.
************************************************************************/
 
#include <stdio.h>
#include <iostream.h>

//#define NMAX 10700
#define NMAX 26000                //maximum number of known nodes 
#define NRMAX 6                   //maximum number number of cells = NR     
//#define LLC (NRMAX*NRMAX)
#define MIN0(a,b) a < b ? a : b   //min. macro
#define MAX0(a,b) a > b ? a : b   //max. macro
#define NFUN 5                    //number of test functions


class Interpolator {
  protected:
         int M,                 //Dimensions of hypervolume
             N,                 //Number of known nodes
            NT,                 //Number of terms in quadratic equation
                                //  with M variables: NT=(M+1)(M+2)/2
            NQ,                 //Number of data points to be used in 
                                //  the least squares fit for 
                                //  coefficients defining the nodal fxns.
                                //  Q(K).  
                                //  NT-1=M(M+3)/2 <= NQ <= MIN0(50,N-1) 
            NW,                 //Number of nodes within (and defining)
                                //  the radii of influence R(K) which enter
                                //  into the weights W(K).
                                //  1 <= NW <= MIN0(50,N-1)
            NR,                 //Number of intervals along each 
                                //  coordinate axis defining the cell grid
                                //  described in function STOREM. 
                                //  Recommended value is nearest integer to
                                //  (N/3)**(1/M) with NR >=1. 
          **IW,                 //Integer work space array of length >= 5*M.
   LCELL[NMAX+1],               //Array of length NR**M containing nodal 
                                //  indices associated with cells. See STOREM.
   LNEXT[NMAX+1];               //Array of length N containing next-node
                                //  indices.  See STOREM.
float      *WS,                 //Work space array of length >= NT**2 for
                                //  NT=(M+1)(M+2)/2.
         *XMIN, *DX,            //Arrays of length M containing minimum nodal
                                //  coordinates and cell dimensions,
                                //  respectively.  See STOREM.
           **A,                 //N by NT-1 array containing the coefficients
                                //  for quadratic nodal function Q(K) in row 
                                //  K.
           **X,                 //M by N array whose Kth column contains the
                                //  cartesian coordinates of node K.  The
                                //  actual row dimension in the calling
                                //  program must be M.
             Q,                 //Weighted sum of quadratic nodal fxns. 
                                //  defined in QSHEPM.
          RMAX,                 //Square root of the largest element in RSQ - 
                                //  maximum radius R(K).
       W[NMAX+1],               //stores values at known nodes
     RSQ[NMAX+1];               //Array of length N containing the squares of
                                //  the radii R(K) which enter into the
                                //  weights W(K).
//            DQ;                 


  public:
        Interpolator (int);
       ~Interpolator ();

        float **QSMTST3(int*, float[][3], char);
        void QSHEPM(int&);
        double QSMVAL(float*);
        void QSMGRD();
        void GETNPM(int,int&,double&);
        void GIVENS(float*,float*,double&,double&);
        void ROTATE(int,double,double,int,int);
        void SETUPM(int,double,int,double,double,double,int);
        void STOREM(int&);
        void get_node_coords(float*,int);
        float **getX(){ return X;};     // returns array of coordinates
        float *getW(){ return W;};      // returns array of weights
        int getN(){ return N; };        // returns number of known nodes
        void TFUN3(int,int,float*,int,int,float*,float,float,float);

 };
