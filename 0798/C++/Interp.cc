/***********************************************************************
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87

                                        Translated into C++ 
                                            by Karen Minser
                                         Univ. of Tennessee
                                                     7/6/98


                             CLASS Interpolator
The following functions contain the C++ code translated from Renka's Fortran
interpolation code based on the Modified Shepard's Method.  Many of the
variables and parameters used and passed around through the Renka's fortran
subroutines have been changed to class variables so a lot of the parameter
passing has been eliminated.
*************************************************************************/
#include <stdio.h>
#include <iostream.h>
#include <std/cmath.h>
#include "Interp.h"
/************************************************************************
                            Interpolator::Interpolator
Constructor function for Interpolator.  Sets up many of the arrays
that were established as global variables in Renka's original QSMTST3
subroutine. 
   Input Parameters:
      DIM      -dimensions of hypervolume, assigned to M
************************************************************************/
Interpolator::Interpolator(int DIM) : M(DIM) {
     int i,j,                // loop indices
//         LIW,          
         LWS;                // size of WS array = NT**2 
//         LA;                 // size of A array 

   NT = (M+1) * (M+2)/2;

//   LIW = 5 * M;
   LWS = NT * NT;
//   LA = NMAX * (NT-1);

//   IW = new int[LIW];

   // define arrays
   IW = new int*[M+1];        // integer work space array >= 5*M
   if(!IW){
     cerr << "IW array not allocated\n";
     exit(1);
    }
   for(i=1;i<M+1;i++){
     IW[i] = new int[6];
     if(!IW[i]){
       cerr << "IW array not allocated\n";
       exit(1);
      }
    }

   WS = new float[LWS+1];     // float work space array >= NT**2
   if(!WS){
     cerr << "WS array not allocated\n";
     exit(1);
    }

   for(i=0;i<LWS+1;i++){
      WS[i] = 0.0;
    }

   XMIN = new float[M+1];     // array of min. nodal coordinates
   if(!XMIN){
     cerr << "XMIN array not allocated\n";
     exit(1);
    }

   DX = new float[M+1];       // array of min. cell dimensions
   if(!DX){
     cerr << "DX array not allocated\n";
     exit(1);
    }

 } // end of constructor function

/************************************************************************
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 
                                 QSMTST3
This function contains the majority of the original fortran function,
QSMTST3, however, parts of the original are found in main and the 
constructor function, Interpolator.  This function finishes the setup for 
the intepolation that takes place in function QSHEPD.  Input from the user
is obtained either from the command line or an input file.  The file
containing the data referring to the known nodes is opened and read. 
The axes and their coordinate values are determined.
  Input Parameters:
    Axes_Dim          -array of integers containing the dimensions of axes. 
    Axes_Ranges       -array of floats containing the range of each axis.
    ans               -character representing user's choice to Test or not.

Upon completion, this function returns an array of floats, Axes_Ticks,
that contain the actual values found at each tick mark along each axis.
**************************************************************************/
float
**Interpolator::QSMTST3(int *Axes_Dim,float Axes_Ranges[][3],char ans) {

   int i,j,                   // loop indices
       NWMAX;                 // max. value for NW, MIN0(50,N-1)

   float step,                // length between tick marks on axis
  **Axes_Ticks;               // array containing tick mark values on axes
   FILE *fin;                 // file containing data of known nodes

   // open data file corresponding to Testing or real Interpolation 
   if(ans == 'N'){           
      if((fin=fopen("DATA.c","r"))==NULL){
        cerr << "Cannot open DATA file\n";
        exit(1);
       }
     }
   else {
      if((fin=fopen("DATA.test","r"))==NULL){
        cerr << "Cannot open DATA file\n";
        exit(1);
       }
    }



   fscanf(fin,"%d",&N);        // get N from datafile - make sure legit
   if((N < NT)||(N > NMAX)){
     cout << "       *** ERROR -- N = " << N << "  MAXIMUM VALUE = "; 
     cout << NMAX << " ***" << endl;
     exit(1);
    }

   // allocate array for cartesian coordinates of known nodes
   X = new float*[M+1];
   if(!X){
     cerr << "X array not allocated\n";
     exit(1);
    }
   for(i=1;i<M+2;i++){
     X[i] = new float[N+1];
     if(!X[i]){
       cerr << "X array not allocated\n";
       exit(1);
      }
    }

   // read in coordinates from datafile of known nodes - store in X
   for(i=1;i<N+1;i++){
     for(j=1;j<M+1;j++){
      fscanf(fin,"%f",&X[j][i]);
//cout << X[j][i] << " " ;
      }
     fscanf(fin,"%f",&W[i]);  // get value of known nodes - store in W
//cout << W[i] << endl;
    }
     
   // allocate array A, size = N*(NT-1)
   A = new float*[N+1];
   if(!A){
     cerr << "A array not allocated\n";
     exit(1);
    }
   for(i=1;i<N+1;i++){
     A[i] = new float[NT];
     if(!A[i]){
       cerr << "A array not allocated\n";
       exit(1);
      }
    }

   for(i=1;i<N+1;i++){
    for(j=1;j<NT;j++){
       A[i][j] = 0;
     }
   }

   NWMAX = MIN0(50,N-1);
   
   cout << "N = " << N << " NWMAX = " << NWMAX << endl;
   cout << "\n\n";

   // obtain values for input parameters for NQ, NW, and NR - check if legit.
   NQ = 0;
   while((NQ < NT-1)||(NQ > NWMAX)) {
     cout << "QSMTST3:  N = " << N << endl;
     cout << "Specify the number of nodes NQ for the least squares fit\n";
     cout << "NQ = 17 is recommended.  9 < NQ < Min(50,N-1) or NQ = 0 to\n";
     cout << "terminate.\n";
     scanf("%d",&NQ);
     cout << "NQ = " << NQ << endl;
     if(NQ==0) exit(0);
   }   

   NW=0;
   while((1 > NW)||(NW > NWMAX)){
     cout << "Specify the number of nodes NW for the weights\n";
     cout << "NW = 32 is recommended.  1 < NW < Min(50,N-1)\n";
     scanf("%d",&NW);
     cout << "NW = " << NW << endl;
    }
   
   NR=0;
   while((NR < 1)||(NR > NRMAX)){
     cout << "Specify the number of rows, columns, and planes NR in the\n";
     cout << "uniform grid of cells used to locate nearest neighbors.\n";
     cout << "NR = (N/3)**(1/M) is recommended. 1 < NR < NRMAX\n";
     scanf("%d",&NR);
     cout << "NR = " << NR << endl;
    }

   //allocate Axes_Ticks array
   Axes_Ticks = new float*[M+1];
   if(!Axes_Ticks){
     cerr << "Axes_Ticks array not allocated\n";
     exit(1);
    }

   for(i=1;i<M+1;i++){
     Axes_Ticks[i] = new float[Axes_Dim[i]+1];
     if(!Axes_Ticks[i]){
       cerr << "Axes_Ticks array not allocated\n";
       exit(1);
      }
     //determine tick values along axis
     Axes_Ticks[i][1] = Axes_Ranges[i][1];
     step = (Axes_Ranges[i][2] - Axes_Ranges[i][1] )/(Axes_Dim[i]-1);
     for(j=2;j<Axes_Dim[i]+1;j++){
       Axes_Ticks[i][j] = Axes_Ticks[i][j-1] + step;
      }
    }

   for(i=1;i<M+1;i++){
    cout << "Axis: " << i << endl;
    for(j=1;j<Axes_Dim[i]+1;j++){
      cout << Axes_Ticks[i][j] << " ";
     }
    cout << endl;
   }

  return Axes_Ticks;
 }  //end of QSMTST3





/*******************************************************************

                         QSHEPM
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 
   This subroutine computes a set of parameters defining a
 smooth (once continuously differentiable) function Q(P)
 which interpolates data values F at scattered nodes X in
 Euclidean M-space (P and X are vectors of length M).  The
 interpolant Q may be evaluated at an arbitrary point by
 function QSMVAL, and its gradient may be computed by sub-
 routine QSMGRD.
   The interpolation scheme is a modified quadratic Shepard
 method:

 Q = (W(1)*Q(1)+W(2)*Q(2)+..+W(N)*Q(N))/(W(1)+W(2)+..+W(N))

 for functions W(K) and Q(K) of M variables.  The nodal
 functions are given by

      Q(K)(P) = SUM[A(K,L)*(P(I)-X(I,K))*(P(J)-X(J,K))] +
                SUM[A(K,MT+I)*(P(I)-X(I,K)] + F(K)

 where the sums are over all I and J such that 1 .LE. I
 .LE. J .LE. M, with L = (J-1)*J/2 + I and MT = M*(M+1)/2.
 Thus, Q(K) is a quadratic function which interpolates the
 data value at node K.  Its coefficients A(K, ) are ordered
 by ordering the upper triangular array of MT quadratic
 terms by columns, followed by the M linear terms.  These
 coefficients are obtained by a weighted least squares fit
 to the closest NQ data points with weights similar to W(K)
 (defined below).  The radius of influence for the least
 squares fit is fixed for each K, but varies with K.
   The weights are taken to be

      W(K)(P) = ( (R(K)-D(K))+ / (R(K)*D(K)) )**2

 where (R(K)-D(K))+ = 0 if R(K) .LE. D(K), and D(K)(P) is
 the Euclidean distance between P and X( ,K).  The radius
 of influence R(K) varies with K and is chosen so that NW
 nodes are within the radius.  Note that W(K) is not de-
 fined at node K, but Q(P) has limit F(K) as P approaches
 X( ,K).
   A cell-based search method, described in subroutines
 STOREM and GETNPM, is used to improve efficiency in both
 the preprocessing phase (QSHEPM) and evaluation phase
 (QSMVAL and QSMGRD).  On a sequential processor, with a
 reasonably uniform distribution of nodes, operation counts
 are O(N) for preprocessing, and constant for each evalua-
 tion.  However, computation time increases exponentially
 with the dimension M.  Also, preprocessing time increases
 linearly with NQ, and evaluation time is proportional to
 NW.

   Input Parameters:
           IER     -Error indicator:
                           IER = 0 if no errors were encountered.
                           IER = 1 if M, N, NQ, NW, or NR is out of range
                                 on input.
                           IER = 2 if duplicate nodes were encountered.
                           IER = 3 if all nodes lie in an affine subspace
                                 of dimension M-1.
   Functions Called:
            GETNPM, GIVENS, ROTATE, SETUPM, STOREM
        


******************************************************************/

void
Interpolator::QSHEPM(int &IER){
  int   I,i,       //loop and array indices
       IB,         // "    "    "      "
     IERR,         // error indicator for STOREM
   II, IJ,         // Indices for WS -- elements [I,I] and [I,J] of
                   //   the augmented regression matrix
     IRM1,         // IROW - 1
     IROW,         // Number of rows currently included in regression matrix
        J,         // Row index for X, column index for A, and loop index
       JJ,         // Index for WS -- element [J,J]
        K,         // Nodal index for F, LNEXT, columns of X, and rows of A
     LMAX,         // Maximum number of NPTS elements (must be consistent with
                   //   the dimension statement above)
      LNP,         // Current length of NTPS
       MT,         // Number of quadratic terms in Q(K)
      NEQ,         // Number of equations in the least squares fit
    NI,NJ,         // Indices for WS - elements [NT,I] and [NT,J]
       NP,         // NPTS element
 NPTS[51],         // Array containing the indices of a sequence of
                   //  nodes to be used in the least squares fit
                   //  or to compute RSQ.  The nodes are ordered
                   //  by distance from K, and the last element
                   //  (usually indexed by LNP) is used only to
                   //  determine RQ, or RSQ[K] if NW > NQ
   NQWMAX,         // MAX0(NQ,NW)
      NTT,         // Number of terms in each quadratic nodal function Q(K). 
                   //   Note: Changed from original NT so won't conflict 
                   //   with class var NT
     NTM1,         // NT-1 (the number of coefficients which, along with F(K),
                   //   define Q(K))
   NTM1NT,         // (NT-1)*NT
     NTP1,         // NT+1 (offset between diagonal elements in WS) 
 necessary,        // conditional for do-while loop
necessary4,        // conditional for do-while loop  
necessary7,        // conditional for do-while loop
      brk;         // indicates if necessary to break out of loop

 double AV,        // Root-mean-square distance between K and the
                   //  nodes in the least squares system (unless
                   //  additional nodes are introduced for stability).  
                   //  The first MT columns of the matrix (quadratic terms) 
                   //  are scaled by 1/AVSQ, the last M (linear terms) by 1/AV
     AVSQ,         //  AV*AV
        C,         // First component of the plane rotation used to
                   //  zero the lower triangle of the regression
                   //  matrix; computed by GIVENS
     DMIN,         // Minimum of the magnitudes of the diagonal elements of 
                   //  the regression matrix after zeros are introduced 
                   //  below the diagonal
     DTOL = 0.01,  // Tolerance for detecting an ill-conditioned system.  
                   //   The system is accepted when DMIN*RQ >= DTOL.
       FK,         // Data value at node K - F(K)
       RQ,         // Radius of influence which enters into the weights for 
                   //  Q(K) (see subroutine SETUPM)
       RS,         // Squared distance between K and NPTS[LNP]; used to 
                   //  compute RQ and RSQ[K]
     RSMX,         // Maximum RSQ element encountered
    RSOLD,         // Squared distance between K and NPTS[LNP-1]; used to 
                   //  compute a relative change in RS between succeeding 
                   //  NPTS elements.
 RTOL = 1.0E-05,   // Tolerance for detecting a sufficiently large relative 
                   //  change in RS.  If the change is not greater than RTOL, 
                   //  the nodes are treated as being the same distance from K
      RWS,         // Current value of RSQ[K]. 
        S,         // Second component of the plane Givens rotation
       SF = 1.0,   // Marquardt stabilization factor used to damp out the 
                   //  first MT solution components (second partials 
                   //  of the quadratic) when the system is ill-conditioned.
                   //  As SF increases the fitting function approaches a 
                   //  linear
      SUM,         // Sum of squared Euclidean distances between node K and 
                   //  the nodes used in the least squares fit (unless 
                   //  additional nodes are added for stability)
        T;         // Temporary variable for accumulating a scalar product in
                   //  the back substitution loop.

   NQWMAX = MAX0(NQ,NW);
   LMAX = MIN0(50,N-1);
   NTM1 = (M*(M+3))/2;
   NTT = NTM1 + 1;
   NTP1 = NTT + 1;
   MT = NTM1 - M;
   NTM1NT = NTM1*NTT;
   C=0;
   S=0;

   // make sure variables are legit - if not, exit the program
   if((M<1)||(NQ<NTM1)||(NW<1)||(NQWMAX>LMAX)||(NR<1)){
     IER = 1;
     return;
    }

   // Create the cell data structure, and initialize RSMX.
   STOREM(IERR);
   if(IERR!=0){    
     IER = 3;           // All nodes are coplanar
     return;            // Return from QSHEPM now
    }

   RSMX = 0;

   // OUTER LOOP on node K
   for(K=1;K<N+1;K++){
     FK = W[K];         // Note: Fortran code used F as W since W was passed 
                        //  in as a parameter. Here W is class var. so 
                        //  no parameter passing was necessary.

     // Mark node K to exclude it from the search for nearest neighbors.
     LNEXT[K] = -LNEXT[K]; 

     // Get the index of the nearest neighbor, test for duplicate nodes, and 
     //   initialize loop on remaining NPTS elements.
     GETNPM(K,NP,RS);
     if(RS==0){
       IER = 2;        // Duplicate nodes encountered
       return;
      } 

     LNP = 1;
     NPTS[1] = NP;
     SUM = RS;
     RWS = 0.0;
     RQ = 0.0;
  
     necessary = 1;    // initialize conditional varible for do-while loop
     brk = 0;          // initialize break variable
     // compute NPTS, LNP, RWS, RQ, NEQ, and AVSQ
     do { // do-while #1  -  similar to loop #1 in fortran code
       LNP++;
       RSOLD = RS;
       GETNPM(K,NP,RS);
       if(NP == 0){
         IER = 2;      // Duplicate nodes encountered
         return;
        }

       NPTS[LNP] = NP; 
       if(((RS - RSOLD)/RSOLD) < RTOL){
         SUM += RS;                    // GO TO 2
        }
       else {
         if((RWS == 0)&&(LNP > NW)){
           RWS = RS;
          }
         if((RQ == 0)&&(LNP > NQ)){     //RQ=0 (not yet computed) and 
                                        //LNP>NQ.  RQ=SQRT(RS) is 
                                        //sufficiently large to (strictly)
                                        //include NQ nodes. The least squares
                                        //fit will include NEQ=LNP-1 
                                        //equations for
                                        //NT-1 <= NQ < LMAX <= N-1 
           RQ = (float)sqrt((double)RS);
           NEQ = LNP - 1;
           AVSQ = SUM/(float)NEQ;
          }

         // Test for termination on LNP > NQ and LNP > NW
         if(LNP > NQWMAX){
           brk = 1;
           necessary = 0;  // GO TO 3
          }
         else{
           SUM = SUM + RS;
          }
        }
       if(LNP == LMAX){    // #2 - loop conditional
         necessary = 0; 
        }
      }
     while(necessary);    // end of do-while #1

     // If brk==0 then no GO TO's (or breaks) encountered and
     // all LMAX nodes are included in NPTS.  RWS and/or RQ**2 is
     // (arbitrarily) taken to be 10% larger than the distance RS to the
     // last node included.

     if(!brk){
       if(RWS == 0){
         RWS = 1.1*RS;
        }
       if(RQ == 0){
         RQ = (float)sqrt((double)(1.1*RS));
         NEQ = LMAX;
         AVSQ = SUM/(float)NEQ;
       }
      }

     // Store RSQ[K], update RSMX if necessary, and compute AV.
     RSQ[K] = RWS; // #3
     if(RWS > RSMX){
       RSMX = RWS;
      }
     AV = (float)sqrt((double)AVSQ);
//cout << "QSHEPM: AV = " << AV << endl; 

     // Set up the augmented regression matrix (stored by rows in WS, and
     //  zero out the lower triangle and Givens rotations - QR decomposition
     //  with orthogonal matrix Q not stored.
     brk = 0;
     necessary4 = 1;
     I = 0;
     do {   //do-while #4
       I++;
       NP = NPTS[I];
       IROW = MIN0(I,NTT);
       IRM1 = IROW-1;
       IJ = IRM1*NTT + 1;
       SETUPM(K,FK,NP,AV,AVSQ,RQ,IJ); 
       if(I != 1){            // else go back to top of do-while #4. Same
                              //   as GO TO 4
        JJ = 1;
        for(J=1;J<IRM1+1;J++){
          GIVENS(&WS[JJ],&WS[IJ],C,S);
          IJ++;
          JJ++;
          ROTATE(NTT-J,C,S,JJ,IJ);
          JJ += NTT;
         }
        if(I >= NEQ){         // else go back to top of do-while #4. Same
                              //   as GO TO 4
          //Test the system for ill-conditioning
          DMIN = abs(WS[1]);
          JJ = 1;
          for(J=2;J<NTM1+1;J++){
            JJ += NTP1;
            DMIN = MIN0(DMIN,abs(WS[JJ]));
           }
          if(DMIN*RQ < DTOL){  // else will GO TO #13
            if(NEQ != LMAX){  // else will GO TO #9
              // Increase RQ and add another equation to the system to
              //   improve the conditioning. The number of NPTS elements
              //   is also increased if necessary.
              necessary7 = 1;
              do{            //do-while #7
                RSOLD = RS;
                NEQ++;
                if(NEQ != LMAX){  //else break out of #7 and GO TO top of #4
                  if(NEQ < LNP){
                    NP = NPTS[NEQ+1];
                    RS = 0.0;
                    for(J=1;J<M+1;J++){ //#8
//                      RS += (float)pow((float)(X[J][NP]-X[J][K]),(float)2);
                    RS += (float)(X[J][NP]-X[J][K])*(float)(X[J][NP]-X[J][K]);
                     }
                   } 
                  else{
                    NEQ = LNP;
                    LNP++;
                    GETNPM(K,NP,RS);
                    if(NP == 0){
                      IER = 2;    //duplicate nodes encountered by GETNPM
                      return;
                     }
                    NPTS[LNP] = NP;
                   }
                  if((RS - RSOLD)/RSOLD >= RTOL){
                    RQ = (float)sqrt((double)RS);
                    necessary7 = 0;   //break out of #7 and GO TO top of #4 
                   }
                  }
                else{             // NEQ==LMAX
                  RQ = (float)sqrt((double)(1.1*RS));
                  necessary7 = 0;
                 } 
                } // end of do-while #7
               while(necessary7);
             }
            else{
             necessary4 = 0;
             brk = 9; //GOTO #9
             }
            }
           else{
            necessary4 = 0; //GOTO #13
            }
         } // end of if I>=NEQ
        } // end of if I!=1
      } // end of do-while #4
     while(necessary4); 

     if(brk == 9){          //#9  
                            //Stabilize the system by damping second
                            // partials - add multiples of the first MT
                            // unit vectors to the first MT equations,
                            // where MT=M*(M+1)/2 is the number of 
                            // quadratic terms. 
       NI = NTM1NT;
       for(I=1;I<MT+1;I++){ //#11
         NI++;
         WS[NI] = SF;
         NJ = NI;
         for(J=I+1;J<NT+1;J++){ //#10
           NJ++;
           WS[NJ] = 0.0;
          }
         JJ = NTP1*I - NTT;
         NJ = NI;
         for(J=I;J<NTM1+1;J++){
           GIVENS(&WS[JJ],&WS[NJ],C,S);
           JJ++;
           NJ++;
           ROTATE(NTT-J,C,S,JJ,NJ);
           JJ += NTT;
          }
        } //end of #11

       //Test the stabilized system for ill-conditioning
       DMIN = abs(WS[1]);
       JJ = 1;
       for(J=2;J<NTM1+1;J++){   // #12
         JJ += NTP1;
         DMIN = MIN0(DMIN,abs(WS[JJ]));
        } //end of #12

       if(DMIN*RQ < DTOL){
         IER = 3;              //No unique solution due to a nearly
                               // singular system.
         return;
        }
      }
     
     // Solve the order NT-1 upper triangular system for the coefficients

     II = NTM1NT - 1;          //   #13
     for(IB=1;IB<NTM1+1;IB++){ //  #15
       I = NT-IB;
       T = 0.0;
       if(I != NTM1){
         IJ = II;
         for(J=I+1;J<NTM1+1;J++){ //#14
           IJ++;
           T += WS[IJ]*A[K][J];
          } // end of #14
        }
       A[K][I] = (WS[(II+IB)]-T)/WS[II];
       II -= NTP1;
      } // end of #15

     // Scale the coefficients to adjust for the column scaling.
     for(I=1;I<MT+1;I++){  // #16
       A[K][I] /= AVSQ;
      } // end of #16

     for(I=MT+1;I<NTM1+1;I++){  // #17
       A[K][I] /= AV;
      } // end of #17

     // Unmark K and the elements of NPTS.
     LNEXT[K] = -LNEXT[K];
     for(I=1;I<LNP+1;I++){   // #18
       NP = NPTS[I];
       LNEXT[NP] = -LNEXT[NP];
      } // end of #18
    }//end of outer for-loop #19 (OUTER LOOP over all nodes K )

   // No errors encountered
   RMAX = (float)sqrt((double)RSMX);
   IER = 0;

} // end of QSHEPM




/******************************************************************
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 

                              QSMVAL
  This function returns the value Q(P) where Q is the weighted sum
of the quadratic nodal functions defined in function QSHEPM. 
  Input Parameters:
      P      -Vector of length M containing the cartesian coordinates
              of the point at which Q is to be evaluated.
Upon completion, QSMVAL returns the function value Q(P) unless M, N,
NR, an element of DX, or RMAX are invalid.  

**********************************************************************/
                         

double
Interpolator::QSMVAL(float *P){
  int I,J,K,         // Loop and array indices
         KP,         // Temp variable for K and index into LNEXT
       IMIN,         // Minimum cell index at I
       IMAX,         // Maximum cell index at I
          L,         // Counter and index
        LCI,         // Index into LCELL
         MT,         // Number of quadratic terms in the nodal fxns. Q(K)
        NNR,         // Local variable for NR
        NTT,         // Total number of terms in Q(K)
        NRP,         // Affects LCI index in proportion to NR
        brk;         // Signals if there was a break out of a do-while loop

  float  DS,         // Used in determining weights
         QK,         // Value of Q(K)(P)
         RD,         // Used in determining weights
        RDS,         //  "    "      "        "
         RM,         // Local variable for RMAX 
         RS,         // RSQ[K]
         SW,         // Accumulated weight values
        SWQ,         // Accumulated weighted nodal function values
          T,         // Variable used in determining QK
         WW,         // W(K)
     qsmval,         // Final value Q(P) that is returned
      cells,         // Conditional for do-while that loops over cells 
      nodes;         // Conditional for do-while that loops over nodes

   NNR = NR;
   RM = RMAX;
      
   MT = (M*(M+1))/2;
   NTT = MT + M + 1;
   if((M < 1)||(N < NTT)||(NNR < 1)||(RM < 0)){
     return -0;
    }

// Cell (L[1],L[2],...,L[M]) has index LCI = SUM[(L[I]-1)*
//   NR**(I-1)]+1 where the sum is over I = 1 to M.  Set
//   IW[][1] and IW[][2] to minimum and maximum cell indices
//   (M-tuples) defining the range of the search for nodes
//   whose radii include P.  The cells which must be searched
//   are those intersected by an M-ball of radius RMAX cen-
//   tered at P.  No cells are within RMAX of P if IMIN =
//   IW[I][1] > IMAX = IW[I][2] for any I.  IW[][3] is initial-
//   ized to IW[][1], and LCI is initialized to the index
//   associated with IW[][3] by Horner's method.

   LCI = 0;
   for(I=M;I>0;I--){  //Loop #1
     if(DX[I] <= 0){
       return -0;
      }
     IMIN = (int)((P[I]-XMIN[I]-RM)/DX[I]) + 1;
     if(IMIN < 1){
       IMIN = 1;
      }
     IMAX = (int)((P[I]-XMIN[I]+RM)/DX[I]) + 1;
     if(IMAX > NNR){
       IMAX = NNR;
      }
     if(IMIN > IMAX){
       qsmval = 0;
       return qsmval;
      }
     LCI = LCI*NNR + IMIN - 1;
     IW[I][1] = IMIN;
     IW[I][2] = IMAX;
     IW[I][3] = IMIN;
    }                //end of Loop #1

   LCI++;

// Accumulate weight values in SW and weighted nodal function
//   values in SWQ.  The weights are W(K) = ((R-D)+/(R*D))**2
//   for R**2 = RSQ[K] and D = distance between P and node K.

   SW = 0.0;
   SWQ = 0.0;

   // outer loop on on cells LCI
   cells = 1;
   do{                   //do-while #2
     brk = 0;            // break out variable to allow for GO TO's
     K = LCELL[LCI];
     if(K != 0){         // else break out to #8
       nodes = 1;
       // Inner loop on nodes K.  Compute WW = W(K) and update SW.
       do{ //#3
         DS = 0.0;
         for(I=1;I<M+1;I++){
//           DS += (float)pow((float)(P[I]-X[I][K]),(float)2);
           DS += (float)(P[I]-X[I][K])*(float)(P[I]-X[I][K]);
          }
         RS = RSQ[K];
         if(RS > DS){  // GOTO #7
          if(DS != 0){  // GOTO #10
            RDS = RS * DS;
            RD = (float)sqrt((double)RDS);
            WW = (RS+DS-RD-RD)/RDS;
            SW += WW;

            // Compute QK = Q(K) (P) and update SWQ
            QK = W[K];
            L = 0;
            for(J=1;J<M+1;J++){
              T = A[K][MT+J];
              for(I=1;I<J+1;I++){
                L++;
                T += A[K][L]*(P[I]-X[I][K]);
               }
              QK += T*(P[J]-X[J][K]);
             }
            SWQ += WW*QK;
           }
          else{  // GOTO #10
            qsmval = W[K];   // P coincides with node K
            return qsmval;
           }
         }  // GOTO #7
         KP = K;
         K = LNEXT[KP];
         if(K == KP){
           nodes = 0;
          }
        }
       while(nodes);       // end of do-while #3 on nodes in cell LCI
      }  // end of if(K!=0)

     NRP = 1;              //statement #8
     for(I=1;I<M+1;I++){   // for-loop #9
       if(IW[I][3] < IW[I][2]){
         IW[I][3]++;
         LCI += NRP;
         brk = 1;
         break;          // break out of for-loop #9 GO TO top of do-while #2
        }
       else{
         if(I < M){
           IW[I][3] = IW[I][1];
           LCI -= (IW[I][2]-IW[I][1])*NRP;
           NRP *= NNR;
          }
        }
      }     // end of for-loop #9
     if(!brk){  // no break in for-loop #9 so do-while #2 terminiates
       cells = 0;
      } 
    }
   while(cells);        // end of do-while #2 on cells. NRP = NR**(I-1)
   
   // SW==0 iff P is not within the radius R(K) for any node K
   if(SW == 0){
     qsmval = 0;  //All weights are 0 at P
    }
   else{
     qsmval = SWQ/SW;
    }
  
   return qsmval;
   
 } // End of QSMVAL




/********************************************************************

 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 
                                GETNPM
    Given a set of N nodes and the data structure defined in
  subroutine STOREM, this subroutine uses the cell method to
  find the closest unmarked node NP to a specified point P.
  NP is then marked by setting LNEXT[NP] to -LNEXT[NP].  (A
  node is marked if and only if the corresponding LNEXT ele-
  ment is negative.  The absolute values of LNEXT elements,
  however, must be preserved.)  Thus, the closest j nodes to
  P may be determined by a sequence of j calls to this rou-
  tine.  Note that if the nearest neighbor to node K is to
  be determined (P = X[][K]), then K should be marked before
  the call to this routine.
    The search is begun in the cell containing (or closest
  to) P and proceeds outward in (M-dimensional) rectangular
  layers until all cells which contain points within dis-
  tance R of P have been searched, where R is the distance
  from P to the first unmarked node encountered (infinite if
  no unmarked nodes are present).
    Input Parameters:
       K       -current node whose nearest unmarked neighbor is
                to be found.
    Output Parameters: NP      -Nodel index (column index for X) of the nearest
                unmarked node to P, or 0 if all nodes are marked,
                M < 1, NR < 1, or an element of DX is not
                positive.  LNEXT[NP] < 0 unless NP = 0.
       DSQ     -Squared Euclidean distance between P and node
                NP, or 0 if NP = 0.

***********************************************************************/

void
Interpolator::GETNPM(int K,int &NP,double &DSQ){

  int I,J,          // Loop indices LCI,          
      LCI,          // Index for LCELL,
                    //   Sum[(IW[I][3]-1)*NR**(I-1)]+1 where the sum is over
                    //   I = 1 to M
      L=0,          // Values of LCELL used as indices as well
      LL,           // Temp var. for determining min. and max. values
      LMIN,         // Current nodel index of the nearest unmarked node 
      LN=0,         // Values of LNEXT used as indices as well
      nodes,        // Conditional for do-while loop over nodes
      cells,        // Conditional for do-while loop over cells
      layers,       // Conditional for do-while loop over layers
      brk;          // Signal for breaking out of loop

   double RRSQ,     // Local RSQ (see Interp.h)
     RSMIN = 0.0,   // RSQ of LMIN
             R;     // Square root of RSMIN used to find min and max
                    //  nodal indices of search range of nearest 
                    //  neighbor

   bool FIRST;      // Signals first unmarked neighbor of K. Is equal to 
                    // TRUE iff the first unmarked node has yet to be
                    // encountered.
  /* IW initialization:  
       IW[I][4],IW[I][5] = Minimum and maximum cell indices, re-
                           spectively, defining the range of the
                           search.
       IW[I][1],IW[I][2] = Minimum and maximum cell indices, re-
                           spectively, defining the layer whose
                           intersection with the range is cur-
                           rently being searched.  These are in-
                           itialized to the cell containing or
                           closest to P.
                IW[I][3] = Indices of the current cell in the range
                           Max(IW(I,1),IW(I,4)) to Min(IW(I,2),IW(I,5)).
  */

  if((M<1)||(NR<1)){    // Check M and NR to make sure legit
    NP = 0;
    DSQ = 0.0;
    return;
   }

  FIRST = true;

  for(I=1;I<M+1;I++){ // DO-LOOP #1
    IW[I][4] = 1;
    IW[I][5] = NR;
    if(DX[I] <= 0){    // Invalid value - "GO TO 18"
      NP = 0;
      DSQ = 0.0;
      return;
     }
    LL =  (int)MIN0(NR,(int)((X[I][K]-XMIN[I])/DX[I])+1);
    IW[I][1] = MAX0(1,LL);
    IW[I][2] = IW[I][1];
   } // end of DO-LOOP #1

  layers = 1;  
  // Top of outer loop on layers: initialize IW[I][3]
  do {  
    brk = 0;
    for(I=1;I<M+1;I++){      // DO-LOOP at #2
      IW[I][3] = MAX0(IW[I][1],IW[I][4]);
     } // end of DO-LOOP at #2

    cells = 1;
    // Top of loop on cells: bypass cells interior to the layer.
    do {
      brk = 0;
      for(I=1;I<M+1;I++){   // DO-LOOP at #4
        if((IW[I][3]==IW[I][1])||(IW[I][3]==IW[I][2])){  // #5
          brk=1;     // if #5==true break out and GO TO #6
          break;     //   else GO TO #12
         }
        } // end of DO-LOOP at #4

        if(brk){   // #6  
          // Compute LCI by Horner's method and test for an empty cell
          LCI = 0;                
          for(J=M;J>0;J--){  // DO-LOOP #7
            LCI = LCI*NR + IW[J][3] - 1;
           } // end of DO-LOOP #7

          LCI++; 
          L = LCELL[LCI];
          if(L != 0){    // Else GO TO #12
            nodes = 1;

            // Loop on nodes in cell LCI
            do{   // Do-While #8
              LN = LNEXT[L]; 
              if(LN >= 0) {   // Else GO TO #11

                // Node L is not marked. Set RSQ to its squared distance from P
                RRSQ = 0.0;
                for(J=1;J<M+1;J++){    // DO-LOOP #9
//                  RRSQ += (double)pow((float)(X[J][K]-X[J][L]),(float)2);
                  RRSQ += (double)(X[J][K]-X[J][L])*(double)(X[J][K]-X[J][L]);
                 }  // end of DO-LOOP #8
                if(FIRST){
                 /*  Node L is the first unmarked neighbor of P encountered,
                      and hence the first candidate for NP.  Initialize LMIN
                      and RSMIN to L and RSQ, and update the search range
                      IW[I][4] and IW[I][5] to the smallest M-box containing
                      a hypersphere of radius R = SQRT(RSMIN) centered at P
                      and contained in [1][NR]**M.  FIRST is reset to FALSE.
                  */

                  LMIN = L;
                  RSMIN = RRSQ;
                  R = sqrt((double)RSMIN);
                  for(J=1;J<M+1;J++){  // DO-LOOP #10
                    IW[J][4] = MAX0(1,(int)((X[J][K]-R-XMIN[J])/DX[J])+1);
                    IW[J][5] = MIN0(NR,(int)((X[J][K]+R-XMIN[J])/DX[J])+1);
                   } // end of DO-LOOP #10
                  FIRST = false;
                 }
                else{
                  if(RRSQ < RSMIN){
                    // Update LMIN and RSMIN for node L closer than LMIN to P.
                    LMIN = L;
                    RSMIN = RRSQ;
                   }
                 }
               } // end if(LN >=0 )

              // Test for termination of loop on nodes in cell LCI.
              if((int)abs((float)LN) != L){
                L = (int)abs((float)LN);
               }
              else{
                nodes = 0;  // terminate the loop
               } 
             }
            while(nodes);

           } // end if(L != 0)
  
          // #12 - Bottom of loop on cells. Update the current cell indices.
          brk = 0;
          for(J=1;J<M+1;J++){         // DO-LOOP #13
            int tmp = MIN0(IW[J][2],IW[J][5]);
            if(IW[J][3] < tmp){
              IW[J][3] = IW[J][3] + 1;
              brk = 1;
              break;
             }
            else {
              IW[J][3] = MAX0(IW[J][1],IW[J][4]);
             }
           } // end of #12
          }
         else {  // #12 from #5 
           // Bottom of loop on cells.  Update the current cell indices.
           for(J=1;J<M+1;J++){    // DO-LOOP #13
             int tmp = MIN0(IW[J][2],IW[J][5]);
            if(IW[J][3] < tmp){
              IW[J][3] = IW[J][3] + 1;
              brk = 1;
              break;
             }
            else {
              IW[J][3] = MAX0(IW[J][1],IW[J][4]);
             }
           }  // end of DO-LOOP #13
         }
       if(!brk){    // If no break-out above, then drop out of do-while on cells
         cells = 0;
        }
      } // end of inner do-while 
     while(cells);

     // Test for termination of loop on cell layers.
     brk = 0;
     for(I=1;I<M+1;I++){   // DO-LOOP #14
       if((IW[I][1] > IW[I][4])||(IW[I][2] < IW[I][5])) {

          // Update IW[I][1] and IW[I][2] to the next layer out
          for(J=1;J<M+1;J++){  // #15
            IW[J][1] = IW[J][1] - 1;
            IW[J][2] = IW[J][2] + 1;
           }
          brk = 1;
          break;   // break out of for-loop and go back to outer do-while
                   //  "GO TO 2"
        }
      }  // end of DO-LOOP #14

     if(!brk){     // Terminate loop on layers
       layers = 0; 
      }
    } // end of outer do-while
   while(layers); 


   // Unless no unmarked nodes were encountered, LMIN is the closest
   //   unmarked node to P.
   if(FIRST) {  // no unmarked nodes encountered
     NP = 0;
     DSQ = 0.0;
     return;
    }
   else {       
     NP = LMIN;
     DSQ = RSMIN;
     LNEXT[LMIN] = -LNEXT[LMIN];
     return;
   } 
  
 }  // end of GETUPM




/************************************************************************
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 

                                   GIVENS
    This routine constructs the Givens plane rotation
      ( C  S)
  G = (     ), where C*C + S*S = 1, which zeros the second
      (-S  C)
  entry of the 2-vector (A B)-Transpose.  A call to GIVENS
  is normally followed by a call to ROTATE, which applies
  the transformation to a 2 by N matrix.  This routine was
  taken from LINPACK.
 
  On input:
 
        A,B = Components of the 2-vector to be rotated..
 
  On output:
 
        A = Value overwritten by R = +/-SQRT(A*A + B*B).
 
        B = Value overwritten by a value Z which allows C
            and S to be recovered as follows:
              C = SQRT(1-Z*Z), S=Z     if ABS(Z) .LE. 1.
              C = 1/Z, S = SQRT(1-C*C) if ABS(Z) .GT. 1.
 
        C = +/-(A/R).
 
        S = +/-(B/R).
 
**************************************************************************/
void
Interpolator::GIVENS(float *A, float *B, double &C, double &S){
   float AA,       // Local copy of A 
         BB,       // Local copy of B
          R,       // C*A + S*B = +/-SQRT(A*A+B*B)
        U,V;       // Variables used to scale A and B for computing R


   AA = *A;
   BB = *B;

   if(abs(AA) > abs(BB)){
     U = AA + AA;
     V = BB/U;
     R = (float)sqrt((double)(0.25+V*V)) * U;
     C = AA/R;
     S = V * (C+C);
     // Note: R has the sign of A, C > 0, and S has SIGN(A)*SIGN(B)
     *B = S;
     *A = R;
     return;
    }
   else{    // abs(AA) <= abs(BB)
     if(BB == 0.0){ // #1
       C = 1;       // #2 -  A = B = 0
       S = 0;
       return;
      }
     else{   // BB != 0
       U = BB + BB;
       V = AA/U;

       // Store R in A
       *A = (float)sqrt((double)(0.25+V*V)) * U;
       S = BB/(*A);
       C = V * (S+S);
      
       // Note: R has the sign of B, S > 0, and C has SIGN(A)*SIGN(B).
       *B = 1.0;
       if(C != 0.0)
         *B = 1.0/C;
       return;
      }
    }
}  // end of GIVENS




/*****************************************************************************
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 

                                 ROTATE
                                             ( C  S)
    This routine applies the Givens rotation (     ) to the
                                             (-S  C)
                (X[1] ... X[N])
  2 by N matrix (             ).
                (Y[1] ... Y[N])
 
  On input:
 
        N = Number of columns to be rotated.
 
        C,S = Elements of the Givens rotation.  These may be
              determined by subroutine GIVENS.
 
        II,JJ = Indices into arrays of length .GE. N containing the vectors
              to be rotated.
 
  Parameters N, C, and S are not altered by this routine.
 
  On output:
 
        IJ,JJ = Rotated vectors at indices IJ and JJ. 
                Note: instead of passing vectors which are part of a larger
                      vector as was done in the original fortran code, 
                      the indices marking the smaller vectors were passed
                      instead.  The larger vector, WS, is a class variable.

****************************************************************************/
void
Interpolator::ROTATE(int NN,double C,double S,int IJ,int JJ){
   int I;           // Loop index

 float XI,          // X_PTR[I]
       YI,          // Y_PTR[I]
   *X_PTR,          // Points to vector in WS at IJ
   *Y_PTR;          // Points to vector in WS at JJ

   X_PTR = &WS[IJ-1];
   Y_PTR = &WS[JJ-1];

   if((NN <= 0)||((C == 1)&&(S == 0))){
     return;
    }
   else{
     for(I=1;I<NN+1;I++){
       XI = X_PTR[I];
       YI = Y_PTR[I];
       X_PTR[I] = C*XI + S*YI;
       Y_PTR[I] = -S*XI + C*YI;
      }
    }
 }  // end of ROTATE


/*************************************************************************
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 

                                SETUPM
    This routine sets up the L-th row of an augmented re-
  gression matrix for a weighted least squares fit of a
  quadratic function Q(X) to a set of data values F, where
  Q(XK) = FK and Q is a function of M variables (X is an
  element of Euclidean M-space).  The upper triangular
  array of M(M+1)/2 quadratic terms (XL[I]-XK[I])*
  (XL[J]-XK[J]), I <= J, are ordered by columns and scaled
  by W/S2.  These are followed by the M linear terms, scaled
  by W/S1, and the last element of ROW is the right hand
  side W*(FL-FK).  The weight is W = (R-D)/(R*D) if R > D,
  or W = 0 if R <= D, where D is the distance between
  nodes L and K.
 
  On input:
            K = current node and column of X which contains the
                Cartesian coordinates of node K 

           FK = interpolated value at node K
 
           NP = index marking column of X which contains the
                coordinates of node L.  Also, the index in
                array W which is the data value at node L.
 
        S1,S2 = Reciprocals of the scale factors for the
                linear and quadratic terms, respectively.
 
            R = Radius of influence about node K defining the
                weight W.

          ROW = Array of length NT = (M+1)(M+2)/2.

  On output:
 
         ROW = Array containing row L of the augmented re-
               gression matrix as described above, or zero
               vector if M < 1 or nodes K and L coincide.

*********************************************************************/ 
void
Interpolator::SETUPM(int K,double FK,int NP,double S1,double S2,double R,int ROW){
   int I,J,          // Loop indices
        IJ,          // Counter and index
        MT,          // Number of quadratic terms in the nodal functions Q(K) 
       NNT;          // Total number of terms in Q(K)
      

  float   D,         // Distance between nodes L and K
   W0,W1,W2,         // Weight values
        DXJ,         // Difference between respective coordinates of L and K 
    *WS_PTR;         // Pointer into array WS


    WS_PTR = &WS[ROW-1];

    MT = (M*(M+1))/2;
    NNT = MT + M + 1;
    if(M >= 1){   // else "GO TO 4"

      // Compute distance D , weights W=(R-D)/(R*D), W1=W/S1
      D = 0.0;
      for(I=1;I<M+1;I++){
//        D += (float)pow((float)(X[I][NP]-X[I][K]),(float)2);
        D += (float)(X[I][NP]-X[I][K])*(float)(X[I][NP]-X[I][K]);
       }  
      D = (float)sqrt((double)D);
      if((D > 0)&&(D < R)){  // else "GO TO 4"
        W0 = (R-D)/R/D;
        W1 = W0/S1;
        W2 = W0/S2;

        // Store the row
        IJ = 0;
        for(J=1;J<M+1;J++){
          DXJ = X[J][NP] - X[J][K];
          for(I=1;I<J+1;I++){
            IJ++;
            WS_PTR[IJ] = (X[I][NP] - X[I][K])*DXJ*W2;
           }
          WS_PTR[MT+J] = DXJ*W1;
         }

        WS_PTR[NNT] = (W[NP]-FK)*W0;
        return;
       } 
    }

   // M < 1, D = 0 (nodes K and L coincide), or node L is outside the
   //   radius of influence.  Set WS to the zero vector.

   for(I=1;I<NNT+1;I++){   //#4
     WS_PTR[I] = 0.0;
    }
  
 } // end of SETUPM
     



/****************************************************************************
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 

                                 STOREM
   Given a set of N nodes arbitrarily distributed in
 Euclidean M-space, this function creates a data struc-
 ture for a cell-based method for solving closest-point
 problems.  The smallest M-dimensional box containing the
 nodes is partitioned into an NR**M uniform grid of cells,
 and the indices of the nodes contained in each cell are
 stored in the data structure.  For a uniform random dis-
 tribution of nodes, the nearest node to an arbitrary point
 can be determined in constant expected time.  Refer to
 function GETNPM.

   Input Parameter:
        IER         -returns results of STOREM
                     IER = 0 if no errors were encountered.
                     IER = 1 if M < 1, N < 2, or NR < 1.
                     IER = 2 if DX[I] = 0 for some I (the nodes lie
                           in an affine subspace of dimension
                           less than M).

Upon completion:
        LCELL - Array of length NR**M containing nodal indi-
                ces of the first node in each cell, with 0
                entries corresponding to empty cells.  The
                nodes in each cell are ordered by their ind-
                ices, so that the first node refers to the
                one with smallest index.  Cells are ordered
                in the manner of FORTRAN array storage, with
                the first subscript varying fastest.  Thus,
                LCELL[LCI] contains the index of the first
                node in cell (L[1],L[2],...,L[M]) for LCI =
                SUM[(L[I]-1)*NR**(I-1)]+1 (with the sum over
                I = 1 to M), where the cell is defined by
                intervals [XMIN[I]+(L[I]-1)*DX[I],XMIN[I]+
                L[I]*DX[I]], for each L[I] in the range 1 to
                NR.  LCELL is not defined if IER .NE. 0 on
                output.
 
        LNEXT - Array of length N containing next-node indi-
                ces such that LNEXT[K] is the index of the
                next node in the cell which contains node K,
                or LNEXT[K] = K if K is the last node in the
                cell.  If, for example, cell LCI contains
                nodes 2, 3, and 5 (and no others), then
                LCELL[LCI] = 2, LNEXT[2] = 3, LNEXT[3] = 5,
                and LNEXT[5] = 5.  LNEXT is not defined if
                IER .NE. 0 on output.
 
        XMIN - Array of length M containing the minimum
               nodal coordinates unless IER = 1.
 
        DX - Array of length M containing the cell dimen-
             sions (interval lengths) unless IER = 1.  DX[I]
             = (XMAX[I]-XMIN[I])/NR where XMIN and XMAX con-
             tain the extrema defining the smallest M-box
             which contains the nodes.

**********************************************************************/

void
Interpolator::STOREM(int &IER){

   int I,K,J,       // Loop indices
         L,         // Index for LNEXT
       LCI,         // Index into LCELL
        LI,         // Used to help determine LCI
        LN,         // Used with L
        NC;         // Loop index = total number of cells

 float RNR;         // Local variable for NR

   if((M<1)||(N<2)||(NR<1)){  // Check for legitimate variables
     IER = 1;   // #8
     return;
    }

   // Compute the corner coordinates XMIN and XMAX (stored in DX) of
   //   the smallest box containing the nodes.
   for(I=1;I<M+1;I++){ 
     XMIN[I] = X[I][1];
     DX[I] = X[I][1];
    }
   
   for(K=2;K<N+1;K++){
     for(I=1;I<M+1;I++){
       if(X[I][K] < XMIN[I]) XMIN[I] = X[I][K];
       if(X[I][K] > DX[I])  DX[I] = X[I][K];
      }
    }

   // Compute the interval widths, and test for a zero width.
   RNR = (float)NR;
   for(I=1;I<M+1;I++){
     DX[I] = (DX[I] - XMIN[I])/RNR;
     if(DX[I]==0){
       IER = 2;   // #9
       return;
      }
    }

   // Initialize LCELL to zeros.
   NC = (int)pow((float)NR,(float)M);
   for(LCI=1;LCI<NC+1;LCI++){ 
     LCELL[LCI] = 0; 
    }

   // Outer loop on nodes
   for(K=1;K<N+1;K++){
     // Compute the index LCI of the cell containing node K 
     //  (using Horner's method).
     LCI = 0;
     for(I=M;I>0;I--){
       LI = (int)( (X[I][K]-XMIN[I])/DX[I] ) + 1;
       if(LI > NR) LI=NR;
       LCI = LCI*NR + LI-1;
      }
     LCI++;

     // Add K to the data structure as the last node in cell LCI.
     LN = LCELL[LCI];

     if(LN==0){
       // K is the first node in cell LCI.
       LCELL[LCI] = K;
      }
     else {
      // Find the index L such that LNEXT[L] = L
      do {
        L = LN;
        LN = LNEXT[L];
       }
      while(LN!=L);
      LNEXT[L] = K;
     } 
    // Bottom of loop:  set LNEXT[K] = K to indicate that K is
    //                  the last node in the cell
    LNEXT[K] = K;
   }
  
  IER = 0;  // No errors encountered
 }  // end of STOREM





/*************************************************************************
                              get_node_coords
This function gets the coordinates of node K from matrix X and stores 
them in vector PP.
    Input Parameters:
            PP      -vector of length M which will obtain the 
                     coordinates of K
             K      -current node - column in X where coordinates are found
****************************************************************************/
void
Interpolator::get_node_coords(float *PP, int K){

  int i;     // loop index

  for(i=1;i<M+1;i++){
    PP[i] = X[i][K];
   }
 } // end of get_node_coords




/***************************************************************************

 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 

                                 TFUN3
This function contains the test functions used to test the interpolation
code.  This works only for 3 dimensions.  It computes the value of the 
trivariate test functions.

  Input Parameters:
       KF      -designates which function to use
                KF = 1 - EXPONENTIAL
                KF = 2 - CLIFF
                KF = 3 - SADDLE
                KF = 4 - GENTLE
                KF = 5 - STEEP

        K      -index into array W which will contain the value computed
                  value at node K.

       PP      -array which contains the X, Y, Z coordinates of node K

    iflag      -flag used in original fortran code to signal whether to
                evaluate the first partial derivatives or not.  NOT used
                here, but could easily be added.

    Fflag      -designates whether to store computed value in array F or W.

        F      -array F to hold computed value.
     
      DUM      -dummy variables not used at this point.

****************************************************************************/
void
Interpolator::TFUN3(int KF,int K,float *PP,int iflag, int Fflag, float *F,float DUM,float DUM,float DUM){

  float *tmp,     // temp variable that holds function values
       X,Y,Z;     // coordinates of node K
   
//  cout << "In TFUN3  - KF = " << KF << " K = " << K << endl;
  if((KF < 1)||(KF > 6)){
    return;
   }

  if(Fflag){
    tmp = F;
   }
  else{
    tmp = &W[K];
  } 

  X = PP[1];  Y = PP[2];  Z = PP[3];
//  cout << X << " " << Y << " " << Z << endl;

  switch (KF){

   case 1:  // Exponential
    *tmp = 0.75*(float)exp(-(pow((double)(9*X-2),2) + pow((double)(9*Y-2),2) +                 pow((double)(9*Z-2),2))/4) +
           0.75*(float)exp(-(pow((double)(9*X+1),2))/49 - (9*Y+1)/10 - 
                 (9*Z+1)/10) +
           0.50*(float)exp(-(pow((double)(9*X-7),2) + pow((double)(9*Y-3),2) +
                 (9*Z+1))/4) - 
           0.20*(float)exp(-pow((double)(9*X-4),2) - pow((double)(9*Y-7),2) - 
                 pow((double)(9*Z-5),2));
    if(iflag != 1) {
      return;
     }

    break;

   case 2:   // Cliff 
    *tmp = ((float)tanh((double)(9*(Z-X-Y)))+1)/9;
    if(iflag != 1) {
      return;
     }
    break;

   case 3:  // Saddle
    *tmp = (1.25 + (float)cos((double)(5.4*Y)))*(float)cos((double)(6*Z)) /
           (6 + 6*(float)pow((double)(3*X-1),2));
    if(iflag != 1) {
      return;
     }
    break;
 
   case 4:  // Gentle
    *tmp = (float)exp(-5.0625*(pow((double)(X-0.5),2)+pow((double)(Y-0.5),2)
           + pow((double)(Z-0.5),2)))/3;
    if(iflag != 1){
      return;
     }
    break;

   case 5:  // Steep
    *tmp = (float)exp(-20.25*(pow((double)(X-0.5),2)+pow((double)(Y-0.5),2)
           + pow((double)(Z-0.5),2)))/3;
    if(iflag != 1){
      return;
     }
    break;
  }
 }  // end of TFUN3


// Destructor function for Interpolator: Deallocates array space.
Interpolator::~Interpolator() {
  delete []X;
  delete []IW;
  delete []WS;
  delete []XMIN;
  delete []DX;
  delete []A;
 }
