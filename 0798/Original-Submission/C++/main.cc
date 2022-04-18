/***********************************************************************
 
                                               ROBERT RENKA
                                       UNIV. OF NORTH TEXAS
                                             (817) 565-2767
                                                   10/25/87
 
                                        Translated into C++ 
                                            by Karen Minser  
                                         Univ. of Tennessee
                                                     7/6/98
 
                                  QSHEPD5
 This is the main C++ driver program for the interpolation code, QSHEPD5.  
Here, the user is allowed to set up the hypervolume dimension and size, 
select whether this is a "test" which will work with predetermined
test functions or a real interpolation run.  For the test run, only 3
dimensions are allowed.  If it is a true interpolation run, dimensions can
range from 2 to 5, inclusive.  All the above input can be included in an 
input file along with the parameters:  
  NQ = Number of data points to be used in the least squares fit for
       coefficients defining the nodal functions Q(K).
  NW = Number of nodes within (and defining) the radii of influence R(K)
       which enter into the weights W(K).
  NR = Number of intervals along each coordinate axis defining the cell
       grid described in the class function Interpolator::STOREM.
There are two classes defined: Interpolator and Net.  The Interpolator
class contains the actual interpolation code, and the Net class contains
the code that allows the interpolation results to be stored in Netcdf
format and saved as a Netcdf file.

  Functions Called:
   Interpolator::QSMTST3 - Reads in the input parameters and sets up
                           the hypervolume axes.

**************************************************************************/

#include <stdio.h>
#include <iostream.h>
#include <stdlib.h>
#include <std/cmath.h>
#include "Interp.h"
#include "Net.h"
#include "Time.h"
#define DIM_MIN 2    // minimum dimensions allowed
#define DIM_MAX 5    // maximum dimensions allowed

// Tests the interpolator with predetermined functions.
void TEST_IT(Interpolator*,int*,float**,Net**);

// Performs actual interpolation
void FOR_REAL(int,Interpolator*,int*,float**,Net*);

main()
 {
   int DIM;         // dimensions of hypervolume
   int *Axes_Dim;   // length of each axis
   int i,j,k,l,m;   // loop indices
   int IER=0;       // returns error value
 
   float **Axes_Ticks,  // actual values along axes
         interp_value;  // interpolated value

   char ans;            // answer to test or not 

// Set up hypervolume
   cout << "Enter dimensions of interpolator: ";
   cin >> DIM;
   cout << "\n";
   if((DIM < DIM_MIN)||(DIM > DIM_MAX)){
     cout << DIM << " dimensions are not available. Only dimensions";
     cout << " 2 - 5 are working at this time.  Sorry!\n";
     exit(1);
    }

   Axes_Dim = new int[DIM+1];
   if(!Axes_Dim){
     cerr << "Cannot allocate space for Axes_Dim" << endl;
     exit(1);
    }

   float Axes_Ranges[DIM+1][3];   // range of axes

   if(!Axes_Ranges){
     cerr << "Cannot allocate space for Axes_Ranges" << endl;
     exit(1);
    }

   for(i=1;i<DIM+1;i++){
     cout << "Enter length of Axis " << i << ": ";
     cin >> Axes_Dim[i];
     cout << "\nEnter min and max of Axis " << i << ": "; 
     cin >> Axes_Ranges[i][1];
     cin >> Axes_Ranges[i][2];
     cout << "\n";
   }

   cout << "Dimension is " << DIM << endl;

   for (i=1;i<DIM+1;i++){
     cout << "Axis " << i << ": range = " << Axes_Dim[i] << "  min = "<< Axes_Ranges[i][1] << "  max = " << Axes_Ranges[i][2] << endl;
    }

// Instantiate an Interpolator object
  Interpolator *interp = new Interpolator(DIM);

// Determine if testing or not
  do{
    cout << "Testing? [Y][N] ";
    cin >> ans;
    cout << ans << endl;
   }
  while((ans!='Y')&&(ans!='N'));

  if(ans=='Y'){   // Testing interpolator only
     //Set up axes' values according to dimensions and ranges
     Axes_Ticks = interp->QSMTST3(Axes_Dim,Axes_Ranges,ans);
     Net **netspace = new Net*[NFUN*2]; //create array for Net objects
     int counter = 0;
     for(i=1;i<NFUN+1;i++){  //test various functions
       for(j=1;j<3;j++){
         // instantiate Net objects for all test functions 
         netspace[counter++] = new Net(DIM,Axes_Dim,Axes_Ranges,Axes_Ticks,i,j);
        }
      }
     // execute test interpolations
     TEST_IT(interp, Axes_Dim, Axes_Ticks,netspace);
     for(i=0;i<NFUN*2;i++){
       netspace[i]->~Net();  // destructor for Net objects
      }
   }
  else{ // Do real interpolation on unknown space
     //Set up axes' values according to dimensions and ranges
     Axes_Ticks = interp->QSMTST3(Axes_Dim,Axes_Ranges,ans);
     i=0; j=0;  //indices for file name created
     // instantiate Net object for interpolation results
     Net *netspace = new Net(DIM,Axes_Dim,Axes_Ranges,Axes_Ticks,i,j);
     //execute interpolation
     FOR_REAL(DIM,interp,Axes_Dim,Axes_Ticks,netspace);
     netspace->~Net();      // destructor for Net object
   }
     
  interp->~Interpolator();  // destructor for Interpolator object


  // clean up
  delete [] Axes_Dim;
  delete [] Axes_Ranges;
  delete [] Axes_Ticks;

 } // end of main




/************************************************************************

                                 FOR_REAL

This function performs interpolation on a real data set.  It allows for
interpolation of 2D-5D hypervolumes.  The interpolated results are
stored in netcdf format in a netcdf file.
  Input Parameters:
   DIM          -dimensions of hypervolume
   interp       -pointer to an Interpolator class object where the 
                 interpolation variables and functions are defined.
   Axes_Dim     -array of integers containing the dimensions of the axes.
   Axes_Ticks   -2D array of floats containing the actual values along
                 each axis.
   netspace     -pointer to a Net object which contains the variables and
                 functions needed to create the netcdf file. 

 Functions Called:
   Interpolator::QSHEPM - Builds the interpolator.
   Interpolator::QSMVAL - Performs interpolation and fills remaining
                          hypervolume.

*************************************************************************/
void
FOR_REAL(int DIM,Interpolator *interp,int *Axes_Dim,float **Axes_Ticks,Net *netspace) {
  
  int i,j,k,l,m,                // loop indices
      IER;                      // return value from QSHEPM  

  float *PP,                    // array of axes' values defining point for
                                // interpolation
        interp_value,           // actual interpolated value
        *interp_vector;         // vector of interpolation values to be
                                // stored in a netcdf file

  double time1,time2,QSMVAL_time=0; // timing variables

  long *coords,                 // coordinate values defining locations
                                // in netcdf file 
       *edges;                  // edge values defining length of vectors
                                // to be stored in netcdf file

  // set up arrays 
  PP = new float[DIM+1];
   if(!PP){
     cerr << "Cannot allocate space for PP" << endl;
     exit(1);
    }

  coords = new long[DIM];
   if(!coords){
     cerr << "Cannot allocate space for coords" << endl;
     exit(1);
    }

  interp_vector = new float[Axes_Dim[DIM]];
   if(!interp_vector){
     cerr << "Cannot allocate space for int_vector" << endl;
     exit(1);
    }

  // set up timers and build the interpolator based on scattered data 
  time1 = timer();
  interp->QSHEPM(IER);
  time2 = timer();
  cout << "QSHEPM CPU time = " << time2-time1 << endl;
  
  // check IER to determine if QSHEPM completed OK.  If not, print out
  // reasons and exit the program.
  if(IER != 0 ){
     if( IER == 2){
       cout << "   *** ERROR IN QSHEPM -- DUPLICATE NODES ENDOUNTERED ***\n";
      }
     else{
      if( IER == 3){
       cout << "   *** ERROR IN QSHEPM -- ALL NODES ARE COPLANAR ***\n";
       }
      }
     interp->~Interpolator();
     exit(1);
   }

  // QSHEPM finished OK, so continue with setting up netcdf edges variable
  edges = new long[DIM];
  for(i=0;i<DIM-1;i++){
    edges[i] = 1;
   }
  edges[DIM-1] = Axes_Dim[DIM];

// Determine which dimension is needed, set up values for interpolation in
//   array PP, perform interpolation, fill up vector of interpolated values, 
//   and store vector in netcdf file.

  switch(DIM){
    case 2:
      coords[1] = 0;                // starting point for vector
      for(i=1;i<Axes_Dim[1]+1;i++){
        PP[1] = Axes_Ticks[1][i];
        coords[0] = i-1;
        for(j=1;j<Axes_Dim[2]+1;j++){
          PP[2] = Axes_Ticks[2][j];
          interp_vector[j-1] = interp->QSMVAL(PP);  //store value
//          cout << PP[1] << " " << PP[2] << " " << interp_vector[j-1] << endl;
         }
        netspace->All_Knodes(coords,edges,interp_vector); //put in netcdf file
       }
      break;

    case 3:
      coords[2] = 0;
      for(i=1;i<Axes_Dim[1]+1;i++){
        PP[1] = Axes_Ticks[1][i];
        coords[0] = i-1;
        for(j=1;j<Axes_Dim[2]+1;j++){
          PP[2] = Axes_Ticks[2][j];
          coords[1] = j-1;
          for(k=1;k<Axes_Dim[3]+1;k++){
            PP[3] = Axes_Ticks[3][k];
            time1 = timer();
            interp_vector[k-1] = interp->QSMVAL(PP);
            time2 = timer();
            QSMVAL_time += time2 - time1;
//            cout << "vect = " << interp_vector[k-1] << endl;
           }
          netspace->All_Knodes(coords,edges,interp_vector);
         }
       }
      cout << "QSMVAL CPU time = " << QSMVAL_time << endl;
      break;

    case 4:
      coords[3] = 0;
      for(i=1;i<Axes_Dim[1]+1;i++){
        PP[1] = Axes_Ticks[1][i];
        coords[0] = i-1;
        for(j=1;j<Axes_Dim[2]+1;j++){
          PP[2] = Axes_Ticks[2][j];
          coords[1] = j-1;
          for(k=1;k<Axes_Dim[3]+1;k++){
            PP[3] = Axes_Ticks[3][k];
            coords[2] = k-1;
            for(l=1;l<Axes_Dim[4]+1;l++){
              PP[4] = Axes_Ticks[5][l];
              interp_vector[l-1] = interp->QSMVAL(PP);
              cout << "vect = " << interp_vector[l-1] << endl;
             }
            netspace->All_Knodes(coords,edges,interp_vector);
           }
         }
       }
      break;

    case 5:
      coords[4] = 0;
      for(i=1;i<Axes_Dim[1]+1;i++){
        PP[1] = Axes_Ticks[1][i];
        coords[0] = i-1;
        for(j=1;j<Axes_Dim[2]+1;j++){
          PP[2] = Axes_Ticks[2][j];
          coords[1] = j-1;
          for(k=1;k<Axes_Dim[3]+1;k++){
            PP[3] = Axes_Ticks[3][k];
            coords[2] = k-1;
            for(l=1;l<Axes_Dim[4]+1;l++){
              PP[4] = Axes_Ticks[5][l];
              coords[3] = l-1;
              for(m=1;m<Axes_Dim[5]+1;m++){
                PP[5] = Axes_Ticks[5][m];
                interp_vector[m-1] = interp->QSMVAL(PP);
             //   cout << i << " "<<j<<" "<<k<<" "<<l<<" "<<m<<"  vect = " << interp_vector[m-1] << endl;
               }
              netspace->All_Knodes(coords,edges,interp_vector);
             }
           }
         }
       }
      break;

    default:
      cout << DIM << " dimensions is not available.\n";
      cout << "You have big problems if you got to this point!\n";
      exit(1);
   }

   // clean up
   delete [] PP;
   delete [] coords;
   delete [] interp_vector;
   delete [] edges;
 } // end of FOR_REAL


/*************************************************************************
                                TEST_IT
This function performs interpolation of known functions using 3D space
only.  The interpolator is built based on random known nodes created by the
test functions.  The interpolation is then completed for all other nodes
and compared to the true data set defined by the known function.
  Input Parameters:
    interp        -pointer to Interpolator object containing class variables
                   and functions.
    Axes_Dim      -array of integers defining the dimensions of each axis.
    Axes_Ticks    -contains actual axes' values.
    netspace      -pointers to Net objects for storing interpolated results
                   in netcdf files.

  Functions Called:
    TFUN3         -contains the test functions and determines the true value
                   of the current test function for the current coordinates
    QSHEPM        -builds the interpolator from scattered data
    QSMVAL        -completes interpolation for all points.

**************************************************************************/ 

void
TEST_IT(Interpolator *interp,int *Axes_Dim,float **Axes_Ticks,Net **netspace){
  
  int I,J,K,       // loop indices
      KF,          // current test function 
      DIM=3,       // set dimensions to 3
      NP,          // total points in 3D space 
      N,           // total known points to interpolate over in 3D space
      IER,         // result of QSHEPM
      i,j,k;       // loop indices

  float *PP,       // array of coordinate values for interpolation
        *interp_vector; // contains interpolated value

  float FT[Axes_Dim[1]+1][Axes_Dim[2]+1][Axes_Dim[3]+1], //3D array which
                   //  contains true test function values
        RMSER,     //error variable
        ERMEAN,    // mean error of interpolation
        ERMAX,     // max error of interpolation
//        PW,   
        DUM,       // dummy input parameters for function TFUN3
  interp_diff,     // difference between interpolated value and true value
  interp_value;    // interpolated value

  double time1,time2,  // timing variables
         time3,time4;

  long *edges,     //coordinate and edge length arrays for netcdf format
       *coords;

   // initialize vectors
   interp_vector = new float[Axes_Dim[DIM]];
   if(!interp_vector){
     cerr << "Cannot allocate space for int_vector" << endl;
     exit(1);
    }


   N = interp->getN();     // get value of N = known nodes
   PP = new float[DIM+1];
   if(!PP){
     cerr << "Cannot allocate space for PP" << endl;
     exit(1);
    }

  coords = new long[DIM];
   if(!coords){
     cerr << "Cannot allocate space for coords" << endl;
     exit(1);
    }

  edges = new long[DIM];
  for(i=0;i<DIM-1;i++){
    edges[i] = 1;
   }
  edges[DIM-1] = Axes_Dim[DIM];

  for(I=1;I<Axes_Dim[1]+1;I++){
    for(J=1;J<Axes_Dim[2]+1;J++){
      for(K=1;K<Axes_Dim[3]+1;K++){
         FT[I][J][K] = 0;
       }
     }
   }

  int index = 0;
  for(KF=1;KF<NFUN+1;KF++){  // loop over total number of test functions
    cout << endl;
    cout << "                FUNCTION NO. " << KF << endl;
    for(K=1;K<N+1;K++){
      interp->get_node_coords(PP,K);  // get coordinates of scattered nodes
      // get test function value at coordinates in PP for interpolation
      interp->TFUN3(KF,K,PP,0,0,&DUM,DUM,DUM,DUM);
     }
   
    // determine true function values at all coordinates for comparison of
    //   interpolation results
    coords[2] = 0;
    for(I=1;I<Axes_Dim[1]+1;I++){
      coords[0] = I-1;
      PP[1] = Axes_Ticks[1][I];
      for(J=1;J<Axes_Dim[2]+1;J++){
        coords[1] = J-1;
        PP[2] = Axes_Ticks[2][J];
        for(K=1;K<Axes_Dim[3]+1;K++){
           PP[3] = Axes_Ticks[3][K];
           interp->TFUN3(KF,K,PP,0,1,&FT[I][J][K],DUM,DUM,DUM);
           // FT contains true function value
           interp_vector[K-1] = FT[I][J][K];
           // store true value in interp_vector
         }
        // store interp_vector in netcdf file
        netspace[index]->All_Knodes(coords,edges,interp_vector);
       }
     }

    // Now, do actual interpolation over test function's scattered data points
    time1 = timer();
    interp->QSHEPM(IER);  // build interpolator
    time2 = timer();

    // check that QSHEPM ran without errors, exit if there were errors
    if(IER != 0 ){
       if( IER == 2){
         cout << "   *** ERROR IN QSHEPM -- DUPLICATE NODES ENDOUNTERED ***\n";
        }
       else{
        if( IER == 3){
         cout << "   *** ERROR IN QSHEPM -- ALL NODES ARE COPLANAR ***\n";
         }
        }
       interp->~Interpolator();
       exit(1);
     }

    cout << "  PREPROCESSING (LEAST SQUARES)  Time = " << time2-time1 << endl;
    RMSER = 0.0;
    ERMEAN = 0.0;
    ERMAX = 0.0;

    
    time1 = timer();
    
    // perform interpolation, compare to true values, calculate errors
    index++;
    coords[2] = 0;
    for(i=1;i<Axes_Dim[1]+1;i++){
      coords[0] = i-1;
      PP[1] = Axes_Ticks[1][i];
      for(j=1;j<Axes_Dim[2]+1;j++){
        coords[1] = j-1;
        PP[2] = Axes_Ticks[2][j];
        for(k=1;k<Axes_Dim[3]+1;k++){
          PP[3] = Axes_Ticks[3][k];
          interp_vector[k-1] = interp->QSMVAL(PP);  // get interpolated value
          // determine difference between the true value and interpolation
          interp_diff = interp_vector[k-1] - FT[i][j][k];

//          cout << PP[1] << " " << PP[2] << " " << PP[3] << " " << FT[i][j][k] << interp_value << endl;
          // calculate errors
          RMSER += interp_diff * interp_diff;
          ERMEAN += abs(interp_diff);
          ERMAX = MAX0(ERMAX,abs(interp_diff));
         }
         // store interpolated values in netcdf file
         netspace[index]->All_Knodes(coords,edges,interp_vector);
       }
     }
 
    time2 = timer();

    // finish calculating errors
    NP = 8000;  // total 3D space - 20x20x20
    RMSER = (float)sqrt((double)(RMSER/(float)NP));
    ERMEAN /= (float)NP;
    cout << "Interpolation:\n"; 
    cout << "ERMAX = " << ERMAX << endl;
    cout << "ERMEAN =  " << ERMEAN << endl;
    cout << "RMSER = " << RMSER << endl;
    cout << "Interpolation Time:  " << time2-time1 << endl; 

    index++ ;  //index into array of netcdf files

   }

  // clean up
  delete [] PP;
  delete [] interp_vector;
  delete [] edges;
 }  // end of TEST_IT
