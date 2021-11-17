//file: demoMP.cc
//it will run both implementations of MP up
// a max_iter number of components or quiting
//by epsilon. The functionality of the code
//below should be clear once the MP is understood. 

#include <time.h>
#include <fstream>
#include "input_maker.h"
#include <vector.h>
#include "Interval.h"
#include "Gabor.h"
#include "Partition.h"
#include "ShiftGaborMP.h"
#include "FFTGaborMP.h"


int main( int argc, char* argv[])
{
  if( argc != 2 )
    { 
      cout << "Must give parameter m to MPDemo. E.g. type: MPDemo 7\n" 
	   << "Alternatively, if you want results to go to a file with \n "
	   << "a different suffix, say .temp, then type : main 7.temp \n"
	   << "In this case results will go to file results_for_m=7.temp"
	   << endl;
      exit(1);
    }
  
  int m = atoi( argv[1] );
  cout << "m = " << m << endl;
  int dim = 1<<m;
  real epsilon = 1 - 0.99;
  int max_iter = 1000;  //set to different value if desired
  real error, error1;
  char filename[60];
  strcpy( filename, "results");
  strcat(filename, argv[1]);
  ofstream outfile( filename, ios::out);
  assert(outfile);
  
  Interval I(0, dim-1);              // The interval I is used to evaluate the Gabor functions
  for(int i=0; i<dim; i++) I[i] = i; // when making a reconstruction of the input signal
  
  Interval CleanSignal(0, dim-1); // The input signal is loaded here; it has l^2 norm =1.
  real coef, coef1;
  for( int inp = 4; inp < 5; inp++)//may be used to loop over many input
    {                              //functions if desired. The allowable values are 1, 2, ... 10.
      Input_maker(inp, CleanSignal);  // Load in the input signal
      Interval Rf, RecSignal, RecSignal1, Rf1;
      vector <RealGabor> G, G1;
      vector <real> Gcoef, Gcoef1;
      Interval f_approx(0, dim-1), f_approx1(0, dim-1); //This will be the approximation to the
                                                        //original given by MP
      RealGabor Gtemp, Gtemp1;
      real scale = 2.0; // Choose the scale variable scales = a^j
      Partition Part(dim, scale);

      cout << endl << "////////////// Analysis of function number " << inp << "///////////" 
	   << endl << endl;

      cout << "******************  For Shift Gabor ***********************" << endl << endl;

      clock_t timer = clock();
      error = RunShiftGaborMP( max_iter, epsilon, CleanSignal, Part, RecSignal, Rf, G, Gcoef);
      timer = clock() - timer;
      timer /= CLOCKS_PER_SEC;
      cout << endl << "time to execute Shift Matching Pursuit = " 
	   << timer << " seconds (approximately) " << endl;
      cout <<"for " << G.size()
	   << " iteration(s) of matching pursuit" << endl << endl;


      cout << "error for Shift, scale = "<< scale << ": " << error << endl << endl;


      cout << "******************  For FFT Gabor ***********************" << endl << endl;
      timer = clock();
      error1 = RunFFTGaborMP(max_iter, epsilon, CleanSignal, RecSignal1,Rf1, G1, Gcoef1);
      timer = clock() - timer;
      timer /= CLOCKS_PER_SEC;
      cout << endl << "time to execute FFT Matching Pursuit= " 
	   << timer << " seconds (approximately) " << endl;
      cout <<"for " << G1.size()
	   << " iteration(s) of matching pursuit" << endl << endl;
      cout << "error for FFT" << " : " << error1 << endl << endl;      
      
      outfile << "********* INPUT_MAKER " << inp << " *********" << endl;
      outfile << "number of iterations performed for Shift " << G.size() << "  error = " << error << endl;
      outfile << "number of iterations performed for FFT " << G1.size() << "  error = " << error << endl;
      
            
      for(unsigned int z = 0; z < G.size(); z++)
	//Here the synthesis is done one Gabor at a time to facilitate the 
	//step by step viewing of the reconstruction for Shift MP.
	{
	  Gtemp = G[z];
	  coef = Gcoef[z];
	  outfile << "iteration " << z << endl;
	  outfile << "\t\tcoeficient for shift: "<< coef <<endl;
	  outfile << "\t(" << Gtemp.s << ", " << Gtemp.u << ", "
		  << Gtemp.v << ", " << Gtemp.w << ")" << endl;
	  Gtemp.createSample(I);
	  f_approx += Gtemp.Sample * coef;
	  
	}  // end of z-loop

      for(unsigned int z = 0; z < G1.size(); z++)
	//Here the synthesis is done one Gabor at a time to facilitate the 
	//step by step viewing of the reconstruction for FFT MP.

	{
	  Gtemp1 = G1[z];
	  coef1 = Gcoef1[z];
	  outfile << "\t\tcoeficient for FFT: "<< coef1 <<endl;
	  outfile << "\t(" << Gtemp1.s << ", " << Gtemp1.u << ", "
		  << Gtemp1.v << ", " << Gtemp1.w << ")" << endl;
	  Gtemp1.createSample(I);
	  f_approx1 += Gtemp1.Sample * coef1;
	  
	}  // end of z-loop
      // Here we may output the point values of the reconstructed function to the screen
      //or to a file. See "Interval.h" for overloading of << and >> operators.
      cout << "Point values of Shift MP reconstruction " << endl;
      cout << Rf << endl;
      cout << "Point values of FFT MP reconstruction " << endl;
      cout << Rf1 << endl;
    }
      
  outfile << "\n\n\n";   
  return 0;
}

