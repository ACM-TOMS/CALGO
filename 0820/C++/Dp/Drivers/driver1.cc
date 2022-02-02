/*  This main tests RunShiftGabor and RunFFTGabor
    Last Modified December 10, 2000
    Partition for RunShiftGabor not exactly the same as for RunFFTGabor     
*/
#include <vector.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include "Gabor.h"
#include "ShiftGaborMP.h"
#include "Partition.h"
#include "FFTGaborMP.h"
#include "input_maker.h"
#include <cctype>
#include <cstdlib>

#define STRSIZE 80
char pauseABit[1];  /* variable to hold user input */

int main()
{
  int m = 6;
  int dim = 1<<m;
  real epsilon = 0.0000000000000001;
  int max_iter = 4;
  real error, error1;
  double biggest, biggest_in_case, parameter_biggest;
  double FFT_largest[2], Run_largest[2];
  double threshold = 0.00000001;
 
  ifstream infile("data1");
  
  if(!infile) 
    { 
      cout << "Failed to open data1 file.\n"; 
      return 1; 
    }
  
  Interval I(0, dim-1);
 
  for(int i=0; i<dim; i++) I[i] = i;
  real coef, coef1;
  char str[255];    /* buffer to walk over irrelevant text in the data file */

  biggest = -1;
  parameter_biggest = -1;
  biggest_in_case = -1;
  FFT_largest[1] = -1;
  Run_largest[1] = -1;

  cout << endl << "*************** Start of RunShiftGaborMP tests ***************" << endl << endl;

  /* the RunShiftGaborMP test cases */
  for( int inp = 1 ; inp <=9 ; inp++)
    {   
      Interval Rf, RecSignal;
      vector <RealGabor> G;
      vector <real> Gcoef;
      Interval f_approx(0, dim-1);
      RealGabor Gtemp;
      Partition Part( dim, (2.0) );
      real *Ip1, *Ip2;
      real input1, input2, input3, input4;   /* holds the input values from the data file */
      double difference1, difference2, difference3, difference4;  /* holds the absolute difference */
     
      Interval CleanSignal(0, dim-1);
      Input_maker(inp, CleanSignal);   /* select a case to test */
     
      
      cout << "Case " << inp << " for RunShiftGaborMP" << endl; 
      error = RunShiftGaborMP(max_iter, epsilon, CleanSignal, Part,
			       RecSignal, Rf, G, Gcoef);
      
      infile.precision(16);
      cout.precision(16);
    
      /* passover irrelevant text in data file */
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE, '\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE, '\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE, ':'); 
      
      infile >> input1;  /* get input value */
 
      difference1 = fabs(input1-error); 
      
      if (difference1 > threshold)
	{
	  cout << "Case " << inp << " of error tolerance for RunShiftGaborMP failed." << endl;
	  return 1;
	}
      if (difference1 > biggest)
	biggest = difference1;
      if (difference1 > biggest_in_case)
	biggest_in_case = difference1;

      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE, '\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE, '\n');
    
      Ip1 = RecSignal.origin;
      Ip2 = Rf.origin;
      for (int i=0; i<RecSignal.length; i++)
	{
	  infile >> input1 >> input2;
	  difference1 = fabs(input1-Ip1[i]);
	  difference2 = fabs(input2-Ip2[i]);

	  if (difference1 > threshold || difference2 > threshold)
	    {
	      cout << "Case " << inp << " of RecSignal and Rf for RunShiftGaborMP failed." << endl;
	      return 1;
	    }
	  if (difference1 > biggest)
	    biggest = difference1;
	  if (difference1 > biggest_in_case)
	    biggest_in_case = difference1;
	  if (difference2 > biggest)
	    biggest = difference2;
	  if (difference2 > biggest_in_case)
	    biggest_in_case = difference2;
	 }
	 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n');
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
  
      for (unsigned int i=0; i<G.size(); i++)
	{
	  infile >> input1 >> input2 >> input3 >> input4;
	  difference1 = fabs(input1-G[i].s);
	  difference2 = fabs(input2-G[i].u);
	  difference3 = fabs(input3-G[i].v);
	  difference4 = fabs(input4-G[i].w);

	  if (difference1 > threshold || difference2 > threshold || difference3 > threshold || difference4 > threshold)
	    {
	      cout << "Case " << inp << " of Gabor G1 for RunFFTGaborMP, minor parameter variation." << endl;
	    }
	}
    
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
 
      for (unsigned int i=0; i<G.size(); i++)
	{
	  memset(str, '\0', STRSIZE); 
	  infile.getline(str, STRSIZE,'\n'); 
	  memset(str, '\0', STRSIZE); 
	  infile.getline(str, STRSIZE,'\n'); 
	  memset(str, '\0', STRSIZE); 
	  infile.getline(str, STRSIZE,'\n'); 
	  Ip1 = G[i].Sample.origin;
	  
	  for (int j=0; j<G[i].Sample.length; j++)
	    {
	      infile >> input1;
	      difference1 = fabs(input1-Ip1[j]);

	      if (difference1 > threshold)
		{
		  cout << "Case " << inp << " of Gabor G Sample member for RunShiftGaborMP failed." << endl;
		  return 1;
		}
	      if (difference1 > biggest)
		biggest = difference1;
	      if (difference1 > biggest_in_case)
		biggest_in_case = difference1;
	    }
	 } 
   
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 

      for (unsigned int i=0; i<Gcoef.size(); i++)
	{
	  infile >> input1;
	  difference1 = fabs(Gcoef[i]-input1);

	  if (difference1 > threshold)
	    {
	      cout << "Case " << inp << " of coefficient for RunShiftGaborMP failed." << endl;
	      return 1;
	    }
	  if (difference1 > biggest)
	    biggest = difference1;
	  if (difference1 > biggest_in_case)
	    biggest_in_case = difference1;
	}
    
      for(unsigned int z=0; z < G.size(); z++)
	{
	  Gtemp = G[z];
	  coef = Gcoef[z];
	  f_approx += Gtemp.Sample * coef;
	} 

      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
  
      Ip1 = f_approx.origin;
      for (int i=0; i<f_approx.length; i++)
	{
	  infile >> input1;
	  difference1 = fabs(input1-Ip1[i]);

	  if (difference1 > threshold)
	    {
	      cout << "Case " << inp << " of Interval f_approx for RunShiftGaborMP failed." << endl;
	      return 1;
	    }
	  if (difference1 > biggest)
	    biggest = difference1;
	  if (difference1 > biggest_in_case)
	    biggest_in_case = difference1;
	}
    
      if (biggest > Run_largest[1])
	{
	  Run_largest[0] = inp;
	  Run_largest[1] = biggest;
	}
      cout << "Case " << inp << " passed." << endl;
      cout << "The largest deviated value for this case is: " << biggest_in_case << endl;
      biggest_in_case = -1;
//      cout << endl << "Press 'y' to continue, Control-C to quit: ";
//      cin >> pauseABit;
        memset(str, '\0', STRSIZE); 
        infile.getline(str, STRSIZE,'\n');
        memset(str, '\0', STRSIZE); 
        infile.getline(str, STRSIZE,'\n');
//      cout << endl;
    } /* end of RunShiftGaborMP loop */
   
  cout << "All cases for RunShiftGaborMP passed" << endl;
  cout << "The largest deviated value came from case " << Run_largest[0] << ": " << Run_largest[1] << endl;  

//  cout << endl << "Press 'y' to continue, Control-C to quit: ";
//  cin >> pauseABit;

  /* the RunFFTGaborMP test cases */

  biggest_in_case = -1;
  cout << endl << "*************** Start of RunFFTGaborMP tests ***************" << endl << endl;

  for( int inp = 1 ; inp <= 9; inp++)
    {      
      Interval RecSignal1, Rf1;
      vector <RealGabor> G1;
      vector <real> Gcoef1;
      Interval f_approx1(0, dim-1);
      RealGabor Gtemp1;
      Partition Part( dim, (2.0) );
      real *Ip1, *Ip2;
      real input1, input2, input3, input4;
      double difference1, difference2, difference3, difference4;

      Interval CleanSignal(0, dim-1);
      Input_maker(inp, CleanSignal);  

      cout << "Case " << inp << " for RunFFTGaborMP" << endl; 
      error1 = RunFFTGaborMP(max_iter, epsilon, CleanSignal,
			     RecSignal1,Rf1, G1, Gcoef1);
     
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,':'); 
      infile >> input1;

      difference1 = fabs(input1-error1);
      if (difference1 > threshold)
	{
	  cout << "Case " << inp << " of error tolerance for RunFFTGaborMP failed." << endl;
	  return 1;
	}
      if (difference1 > biggest)
	biggest = difference1;
      if (difference1 > biggest_in_case)
	biggest_in_case = difference1;
     
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
    
      Ip1 = RecSignal1.origin;
      Ip2 = Rf1.origin;
      for (int i=0; i<RecSignal1.length; i++)
	{
	  
	  infile >> input1 >> input2;
	  difference1 = fabs(input1-Ip1[i]);
	  difference2 = fabs(input2-Ip2[i]);
	   if (difference1 > threshold || difference2 > threshold)
	    {
	      cout << "Case " << inp << " of RecSignal1 and Rf1 for RunFFTGaborMP failed." << endl;
	      return 1;
	    }
	  if (difference1 > biggest)
	    biggest = difference1;
	  if (difference1 > biggest_in_case)
	    biggest_in_case = difference1;
	  if (difference2 > biggest)
	    biggest = difference2;
	  if (difference2 > biggest_in_case)
	    biggest_in_case = difference2;
	  }
	
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      
      for (unsigned int i=0; i<G1.size(); i++) // Check for possible parameter variations
	{
	  infile >> input1 >> input2 >> input3 >> input4;
	  difference1 = fabs(input1-G1[i].s);
	  difference2 = fabs(input2-G1[i].u);
	  difference3 = fabs(input3-G1[i].v);
	  difference4 = fabs(input4-G1[i].w);
	  if (difference1 > threshold || difference2 > threshold || difference3 > threshold || difference4 > threshold)
	    {
	      cout << "Case " << inp << " of Gabor G1 for RunFFTGaborMP, minor parameter variation." << endl;
	    }
	}
    
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      
      for (unsigned int i=0; i<G1.size(); i++)
	{
	  memset(str, '\0', STRSIZE); 
	  infile.getline(str, STRSIZE,'\n'); 
	  memset(str, '\0', STRSIZE); 
	  infile.getline(str, STRSIZE,'\n'); 
	  memset(str, '\0', STRSIZE); 
	  infile.getline(str, STRSIZE,'\n'); 
	  Ip1 = G1[i].Sample.origin;
	 
	  for (int j=0; j<G1[i].Sample.length; j++)
	    {
	      infile >> input1;
	      difference1 = fabs(input1-Ip1[j]);
	      if (difference1 > threshold)
		{
		  cout << "Case " << inp << " of Gabor G1 Sample member for RunFFTGaborMP failed." << endl;
		  return 1;
		}
	      if (difference1 > biggest)
		biggest = difference1;
	      if (difference1 > biggest_in_case)
		biggest_in_case = difference1;
	    }
	}
 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      for (unsigned int i=0; i<Gcoef1.size(); i++)
	{
	  infile >> input1; 
	  difference1 = fabs(input1-Gcoef1[i]);
	  if (difference1 > threshold)
	    {
	      cout << "Case " << inp << " of coefficient for RunFFTGaborMP failed." << endl;
	      return 1;
	    }
	  if (difference1 > biggest)
	    biggest = difference1;
	  if (difference1 > biggest_in_case)
	    biggest_in_case = difference1;
	}
      
      for(unsigned int z=0; z < G1.size(); z++)
	{
	  Gtemp1 = G1[z];
	  coef1 = Gcoef1[z];
	  f_approx1 += Gtemp1.Sample * coef1;
	}
      
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 
      memset(str, '\0', STRSIZE); 
      infile.getline(str, STRSIZE,'\n'); 

      Ip1 = f_approx1.origin;
      for (int i=0; i<f_approx1.length; i++)
	{
	  infile >> input1;
	  difference1 = fabs(input1-Ip1[i]);
	  if (difference1 > threshold)
	    {
	      cout << "Case " << inp << " of Interval f_approx for RunShiftGaborMP failed." << endl;
	      return 1;
	    }
	  if (difference1 > biggest)
	    biggest = difference1;
	  if (difference1 > biggest_in_case)
	    biggest_in_case = difference1;
	}

      if (biggest > FFT_largest[1])
	{
	  FFT_largest[0] = inp;
	  FFT_largest[1] = biggest;
	}
      cout << "Case " << inp << " passed." << endl;
      cout << "The largest deviated value for this case is: " << biggest_in_case << endl;
      biggest_in_case = -1;
///    cout << endl << "Press 'y' to continue, Control-C to quit: ";
///    cin >> pauseABit;
///    cout << endl;

    } /* end of RunFFTGaborMP loop */

  cout << "All cases for RunFFTGaborMP passed" << endl;
  cout << "The largest deviated value came from case " << FFT_largest[0] << ": " << FFT_largest[1] << endl;

  cout << endl << "**************************** SUMMARY ****************************" << endl << endl;
  cout << "The largest deviated value for all cases from RunShiftGaborMP is: " << endl;
  cout << Run_largest[1] << " from case: " << Run_largest[0] << endl << endl;
  cout << "The largest deviated value for all cases from RunFFTGaborMP is: " << endl;
  cout << FFT_largest[1] << " from case: " << FFT_largest[0] << endl;
  cout << endl << "All test cases passed." << endl << endl;

}
   

