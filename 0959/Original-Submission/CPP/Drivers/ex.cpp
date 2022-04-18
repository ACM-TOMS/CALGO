/*   ex		          */

#include <iostream>
#include <fstream>
#include "VBF.h"
  

/******************************************************************/
int main(int argc, char *argv[]) 
{
   using namespace VBFNS;
   
   VBF		F;
   NTL::vec_long vec_F;
   NTL::vec_ZZ	 c;
   NTL::mat_GF2 A, T;
   NTL::mat_ZZ  W, LP, DP;
   NTL::mat_ZZ  Ac;
   long a;
   int  i, n;
   char file[33];

   // Load VBF definitions
   
   sprintf(file,"%s.dec",argv[1]); 
   ifstream input(file);
   if(!input)
   {
      cerr << "Error opening " << file << endl;
      return 0;
   }
   input >> vec_F;
   n = atoi(argv[2]);
   F.putDecTT(vec_F,n);
   input.close();

   sprintf(file,"%s.anf",argv[1]);
   ofstream output(file);
   if(!output)
   {
      cerr << "Error opening " << file << endl;
      return 0;
   }

   A = ANF(F); 
   cout << "Argument Dimension = " << F.n() << endl;
   cout << "Argument space has " << F.spacen() << " elements."<< endl;
   cout << "Image Dimension = " << F.m() << endl;
   cout << "Image space has " << F.spacem() << " elements." << endl << endl;
   cout << "Writing Algebraic Normal Form to file: " << file << endl;
   cout << "[Columns = Image components]" << endl;
   output << A << endl;
   output.close();

   sprintf(file,"%s.tt",argv[1]);
   ofstream output1(file);
   if(!output1)
   {
      cerr << "Error opening " << file << endl;
      return 0;
   }

   T = TT(F);   
   cout << endl << "Writing Truth Table to file: " << file << endl;
   cout << "[Columns = Image components]" << endl;
   output1 << T << endl;
   output1.close();

   sprintf(file,"%s.wal",argv[1]);
   ofstream output2(file);
   if(!output2)
   {
      cerr << "Error opening " << file << endl;
      return 0;
   }

   W = Walsh(F); 
   cout << endl << "Writing Walsh Spectrum to file: " << file <<endl;
   output2 << W << endl;
   output2.close();

   sprintf(file,"%s.lp",argv[1]);
   ofstream output3(file);
   if(!output3)
   {
      cerr << "Error opening " << file << endl; 
      return 0; 
   }
   
   LP = LAT(F);
   cout << endl << "Writing Linear Profile to file: " << file << endl;
   cout << "[To normalize divide by " << LP[0][0] << "]" << endl;
   output3 << LP << endl;
   output3.close();

   sprintf(file,"%s.dp",argv[1]);
   ofstream output4(file);
   if(!output4)
   {
      cerr << "Error opening " << file << endl;
      return 0;
   }
   
   DP = DAT(F);
   cout << endl << "Writing Differential Profile to file: " << file << endl;
   cout << "[To normalize divide by " << DP[0][0] << "]" << endl;
   output4 << DP << endl;
   output4.close();

   sprintf(file,"%s.pol",argv[1]);
   ofstream output5(file);
   if(!output5)
   {
      cerr << "Error opening " << file << endl;
      return 0;
   }

   cout << endl << "Writing the polynomials in ANF to file: " << file << endl;
   Pol(output5,F);
   output5.close();

   sprintf(file,"%s.ls",argv[1]);
   ofstream output6(file);
   if(!output6)
   {
      cerr << "Error opening " << file << endl;
      return 0;   
   }

   A = LS(F);
   cout << endl << "Writing Linear structures to file: " << file << endl;
   output6 << A << endl; 
   output6.close();

   sprintf(file,"%s.ac",argv[1]);
   ofstream output7(file);
   if(!output7)
   {
      cerr << "Error opening " << file << endl;
      return 0;   
   }

   Ac = AC(F);
   cout << endl << "Writing Autocorrelation Spectrum to file: " << file << endl;
   output7 << Ac << endl;
   output7.close();

   sprintf(file,"%s.cy",argv[1]);
   ofstream output8(file);
   if(!output8)
   {
      cerr << "Error opening " << file << endl;
      return 0;   
   }

   c = Cycle(F);
   cout << endl << "Writing Cycle Structure to file: " << file << endl;
   for (i = 0; i < c.length(); i++)
   {
      if (c[i] > 0)
      {
	output8 << i << "," << c[i] << endl;
      }
   }
   output8.close();

   cout << endl <<  "Nonlinearity: " << nl(F) << endl;
   nlr(a,F,2);
   cout << "Second order Nonlinearity: " << a << endl;
   cout << "Linearity distance: " << ld(F) << endl;
   cout << "Algebraic degree: " << deg(F) << endl;
   cout << "Algebraic immunity: " << AI(F) << endl;
   cout << "Absolute indicator: " << maxAC(F) << endl;
   cout << "Sum-of-squares indicator: " << sigma(F) << endl;
   cout << "Linear potential: " << lp(F) << endl;
   cout << "Differential potential: " << dp(F) << endl;
   cout << "Maximum Nonlinearity (if n is even): " << nlmax(F) << endl;
   cout << "Maximum Linearity distance: " << ldmax(F) << endl;

   int type;
   typenl(type, F);

   if (type == BENT)
   {
     cout << "It is a bent function" << endl;
   } else if (type == ALMOST_OPTIMAL)
   {
     cout << "It is an almost optimal function" << endl;
   } else if (type == LINEAR)
   {
     cout << "It is a linear function" << endl;
   }

   cout << "The fixed points are: " << endl;
   cout << fixedpoints(F) << endl;
   cout << "The negated fixed points are: " << endl;
   cout << negatedfixedpoints(F) << endl;
   cout << "Correlation immunity: " << CI(F) << endl;
   if (F.getbal())
   {
     cout << "It is a balanced function" << endl;
   } else
   {
     cout << "It is a non-balanced function" << endl;
   }
   cout << "The function is PC of degree " << PC(F) << endl;
 
/* Finish **********************************************************/

  return 0;
}
//


