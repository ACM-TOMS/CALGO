#ifndef VBF_GF2EX__H
#define VBF_GF2EX__H

#include <NTL/GF2EX.h>
#include "vbf_vec_GF2E.h"

NTL_CLIENT

class vbf_GF2EX: public GF2EX {

  GF2EX str2GF2EX(string& str, const long& n);
  // convert string str to GF2EX

}; // end class vbf_GF2EX

GF2EX str2GF2EX(string& str, const long& n)
{
   vector<string> monoms;
   unsigned long i, spacen;
   mat_GF2 A,B;
   vec_GF2E x;
   GF2EX f;

   Tokenize(str, monoms,"+");
   spacen = (1 << n);
   A.SetDims(spacen,n);

   for (i = 0; i < monoms.size(); i++)
   {
      vector<string> values;
      Tokenize(monoms[i], values, "x^");

      if (values.size() == 1)
      {
         NTL_SNS string::size_type pos = values[0].find_first_of("x",0); 
         if (NTL_SNS string::npos != pos )
	 {
            istringstream Scoef(values[0].substr(0,values[0].length()-2));
            long coef;

            Scoef >> std::hex >> coef >> std::dec;
            A[1] = to_vecGF2(coef,n);
	 } else {
            istringstream Scoef(values[0]);
            long coef;

            Scoef >> std::hex >> coef >> std::dec;
	    A[0] = to_vecGF2(coef,n);
         }
      } else 
      {
         istringstream Scoef(values[0]);
         istringstream Sexp(values[1]);
         long coef, exp;

	 Sexp >> exp;
         Scoef >> std::hex >> coef >> std::dec;
         A[exp] = to_vecGF2(coef,n); 
      }
   }

   reverse(B,A);
   conv(x,B);
   f = to_GF2EX(x);

   return f;
}

#endif
