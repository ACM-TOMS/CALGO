#ifndef VBF_GF2X__H
#define VBF_GF2X__H

#include <NTL/GF2X.h>

NTL_CLIENT

class vbf_GF2X: public GF2X {

  GF2X str2GF2X(string& str);

}; // end class vbf_GF2X

GF2X str2GF2X(string& str)
{
   vector<string> monoms;
   long i, n, max_exp = 0;
   vec_GF2 x;
   GF2X f;

   Tokenize(str, monoms,"+");
   n = monoms.size();

   for (i = 0; i < n; i++)
   {
      vector<string> values;

      Tokenize(monoms[i], values, "^"); 
      if (values.size() > 1)
      {
         istringstream Sexp(values[1]);
         long exp;

         Sexp >> exp;
         if (exp > max_exp) max_exp = exp; 
      } else {
         if (monoms[i] == "x")
         {
            if (max_exp < 1) max_exp = 1;  
         }
      }
   }

   x.SetLength(max_exp+1);
   clear(x);

   for (i = 0; i < n; i++)
   {
      vector<string> values;

      Tokenize(monoms[i], values, "^");
      if (values.size() > 1)
      {
         istringstream Sexp(values[1]);
         long exp;

         Sexp >> exp;
         x[exp] = 1;
      } else {
         if (monoms[i] == "x")
         {
	    x[1] = 1;
         } else {
	    x[0] = 1;
	 } 
      }
   }  

   f = to_GF2X(x);

   return f;
}

#endif
