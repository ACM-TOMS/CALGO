#ifndef VBF_vec_RR__H
#define VBF_vec_RR__H

#include <NTL/vec_RR.h>

NTL_CLIENT

class vbf_vec_RR: public vec_RR {

void convolution(vec_RR& x, const vec_RR& a, const vec_RR& b);

}; // end class vbf_vec_RR

void convolution(vec_RR& x, const vec_RR& a, const vec_RR& b)
{
   long n = a.length();
   if (b.length() != n) Error("vector convolution: dimension mismatch");

   x.SetLength(n);

   long i,j,k;
   RR temp;
   vec_GF2 vi, vj, v;

   for (i = 0; i < n; i++)
   {
      temp = 0;
      vi = to_vecGF2(i,n);
      for (j = 0; j < n; j++)
      {
         vj = to_vecGF2(j,n);
         v = vi+vj;
         k = conv_long(v);
         temp += a[j] * b[k];
      }
      div(x[i],temp,n);
   }
}

#endif
