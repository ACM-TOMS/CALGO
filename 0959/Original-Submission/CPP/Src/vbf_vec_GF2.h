#ifndef VBF_vec_GF2__H
#define VBF_vec_GF2__H

#include <NTL/vec_GF2.h>
#include "vbf_ZZ.h"

NTL_CLIENT

class vbf_vec_GF2: public vec_GF2 {

  int degbf(const vec_GF2& X, int n);
  // Calculate the degree of the Boolean function
  // associated with the vector vec_GF2

  long conv_long(const vec_GF2& a);
  // Calculate the long with binary representation a

  vec_GF2 to_vecGF2(long x, int n);
  //binary representation of x with a n-bit vector

  void opposite(vec_GF2&x, const vec_GF2& a);
  // x = negated of a

}; // end class vbf_vec_GF2

int degbf(const vec_GF2& X, int n)
{
   int 	i, j, num;
   long *v = NULL;
      	
   for (i = n; i > 0; i--)
   {
      num = to_long(numofweight(n,i));
      v = (long *) malloc(num * sizeof(long));
      vectors_weight(v, n, i);

      for (j = 0; j < num; j++) 
      { 
        if (!IsZero(X[v[j]])) return i;
      }

   }  	
   return 0;	
}

long conv_long(const vec_GF2& a)
{
   long i, temp, total = 0, n = a.length();

   for (i = 0; i < n; i++)
   {
      if (a[i] == 1)
      { 
         temp = 1;
      }
      else
      {
	 temp = 0;
      }
      total += temp << (n-i-1);
   }
  
   return total;
}

vec_GF2 to_vecGF2(long x, int n)
{
   int 		i;
   vec_GF2 	bin;

   bin.SetLength(n);   	
   for (i = 0; i < n; i++)  
   {
      bin[n-1-i] = bit(x,i);
   }
   return bin;
} 

void opposite(vec_GF2&x, const vec_GF2& a)
{
   long i, n = a.length();

   x.SetLength(n);
   for (i = 0; i < n; i++)
   {
      x[i] = a[i]+1;
   }
  
}

#endif
