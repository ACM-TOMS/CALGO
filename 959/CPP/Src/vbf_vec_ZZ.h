#ifndef VBF_vec_ZZ__H
#define VBF_vec_ZZ__H

#include <NTL/vec_ZZ.h>

NTL_CLIENT

class vbf_vec_ZZ: public vec_ZZ {

  long IsImpulse(const vec_ZZ& a);

  int IsVecCharFunct(const vec_ZZ& a, long& j);

  void convolution(vec_ZZ& x, const vec_ZZ& a, const vec_ZZ& b);

  int IsConstant(const vec_ZZ& a);

  ZZ maxvalue(const vec_ZZ& a);

  ZZ maxvalue_abs(const vec_ZZ& a);

  vec_ZZ to_vecZZ(long x, int n);

  ZZ conv_ZZ(const vec_GF2& a);

  void swap(vec_ZZ& x, const vec_ZZ& a, long i, long j);

  void reverse(vec_ZZ& x, const vec_ZZ& a);

  void exchange(vec_ZZ& x, const vec_ZZ& a, long size);

}; // end class vbf_vec_ZZ

long IsImpulse(const vec_ZZ& a)
{
   long n = a.length();
   long i, cont = 0;

   for (i = 0; i < n; i++)
      if (!IsZero(a[i]) && (cont == 0)) {
         cont = 1;
      } else if (!IsZero(a[i]) && (cont == 1)) {
         return 0;
      }

   return 1;
}

int IsVecCharFunct(const vec_ZZ& a, long& j)
{
   long len = a.length();
   long i, numones = 0;
   vec_ZZ c;

   for (i = 0; i < len; i++) 
   { 
      if (a[i] == 1)
      {
      	 numones++;
      	 j = i;
      }	 
   }
   
   if (numones == 1)
   {
      vec_ZZ b(INIT_SIZE, len);
      b[j]=1;
      c = a-b;
      if (IsZero(c))
      {
         return 1;
      }
      else
      {
      	 return 0;
      }	    
   }
   else
   {
      return 0;
   }      
}

void convolution(vec_ZZ& x, const vec_ZZ& a, const vec_ZZ& b)
{
   long n = a.length();
   if (b.length() != n) Error("vector convolution: dimension mismatch");

   x.SetLength(n);

   long i,j,k;
   ZZ temp;
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

int IsConstant(const vec_ZZ& a)
{
   long n = a.length();
   long i;
   
   for (i = 0; i < n-1; i++) 
      if (a[i] != a[i+1]) return 0;

   return 1;
}

ZZ maxvalue(const vec_ZZ& a)
{
   ZZ   max;
   long i;
   long n = a.length();

   max = 0;
   for (i = 0; i < n; i++)
   {
       if (a[i] > max) max = a[i];
   }  

   return max;
}

ZZ maxvalue_abs(const vec_ZZ& a)
{
   ZZ   max,temp;
   long i;
   long n = a.length();

   max = 0;
   for (i = 0; i < n; i++)
   {
	abs(temp,a[i]);
        if (temp > max) max = temp;
   }

   return max;
}

vec_ZZ to_vecZZ(long x, int n)
{
   int 		i;
   vec_ZZ 	bin;

   bin.SetLength(n);   	
   for (i = 0; i < n; i++)  
   {
      bin[n-1-i] = bit(x,i);
   }
   return bin;
} 

ZZ conv_ZZ(const vec_GF2& a)
{
   long i, n = a.length();
   ZZ total;

   total = 0;
   for (i = 0; i < n; i++)
   {
      if (a[i] == 1)
      { 
         total += power2_ZZ(n-i-1);
      }      
   }
  
   return total;
}

void swap(vec_ZZ& x, const vec_ZZ& a, long i, long j)
{
   long n = a.length();
   x.SetLength(n);
   VectorCopy(x,a,n); 
   x(i) = a(j);
   x(j) = a(i);
}

void reverse(vec_ZZ& x, const vec_ZZ& a)
{
	 long i;
   long n = a.length();
   x.SetLength(n);
   for (i = 0; i < n; i++)
   {
      x[i] = a[n-1-i];
   }
}

void exchange(vec_ZZ& x, const vec_ZZ& a, long size)
{
	 long i, j, step;
   long n = a.length();
   x.SetLength(n);
   step = n/size;
   
   for (i = 0; i < step; i=i+2)
   {
   	  for (j = 0; j < size; j++)
   	  {
         x[i*size+j] = a[(i+1)*size+j];
      }
   	  for (j = 0; j < size; j++)
   	  {
         x[(i+1)*size+j] = a[i*size+j];
      }
   }
}

#endif
