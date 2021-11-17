#ifndef VBF_ZZ__H
#define VBF_ZZ__H

#include <NTL/ZZ.h>

NTL_CLIENT

class vbf_ZZ: public ZZ {

  int logtwo(long x);

  ZZ numofweight(int n, int w);
  // computes the number of n-bit vectors with weight w

  void vectors_weight(long *x, int n, int w);
  // Return an array of long with a n-bit representation and weight w

}; // end class vbf_ZZ

int logtwo(long x)
{
   int n = 0;
   long tmp;

   while (x > 1)
   {
      tmp = x >> 1;
      n++;
      x = tmp;
   }
   return n;
}

ZZ numofweight(int n, int w)
{
   int 	    i;
   NTL::ZZ  x, y, z;
   
   x = 1;
   for (i = n; i > (n-w); i--)
   {
      x = i*x;
   }		
 
   y = 1;
   for (i = 2; i <= w; i++)
   {
      y = i*y;
   }		
     
   z = x/y;

   return z;
}

void vectors_weight(long *x, int n, int w)
{
   int 	i, j, k, l;
   int  length = n-w+1;
   long cont = 0;
   long	**a = NULL; 
   long	**b = NULL;  
   long	*la = NULL;    
   long	*lb = NULL;
      	
   a = (long **) malloc(length * sizeof(long *));
   la = (long *) malloc(length * sizeof(long));
   b = (long **) malloc(length * sizeof(long *));
   lb = (long *) malloc(length * sizeof(long));

   for (i = 0; i < length; i++) 
   { 
      a[i] = (long *) malloc(sizeof(long));
      la[i] = 1;
      a[i][0] =  (1 << (i+w-1));
   }
 
   for (i = w-1; i > 0; i--) 
   { 
      for (j = 0; j < length; j++) 
      { 
      	 lb[j] = 0;
      	 b[j] = NULL;
      	 for (k = length-1; k >= j; k--)
      	 {
      	    lb[j] += la[k];
      	 }   
      }
      for (j = 0; j < length; j++) 
      { 
      	 long cont = 0;
      	 b[j] = (long *) malloc(lb[j] * sizeof(long));
      	 for (k = length-1; k >= j; k--)
      	 {
      	    for (l = 0; l < la[k]; l++)
      	    {	
      	       b[j][cont] = a[k][l] + (1 << (j+i-1));
      	       cont++;
      	    }  
      	 }   
      }
      for (j = 0; j < length; j++) 
      { 
      	 la[j] = lb[j];
      	 a[j] = NULL;
      	 a[j] = (long *) malloc(la[j] * sizeof(long));
      	 for (k = 0; k < la[j]; k++)
      	 {
      	    a[j][k] = b[j][k];
      	 }   
      }
   }
   
   for (i = 0; i < length; i++) 
   { 
      for (j = 0; j < la[i]; j++)
      {
      	 x[cont] = a[i][j];
      	 cont++;
      }   
   }
    	     
}   

#endif
