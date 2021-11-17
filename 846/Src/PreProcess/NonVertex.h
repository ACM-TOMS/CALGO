#ifndef _Non_Vertex_
#define _Non_Vertex_

double* sb;  //Global variable for RSimplex

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <numeric>
#include <fstream>
#include <iostream>

#include "LowerTriangular.h"
#include "RSimplex.h"

void NonVertex
   (
   int &nVar, int &nSpt, int* Spt1stIdx, int** Spt, int* ynVtx
   )

/*
   Purpose:
   ========

   To detect the non-vertex points of each support
*/

{
   sb = new double [nVar+2];

   int i, i1, info, itmp, j, jStrt, jEnd, k, PtIn, rnk;
   double tmp;

   double** A = new double* [Spt1stIdx[nSpt]+1];
   A[0] = new double [(Spt1stIdx[nSpt]+1)*(nVar+2)];
   for (i=1; i<=Spt1stIdx[nSpt]; i++) A[i] = A[i-1] + nVar+2;

   double* b = new double [Spt1stIdx[nSpt]+1];

   int* ib = new int [nVar+2];
   int* ib0 = new int [nVar+2];
   double* xb = new double [nVar+2];

   double** binv = new double* [nVar+2];
   binv[0] = new double [(nVar+2)*(nVar+2)];
   for (i=1; i<nVar+2; i++) binv[i] = binv[i-1] + nVar+2;

   for (i=0; i<Spt1stIdx[nSpt]; i++) ynVtx[i] = 1;
   // = 0: non-vertex point
 
   // To form the matrix for simplex method

   // To choose a random point
   srand( unsigned(time(0)) );
   for (i=0; i<nVar; i++) b[i] = double(rand())/RAND_MAX;
 
   // To lift the points
   for (j=0; j<Spt1stIdx[nSpt]; j++)
   {
      tmp = 0.0e0;
      for (i1=0; i1<nVar; i1++) tmp=tmp+(Spt[j][i1]-b[i1])*(Spt[j][i1]-b[i1]);
      A[j+1][0] = -sqrt(tmp);    // j+1: shift a slot for the added 1st column
      for (i1=0; i1<nVar; i1++) A[j+1][i1+1] = Spt[j][i1];
      A[j+1][nVar+1] = 1.0e0;
   }
 
   for (i=0; i<nSpt; i++)  // for each support
   {
      jStrt = Spt1stIdx[i];  // The first row is the lifting direction
      jEnd = Spt1stIdx[i+1];  // +1: shift a slot for the added 1st column
 
      A[jStrt][0] = 1.0e0;
      for (i1=1; i1<nVar+2; i1++) A[jStrt][i1] = 0.0e0;
 
      rnk = 0;
      // To lower-trianglize the matrix A
      LowerTriangular(A,jStrt,jEnd,nVar+2,rnk,ib0);

      for (j=0; j<rnk; j++)
      {
	 for (i1=0; i1<rnk; i1++) binv[j][i1] = 0.0e0;
         binv[j][j] = 1.0e0;
      }
 
      for (j=jStrt+1; j<=jEnd; j++) // for each point
      {
         PtIn = j;
         for (i1=0; i1<rnk; i1++) ib[i1] = -1;
         for (i1=0; i1<rnk; i1++)
         {
            if (ib0[i1] == j)
            {
               ib[i1] = PtIn;
               PtIn = -PtIn;
               break;
	    }
         }			
         if (PtIn >= 0)
	 {
            // The j-th row will get into the basis.
            // To decide which column will get out of the basis
            for (k=0; k<rnk; k++) 
            {
	       xb[k] = 0.0e0;
               for (i1=0; i1<rnk; i1++) xb[k] += binv[k][i1]*A[j][i1];
	    }
 
            // To find the biggest component from xb[0] to xb[rnk-1]
            tmp  = 0.0e0;
            itmp = -1;
            for (k=1; k<rnk; k++)
	    {
               if (fabs(xb[k]) > tmp && ib[k] == -1)
	       {
                  tmp = fabs(xb[k]);
                  itmp = k;
               }
            }
            if (itmp == -1)
            {
	       info = -1;
	       abort();  // It is impossible in our problem
            }

            // itmp-th row will get out of the basis

            // To find the new B^{-1} after itmp out and j in

            for (k=0; k<itmp; k++) xb[k] /= xb[itmp];
            for (k=itmp+1; k<rnk; k++) xb[k] /= xb[itmp];
            for (i1=0; i1<rnk; i1++)
            {
	       tmp = binv[itmp][i1];
               for (k=0; k<itmp; k++) binv[k][i1] -= xb[k]*tmp;
               for (k=itmp+1; k<rnk; k++) binv[k][i1] -= xb[k]*tmp;
               binv[itmp][i1] = tmp/xb[itmp];
            }
            ib[itmp] = j;
	 }
 
         for (i1=0; i1<rnk; i1++)
         {
	    if (ib[i1] != -1)
	    {
               xb[i1] = 1.0e0;
            }
	    else
            {
	       xb[i1] = 0.0e0;
               ib[i1] = ib0[i1];
            }
         }
 
         for (i1=0; i1<rnk; i1++) b[i1] = A[j][i1];

 
         RSimplex(A,b,jStrt,jEnd,rnk,xb,ib,binv,info);
         if (info > -1)
	 {
            if (fabs(xb[info]) > 1.0e-10)
            {
	       // It is a non-vertex
               ynVtx[j-1] = 0;  // -1: since the added 1st column
	    }
         }
         for (i1=0; i1<rnk; i1++) ib0[i1] = ib[i1];
      } 
   }

   delete [] A[0];
   delete [] A;
   delete [] b;
   delete [] ib;
   delete [] ib0;
   delete [] xb;
   delete [] binv[0];
   delete [] binv;
   delete [] sb;

   return;
}
#endif
