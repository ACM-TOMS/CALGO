/*
 * This package is to calculate the mixed volume of a system of n polynomials
 * in n variables.
 * 
 * It can take the polynomial system itself  or the support of the polynomial
 * system as its initial input.
 *
 * The output of this code is the mixed volume of the polynomial system.
 * 
 * For the format of the initial input files, see the file ReadMe included,
 * or the examples included in the subdirectories "equations" and "supports".
 *
 * Version 1.0 Copyright (C) 2003  Tangan Gao, T. Y. Li, Xing Li, and Mengnien Wu.
 * 
*/
/*
 * This package can be divided into three main parts:
 * 
 * Part 1: The codes included in the subdirectory "PolyReader"
 *
 *         This part reads the polynomial system from the initial input file
 *         then generates the support of the polynomial system.
 *
 *         The codes in this part are developed by
 *                   Xing Li
 *
 * Part 2: The codes included in the subdirectory "PreProcess"
 *
 *         This part removes the non-vertex points from each support,
 *         forms a smaller support which consists of only the vertex points,
 *         figures out the type of the support of the polynomial system,
 *         then re-arranges the support.
 *
 *         The codes in this part are developed by
 *                   Tangan Gao
 *                   Mengnien Wu
 *                   T.Y. Li
 *
 * Part 3: The codes included in the subdirectory "MixedVol"
 *
 *         This part calculates the mixed volume of the vertex support with
 *         its type given.
 *
 *         The codes in this part are developed by
 *                   Tangan Gao
 *                   Mengnien Wu
 *                   T.Y. Li
 *
 * 
 *
 * The authors greatly appreciate any comments from users. For the bug
 * reports, questions, suggestions for improvement, or applications that
 * have been solved successfully, contact with:  
 *       tgao@csulb.edu, li@math.msu.edu  or  mwu@mail.tku.edu.tw
*/

////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <string>
#include <ctime>
#include <fstream>
#include <iostream>

#include "PolyReader/PolynomialSystemReader.h"
#include "PolyReader/PolynomialException.h"
#include "PreProcess/Pre4MV.h"
#include "MixedVol/MixedVol.h"

using namespace std;

int main(int argc, char** argv)
{
	if ( argc != 3 || strlen(argv[1]) != 2 || argv[1][0] != '-' 
	     || (argv[1][1] != 'p' && argv[1][1] != 's'))
	{
	   cout <<"*******************************" << endl;
	   cout <<"Usage: " << endl;
	   cout << argv[0] << " -p yourPolynomialFile" << endl;
	   cout <<"  or" << endl;
	   cout << argv[0] << " -s yourSupportFile" << endl;
	   cout <<"*******************************" << endl;
	   exit(0);
        }

	int i, j, k;

	int nPts = 0;
	int nVar;
	int nSpt;
	int* Spt1stIdx;
	int** Spt;

	if ( argv[1][1] == 'p' )   // The polynomial system is given
	{
	   fstream IN;
	   IN.open(argv[2],ios::in);
           if ( !(IN.is_open()) )
           {
              cout << "Your polynomial file " << argv[2] << " does not exist" << endl;
              exit(0);
           }

           PolynomialSystemReader ps( IN );
           ps.GetDimensions(nVar,nPts);  // # of variables and # of terms
	   Spt1stIdx = new int [nVar+1];
	   Spt = new int* [nPts];
	   Spt[0] = new int [nPts*nVar];
	   for (i=1; i<nPts; i++) Spt[i] = Spt[i-1] + nVar;

           ps.GetSupport(Spt1stIdx,Spt);
	}
	else // if ( argv[1][1] == 's' )  // The support is given
	{
	   fstream IN;
	   IN.open(argv[2],ios::in);
           if ( !(IN.is_open()) )
           {
              cout << "Your support file " << argv[2] << " does not exist" << endl;
              exit(0);
           }

	   if ( !(IN >> nVar) )   // # of variables in the system
	   {
	      cout << "Your support file " << argv[2] << " has an error!" << endl;
	      exit(0);
	   }

           Spt1stIdx = new int [nVar+1];
           for (i=0; i<nVar; i++)
	   {
              if ( !(IN >> Spt1stIdx[i]) )  // # of points in the i-th support
	      {
	         cout << "Your support file " << argv[2] << " has an error!" << endl;
	         exit(0);
	      }
	      nPts += Spt1stIdx[i];
	   }

           Spt = new int* [nPts];
	   Spt[0] = new int [nVar*nPts];
           for (i=1; i<nPts; i++) Spt[i] = Spt[i-1] + nVar;
           for (i=0; i<nPts; i++)
	   {
	      for (j=0; j<nVar; j++)
	      {
		 if ( !(IN >> Spt[i][j]) )
	         {
	            cout << "Your support file " << argv[2] << " has an error!" << endl;
	            exit(0);
	         }
	      }
           }

	   IN.close();
	}
	Spt1stIdx[nVar] = nPts;
	for (i=nVar-1; i>=0; i--)
	{
	   Spt1stIdx[i] = Spt1stIdx[i+1] - Spt1stIdx[i];
	}

        // Quick return if 1-variable system or less than two terms
	k = -1;
	for (i=0; i<nVar; i++)
	{
	   if ( Spt1stIdx[i+1]-Spt1stIdx[i] < 2 )
           {
	      k = i;
	      break;
	   }
	}
	if ( k >= 0 )
        {
	   if ( argv[1][1] == 'p' )
	   {
              cout << "The " << k+1 << "-th polynomial has less than 2 terms" << endl;
	   }
	   else
           {
              cout << "The " << k+1 << "-th support has less than 2 points" << endl;
	   }
	   exit(0);   // end of the case: too few terms
	}
	
	if ( nVar == 1 )
        {
	   int kmin, kmax;
	   kmin = INT_MAX/2;
	   kmax = -INT_MAX/2;
	   for (i=0; i<Spt1stIdx[1]; i++) 
           {
	      kmax = max(kmax, Spt[i][0]);
	      kmin = min(kmin, Spt[i][0]);
	   }
	   if ( argv[1][1] == 'p' )
	   {
              cout << "A 1-variable polynomial, its mixed volume is " << kmax-kmin << endl;
	   }
	   else
           {
              cout << "A 1-dimensional support, its mixed volume is " << kmax-kmin << endl;
	   }
	   exit(0);   // end of the case: 1-variable
	}
	for (i=0; i<Spt1stIdx[nVar]; i++)
	{
	   int kmin, kmax;
	   kmin = INT_MAX/2;
	   kmax = -INT_MAX/2;
	   k=0;
	   for (j=0; j<nVar; j++)
           {
	      kmax = max(kmax, Spt[i][j]);
	      kmin = min(kmin, Spt[i][j]);
	      k += abs(Spt[i][j]);
	   }
	   if ( kmin != 0 || kmax > 1 || k > 1 )
           {
	      j = -1;
	      break;
	   }
	}
	if ( j != -1 )
        {
	   if ( argv[1][1] == 'p' )
	   {
              cout << "A linear system, its mixed volume <= 1" << endl;
	   }
	   else
           {
              cout << "A linear support, its mixed volume <= 1" << endl;
	   }
	   exit(0);   // end of the case: linear system
	} 
	// end of quick return
	   
	// To preprocess the support

	int* SptType = new int [nVar];
        int* SptVtx1stIdx = new int [nVar+1];
        int** SptVtx = new int* [nPts];
	SptVtx[0] = new int [nVar*nPts];
        for (i=1; i<nPts; i++) SptVtx[i] = SptVtx[i-1] + nVar;
	int* NuIdx2OldIdx = new int [nPts];
	nSpt = nVar;

	Pre4MV(nVar,nSpt,SptType,Spt,Spt1stIdx,SptVtx,SptVtx1stIdx,NuIdx2OldIdx);

	// To calculate the mixed volume of the support

	srand( unsigned(time(0)) );
	double* lft = new double[SptVtx1stIdx[nSpt]];
	for (j=0; j<SptVtx1stIdx[nSpt]; j++)
        {
           lft[j] = 1.0+(3*double(rand())-double(rand()))/RAND_MAX;
        }

	int CellSize = nSpt;
	for (j=0; j<nSpt; j++ )
	   CellSize += SptType[j];
	CellStack  MCells( CellSize );
	int MVol;

        MixedVol(nVar,nSpt,SptType,SptVtx1stIdx,SptVtx,lft,MCells,MVol);

	if ( argv[1][1] == 'p' )
	{
           cout << "The mixed volume of this system is " << MVol << endl;
	}
	else
        {
           cout << "The mixed volume of this support is " << MVol << endl;
	}

	// Clean the memory

	while( ! MCells.IsEmpty() ) MCells.Pop();
	
	delete [] Spt1stIdx;
	delete [] Spt[0];
	delete [] Spt;
	delete [] SptType;
	delete [] SptVtx1stIdx;
	delete [] SptVtx[0];
	delete [] SptVtx;
	delete [] NuIdx2OldIdx;
        delete [] lft;


        return 0;
}
