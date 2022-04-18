/*
   --------------------------------------------------------------
   Testexample  of ADOL-C version 1.6 as of January 1,   1995
   --------------------------------------------------------------
*/

#include <stream.h>
#include <stdio.h>

#ifdef __GNUG__
#include <std.h>
#else
#include <stdlib.h>
#endif
#include "usrparms.h"
extern "C" {
#include "taputil3.h"
}
#include "adouble.h"    // These includes provide the compiler with
#include "adutils.h"    // definitions and utilities for `adoubles'.

void gausselim(int n, adoublem& A, adoublev& bv)
{
  along i;
  adoublev temp(n);
  adouble r,rj,temps;
  int j,k,ik;
  for (k=0; k < n; k++)
  {
    for(j=0;j<n;j++)
      cout << value(A[k][j]) << "  ";
    cout << "             " << value(bv[k]) << "\n";
  }
  cout << "initial state ----------------------\n";
  for (k=0; k < n; k++) /* elimination loop */
  {
    i = k;
    r = fabs(A[k][k]); /* initial pivot element */
    for (j=k+1; j<n; j++)
    {
      rj = fabs(A[j][k]); /* look for a greater element in the same column */
       condassign(i,(rj >r),j);
       condassign(r,(rj >r),rj);
    } // endfor
     cout << i << "index \n";
    temp = A[i];
    A[i] = A[k];
    A[k] = temp; /* exchange  of rows */
    temps = bv[i];
    bv[i]=bv[k];
    bv[k]=temps;
    if (!value(A[k][k]))
    {
      cout << " Matrix does not have full rank!\n";
      exit(-1);
    } // endif
    cout << "changed rows: ---------------------\n";
    for (ik=0; ik < n; ik++)
    {
    for(j=0;j<n;j++)
      cout << value(A[ik][j]) << "  ";
    cout << "             " << value(bv[ik]) << "\n";
    }
    temps= A[k][k];
    A[k] = A[k]/temps;
    bv[k] = bv[k] /temps;
    for (j=k+1; j<n; j++)
    {
      temps= A[j][k];
      A[j] -= temps*A[k];
      bv[j] -= temps*bv[k];
    } // endfor
    cout << "step:---------------------------\n";
    for (ik=0; ik < n; ik++)
    {
    for(j=0;j<n;j++)
      cout << value(A[ik][j]) << "  ";
    cout << "             " << value(bv[ik]) << "\n";
    }
  } // endfor elimination loop
  temp=0.0;
  for(k=n-1; k >= 0; k--)
  {
     temp[k] = (bv[k]-(A[k]*temp))/A[k][k];
      cout << value(temp[k]) << "\n";
  }
  bv=temp;
  return;
} // end gausselim


void main() 
{
  int i,j,k,l=1,ok=1;
  short tag = 1;
  double epsilon=0.0000000001; // max. allowed difference between results
  int dum=1;
  const int max_deg=4; // maximal order of derivation
  const int tayl_num=2; // Number of taylor series
  int pages=1; // for eventually writing the tape to a file
  const int size=5;
  const int indep=size*size+size;
  const int depen=size;
  const int laglength=2;
  double*lagras=new double[depen];
  for(i=0;i<depen;i++)
   lagras[i]=i+1;
  double** lagrav=myalloc(laglength,depen);
  for(j=0;j<laglength;j++)
  {
    for(i=0;i<depen;i++)
      lagrav[j][i]=j+i+1;
  } // endfor
  short** nonzero=new short*[laglength];
  for(i=0;i<laglength;i++)
    nonzero[i]=new short[indep];
  double** resultshos=myalloc(indep,max_deg);
  double*** resultshov=myalloc(laglength,indep,max_deg);
  double* resultsfos=new double[indep];
  double** resultsfov=myalloc(laglength,indep);
  double* basepoint=new double[indep];
  double* valuepoint=new double[depen];
  double*** arguments=myalloc(indep,tayl_num,max_deg);
  double*** scalargs=myalloc(tayl_num,indep,max_deg+1);
  double*** scalres=myalloc(tayl_num,depen,max_deg+1);
  double*** taylors=myalloc(depen,tayl_num,max_deg);

  double yp[size],xp[size*size+size];  // passive variable
  adoublem A(size,size);
  adoublev bv(size);  // active variables
  int N=size*size;
  trace_on(tag,dum);   // Begin taping all calculations with 'adoubles'
  for(i=0;i<size;i++)
  {
    for(j=0;j<size;j++)
    {
      A[i][j]<<=pow(1+j,i); /* indep. vars */
      xp[i*size+j]=pow(1+j,i); /* args for forward */
    } /* endfor */
  } /* endfor */
  for(i=0;i<size;i++)
  {
    bv[i]<<=-i-1; /* indep. vars */
    xp[N+i]=-i-1; /* args for forward */
  }
  gausselim(size,A,bv);
  bv >>= yp;
  trace_off(); 
  int buf_size,maxlive,deaths;
  int tape_stats[11];  /* tape stats */
  tapestats(tag,tape_stats);

  maxlive = tape_stats[2];
  deaths = tape_stats[3];
  buf_size =  tape_stats[4];

  // initialization for Forward - testing +++++++++++++++++++++++++++++++++++++

  basepoint=xp;
  for(i=0;i<tayl_num;i++)
  {
    for(j=0;j<indep;j++)
    {
      scalargs[i][j][0]=xp[j];
    } // endfor 
  } // endfor 
  for(i=0;i<indep;i++)
  {
    for(j=0;j<tayl_num;j++)
    {
      for(k=0;k<max_deg;k++)
      {
        arguments[i][j][k]=i+j+k+3;
        scalargs[j][i][k+1]=i+j+k+3;
      } // endfor
    } // endfor
  } // endfor

  // calculation ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(i=0;i<tayl_num;i++)
    hos_forward(tag,depen,indep,max_deg,1,scalargs[i],scalres[i]);
  hov_forward(tag,depen,indep,max_deg,tayl_num,basepoint,arguments,valuepoint,taylors);
  // test for correctness +++++++++++++++++++++++++++++++++++++++++++++++
  for(i=0;i<depen;i++)
  {
    cout << "dependent variable number "<< i << "\n";
    for(j=0;j<tayl_num;j++)
    {
      cout << "taylor serie number " << j << "\n";
      cout <<"hov_f. valuepoint["<<i<<"]: "<< valuepoint[i] << " =?     hos_f. scalres["<<j<<"]["<<i<<"][0]: " << scalres[j][i][0] << "   =? yp["<<i<<"] : " << yp[i] << "\n";
      if (fabs(valuepoint[i]-scalres[j][i][0])>epsilon)
      {
        cout << "difference is here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
        ok=0;
      }
      for(k=0;k<max_deg;k++)
      {
        cout <<"hov_f. taylors["<<i<<"]["<<j<<"]["<<k<<"]:"<< taylors[i][j][k] <<"  =?   hos_f. scalres["<<j<<"]["<<i<<"]["<<k+1<<"]: " << scalres[j][i][k+1] << "   <- inp. coeff.: " << scalargs[j][i][k+1] << "\n";
        if (fabs(taylors[i][j][k]-scalres[j][i][k+1])>epsilon)
        {
          cout << "difference is here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
          ok=0;
        } // endif
      } // endfor
    } // endfor
  } // endfor

  // some preparation for the 4 different reverse modes: -----------------
  hos_forward(tag,depen,indep,max_deg,max_deg,scalargs[0],scalres[0]);
  cout << "reverse sweeps will be done for the first taylor serie only\n";
  hos_reverse(tag,depen,indep,max_deg-1,lagras,resultshos);
  hov_reverse(tag,depen,indep,max_deg-1,laglength,lagrav,resultshov,nonzero);
  hos_forward(tag,depen,indep,max_deg,1,scalargs[0],scalres[0]);
  fos_reverse(tag,depen,indep,lagras,resultsfos);
  fov_reverse(tag,depen,indep,laglength,lagrav,resultsfov);
  // output
  for (i=0;i<laglength;i++)
  {
    if (i==0)
    {
    for (j=0;j<indep;j++)    
    {
      for (k=0;k<max_deg;k++)
      {
        if (k==0)
        {
          cout << "reshov["<<i<<"]["<<j<<"]["<<k<<"]: "<<  resultshov[i][j][k] << " =? reshos["<<j<<"]["<<k<<"]: "<<  resultshos[j][k] << " =? resfov["<<i<<"]["<<j<<"]: "<<  resultsfov[i][j] << " =? resfos["<<j<<"]: "<< resultsfos[j] << "\n";
          if ((fabs(resultshov[i][j][k]-resultshos[j][k])>epsilon) || (fabs(resultshov[i][j][k]- resultsfov[i][j])>epsilon) || (fabs(resultshov[i][j][k]-resultsfos[j])>epsilon))
          {
            cout << "difference is here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
            ok=0;
          } // endif
        } // endif
        else
        {
          cout << "reshov["<<i<<"]["<<j<<"]["<<k<<"]: "<<  resultshov[i][j][k] << " =? reshos["<<j<<"]["<<k<<"]: "<<  resultshos[j][k] << "\n";
          if (fabs(resultshov[i][j][k]-resultshos[j][k])>epsilon)
          {
            cout << "difference is here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
            ok=0;
          } // endif
        } // endelse
      } // endfor
    } // endfor
    } // endif 
    else
    {
    for (j=0;j<indep;j++)    
    {
      for (k=0;k<max_deg;k++)
      {
        if (k==0)
        {
          cout << "reshov["<<i<<"]["<<j<<"]["<<k<<"]: "<<  resultshov[i][j][k] << " =? resfov["<<i<<"]["<<j<<"]: "<<  resultsfov[i][j]  << "\n";
          if (fabs(resultshov[i][j][k]- resultsfov[i][j])>epsilon)
          {
	    cout << "difference is here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
            ok=0;
          } // endif
        } // endif
        else
          cout << "reshov["<<i<<"]["<<j<<"]["<<k<<"]: "<<  resultshov[i][j][k] << "\n";
      } // endfor
    } // endfor
    } // endelse
  } // endfor 
  for (i=0;i<laglength;i++)
  {
    for (j=0;j<indep;j++)
    cout << "nonzero["<<i<<"]["<<j<<"]: " << nonzero[i][j] << "\n";
  } // endfor 
  if (! ok)
    cout << "calculation  is not ok<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n This message may be caused by very small differences (not necessary \nrecognizable in the output)\n";
} // endmain  

