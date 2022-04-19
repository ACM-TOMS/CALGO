#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <complex.h>
#include "./solvers.h"
extern double ranf(void);
int main(int argc, char **argv)
{
  complex double csol[4]; 
  double c[5];
  int numtrials, its, caso, okHQR, num;
  srand48(4242);
  numtrials=100000000;
  if (argc==1)
    caso=-1; // do nothing
  else
    caso = atoi(argv[1]);
  if (argc > 2)
    numtrials = atoi(argv[2]);
  printf("Timing Test of ");
  switch (caso)
    {
    case 1:
      printf("OQS\n");
      break;
    case 2:
      printf("FLO\n");
      break;
    case 3:
      printf("STR\n");     
      break;
    case 4: 
      printf("ABR\n");
      break;
    case 5: 
      printf("FQS\n");
      break;	
    case 6:
      printf("HQR\n");
      break;
    case 7: 
      printf("SHM\n");
      break;	
    default:
      printf("Nothing\n");
      break;
    }
  for (its=0; its < numtrials; its++)
    {
      c[4]=1.0;
      c[3]=ranf()-0.5;
      c[2]=ranf()-0.5;
      c[1]=ranf()-0.5;
      c[0]=ranf()-0.5;
      switch (caso)
	{
	case 1:
	  oqs_quartic_solver(c, csol);     
	  break;
	case 2:
	  CquarticRoots(c,&num, csol);
	 break;
	case 3:
	  CLDLT_quartic(c, csol);     
	  break;
	case 4: 
	 csolve_quartic_abramovitz_cmplx(c,csol);
	 break;
	case 5: 
	 fast_quartic_solver(c, csol); 
	 break;	
	case 6:
	 solve_numrec(c, csol, &okHQR);
	 break;
	case 7: 
	 csolve_quartic_shmakov(c, csol); 
	 break;	
	default:
	 // do nothing
	 break;
	}
    }
  exit(0);
}
