#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include <complex.h>
#include "./solvers.h"
extern double ranf(void);
//#define TEST_VERY_LARGE_ROOTS
int main(int argc, char **argv)
{
  long double x1, y1;
  complex long double x1c, x2c, x3c, x4c;
  complex double csol[4]; 
  double sig, c[5];
  int numtrials, its, caso, okHQR, num;
#ifdef TEST_VERY_LARGE_ROOTS
  int k, cc=0;
#endif
  srand48(4242);
  numtrials=100000000;
  sig=1.0;
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
      printf("ODM\n");
      break;
    case 2:
      printf("FLO\n");
      break;
    case 3:
      printf("STR\n");     
      break;
    case 4: 
      printf("FER\n");
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
#ifdef TEST_VERY_LARGE_ROOTS
      x1 = 1E140+1E10*sig*(ranf()-0.5);
      y1 = 1E130+1E10*sig*(ranf()-0.5);
      x1c = x1;
      x2c = y1;
      x3c = 100+1E10*sig*(ranf()-0.5);
      x4c = 1+1E10*sig*(ranf()-0.5);
#else
      x1 = sig*(ranf()-0.5);
      y1 = sig*(ranf()-0.5);
      x1c = x1;
      x2c = y1;
      x1 = sig*(ranf()-0.5);
      y1 = sig*(ranf()-0.5);
      x3c = x1 + I*y1;
      x4c = x1 - I*y1;
#endif
      c[4] = 1.0;
      c[3] = creall(-(x1c+x2c+x3c+x4c));
      c[2] = creall(x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c); 
      c[1] = creall(-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c));
      c[0] = creall(x1c*x2c*x3c*x4c);
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
#ifdef TEST_VERY_LARGE_ROOTS
      for (k=0; k < 4; k++)
	if (isnan(creal(csol[k])) || isinf(creal(csol[k]))||
	    isnan(cimag(csol[k])) || isinf(cimag(csol[k])))
	  cc++;
#endif
    }
#ifdef TEST_VERY_LARGE_ROOTS
  printf("#INFs or NANs = %d\n", cc);
#endif
  exit(0);
}
