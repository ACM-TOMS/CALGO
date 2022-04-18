#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#define Sqr(x) ((x)*(x))
#include <complex.h>
#include "./solvers.h"
#ifndef CMPLX
#define CMPLX(x,y) ((x)+I*(y))
#endif
#ifndef CMPLXL
#define CMPLXL(x,y) ((x)+I*(y))
#endif
extern int perm[24][4]; 

extern void sort_sol_opt(complex double *sol, complex double *exsol);
extern void sort_sol_optl(complex long double *sol, complex double *exsol);

#ifdef OQS_LONG_DBL
long double print_accuracy_atl(char *str, complex long double *csol, complex double exsol[4])
{
  /* we follow FLocke here */
  int k1;
  long double relerr, relerrmax;
  for (k1=0; k1 < 4; k1++)
    {
      relerr=cabsl((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=cabsl((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16LG\n", str, relerrmax);
  return relerrmax;
}
#endif

double print_accuracy_at(char *str, double complex *csol, complex double exsol[4])
{
  /* we follow FLocke here */
  int k1;
  double relerr, relerrmax;
  for (k1=0; k1 < 4; k1++)
    {
      relerr=cabs((csol[k1] - exsol[k1])/exsol[k1]); 
      if (k1==0 || relerr > relerrmax)
        {
          relerrmax=cabs((csol[k1] - exsol[k1])/exsol[k1]); 
        }
    }
  printf("[%s] relative accuracy=%.16G\n", str, relerrmax);
  return relerrmax;
}

void print_roots(char *str, complex long double x1c, complex long double x2c,
                 complex long double x3c, complex long double x4c)
{
  printf("%s\n", str);
  printf("exact root #1= %.36LG+I*(%.36LG)\n", creall(x1c), cimagl(x1c));
  printf("exact root #2= %.36LG+I*(%.36LG)\n", creall(x2c), cimagl(x2c));
  printf("exact root #3= %.36LG+I*(%.36LG)\n", creall(x3c), cimagl(x3c));
  printf("exact root #4= %.36LG+I*(%.36LG)\n", creall(x4c), cimagl(x4c));
}
int main(int argc, char **argv)
{
  complex long double x1c, x2c, x3c, x4c; 
  complex double csol[4];
  complex double csolREF[4];
  double c[5], S;
  int num, k1, okHQR, caso;
#ifdef OQS_LONG_DBL
  long double cl[5];
  complex long double csoll[4];
#endif
  if (argc == 2)
    {
      caso = atoi(argv[1]);
    }
  else
    {
      caso = 1;
    }
  if (caso <= 0 || caso > 24)
    {
      printf("caso must be between 1 and 20\n");
      exit(-1);
    }
  x1c=x2c=x3c=x4c=0.0;
  if (caso > 0)
    {
      switch (caso)
        {
          /* Strobach */
        case 14:
          x1c = x2c = x3c = x4c = 1000;
          print_roots("CASE 14", x1c, x2c, x3c, x4c);
          break;
        case 15:
          x1c = x2c = x3c = 1000;
          x4c = 1E-15;
          print_roots("CASE 15", x1c, x2c, x3c, x4c);
          break;
        case 16:
          x3c = 1E16 + I*1E7;
          x4c = 1E16 - I*1E7;
          x1c = 1 + 0.1*I;
          x2c = 1 - 0.1*I;
          print_roots("CASE 16", x1c, x2c, x3c, x4c);
          break;
        case 17:
          x1c=10000;
          x2c=10001;
          x3c=10010;
          x4c=10100;
          print_roots("CASE 17", x1c, x2c, x3c, x4c);
          break;
        case 18:
          x1c=4E5+I*3E2;
          x2c=4E5-I*3E2;
          x3c=3E4+I*7E3;
          x4c=3E4-I*7E3;
          print_roots("CASE 18", x1c, x2c, x3c, x4c);
          break;
        case 1:
          x1c = 1E9;
          x2c = 1E6;
          x3c = 1E3;
          x4c = 1;
          print_roots("CASE 1", x1c, x2c, x3c, x4c);
          break;
        case 2:
          x1c = 2.003;
          x2c = 2.002;
          x3c = 2.001;
          x4c = 2;
          print_roots("CASE 2", x1c, x2c, x3c, x4c);
          break;
        case 3:
          x1c = 1E53;
          x2c = 1E50;
          x3c = 1E49;
          x4c = 1E47;
          print_roots("CASE 3", x1c, x2c, x3c, x4c);
          break;
        case 4:
          x1c = 1E14;
          x2c = 2;
          x3c = 1;
          x4c = -1;
          print_roots("CASE 4", x1c, x2c, x3c, x4c);
          break;
        case 5:
          x1c = -2E7;
          x2c = 1E7;
          x3c = 1;
          x4c = -1;
          print_roots("CASE 5", x1c, x2c, x3c, x4c);
          break;
        case 6:
          x1c = 1E7;
          x2c = -1E6;
          x3c = 1+I;
          x4c = 1-I;
          print_roots("CASE 6", x1c, x2c, x3c, x4c);
          break;
        case 7:
          x1c = -7;
          x2c = -4;
          x3c = -1E6+I*1E5;
          x4c = -1E6-I*1E5;
          print_roots("CASE 7", x1c, x2c, x3c, x4c);
          break;
        case 8:
          x1c = 1E8;
          x2c = 11;
          x3c = 1E3+I;
          x4c = 1E3-I;
          print_roots("CASE 8", x1c, x2c, x3c, x4c);
          break;
        case 9:
          x1c = 1E7+I*1E6;
          x2c = 1E7-I*1E6;
          x3c = 1+2*I;
          x4c = 1-2*I;
          print_roots("CASE 9", x1c, x2c, x3c, x4c);
          break;
        case 10:
          x1c = 1E4+3*I;
          x2c = 1E4-3*I;
          x3c = -7+1E3*I;
          x4c = -7-1E3*I;
          print_roots("CASE 10", x1c, x2c, x3c, x4c);
          break;
        case 11:
          x1c = 1.001+4.998*I;
          x2c = 1.001-4.998*I;
          x3c = 1.000+5.001*I;
          x4c = 1.000-5.001*I;
          print_roots("CASE 11", x1c, x2c, x3c, x4c);
          break;
        case 12:
          x1c = 1E3+3*I;
          x2c = 1E3-3*I;
          x3c = 1E3+I;
          x4c = 1E3-I;
          print_roots("CASE 12", x1c, x2c, x3c, x4c);
          break;
        case 13:
          x1c = 2+1E4*I;
          x2c = 2-1E4*I;
          x3c = 1+1E3*I;
          x4c = 1-1E3*I;
          print_roots("CASE 13", x1c, x2c, x3c, x4c);
          break;
        case 19:
          x1c = 1E44;
          x2c = 1E30;
          x3c = 1E30;
          x4c = 1.0;
          print_roots("CASE 19", x1c, x2c, x3c, x4c);
          break;
        case 20:
          x1c = 1E14;
          x2c = 1E7;
          x3c = 1E7;
          x4c = 1.0;
          print_roots("CASE 20", x1c, x2c, x3c, x4c);
          break;
        case 21:
          x1c = 1E15;
          x2c = 1E7;
          x3c = 1E7;
          x4c = 1.0;
          print_roots("CASE 21", x1c, x2c, x3c, x4c);
          break;
        case 22:
          x1c = 1E154;//1E102;
          x2c = 1E152;//1E100;
          x3c = 10.0;
          x4c = 1.0;
          print_roots("CASE 22", x1c, x2c, x3c, x4c);
          break;
        case 23:
          /* condition */
          c[4] = 1.0;
          c[3] = 1.0;
          c[2] = 1.0;
          c[1] = 3.0/8.0;
          c[0] = 0.001;
          printf("CASE 23\n");
          break;
        case 24:
          S = 1E30;
          c[4] = 1.0;
          c[3] = -(1.0+1.0/S);
          c[2] = 1.0/S - S*S;
          c[1] = S*S + S;
          c[0] = -S;
          printf("CASE 24\n");
          break;
        }
    }
  if (caso <=22)
    {
      csolREF[0]=x1c;
      csolREF[1]=x2c;
      csolREF[2]=x3c;
      csolREF[3]=x4c;
#ifdef OQS_LONG_DBL
      cl[4] = 1.0;
      cl[3] = creall(-(x1c+x2c+x3c+x4c));
      cl[2] = creall(x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c); 
      cl[1] = creall(-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c));
      cl[0] = creall(x1c*x2c*x3c*x4c);
#endif
      c[4] = 1.0;
      c[3] = creall(-(x1c+x2c+x3c+x4c));
      c[2] = creall(x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c); 
      c[1] = creall(-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c));
      c[0] = creall(x1c*x2c*x3c*x4c);
    }
  else
    csolREF[0]=csolREF[1]=csolREF[2]=csolREF[3]=0.0;
  printf("(%.32G)*x^4+(%.32G)*x^3+(%.32G)*x^2+(%.32G)*x+(%.32G)==0\n", c[4], c[3], c[2], c[1], c[0]);

#ifdef OQS_LONG_DBL
  oqs_quartic_solver_dl(cl, csoll);     
#else
  oqs_quartic_solver(c, csol);     
#endif

#ifdef OQS_LONG_DBL
  sort_sol_optl(csoll, csolREF);
  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        printf("[ODM] root #%d=  %.15LG+I*(%.15LG) [%.15LG + I*(%.15LG)]\n", 
               k1, creall(csoll[k1]), cimagl(csoll[k1]), creall(csolREF[k1]), cimagl(csolREF[k1]));
      else
        printf("[ODM] root #%d=  %.15LG+I*(%.15LG)\n", 
               k1, creall(csoll[k1]), cimagl(csoll[k1]));
    }
  if (caso <=22)
    print_accuracy_atl("ODM", csoll, csolREF);
#else
  sort_sol_opt(csol, csolREF);
  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        printf("[ODM] root #%d=  %.15G+I*(%.15G) [%.15G + I*(%.15G)]\n", 
               k1, creal(csol[k1]), cimag(csol[k1]), creal(csolREF[k1]), cimag(csolREF[k1]));
      else
        printf("[ODM] root #%d=  %.15G+I*(%.15G)\n", 
               k1, creal(csol[k1]), cimag(csol[k1]));
    }
  if (caso <=22)
    print_accuracy_at("ODM", csol, csolREF);
#endif
  printf("===============================\n");

  CquarticRoots(c,&num, csol);
  sort_sol_opt(csol, csolREF);
  for (k1=0; k1 < 4; k1++)
    {
      if (caso <=22)
        printf("[FLO] root #%d=  %.15G+I*(%.15G) [%.15G + I*(%.15G)]\n", 
               k1, creal(csol[k1]), cimag(csol[k1]), creal(csolREF[k1]), cimag(csolREF[k1]));
      else
        printf("[FLO] root #%d=  %.15G+I*(%.15G)\n", 
               k1, creal(csol[k1]), cimag(csol[k1]));
    } 
  if (caso <=22)
    print_accuracy_at("FLO", csol, csolREF);

  printf("===============================\n");

  CLDLT_quartic(c, csol);
  sort_sol_opt(csol, csolREF);

  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        printf("[STR] root #%d=  %.15G+I*(%.15G) [%.15G + I*(%.15G)]\n", 
               k1, creal(csol[k1]), cimag(csol[k1]), creal(csolREF[k1]), cimag(csolREF[k1]));
      else
        printf("[STR] root #%d=  %.15G+I*(%.15G)\n", 
               k1, creal(csol[k1]), cimag(csol[k1]));
      
    } 
  if (caso <=22)
    print_accuracy_at("STR", csol, csolREF);

  printf("===============================\n");
  csolve_quartic_abramovitz_cmplx(c, csol);
  sort_sol_opt(csol, csolREF);

  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        printf("[FER] root #%d=  %.15G+I*(%.15G) [%.15G + I*(%.15G)]\n", 
               k1, creal(csol[k1]), cimag(csol[k1]), creal(csolREF[k1]), cimag(csolREF[k1]));
      else
        printf("[FER] root #%d=  %.15G+I*(%.15G)\n", 
               k1, creal(csol[k1]), cimag(csol[k1]));

    } 
  if (caso <=22)
    print_accuracy_at("FER", csol, csolREF);

  printf("===============================\n");

  fast_quartic_solver(c, csol);
  sort_sol_opt(csol, csolREF);
  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        printf("[FQS] root #%d=  %.15G+I*(%.15G) [%.15G + I*(%.15G)]\n", k1, creal(csol[k1]), cimag(csol[k1]),
               creal(csolREF[k1]), cimag(csolREF[k1]));
      else
        printf("[FQS] root #%d=  %.15G+I*(%.15G)\n", k1, creal(csol[k1]), cimag(csol[k1]));
    }
  if (caso <=22)
    print_accuracy_at("FQS", csol, csolREF);

  printf("===============================\n");
  solve_numrec(c, csol, &okHQR);
  sort_sol_opt(csol, csolREF);
  for (k1=0; k1 < 4; k1++)
    {
      if (caso <= 22)
        printf("[HQR] root #%d=  %.15G+I*(%.15G) [%.15G + I*(%.15G)]\n", 
               k1, creal(csol[k1]), cimag(csol[k1]), creal(csolREF[k1]), cimag(csolREF[k1]));
      else
        printf("[HQR] root #%d=  %.15G+I*(%.15G)\n", 
               k1, creal(csol[k1]), cimag(csol[k1]));
    } 
  if (caso <=22)
    print_accuracy_at("HQR", csol, csolREF);
#if 0
  csolve_quartic_shmakov(c, csol);
  sort_sol_opt(csol, csolREF);
  for (k1=0; k1 < 4; k1++)
    {
      printf("[SHM] root #%d=  %.15G+I*(%.15G) [%.15G + I*(%.15G)]\n", 
             k1, creal(csol[k1]), cimag(csol[k1]), creal(csolREF[k1]), cimag(csolREF[k1]));
    } 
  if (caso <= 18)
    print_accuracy_at("SHM", csol, csolREF);
#endif
  exit(-1);
}
