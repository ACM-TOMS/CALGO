#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<unistd.h>
#include <complex.h>
#include "./solvers.h" 
#define MAXSOLV 10
#define PEPTS 150
extern double ranf(void);
extern int perm[24][4];
extern void sort_sol_opt(complex double *sol, complex double *exsol);

int cmplxreal=0, restart, dojust=-1;
long double cumPEHQRL[PEPTS],PEHQRL[PEPTS]; 
double cumPEall[MAXSOLV][PEPTS], PEall[MAXSOLV][PEPTS];
complex double csolall[MAXSOLV][4];
complex double exsol[4];
char algs[8][64] = {"HQR", "FQS", "FER", "STR", "FLO", "ODM", "SHM", "HQRL"};
char *ic2algo(int ic)
{
  if (ic > 8)
    {
      exit(-1);
    }
  return algs[ic];
}
extern void sort_sol_optl(complex long double *sol, complex double *exsol);

void print_legend(FILE *f)
{
  int ic;
  if (dojust < 0)
    {
      for (ic=0; ic < 8; ic++)
	fprintf(f,"@    s%d legend \"%s\"\n", ic, ic2algo(ic));
    }
  else
    fprintf(f,"@    s%d legend \"%s\"\n", dojust, ic2algo(dojust));
}
int maxic=7, icref;
char fname[256];
void save_PE(long long int numtrials, int numpts, double dlogdE, double logdEmin)
{
  FILE *f;
  int k, kk, ic;
  for (ic=0; ic < maxic; ic++)
    {
     if (dojust >= 0 && ic != dojust)
	continue;
      sprintf(fname,"P_of_eps_rel-%s.dat", ic2algo(ic));
      f = fopen(fname, "w+");
       for (k=0; k < numpts; k++)
	{
	  if (PEall[ic][k]==0) 
	    continue;
	  fprintf(f, "%.32G %.32G\n", k*dlogdE+logdEmin, PEall[ic][k]/((double)numtrials)/4.);
	}
      fclose(f);
    }
  if (dojust==-1 || (dojust >=0 && dojust==8))
    {
      f = fopen("P_of_eps_rel-HQRL.dat", "w+");
      for (k=0; k < numpts; k++)
	{
	  if (PEHQRL[k]==0)
	    continue;
	  fprintf(f, "%.32G %.32LG\n", k*dlogdE+logdEmin, PEHQRL[k]/((long double)numtrials)/4.);
	}
       fclose(f);
    }
  for (ic=0; ic < maxic; ic++)
    {
      if (dojust >= 0 && ic != dojust)
	continue;
      for (k=0; k < numpts; k++)
	{
	  cumPEall[ic][k] = 0.0;
	  for (kk=k; kk < numpts; kk++) 
	    {
	      cumPEall[ic][k] += PEall[ic][kk]/((double)numtrials)/4.0;
	    }
	}
    }
  if (dojust==-1 || (dojust >= 0 && dojust == 7))
    {
      for (k=0; k < numpts; k++)
	{
	  cumPEHQRL[k] = 0.0;
	  for (kk=k; kk < numpts; kk++) 
	    {
	      cumPEHQRL[k] += PEHQRL[kk]/((long double)numtrials)/4.0;
	    }
	}
    }
  for (ic=0; ic < maxic; ic++)
    {
      if (dojust >= 0 && ic != dojust)
	continue;
      
      sprintf(fname,"F_of_eps_rel-%s.dat", ic2algo(ic));
      f = fopen(fname, "w+");
      for (k=0; k < numpts; k++)
	{
	  if (cumPEall[ic][k]==0)
	    continue;
	  fprintf(f, "%.32G %.32G\n", k*dlogdE+logdEmin, cumPEall[ic][k]);
	}
      fclose(f);
    }
  if (dojust==-1 || (dojust >= 0 && dojust == 7))
    {
      f = fopen("F_of_eps_rel-HQRL.dat", "w+");
      if (cmplxreal!=5)
	{
	  for (k=0; k < numpts; k++)
	    {
	      if (cumPEHQRL[k]==0)
		continue;
	      fprintf(f, "%.32G %.32LG\n", k*dlogdE+logdEmin, cumPEHQRL[k]);
	    }
	}
      else
	fprintf(f,"-16 1\n&\n");
      fclose(f);
    }
}

int main(int argc, char **argv)
{
  complex long double x1c, x2c, x3c, x4c,csolHQRL[4]; 
  double logdE, dlogdE, logdEmax, logdEmin, sig, sig2, c[5];
  long double dE, x1, cl[5];
  long double y1;
  long long int numtrials, its, numout, itsI;
  int numpts, ilogdE;
  int num, k, k2, ic=0, okHQR, okHQRL, nsample;

  srand48(4242);
  
  for (ic=0; ic < 8; ic++)
    for (k=0; k < PEPTS; k++)
      PEall[ic][k] = 0;
  for (k=0; k < PEPTS; k++)
    PEHQRL[k] = 0;

  sig = 1.0;
  sig2= 1.0;
  logdEmax=10.0;
  logdEmin=-16.0;
  numpts = PEPTS; 
  dlogdE = (logdEmax -logdEmin)/numpts;

  if (argc>=2)
    numtrials=atoll(argv[1]);
  else 
    numtrials=1000000000;

  restart = 0;
  itsI = 0;
  if (argc>=3)
    numout=atoll(argv[2]);
  else
    numout=100;
  nsample=cmplxreal;
  if (argc>=4)
    cmplxreal = nsample = atoi(argv[3]);
  if (cmplxreal < 0 || cmplxreal > 5)
    {
      printf("cmplxreal must be between 0 and 5!\n");
      exit(-1);
    }
  if (cmplxreal==3)
    {
      sig = 1.0;
      sig2= 1E6;
      cmplxreal=1;
    }
  else if (cmplxreal==4)
    {
      sig = 1E6;
      sig2 = 1E6;
      cmplxreal = 2;
    } 
  if (argc  >= 5)
    dojust=atoi(argv[4]);
  if (dojust > 8)
    {
      printf("which test should I have to perform?!?\n Last arg is too big!\n");
      exit(-1);
    }
  if (numtrials < 0)
    {
      printf("number of trials must be a positive integer!\n");
      exit(-1);
    }
  else
    {
      printf("numtrials=%lld numout=%lld cmplxreal=%d dojust=%d\n", 
	     numtrials, numout, cmplxreal, dojust);
    }

  x1c=x2c=x3c=x4c=0;
  for (its=itsI; its < numtrials; its++)
    {
      if (its > 0 && (its % (numtrials/numout) == 0))
	{
          if (nsample == 0)
	    printf("[SAMPLE A sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
          else if (nsample==1)
	    printf("[SAMPLE B sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
	  else if (nsample==2)
	    printf("[SAMPLE C sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
          else if (nsample==3)
	    printf("[SAMPLE D sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
	  else if (nsample==4)
            printf("[SAMPLE E sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
          else if (nsample==5)
            printf("[SAMPLE F sig=%G %G]>>> its=%lld/%lld\n", sig, sig2, its, numtrials);
	  save_PE(its, numpts, dlogdE, logdEmin);
	  sync();  
	}
      /* generate 4 random roots */
      if (cmplxreal==2) /* 4 complex */
	{
	  x1 = sig2*(ranf()-0.5);
	  y1 = sig2*(ranf()-0.5);
	  x1c = x1 + I*y1;
	  x2c = x1 - I*y1;
	  x1 = sig*(ranf()-0.5);
	  y1 = sig*(ranf()-0.5);
	  x3c = x1 + I*y1;
	  x4c = x1 - I*y1;
	}
      else if (cmplxreal==1) /* two complex two real */
	{
	  x1 = sig2*(ranf()-0.5);
	  y1 = sig2*(ranf()-0.5);
	  x1c = x1 + I*y1;
	  x2c = x1 - I*y1;
	  x1 = sig*(ranf()-0.5);
	  y1 = sig*(ranf()-0.5);
	  x3c = x1;
	  x4c = y1;
	}
      else if (cmplxreal==0)/* four real */
	{
	  x1c = sig*(ranf()-0.5);
	  x2c = sig*(ranf()-0.5);
	  x3c = sig*(ranf()-0.5);
	  x4c = sig*(ranf()-0.5);
	}
      
      if (cmplxreal == 5)
	{
	  cl[4]=c[4]=1.0;
	  cl[3]=c[3]=ranf()-0.5;
	  cl[2]=c[2]=ranf()-0.5;
	  cl[1]=c[1]=ranf()-0.5;
	  cl[0]=c[0]=ranf()-0.5;
	}
      else
	{
	  c[4] = 1.0;
	  cl[4] = 1.0;
	  c[3] = creall(-(x1c+x2c+x3c+x4c));
	  c[2] = creall(x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c); 
	  c[1] = creall(-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c));
	  c[0] = creall(x1c*x2c*x3c*x4c);
	  cl[3] = creall(-(x1c+x2c+x3c+x4c));
	  cl[2] = creall(x1c*x2c + (x1c+x2c)*(x3c+x4c) + x3c*x4c); 
	  cl[1] = creall(-x1c*x2c*(x3c+x4c) - x3c*x4c*(x1c+x2c));
	  cl[0] = creall(x1c*x2c*x3c*x4c);
	  exsol[0] = x1c;
	  exsol[1] = x2c;
	  exsol[2] = x3c;
	  exsol[3] = x4c;
	}

      ic = 0;
      if (dojust==-1 || dojust == ic)
	solve_numrec(c, csolall[ic], &okHQR);
      ic++;
      if (dojust==-1 || dojust == ic)
	{
	  fast_quartic_solver(c, csolall[ic]); 
	}
      ic++;
      if (dojust==-1 || dojust == ic)
	csolve_quartic_abramovitz_cmplx(c, csolall[ic]);
      ic++;
       
      if (dojust==-1 || dojust == ic)
	CLDLT_quartic(c, csolall[ic]);     
      ic++;
      if (dojust==-1 || dojust == ic)
	CquarticRoots(c,&num, csolall[ic]);
      ic++;	
      if (dojust==-1 || dojust == ic)
	{
	 oqs_quartic_solver (c, csolall[ic]);
	}
      ic++;  
      if (dojust==-1 || dojust == ic)
	{
	  csolve_quartic_shmakov(c, csolall[ic]); 
	}
      ic++;
      if (dojust==-1 || dojust == ic || cmplxreal==5)
	{
	  solve_numrecl(cl, csolHQRL, 4, &okHQRL);
	}
      if (cmplxreal==5)
	{
	  for (k2=0; k2 < 4; k2++) 
	    exsol[k2] = csolHQRL[k2];
	}
      for (ic = 0; ic < maxic; ic++)
	{
	  if (dojust==-1 || dojust==ic)
	    sort_sol_opt(csolall[ic], exsol);
	  if (ic==0 && !okHQR)
	    continue;
	}
      if ((dojust == -1 || dojust == 7)&&okHQRL)
	sort_sol_optl(csolHQRL, exsol);

      if ((dojust==-1 || dojust==7) && okHQRL)
	{
	  for (k=0; k < 4; k++)
	    {
	      dE = (exsol[k]!=0)?cabsl((csolHQRL[k]- exsol[k])/exsol[k]):cabsl(csolHQRL[k]- exsol[k]); 
	      if (dE > 0.0)
		{
		  logdE=log10l(dE)-logdEmin;
		  ilogdE=(int)(logdE/dlogdE);
		  if (ilogdE >= 0 && ilogdE < numpts)
		    {
		      (PEHQRL[ilogdE])++;
		    }
		}
	    }
	}
      for (ic=0; ic < maxic; ic++)
	{
	  if (dojust >= 0 && dojust!=ic)
	    continue;
	  if (ic==0 && !okHQR)
	    continue;
	  for (k=0; k < 4; k++)
	    {	
	     dE = (exsol[k]!=0)?cabs((csolall[ic][k] - exsol[k])/exsol[k]):
	       cabs(csolall[ic][k]) - exsol[k]; 
	      if (dE > 0.0)
		{
		  logdE=log10(dE)-logdEmin;
		  ilogdE=(int)(logdE/dlogdE);
		  if (ilogdE >= 0 && ilogdE < numpts)
		    {
		      (PEall[ic][ilogdE])++;
		    }
		}
	    }
	}
    }
  save_PE(numtrials, numpts, dlogdE, logdEmin);
  printf("Finished\n");
  sync();
  exit(0);
}
