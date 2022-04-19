#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>

#define rabs(x) ((x)*((x)>=0? 1.0: -1.0))
#define min(x,y) ((x)>(y)? (y):(x))
#define max(x,y) ((x)>(y)? (x):(y))

typedef double real;
typedef real (*function)(real*, int);
typedef unsigned long long opcount;

/* maxops computes and returns the function C(d,m) = binom(d+m-1,m-1).*/
/* d and m must be >= 1.*/
/* By Appendix A of "A Recursive Algorithm for the Infinity-Norm*/
/* Fixed Point Problem" by Shellman and Sikorski, this function is*/
/* equivalent to sigma(d,m), where sigma(d,m) is defined as*/
/* sigma(d,m) = m if d = 1,*/
/* sigma(d,m) = sum(i = 1 to m) of sigma(d-1,i) if d > 1.*/
/* maxops computes C(d,m) by computing sigma(d,m) recursively.*/
opcount maxops(int d, int m)
{
/* A pointer to a table of sigma(d,m) values that are computed*/
/* by recursive calls to maxops.*/
  static opcount *table = NULL;
/* The top-level value of m, used to index into the sigma table.*/
  static int maxm;
  opcount *oldtable = table;
  opcount ops = 0;
  int i, index;
  if (d < 1)
    {
      printf("ERROR in maxops: d must be >= 1\n");
      return 0;
    }
  if (m < 1)
    {
      printf("ERROR in maxops: m must be >= 1\n");
      return 0;
    }
/* Simplest case: sigma(1,m) = m.*/
  if (d==1) return m;
  if (!oldtable)
    {
/* This is the top-level call to maxops.*/
/* Here we allocate and initialize the sigma table.*/
/* A zero value in the table indicates that sigma*/
/* has not yet been computed for this combination*/
/* of d and m values.*/
/* Naturally, there is no need to include a row*/
/* in the table for d = 1.*/
      maxm = m;
      table = (opcount*)malloc((d-1)*maxm*sizeof(opcount));
      for (i = 0; i < (d-1)*maxm; ++i) table[i] = 0;
    }
/* If sigma has already been computed for this combination*/
/* of d and m, then it will be at this index.*/
  index = (d-2)*maxm+(m-1);
  if (table[index] != 0) return table[index];
/* Call maxops recursively for d-1, 1..m.*/
/* Sum the results to get sigma(d,m).*/
  for (i = 1; i <= m; ++i)
    {
      ops += maxops(d-1,i);
    }
/* Add sigma(d,m) to the sigma table.*/
  table[index] = ops;
  if (!oldtable)
    {
/* Deallocate the sigma table if this is the top-level call.*/
      free(table);
      table = NULL;
    }
  return ops;
}

/* Parameters are similar to those of fptest.*/
/* Performs simple iteration.  Calls x^{k+1}=f(x^k)*/
/* log(epsilon)/log(q) times or until ||f(x^k)-x^k)|| <= epsilon(1-q).*/
/* Result satisfies either relative or absolute criterion.*/
/* q is the Lipschitz constant; must be < 1.*/

int simple(function f, real *retval, int d, real epsilon, real q,
	   opcount *ops)
{
  int k;
  opcount i, upper = (opcount)ceil(log(epsilon) / log(q));
  real *nextf;
  real maxdiff;
  if (!retval || !ops)
    {
      printf("ERROR in simple: null pointer\n");
      return 1;
    }
  if (q >= 1)
    {
      printf("ERROR in simple: q must be < 1 for SI\n");
      return 1;
    }
  *ops = 0;
  for (k = 0; k < d; ++k) retval[k] = 0.5;
  if (epsilon < 0) epsilon = -epsilon;
  if (epsilon >= 0.5) return 0;
  nextf = (real*)malloc(d*sizeof(real));
  for (i = 0; i < upper; ++i)
    {
      for (k = 0; k < d; ++k) nextf[k] = (*f)(retval,k);
      *ops += d;
      maxdiff = 0.0;
      for (k = 0; k < d; ++k) maxdiff = max(maxdiff, rabs(nextf[k]-retval[k]));
      if (maxdiff<=epsilon*(1-q)) break;
      for (k = 0; k < d; ++k) retval[k] = nextf[k];
    }
  free(nextf);
  return 0;
}

/* dimension of a test function*/
static int dim;
/* Lipschitz constant of a test function*/
static real lipq;

/* The following test functions correspond to the similarly named functions in the PFix paper.*/

real f1(real *params, int index)
{
  real retval;
  if (index%2==0)
    retval = sin((params[index]+params[(index+1)%dim]+params[(index+2)%dim])/3.0+0.4);
  else
    retval = log((params[index]+1.0)*(params[(index+1)%dim]+1.0)*(params[(index+2)%dim]+1.0))/3.0+0.1;
  return retval;
}

real f2(real *params, int index)
{
  real retval;
  int i,j;
/* saves the old dimension*/
  static int olddim = 0;
  static real *c = NULL;
  if (dim != olddim)
    {
/* The dimension has changed, so reallocate and reinitialize the c table.*/
      real dr = (real)dim;
      olddim = dim;
      free(c);
      c = (real*) malloc(dim*dim*sizeof(real));
      for (i = 0; i < dim; ++i)
	{
	  for (j = 0; j < dim; ++j)
	    {
	      c[i*dim+j] = 0.5-0.5*(i-dr/2.0)*(j-dr/2.0)/((dr/2.0)*(dr/2.0));
	    }
	}
    }
  retval = 0.0;
  for (i = 0; i < dim; ++i)
    {
      retval = max(retval,rabs(params[i]-c[index*dim+i]));
    }
  retval = max(0,1.0-lipq*retval);
  return retval;
}

real f3(real *params, int index)
{
  return 0;
}

real f4(real *params, int index)
{
/* f2(f1(.))*/
  int i;
  real params2[dim];
  for (i = 0; i < dim; ++i)
    {
      params2[i] = f1(params,i);
    }
  return f2(params2,index);
}


/* Runs a PFix test on the domain [0,1]^d with eps < 0.5.*/
/* fname points to a string containing a name for the function.*/
/* f is the test function to use.  PFix computes an approximation of a fixed point of f.*/
/* d is the dimension (number of variables) of f.*/
/* eps (0<eps<0.5) is the error tolerance for the test.*/
/* If residual = 1 then PFix computes a residual criterion approximation.*/
/* If residual = 0 then PFix computes an absolute criterion approximation, and its performance*/
/* is compared to SI for the same problem.*/
int pfixtest(char *fname, function f, int d, real eps, real q, int residual)
{
  real fp[d];
  opcount ops, mops, sops;
  int i, mi, r;
  real v, m;
  time_t start, end;
  real elapsed;

  if (q <= 0.0 || q > 1.0)
    {
      printf("ERROR in pfixtest: q must be in the interval (0,1]\n");
      return -1;
    }
  if (eps <= 0.0 || eps >= 0.5)
    {
      printf("ERROR in pfixtest: eps must be in the interval (0,0.5)\n");
      return -1;
    }
  r = (int) ceil(-log(eps)/log(2.0));
/* set the static variables for use by the test function*/
  dim = d;
  lipq = q;
  printf("\nPFix test: %s criterion\nfunction=%s, dimension=%d, epsilon=%e, r=%d, q=%e\n",
	 residual? "residual": "absolute", fname == NULL? "": fname, d, eps, r, q);
/* compute the upper complexity bound according to the PFix paper*/
  mops = maxops(dim, r) - maxops(dim-1, r) + 2*maxops(dim-1,r+2) - 2*maxops(dim-2,r+2);
  printf("PFix max operations = %llu\n", mops);
  start = time(NULL);

/* set epsilon to epsilon*(1-q) for absolute case with q<1*/
  fixedpoint(f, fp, dim, (!residual && q < 1.0)? eps*(1.0-q): eps, &ops);
  end = time(NULL);
  elapsed = difftime(end, start);
  printf("PFix number of operations = %llu, seconds elapsed = %d\n", ops, (int) ceil(elapsed));
  m = 0.0;
  mi = 0;
  for (i = 0; i < dim; ++i)
    {
      v = f(fp,i);
      if (rabs(fp[i]-v)>m)
	{
	  m = rabs(fp[i]-v);
	  mi = i;
	}
    }
  printf("max difference = %e at component %d\n", m, mi);
  printf("ratio of number of ops to max = %e\n", (real)ops / (real) mops);
  if (m<=eps)
    {
      printf("satisfies residual criterion\n");
    }
  else
    {
      printf("WARNING: does not satisfy residual criterion\n");
    }

  if (!residual)
    {
/* execute SI*/
      printf("SI max operations = %llu\n", (opcount)dim*(opcount)ceil(log(eps) / log(lipq)));
      start = time(NULL);
      simple(f, fp, dim, eps, lipq, &sops);
      end = time(NULL);
      elapsed = difftime(end, start);
      printf("SI number of operations = %llu, seconds elapsed = %d\n", sops, (int) ceil(elapsed));
    }

  return 0;
}

int main()
{
/*The following tests are labeled with the lengths of their execution times (hh:mm:ss) on an*/
/*Intel Xeon 1.7 GHz system running Red Hat Linux.*/
/*These tests correspond to those described in the PFix paper.*/
/*Residual criterion tests:*/
  pfixtest("f1", f1, 6, 0.0000000000001, 1.0, 1);    /*  00:00:01 */
  pfixtest("f1", f1, 35, 0.00001, 1.0, 1);           /*  00:02:52 */
  pfixtest("f1", f1, 100, 0.001, 1.0, 1);            /*  00:06:44 */
  pfixtest("f1", f1, 1000, 0.025, 1.0, 1);           /*  00:02:02 */
  pfixtest("f2", f2, 6, 0.0000000000001, 1.0, 1);    /*  00:00:01 */
  pfixtest("f2", f2, 35, 0.00001, 1.0, 1);           /*  00:00:19 */
  pfixtest("f2", f2, 100, 0.001, 1.0, 1);            /*  00:05:55 */
  pfixtest("f2", f2, 1000, 0.025, 1.0, 1);           /*  15:20:10 */
  pfixtest("f3", f3, 6, 0.0000000000001, 1.0, 1);    /*  00:00:01 */
  pfixtest("f3", f3, 35, 0.00001, 1.0, 1);           /*  00:41:19 */
  pfixtest("f3", f3, 100, 0.001, 1.0, 1);            /*  00:11:12 */
  pfixtest("f3", f3, 1000, 0.025, 1.0, 1);           /*  02:40:02 */
  pfixtest("f4", f4, 6, 0.0000000000001, 1.0, 1);    /*  00:00:02 */
  pfixtest("f4", f4, 35, 0.00001, 1.0, 1);           /*  00:03:43 */
  pfixtest("f4", f4, 100, 0.001, 1.0, 1);            /*  00:12:38 */
  pfixtest("f4", f4, 1000, 0.025, 1.0, 1);           /*  00:24:24 */
/*Absolute criterion tests (compared with SI):*/
  pfixtest("f5", f2, 5, 0.001, 0.99, 0);             /*  PFix: 00:00:01  SI: 00:00:01 */
  pfixtest("f5", f2, 5, 0.00001, 0.9999, 0);         /*  PFix: 00:00:01  SI: 00:00:01 */
  pfixtest("f5", f2, 5, 0.0000001, 0.999999, 0);     /*  PFix: 00:00:01  SI: 00:00:38 */
  pfixtest("f5", f2, 10, 0.001, 0.99, 0);            /*  PFix: 00:00:01  SI: 00:00:01 */
  pfixtest("f5", f2, 10, 0.00001, 0.9999, 0);        /*  PFix: 00:00:01  SI: 00:00:01 */
  pfixtest("f5", f2, 10, 0.0000001, 0.999999, 0);    /*  PFix: 00:00:01  SI: 00:01:47 */
  pfixtest("f5", f2, 25, 0.001, 0.99, 0);            /*  PFix: 00:00:01  SI: 00:00:01 */
  pfixtest("f5", f2, 25, 0.00001, 0.9999, 0);        /*  PFix: 00:00:07  SI: 00:00:04 */
  pfixtest("f5", f2, 25, 0.0000001, 0.999999, 0);    /*  PFix: 00:01:13  SI: 00:09:20 */
}
