#include <stdio.h>
#include <memory.h>

/*type-independent absolute value, minimum, and maximum*/
#define rabs(x) ((x)*((x)>=0? 1.0: -1.0))
#define min(x,y) ((x)>(y)? (y):(x))
#define max(x,y) ((x)>(y)? (x):(y))

/*the floating-point type we use*/
typedef double real;
/*function type*/
typedef real (*function)(real*, int);
/*this type is for function evaluation counts so it has to be LARGE*/
typedef unsigned long long opcount;

/*In the descriptions of the following functions, if a and b are d-vectors*/
/*  and a_i <= b_i for i = 1,..d,*/
/*  then we define [a,b] as the domain [a_1,b_1]X[a_2,b_2]X...X[a_d,b_d]*/
/*  and [a-epsilon,b+epsilon] as the domain*/
/*  [a_1-epsilon,b_1+epsilon]X[a_2-epsilon,b_2+epsilon]X...X[a_d-epsilon,b_d+epsilon].*/

/*The parameters to fixedpoint are similar to those of fixedpoint_ab below, but with the absence*/
/*of a and b.  The unit hypercube [0,1]^d is used as the domain.*/

int fixedpoint(function f, real *retval, int d, real epsilon, opcount *ops)
{
  int i, result;
  real *a, *b;
  if (epsilon <= 0 || epsilon >= 0.5)
    {
      printf("ERROR in fixedpoint: epsilon must be in the interval (0,0.5) \n");
      return -1;
    }
/*unit cube*/
  a = (real*)malloc(d*sizeof(real));
  b = (real*)malloc(d*sizeof(real));
  for (i = 0; i < d; ++i)
    {
      a[i] = 0;
      b[i] = 1;
    }
  result = fixedpoint_ab(f,retval,d,a,b,epsilon,ops);
  free(a);
  free(b);
  return result;
}

/*The parameters to fixedpoint_ab are similar to those of fpoint below, but with the absence*/
/*of the internal variables ad, bd, xpl, xmi, and firstlarge, which are set by fixedpoint_ab.*/

int fixedpoint_ab(function f, real *retval, int d, real *a, real *b, real epsilon, opcount *ops)
{
  int i, result;
  real *ad, *bd, *xpl, *xmi;
  int firstlarge = d;
  if (epsilon <= 0)
    {
      printf("ERROR in fixedpoint_ab: epsilon must be >= 0\n");
      return -1;
    }
/*Allocate large arrays.  This avoids the need to allocate new*/
/*internal arrays at every recursion level.*/
/*Not needed if d=1.*/
/*Each recursive level needs the current dimension - 1 entries,*/
/*so the total number of entries is the sum of 1 through d-1.*/
  if (d>1)
    {
      ad = (real*)malloc(d*(d-1)/2*sizeof(real));
      bd = (real*)malloc(d*(d-1)/2*sizeof(real));
      xpl = (real*)malloc(d*(d-1)/2*sizeof(real));
      xmi = (real*)malloc(d*(d-1)/2*sizeof(real));
    }
/*Find the first component such that b_i-a_i > 2*epsilon.*/
  for (i = 0; i < d-1; ++i)
    {
      if (rabs(b[i]-a[i]) > 2.0*epsilon)
	{
	  firstlarge = i;
	  break;
	}
    }
  result = fpoint(f,retval,d,a,b,epsilon,ops,ad,bd,xpl,xmi,firstlarge);
  if (d>1)
    {
      free(ad);
      free(bd);
      free(xpl);
      free(xmi);
    }
  return result;
}

/*Given a function, fpoint returns a residual fixed point approximation.*/
/*f has type function(real*, int).  A call to f with a pointer to an array of d reals*/
/*  (representing the point x) and an index i in [0,..,d-1] should evaluate and return*/
/*  the index+1'th component of the function f at x.*/
/*  f must be Lipschitz continuous with constant 1, [a,b] must be its domain, and*/
/*  [a-epsilon,b+epsilon] must contain its range.*/
/*retval is a pointer to an array of d reals.  On output retval will be filled in with*/
/*  a point x* in [a,b] which satisfies the residual criterion ||f(x*)-x*|| <= epsilon.*/
/*d is the dimension (number of variables) of the function f.*/
/*a and b are points to arrays of d reals.  They define the domain of f which will be*/
/*  searched for a residual solution.  fpoint expects that a_i <= b_i for i = 1,..,d.*/
/*  NOTE: fpoint does not check that a_i <= b_i for i = 1,..,d, as this would be too*/
/*  expensive.  If this assertion does not hold then the results may be inaccurate.*/
/*epsilon is the error tolerance.*/
/*ops is a pointer to a large integer that will be set to the number of function component*/
/*  evaluations that are performed.*/
/*ad, bd, xpl, and xmi are pointers to arrays for internal use.  These are allocated*/
/*  by functions that call fpoint to avoid having to allocate them at every recursion level.*/
/*  If d = 1 then these arrays are not needed.  Otherwise they must each contain at least*/
/*  d*(d-1)/2 entries.  During each call to fpoint, d-1 elements of each of these arrays*/
/*  are available for use.*/
/*firstlarge is one less than the lowest nonnegative value of i such that b_i-a_i > 2*epsilon,*/
/*  or some integer >= d if ||b-a|| <= 2*epsilon.  Whenever ||b-a|| <= 2*epsilon the algorithm*/
/*  can take a "shortcut" by evaluating only the center point of the domain.*/
/*fpoint returns 0 on success, -1 on error.*/

int fpoint(function f, real *retval, int d, real *a, real *b, real epsilon, opcount *ops,
	   real *ad, real *bd, real *xpl, real *xmi, int firstlarge)
{
  int i,j;
  real val, diff, xa, xb, xm;
  opcount ops2;
  int result = 0;

/*parameter check*/
  if (f==NULL)
    {
      printf("ERROR in fpoint: f is null\n");
      return -1;
    }
  if (ops==NULL)
    {
      printf("ERROR in fpoint: ops is null\n");
      return -1;
    }
  if (retval==NULL)
    {
      printf("ERROR in fpoint: retval is null\n");
      return -1;
    }
  if (d<1)
    {
      printf("ERROR in fpoint: must have d > 0\n");
      return -1;
    }
  if (a==NULL)
    {
      printf("ERROR in fpoint: a is null\n");
      return -1;
    }
  if (b==NULL)
    {
      printf("ERROR in fpoint: b is null\n");
      return -1;
    }
  if (epsilon<=0.0)
    {
      printf("ERROR in fpoint: must have epsilon > 0\n");
      return -1;
    }

/*initialize ops*/
  *ops = 0;

  if (d==1)
    {
/*Univariate bisection (d=1)*/
      retval[0] = (a[0]+b[0])/2.0;
      xa = a[0];
      xb = b[0];
/*If a[0]=b[0] then take the center point.*/
/*If a[0]>b[0] do the same, but the result may be inaccurate.*/
/*Numerical roundoff may cause a[0] to be a tiny amount greater than b[0].*/
      if (xa>=xb) return 0;
      while (1)
	{
	  /* The univariate bisection loop follows the steps described in */
	  /* Shellman and Sikorski's PFix paper. */
	  /* Important note: If fpoint has been called recursively then */
	  /* f may be a function of more than one dimension.  In this case */
	  /* the elements of retval after the first contain the other */
	  /* components of the current evaluation point. */
	  val = (*f)(retval,0);
	  ++*ops;
	  diff = val-retval[0];
	  if (rabs(diff)<=epsilon) return 0;
	  if (xa==a[0] && retval[0]-xa<=epsilon && diff<0)
	    {
	      retval[0] = xa;
	      return 0;
	    }
	  if (xb==b[0] && xb-retval[0]<=epsilon && diff>0)
	    {
	      retval[0] = xb;
	      return 0;
	    }
	  if (diff>0)
	    {
	      xa = min(xb,retval[0]+diff/2.0);
	    }
	  else
	    {
	      xb = max(xa,retval[0]+diff/2.0);
	    }
	  if (xa>=xb)
	    {
	      retval[0] = xa;
	      return 0;
	    }
	  if (xa==a[0] && xb-xa<=epsilon/2.0 && diff<0)
	    {
	      retval[0] = xa;
	      return 0;
	    }
	  if (xb==b[0] && xb-xa<=epsilon/2.0 && diff>0)
	    {
	      retval[0] = xb;
	      return 0;
	    }
	  retval[0] = (xa+xb)/2.0;
	  if (xa!=a[0] && xb!=b[0] && xb-xa<=epsilon) return 0;
	}
    }

/*The program will not reach this point if d=1.*/
  if (firstlarge >= d)
    {
/*||b-a|| is smaller than 2*epsilon.*/
/*As described in the PFix paper, we can take a shortcut by*/
/*evaluating the function at the center point of the domain.*/
/*Since d-1 elements of the ad array are available for use,*/
/*we make use of it as a temporary array, along with an*/
/*additional variable.  This avoids the need for an allocation.*/
      real dth;
      for (i = 0; i < d; ++i)
	{
	  retval[i] = (a[i]+b[i])/2.0;
	}
      for (i = 0; i < d-1; ++i)
	{
	  if (a[i] >= b[i])
	    {
	      ad[i] = a[i];
	    }
	  else
	    {
	      val = (*f)(retval,i);
	      ++*ops;
	      ad[i] = max(a[i],min(b[i],val));
	    }
	}
      if (a[d-1] >= b[d-1])
	{
	  dth = a[d-1];
	}
      else
	{
	  val = (*f)(retval,d-1);
	  ++*ops;
	  dth = max(a[d-1],min(b[d-1],val));
	}
      for (i = 0; i < d-1; ++i)
	{
	  retval[i] = ad[i];
	}
      retval[d-1] = dth;
      return 0;
    }

/*Note that in the following recursive calls, we are passing a d-dimensional*/
/*function and a dimension of d-1.  This is acceptable because retval has d*/
/*components.  The d'th component of retval will remain fixed throughout the*/
/*recursive call, and will be passed to f along with the first d-1 components.*/
/*Note also that we pass the unused portions of the ad, bd, xpl, and xmi*/
/*arrays into the recursive calls, so that they do not need to allocate*/
/*extra storage.*/

/*We recursively search the middle plane in the d'th component for a point*/
/*that is residual in the first d-1 components of f.*/
  retval[d-1] = (a[d-1]+b[d-1])/2.0;
  result = fpoint(f, retval, d-1, a, b, epsilon, &ops2,
		  ad+d-1, bd+d-1, xpl+d-1, xmi+d-1, firstlarge);
  *ops += ops2;
  if (result) goto ending;

/*Multivariate bisection (d>1)*/
/*The multivariate bisection loop follows the steps described in*/
/*Shellman and Sikorski's PFix paper.*/
/*At this point retval contains a point that is residual in the first d-1*/
/*components of f, and whose d'th component is halfway between a_d and b_d.*/
  xa = a[d-1];
  xb = b[d-1];
  if (xa>=xb) goto ending;
  while (1)
    {
      val = (*f)(retval,d-1);
      ++*ops;
      diff = val-retval[d-1];
      if (rabs(diff)<=epsilon) break;
/*The xmi and xpl arrays correspond to the x^- and x^+ variables in the*/
/*PFix paper.  xa and xb are their d-1'th components, respectively.*/
      if (diff>0)
	{
	  xa = retval[d-1];
	  for (j = 0; j < d-1; ++j)
	    {
	      xmi[j] = retval[j];
	    }
	}
      else
	{
	  xb = retval[d-1];
	  for (j = 0; j < d-1; ++j)
	    {
	      xpl[j] = retval[j];
	    }
	}
      if (xa==a[d-1] && xb-xa<=epsilon)
	{
	  for (j = 0; j < d-1; ++j)
	    {
	      ad[j] = max(a[j],retval[j]-(xb-xa));
	      bd[j] = min(b[j],retval[j]+(xb-xa));
	    }
	  retval[d-1] = xa;
	  /* We pass firstlarge=d-1 because we know that ||bd-ad|| <= 2*epsilon. */
	  result = fpoint(f, retval, d-1, ad, bd, epsilon, &ops2,
			  ad+d-1, bd+d-1, xpl+d-1, xmi+d-1, d-1);
	  *ops += ops2;
	  break;
	}
      if (xb==b[d-1] && xb-xa<=epsilon)
	{
	  for (j = 0; j < d-1; ++j)
	    {
	      ad[j] = max(a[j],retval[j]-(xb-xa));
	      bd[j] = min(b[j],retval[j]+(xb-xa));
	    }
	  retval[d-1] = xb;
	  /* We pass firstlarge=d-1 because we know that ||bd-ad|| <= 2*epsilon. */
	  result = fpoint(f, retval, d-1, ad, bd, epsilon, &ops2,
			  ad+d-1, bd+d-1, xpl+d-1, xmi+d-1, d-1);
	  *ops += ops2;
	  break;
	}
      xm = (xa+xb)/2.0;
      if (xa==a[d-1] || xb==b[d-1])
	{
	  /* While computing the recursive domain, also find the first */
	  /* component such that bd[j]-ad[j] > 2*epsilon. */
	  firstlarge = d-1;
	  for (j = 0; j < d-1; ++j)
	    {
	      ad[j] = max(a[j],retval[j]-(xb-xm));
	      bd[j] = min(b[j],retval[j]+(xb-xm));
	      if (firstlarge == d-1 && bd[j]-ad[j] > 2.0*epsilon)
		{
		  firstlarge = j;
		}
	    }
	  retval[d-1] = xm;
	  result = fpoint(f, retval, d-1, ad, bd, epsilon, &ops2,
			  ad+d-1, bd+d-1, xpl+d-1, xmi+d-1, firstlarge);
	  *ops += ops2;
	  if (result) goto ending;
	}
      else
	{
	  /* While computing the recursive domain, also find the first */
	  /* component such that bd[j]-ad[j] > 2*epsilon. */
	  firstlarge = d-1;
	  for (j = 0; j < d-1; ++j)
	    {
	      ad[j] = max(a[j],max(xpl[j],xmi[j])-(xb-xm));
	      bd[j] = min(b[j],min(xpl[j],xmi[j])+(xb-xm));
	      if (firstlarge == d-1 && bd[j]-ad[j] > 2.0*epsilon)
		{
		  firstlarge = j;
		}
	    }
	  retval[d-1] = xm;
	  result = fpoint(f, retval, d-1, ad, bd, epsilon, &ops2,
			  ad+d-1, bd+d-1, xpl+d-1, xmi+d-1, firstlarge);
	  *ops += ops2;
	  if (result || xb-xa<=2.0*epsilon)
	    {
	      break;
	    }
	}
    }

 ending:
  return result;
}
