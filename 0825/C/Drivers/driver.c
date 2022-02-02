#include <stdio.h>
#include <math.h>
#include <limits.h>

typedef unsigned long long ulong;
typedef double (*bivariate_func)(double, double);

//returns minimum of a and b
static double min(double a, double b)
{
  return a < b? a: b;
}

//returns maximum of a and b
static double max(double a, double b)
{
  return a > b? a: b;
}


// THE CODE BELOW THIS LINE IS FOR TESTING PURPOSES.

// The function fptest is for testing purposes.
// It is the same as the previous function fp, but has a 'deep'
// parameter for purposes of comparison between BEDFix and BEFix.
// The previous fp function assumes that deep=1.

// Given tolerance epsilon and bivariate functions f1 and f2,
// returns in *fx, *fy a fixed point approximation satisfying
// the residual criterion
// |f1(*fx,*fy)-*fx| <= epsilon and |f2(*fx,*fy)-*fy| <= epsilon.
// f1, f2 have domain [0, 1] X [0, 1], range [0, 1], and are
// Lipschitz continuous with constant q<=1.
// Sets *total to the number of evaluations of f = (f1, f2).
// Sets *absolute to true if the method is able to determine
// that the solution also satisfies the absolute error
// criterion.  (There may be cases where the solution
// satisfies the absolute criterion but *absolute is
// returned false.  Note that if q=1 there is no way to
// computationally verify whether the solution satisfies
// the absolute criterion, unless the location of fixed
// points is known beforehand.)
// Input deep = 1 to turn on deep cuts (BEDFix), 0 to turn them off
// (BEFix).
// Return code: 0 normal, otherwise error
// To guarantee that the solution satisfies the absolute
// criterion |*fx-x| <= epsilon and |*fy-y| <= epsilon,
// pass in epsilon(1-q) for the epsilon parameter,
// where q<1 is the Lipschitz constant of f.

int fptest(double epsilon, bivariate_func f1, bivariate_func f2,
       double *fx, double *fy, ulong *total, int *absolute, int deep)
{
  //midpoint and diagonal extents of the feasible domain
  double xm = 0.5, ym = 0.5, ext1 = 0.5, ext2 = 0.5;
  if (!fx || !fy || !total || !absolute)
    {
      printf("ERROR: null pointer\n");
      return 1;
    }
  *total = 0;
  *absolute = 0;
  if (epsilon < 0) epsilon = -epsilon;
  if (epsilon >= 0.5)
    {
      *fx = xm;
      *fy = ym;
      *absolute = 1;
      return 0;
    }
  do
    {
      double v1 = (*f1)(xm, ym) - xm, v2 = (*f2)(xm, ym) - ym,
        fv1 = fabs(v1), fv2 = fabs(v2);
      double m_ru = 0.0, m_rd = 0.0, m_lu = 0.0, m_ld = 0.0;
      int goleft = 1, goright = 1, goup = 1, godown = 1;
      ++*total;
      if (fv1 <= epsilon && fv2 <= epsilon)
        {
          //Found a residual solution.
          *fx = xm;
          *fy = ym;
          if (v1 == 0 && v2 == 0) *absolute = 1;
          return 0;
        }
      if (ext1 > epsilon / 4.0 && ext2 > epsilon / 4.0)
        {
          double m = min(fv1,fv2) / 2.0;
          if (v1 > 0.0 && v2 == 0.0)
            {
              goleft = goup = godown = 0;
            }
          else if (v1 < 0.0 && v2 == 0.0)
            {
              goright = goup = godown = 0;
            }
          else if (v1 == 0.0 && v2 > 0.0)
            {
              goleft = goright = godown = 0;
            }
          else if (v1 == 0.0 && v2 < 0.0)
            {
              goleft = goright = goup = 0;
            }
          else if (v1 > 0.0 && v2 > 0.0)
            {
              goleft = godown = 0;
              m_ru = m;
            }
          else if (v1 < 0.0 && v2 > 0.0)
            {
              goright = godown = 0;
              m_lu = m;
            }
          else if (v1 < 0.0 && v2 < 0.0)
            {
              goright = goup = 0;
              m_ld = m;
            }
          else if (v1 > 0.0 && v2 < 0.0)
            {
              goleft = goup = 0;
              m_rd = m;
            }
        }
      if (ext1 <= epsilon / 2.0 && ext2 > epsilon / 4.0)
        {
          if (v1 < -epsilon || v2 > epsilon)
            {
              goright = godown = 0;
              m_lu = max(m_lu,max(-v1/2.0-ext1,v2/2.0-ext1));
            }
          else
            {
              goleft = goup = 0;
              m_rd = max(m_rd,max(v1/2.0-ext1,-v2/2.0-ext1));
            }
        }
      if (ext2 <= epsilon / 2.0 && ext1 > epsilon / 4.0)
        {
          if (v1 > epsilon || v2 > epsilon)
            {
              goleft = godown = 0;
              m_ru = max(m_ru,max(v1/2.0-ext2,v2/2.0-ext2));
            }
          else
            {
              goright = goup = 0;
              m_ld = max(m_ld,max(-v1/2.0-ext2,-v2/2.0-ext2));
            }
        }
      if (!goleft && !goright && !goup && !godown)
        {
          printf("ERROR: entire domain eliminated\n");
          return -1;
        }
      if (!goleft && !goright && !goup)
        {
          //down
          ext1 /= 2.0;
          ext2 /= 2.0;
          xm -= ext1 - ext2;
          ym -= ext1 + ext2;
        }
      else if (!goleft && !goright && !godown)
        {
          //up
          ext1 /= 2.0;
          ext2 /= 2.0;
          xm += ext1 - ext2;
          ym += ext1 + ext2;
        }
      else if (!goleft && !goup && !godown)
        {
          //right
          ext1 /= 2.0;
          ext2 /= 2.0;
          xm += ext1 + ext2;
          ym += ext1 - ext2;
        }
      else if (!goright && !goup && !godown)
        {
          //left
          ext1 /= 2.0;
          ext2 /= 2.0;
          xm -= ext1 + ext2;
          ym -= ext1 - ext2;
        }
      else if (!goleft && !goup)
        {
          //right, down
          ext2 /= 2.0;
          xm += ext2;
          ym -= ext2;
        }
      else if (!goleft && !godown)
        {
          //right, up
          ext1 /= 2.0;
          xm += ext1;
          ym += ext1;
        }
      else if (!goright && !goup)
        {
          //left, down
          ext1 /= 2.0;
          xm -= ext1;
          ym -= ext1;
        }
      else if (!goright && !godown)
        {
          //left, up
          ext2 /= 2.0;
          xm -= ext2;
          ym += ext2;
        }
      else
        {
          printf("ERROR: incorrect domain reduction\n");
          return -1;
        }
      if (deep)
        {
          ext2 -= m_rd / 2.0;
          xm += m_rd / 2.0;
          ym -= m_rd / 2.0;
          ext1 -= m_ru / 2.0;
          xm += m_ru / 2.0;
          ym += m_ru / 2.0;
          ext1 -= m_ld / 2.0;
          xm -= m_ld / 2.0;
          ym -= m_ld / 2.0;
          ext2 -= m_lu / 2.0;
          xm -= m_lu / 2.0;
          ym += m_lu / 2.0;
        }
    }
  while (ext1 + ext2 > epsilon / 2.0);
  *fx = xm;
  *fy = ym;
  *absolute = 1;
  return 0;
}

// Parameters are similar to those of fptest.
// Performs simple iteration.  Calls x^{k+1}=f(x^k)
// log(epsilon)/log(q) times or until ||f(x^k)-x^k)|| <= epsilon(1-q).
// Result satisfies either relative or absolute criterion.
// q is the Lipschitz constant; must be < 1.

int simple(double epsilon, double q, bivariate_func f1, bivariate_func f2,
       double *fx, double *fy, ulong *total)
{
  double xm = 0.5, ym = 0.5;
  ulong i, upper = (ulong)ceil(log(epsilon) / log(q));
  if (!fx || !fy || !total)
    {
      printf("ERROR: null pointer\n");
      return 1;
    }
  if (q >= 1)
    {
      printf("ERROR: q must be < 1 for SI\n");
      return 1;
    }
  *total = 0;
  if (epsilon < 0) epsilon = -epsilon;
  *fx = xm;
  *fy = ym;
  if (epsilon >= 0.5) return 0;
  for (i = 0; i < upper; ++i)
    {
      xm = (*f1)(*fx,*fy);
      ym = (*f2)(*fx,*fy);
      ++*total;
      if (max(fabs(*fx-xm),fabs(*fy-ym))<=epsilon*(1-q)) break;
      *fx = xm;
      *fy = ym;
    }
  return 0;
}





/* parameters for pyramid functions */
typedef struct {
  double cx, dx, hx;
} testcase;

static testcase tests[] = {
  {0.5, 0.5, 0.8},
  {0.6, 0.4, 1.2},
  {0.4, 0.6, 0.9},
  {0.6, 0.98, 0.99},
  {0.98, 0.3, 0.99},
  {0.27, 0.64, 1.01},
  {0.64, 0.27, 0.99},
  {0.0, 0.0, 0.1}
};

#define TESTCASE_COUNT (sizeof(tests) / sizeof(testcase))

/* Lipschitz constant */
static double lc = 1.0;

/* number of test cases to apply, and array of test indices */
static unsigned int select1, select2;

double f(double x, double y, unsigned int select)
{
  int i;
  unsigned int b = 1 << (TESTCASE_COUNT - 1);
  double result = 0;
  x = min(1.0, max(0.0, x));
  y = min(1.0, max(0.0, y));
  for (i = 0; b; b >>= 1, ++i)
    {
      if (select & b)
        {
          result = max(result,
                       min(1, max(tests[i].hx
                                  - max(fabs(x - tests[i].cx),
                                        fabs(y - tests[i].dx))*lc, 0)));
        }
    }
  if (result < 0.0 || result > 1.0)
    {
      printf("ERROR: result out of range\n");
      exit(-1);
    }
  return result;
}

double f1(double x, double y)
{
  return f(x, y, select1);
}

double f2(double x, double y)
{
  return f(x, y, select2);
}

int report_test(ulong l2, double epsilon, bivariate_func f1,
                bivariate_func f2, ulong *average, ulong *absaverage,
                ulong *minimum, ulong *maximum, ulong *testcount,
                ulong *wrongcount, ulong *exactcount,ulong *abscount,
                int output, int deep, int forceabs, double q, int si)
{
  double fx, fy, f1xy, f2xy;
  ulong total;
  int absolute, result;
  // modify the tolerance if computing absolute solution
  if (!si && forceabs) epsilon *= 1-q;
  if (si)
    result = simple(epsilon,q,f1,f2,&fx,&fy,&total);
  else
    result = fptest(epsilon,f1,f2,&fx,&fy,&total,&absolute,deep);
  if (si || forceabs) absolute = 0;
  if (!result)
    {
      double err1, err2, ratio;
      fx = min(1.0, max(0.0, fx));
      fy = min(1.0, max(0.0, fy));
      if (output) printf("\nx=%f y=%f\n",fx,fy);
      if (!si)
	{
	  // verify that residual criterion is satisfied
	  f1xy = f1(fx,fy);
	  f2xy = f2(fx,fy);
	  if (output) printf("f1(x,y)=%f f2(x,y)=%f\n", f1xy, f2xy);
	  err1 = fabs(fx-f1xy)/2.0;
	  err2 = fabs(fy-f2xy)/2.0;
	  if (output) printf("err1=%f err2=%f\n", err1, err2);
	  if (err1 > epsilon || err2 > epsilon)
	    {
	      ++*wrongcount;
	      if (output) printf("ERROR: criterion not satisfied\n");
	    }
	}
      if (output) printf("total number of evaluations = %d\n", total);
      *minimum = min(*minimum, total);
      *maximum = max(*maximum, total);
      *average += total;
      if (absolute) *absaverage += total;
      if (!si && !forceabs)
	{
	  ratio = (double) total / (double) l2;
	  if (output && !forceabs)
	    printf("ratio of evals to 2*log_2(1/epsilon)+1 = %f\n", ratio);
	  if (output && absolute)
	    printf("satisfies absolute criterion\n");
	  if (ratio > 1.0)
	    {
	      if (output) printf("ERROR: ratio > 1\n");
	      ++*wrongcount;
	    }
	  if (ratio == 1.0) ++*exactcount;
	}
      if (output) printf("\n\n\n");
      ++*testcount;
      if (absolute) ++*abscount;
    }
  else
    {
      if (output) printf("\nERROR: failed to find fixed point\n\n\n\n");
      ++*wrongcount;
    }
  return result;
}

// Run all tests and display statistics.
// epsilon = error tolerance
// output = 0 if only final statistics are desired, 1 if each test should
//   output its statistics as well
// deep = 1 if deep cuts (BEDFix), 0 if no deep cuts (BEFix)
// q = Lipschitz constant (must be < 1 if forceabs = 1 or si = 1)
// forceabs = 1 if the BE(D)Fix result must satisfy the absolute criterion,
//   0 if the residual criterion
// si = 1 if simple iteration, 0 if BE(D)Fix
void test(double epsilon, int output, int deep, double q, int forceabs, int si)
{
  ulong testno = 1, testcount = 0, wrongcount = 0, exactcount = 0,
    abscount = 0, average = 0, absaverage = 0,
    maximum = 0, minimum = ULONG_MAX, l2;
  // make sure that q < 1 if absolute solution is sought
  if (q >= 1 && (si || forceabs))
    {
      printf("ERROR: q must be < 1 for absolute criterion or SI\n");
      return;
    }
  // set the global Lipschitz constant
  lc = q;
  l2 = 2*(ulong)ceil(log(1/epsilon) / log(2)) + 1;
  // Runs tests using all possible combinations of tests from the list of
  // test cases.
  for (select1 = 1; select1 < 1 << TESTCASE_COUNT; ++select1)
    {
      for (select2 = 1; select2 < 1 << TESTCASE_COUNT; ++select2)
        {
          int c;
          if (output) printf("pyramid test %llu: %s%s, epsilon = %f\n",
			     testno, si? "SI": (deep? "BEDFix": "BEFix"),
			     si? "": (forceabs? " absolute criterion":
			     " residual criterion"), epsilon);
          ++testno;
          if (output) for (c = 1; c <= 2; ++c)
            {
              unsigned int select = c == 1? select1: select2;
              unsigned int b = 1;
              int i = 0;
              printf("f%d:\n",c);
              for (; b < 1 << TESTCASE_COUNT; b <<= 1, ++i)
                {
                  if (select & b)
                    {
                      printf("cx%d=%f dx%d=%f hx%d=%f\n",
                             i+1, tests[i].cx, i+1,
                             tests[i].dx, i+1, tests[i].hx);
                    }
                }
            }
          report_test(l2,epsilon,f1,f2,&average,&absaverage,&minimum,
                      &maximum,&testcount,&wrongcount,&exactcount,
                      &abscount,output,deep,forceabs,q,si);
        }
    }

  printf("\n\n");
  printf("%s%s method, epsilon = %f, q = %f\n", si? "SI":
	 (deep? "BEDFix": "BEFix"), si? "":
	 (forceabs? " absolute criterion": " residual criterion"),
	 epsilon, q);
  printf("number of tests performed = %llu\n", testcount);
  if (!si)
    printf("number of tests with errors or results not " \
	   "satisfying criterion = %llu\n", wrongcount);
  printf("average number of evaluations = %f\n",
         (double) average / (double) testcount);
  if (!forceabs && !si)
    {
      printf("2*log_2(1/epsilon)+1=%llu\n", l2);
      printf("number of tests with 2*log_2(1/epsilon)+1 evaluations = %llu\n",
             exactcount);
      printf("number of tests satisfying absolute error criterion = %llu\n",
	     abscount);
      printf("average ratio of evaluations to 2*log_2(1/epsilon)+1 = %f\n",
             (double) average / ((double) (testcount * l2)));
      printf("average ratio of evaluations to 2*log_2(1/epsilon)+1\n" \
	     "when absolute criterion satisfied = %f\n",
	     (double) absaverage / ((double) (abscount * l2)));
      printf("minimum and maximum function evaluations = %llu, %llu\n",
             minimum, maximum);
    }
  if (!si && wrongcount>0) printf("ERROR: some tests did not run correctly\n" \
				  "turn on per-test output to see details\n");
  printf("\n\n\n\n");
}

// Runs large test if whichtest=1, small test if whichtest=2.
void testall(int whichtest)
{
  // no per-test output
  int pertestoutput = 0;
  if (whichtest==1)
    {
      // large test
      // test with epsilon = 10^-4, BEFix, q = 1, residual
      test(0.0001, pertestoutput, 0, 1.0, 0, 0);
      // test with epsilon = 10^-4, BEDFix, q = 1, residual
      test(0.0001, pertestoutput, 1, 1.0, 0, 0);
      // test with epsilon = 10^-4, BEDFix, q = 0.9, absolute
      test(0.0001, pertestoutput, 1, 0.9, 1, 0);
      // test with epsilon = 10^-4, SI, q = 0.9
      test(0.0001, pertestoutput, 1, 0.9, 1, 1);
      // test with epsilon = 10^-4, BEDFix, q = 0.99, absolute
      test(0.0001, pertestoutput, 1, 0.99, 1, 0);
      // test with epsilon = 10^-4, SI, q = 0.99
      test(0.0001, pertestoutput, 1, 0.99, 1, 1);
      // test with epsilon = 10^-4, BEDFix, q = 0.999, absolute
      test(0.0001, pertestoutput, 1, 0.999, 1, 0);
      // test with epsilon = 10^-4, SI, q = 0.999
      test(0.0001, pertestoutput, 1, 0.999, 1, 1);
    }
  else
    {
      // small test
      // test with epsilon = 10^-4, BEDFix, q = 1, residual
      test(0.0001, pertestoutput, 1, 1.0, 0, 0);
    }
}

int main(int argc, char *argv[])
{
  int whichtest;
  if (argc<2 || argv[1] == NULL || (argv[1][0] != '1' && argv[1][0] != '2'))
    {
      while (1)
	{
	  char answer;
	  printf("Enter 1 to run the large test, 2 to run the small test:\n");
	  scanf("%c", &answer);
	  if (answer=='1')
	    {
	      whichtest = 1;
	      break;
	    }
	  if (answer=='2')
	    {
	      whichtest = 2;
	      break;
	    }
	}
    }
  else
    {
      whichtest = argv[1][0] == '1'? 1: 2;
    }
  printf("Running %s BEDFix test...\n", whichtest==1? "large": "small");
  testall(whichtest);
  return 0;
}
