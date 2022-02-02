// Here is the C implementation of the BEDFix algorithm.
// The function "fp" executes the algorithm and returns a fixed point approximation
// given two bivariate functions.  It and the functions and definitions preceding
// it may be copied into other software modules; the remainder of the code is for
// testing purposes.
// The function "testall" performs all of the tests described in the BEDFix paper.
// The commands used to compile this module and run the tests are (using GNU C):

// gcc -lm fp.c -o fp
// ./fp 1
//   for the large test
// ./fp 2
//   for the small (quick) test
// ./fp
//   to be prompted for which test to run

// where fp.c is the name of the source file.

// The small test is the same as the large test, but without the most
// time-consuming SI comparison tests.

// The tests were executed successfully on the following platforms:

// Pentium i686 running Red Hat Linux 7.2




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
// Return code: 0 normal, otherwise error
// To guarantee that the solution satisfies the absolute
// criterion |*fx-x| <= epsilon and |*fy-y| <= epsilon,
// pass in epsilon(1-q) for the epsilon parameter,
// where q<1 is the Lipschitz constant of f.

int fp(double epsilon, bivariate_func f1, bivariate_func f2,
       double *fx, double *fy, ulong *total, int *absolute)
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
      // deep cut
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
  while (ext1 + ext2 > epsilon / 2.0);
  *fx = xm;
  *fy = ym;
  *absolute = 1;
  return 0;
}


