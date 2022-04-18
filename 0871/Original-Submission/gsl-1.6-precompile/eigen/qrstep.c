#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* remove off-diagonal elements which are neglegible compared with the
   neighboring diagonal elements */

inline static void
chop_small_elements (const size_t N, const MpIeee d[], MpIeee sd[])
{
  MpIeee d_i=  d[0];

  size_t i;

  for (i = 0; i < N - 1; i++)
    {
      MpIeee sd_i=  sd[i];
      MpIeee d_ip1=  d[i + 1];

      if (fabs (sd_i) < GSL_DBL_EPSILON * (fabs (d_i) + fabs (d_ip1)))
        {
          sd[i] = MpIeee( "0.0" );
        }
      d_i = d_ip1;
    }
}

/* Generate a Givens rotation (cos,sin) which takes v=(x,y) to (|v|,0) 

   From Golub and Van Loan, "Matrix Computations", Section 5.1.8 */

inline static void
create_givens (const MpIeee a, const MpIeee b, MpIeee *c, MpIeee *s)
{
  if (b == 0)
    {
      *c = MpIeee( "1" );
      *s = MpIeee( "0" );
    }
  else if (fabs (b) > fabs (a))
    {
      MpIeee t=  -a / b;
      MpIeee s1=  MpIeee( "1.0" ) / sqrt (MpIeee( "1" ) + t * t);
      *s = s1;
      *c = s1 * t;
    }
  else
    {
      MpIeee t=  -b / a;
      MpIeee c1=  MpIeee( "1.0" ) / sqrt (MpIeee( "1" ) + t * t);
      *c = c1;
      *s = c1 * t;
    }
}

inline static MpIeee trailing_eigenvalue(const size_t n, const MpIeee d[], const MpIeee sd[])
{
  MpIeee ta=  d[n - 2];
  MpIeee tb=  d[n - 1];
  MpIeee tab=  sd[n - 2];

  MpIeee dt=  (ta - tb) / MpIeee( "2.0" );

  MpIeee mu;

  if (dt >= MpIeee( "0" ))
    {
      mu = tb - (tab * tab) / (dt + hypot (dt, tab));
    }
  else
    {
      mu = tb + (tab * tab) / ((-dt) + hypot (dt, tab));
    }

  return mu;
}

static void
qrstep (const size_t n, MpIeee d[], MpIeee sd[], MpIeee gc[], MpIeee gs[])
{
  MpIeee x;MpIeee  z;
  MpIeee ak;MpIeee  bk;MpIeee  zk;MpIeee  ap;MpIeee  bp;MpIeee  aq;MpIeee  bq;
  size_t k;

  MpIeee mu=  trailing_eigenvalue (n, d, sd);

  x = d[0] - mu;
  z = sd[0];

  ak = MpIeee( "0" );
  bk = MpIeee( "0" );
  zk = MpIeee( "0" );

  ap = d[0];
  bp = sd[0];

  aq = d[1];

  if (n == 2)
    {
      MpIeee c;MpIeee  s;
      create_givens (x, z, &c, &s);

      if (gc != NULL)
        gc[0] = c; 
      if (gs != NULL)
        gs[0] = s;

      {
        MpIeee ap1=  c * (c * ap - s * bp) + s * (s * aq - c * bp);
        MpIeee bp1=  c * (s * ap + c * bp) - s * (s * bp + c * aq);

        MpIeee aq1=  s * (s * ap + c * bp) + c * (s * bp + c * aq);

        ak = ap1;
        bk = bp1;

        ap = aq1;
      }

      d[0] = ak;
      sd[0] = bk;
      d[1] = ap;

      return;
    }

  bq = sd[1];

  for (k = 0; k < n - 1; k++)
    {
      MpIeee c;MpIeee  s;
      create_givens (x, z, &c, &s);

      /* store Givens rotation */
      if (gc != NULL)
        gc[k] = c; 
      if (gs != NULL)
        gs[k] = s;

      /* compute G' T G */

      {
        MpIeee bk1=  c * bk - s * zk;

        MpIeee ap1=  c * (c * ap - s * bp) + s * (s * aq - c * bp);
        MpIeee bp1=  c * (s * ap + c * bp) - s * (s * bp + c * aq);
        MpIeee zp1=  -s * bq;

        MpIeee aq1=  s * (s * ap + c * bp) + c * (s * bp + c * aq);
        MpIeee bq1=  c * bq;

        ak = ap1;
        bk = bp1;
        zk = zp1;

        ap = aq1;
        bp = bq1;

        if (k < n - 2)
          aq = d[k + 2];
        if (k < n - 3)
          bq = sd[k + 2];

        d[k] = ak;

        if (k > 0)
          sd[k - 1] = bk1;

        if (k < n - 2)
          sd[k + 1] = bp;

        x = bk;
        z = zk;
      }
    }

  /* k = n - 1 */
  d[k] = ap;
  sd[k - 1] = bk;
}
