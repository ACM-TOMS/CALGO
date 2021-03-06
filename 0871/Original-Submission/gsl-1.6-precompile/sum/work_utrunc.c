#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sum.h>

gsl_sum_levin_utrunc_workspace * 
gsl_sum_levin_utrunc_alloc (size_t n)
{
  gsl_sum_levin_utrunc_workspace * w;

  if (n == 0)
    {
      GSL_ERROR_VAL ("length n must be positive integer", GSL_EDOM, 0);
    }

  w = (gsl_sum_levin_utrunc_workspace *) malloc(sizeof(gsl_sum_levin_utrunc_workspace));

  if (w == NULL)
    {
      GSL_ERROR_VAL ("failed to allocate struct", GSL_ENOMEM, 0);
    }

  w->q_num = (MpIeee*) malloc (n * sizeof (MpIeee));

  if (w->q_num == NULL)
    {
      free(w) ; /* error in constructor, prevent memory leak */

      GSL_ERROR_VAL ("failed to allocate space for q_num", GSL_ENOMEM, 0);
    }

  w->q_den = (MpIeee*) malloc (n * sizeof (MpIeee));

  if (w->q_den == NULL)
    {
      free (w->q_num);
      free (w) ; /* error in constructor, prevent memory leak */

      GSL_ERROR_VAL ("failed to allocate space for q_den", GSL_ENOMEM, 0);
    }

  w->dsum = (MpIeee*) malloc (n * sizeof (MpIeee));

  if (w->dsum == NULL)
    {
      free (w->q_den);
      free (w->q_num);
      free (w) ; /* error in constructor, prevent memory leak */

      GSL_ERROR_VAL ("failed to allocate space for dsum", GSL_ENOMEM, 0);
    }

  w->size = n;
  w->terms_used = 0;
  w->sum_plain = 0;

  return w;
}

void
gsl_sum_levin_utrunc_free (gsl_sum_levin_utrunc_workspace * w)
{
  free (w->dsum);
  free (w->q_den);
  free (w->q_num);
  free (w);
}
