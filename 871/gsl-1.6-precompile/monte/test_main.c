#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

for (I = problems ; I->f != 0; I++) 
{
  size_t i;
  MpIeee sum=  MpIeee( "0" );MpIeee  mean;MpIeee  sumd2=  MpIeee( "0" );MpIeee  sd;MpIeee  res;MpIeee  err; 
  
  gsl_rng * r;

  if (I->dim > 3)
    {
      continue ;
    }

  r = gsl_rng_alloc (gsl_rng_default);

  for (i = 0; i < TRIALS ; i++)
    {
      MONTE_STATE *s = MONTE_ALLOC (I->dim);
      
      I->f->dim = I->dim;
      
      MONTE_INTEGRATE (I->f, I->xl, I->xu, 
                       I->dim, I->calls / MONTE_SPEEDUP, r, s,
                       &res, &err);
      
      gsl_test_abs (res, I->expected_result, 
                    5 * GSL_MAX(err, 1024*GSL_DBL_EPSILON), 
                    NAME ", %s, result[%d]", I->description, i);

      MONTE_ERROR_TEST (err, I->expected_error);

      result[i] = res;
      error[i] = err;
      
      MONTE_FREE (s);
    }

 for (i = 0; i < TRIALS; i++)
   {
     sum += result[i];
   }

 mean = sum / TRIALS ;

 for (i = 0; i < TRIALS; i++) 
   {
     sumd2 += pow(result[i] - mean, MpIeee( "2.0" ));
   }

 sd = sqrt(sumd2 / (TRIALS-MpIeee( "1.0" ))) ;
 
 if (sd < TRIALS * GSL_DBL_EPSILON * fabs (mean))
   {
     sd = MpIeee( "0" );
   }

 for (i = 0; i < TRIALS; i++)
   {
     if (sd == MpIeee( "0" ) && fabs(error[i]) < GSL_DBL_EPSILON * fabs(result[i]))
       {
         error[i] = 0.0;
       }

     gsl_test_factor (error[i], sd, 5.0,
                      NAME ", %s, abserr[%d] vs sd", I->description, i);
   }


  gsl_rng_free (r);
}

