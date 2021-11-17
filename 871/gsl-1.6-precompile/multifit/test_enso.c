#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

const size_t enso_N = 168;
const size_t enso_P = 9;

MpIeee enso_x0[9] =  { MpIeee( "10.0" ), MpIeee( "3.0" ), MpIeee( "0.5" ), MpIeee( "44.0" ), -MpIeee( "1.5" ), MpIeee( "0.5" ), MpIeee( "26.0" ), MpIeee( "0.1" ), MpIeee( "1.5" ) };

MpIeee enso_x[9] =  {
  MpIeee( "1.0510749193E+01" ), 
  MpIeee( "3.0762128085E+00" ),
  MpIeee( "5.3280138227E-01" ),
  MpIeee( "4.4311088700E+01" ),
 -MpIeee( "1.6231428586E+00" ),
  MpIeee( "5.2554493756E-01" ),
  MpIeee( "2.6887614440E+01" ),
  MpIeee( "2.1232288488E-01" ),
  MpIeee( "1.4966870418E+00" )
};

MpIeee enso_sumsq=  MpIeee( "7.8853978668E+02" );

MpIeee enso_sigma[9] =  {
 MpIeee( "1.7488832467E-01" ),
 MpIeee( "2.4310052139E-01" ),
 MpIeee( "2.4354686618E-01" ),
 MpIeee( "9.4408025976E-01" ),
 MpIeee( "2.8078369611E-01" ),
 MpIeee( "4.8073701119E-01" ),
 MpIeee( "4.1612939130E-01" ),
 MpIeee( "5.1460022911E-01" ),
 MpIeee( "2.5434468893E-01" )
};

MpIeee enso_F[168] =  {
    MpIeee( "12.90000" ), 
    MpIeee( "11.30000" ), 
    MpIeee( "10.60000" ), 
    MpIeee( "11.20000" ), 
    MpIeee( "10.90000" ), 
    MpIeee( "7.500000" ), 
    MpIeee( "7.700000" ), 
    MpIeee( "11.70000" ), 
    MpIeee( "12.90000" ), 
    MpIeee( "14.30000" ), 
    MpIeee( "10.90000" ), 
    MpIeee( "13.70000" ), 
    MpIeee( "17.10000" ), 
    MpIeee( "14.00000" ), 
    MpIeee( "15.30000" ), 
    MpIeee( "8.500000" ), 
    MpIeee( "5.700000" ), 
    MpIeee( "5.500000" ), 
    MpIeee( "7.600000" ), 
    MpIeee( "8.600000" ), 
    MpIeee( "7.300000" ), 
    MpIeee( "7.600000" ), 
    MpIeee( "12.70000" ), 
    MpIeee( "11.00000" ), 
    MpIeee( "12.70000" ), 
    MpIeee( "12.90000" ), 
    MpIeee( "13.00000" ), 
    MpIeee( "10.90000" ), 
   MpIeee( "10.400000" ), 
   MpIeee( "10.200000" ), 
    MpIeee( "8.000000" ), 
    MpIeee( "10.90000" ), 
    MpIeee( "13.60000" ), 
   MpIeee( "10.500000" ), 
    MpIeee( "9.200000" ), 
    MpIeee( "12.40000" ), 
    MpIeee( "12.70000" ), 
    MpIeee( "13.30000" ), 
   MpIeee( "10.100000" ), 
    MpIeee( "7.800000" ), 
    MpIeee( "4.800000" ), 
    MpIeee( "3.000000" ), 
    MpIeee( "2.500000" ), 
    MpIeee( "6.300000" ), 
    MpIeee( "9.700000" ), 
    MpIeee( "11.60000" ), 
    MpIeee( "8.600000" ), 
    MpIeee( "12.40000" ), 
   MpIeee( "10.500000" ), 
    MpIeee( "13.30000" ), 
   MpIeee( "10.400000" ), 
    MpIeee( "8.100000" ), 
    MpIeee( "3.700000" ), 
    MpIeee( "10.70000" ), 
    MpIeee( "5.100000" ), 
   MpIeee( "10.400000" ), 
    MpIeee( "10.90000" ), 
    MpIeee( "11.70000" ), 
    MpIeee( "11.40000" ), 
    MpIeee( "13.70000" ), 
    MpIeee( "14.10000" ), 
    MpIeee( "14.00000" ), 
    MpIeee( "12.50000" ), 
    MpIeee( "6.300000" ), 
    MpIeee( "9.600000" ), 
    MpIeee( "11.70000" ), 
    MpIeee( "5.000000" ), 
    MpIeee( "10.80000" ), 
    MpIeee( "12.70000" ), 
    MpIeee( "10.80000" ), 
    MpIeee( "11.80000" ), 
    MpIeee( "12.60000" ), 
    MpIeee( "15.70000" ), 
    MpIeee( "12.60000" ), 
    MpIeee( "14.80000" ), 
    MpIeee( "7.800000" ), 
    MpIeee( "7.100000" ), 
    MpIeee( "11.20000" ), 
    MpIeee( "8.100000" ), 
    MpIeee( "6.400000" ), 
    MpIeee( "5.200000" ), 
    MpIeee( "12.00000" ), 
   MpIeee( "10.200000" ), 
    MpIeee( "12.70000" ), 
   MpIeee( "10.200000" ), 
    MpIeee( "14.70000" ), 
    MpIeee( "12.20000" ), 
    MpIeee( "7.100000" ), 
    MpIeee( "5.700000" ), 
    MpIeee( "6.700000" ), 
    MpIeee( "3.900000" ), 
    MpIeee( "8.500000" ), 
    MpIeee( "8.300000" ), 
    MpIeee( "10.80000" ), 
    MpIeee( "16.70000" ), 
    MpIeee( "12.60000" ), 
    MpIeee( "12.50000" ), 
    MpIeee( "12.50000" ), 
    MpIeee( "9.800000" ), 
    MpIeee( "7.200000" ), 
    MpIeee( "4.100000" ), 
    MpIeee( "10.60000" ), 
   MpIeee( "10.100000" ), 
   MpIeee( "10.100000" ), 
    MpIeee( "11.90000" ), 
    MpIeee( "13.60000" ), 
    MpIeee( "16.30000" ), 
    MpIeee( "17.60000" ), 
    MpIeee( "15.50000" ), 
    MpIeee( "16.00000" ), 
    MpIeee( "15.20000" ), 
    MpIeee( "11.20000" ), 
    MpIeee( "14.30000" ), 
    MpIeee( "14.50000" ), 
    MpIeee( "8.500000" ), 
    MpIeee( "12.00000" ), 
    MpIeee( "12.70000" ), 
    MpIeee( "11.30000" ), 
    MpIeee( "14.50000" ), 
    MpIeee( "15.10000" ), 
   MpIeee( "10.400000" ), 
    MpIeee( "11.50000" ), 
    MpIeee( "13.40000" ), 
    MpIeee( "7.500000" ), 
   MpIeee( "0.6000000" ), 
   MpIeee( "0.3000000" ), 
    MpIeee( "5.500000" ), 
    MpIeee( "5.000000" ), 
    MpIeee( "4.600000" ), 
    MpIeee( "8.200000" ), 
    MpIeee( "9.900000" ), 
    MpIeee( "9.200000" ), 
    MpIeee( "12.50000" ), 
    MpIeee( "10.90000" ), 
    MpIeee( "9.900000" ), 
    MpIeee( "8.900000" ), 
    MpIeee( "7.600000" ), 
    MpIeee( "9.500000" ), 
    MpIeee( "8.400000" ), 
    MpIeee( "10.70000" ), 
    MpIeee( "13.60000" ), 
    MpIeee( "13.70000" ), 
    MpIeee( "13.70000" ), 
    MpIeee( "16.50000" ), 
    MpIeee( "16.80000" ), 
    MpIeee( "17.10000" ), 
    MpIeee( "15.40000" ), 
    MpIeee( "9.500000" ), 
    MpIeee( "6.100000" ), 
   MpIeee( "10.100000" ), 
    MpIeee( "9.300000" ), 
    MpIeee( "5.300000" ), 
    MpIeee( "11.20000" ), 
    MpIeee( "16.60000" ), 
    MpIeee( "15.60000" ), 
    MpIeee( "12.00000" ), 
    MpIeee( "11.50000" ), 
    MpIeee( "8.600000" ), 
    MpIeee( "13.80000" ), 
    MpIeee( "8.700000" ), 
    MpIeee( "8.600000" ), 
    MpIeee( "8.600000" ), 
    MpIeee( "8.700000" ), 
    MpIeee( "12.80000" ), 
    MpIeee( "13.20000" ), 
    MpIeee( "14.00000" ), 
    MpIeee( "13.40000" ), 
    MpIeee( "14.80000" )
};


int
 enso_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee b[9];
  size_t i;

  for (i = 0; i < 9; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < 168; i++)
    {
      MpIeee t=  (i + MpIeee( "1.0" ));
      MpIeee y;
      y = b[0];
      y += b[1] * cos(MpIeee( "2" )*M_PI*t/MpIeee( "12" ));
      y += b[2] * sin(MpIeee( "2" )*M_PI*t/MpIeee( "12" ));
      y += b[4] * cos(MpIeee( "2" )*M_PI*t/b[3]);
      y += b[5] * sin(MpIeee( "2" )*M_PI*t/b[3]);
      y += b[7] * cos(MpIeee( "2" )*M_PI*t/b[6]);
      y += b[8] * sin(MpIeee( "2" )*M_PI*t/b[6]);

      gsl_vector_set (f, i, enso_F[i] - y);
    }

  return GSL_SUCCESS;
}

int
 enso_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee b[9];
  size_t i;

  for (i = 0; i < 9; i++)
    {
      b[i] = gsl_vector_get(x, i);
    }

  for (i = 0; i < 168; i++)
    {
      MpIeee t=  (i + MpIeee( "1.0" ));

      gsl_matrix_set (df, i, 0, -1.0);
      gsl_matrix_set (df, i, 1, -cos(2*M_PI*t/12));
      gsl_matrix_set (df, i, 2, -sin(2*M_PI*t/12));
      gsl_matrix_set (df, i, 3, 
                      -b[4]*(2*M_PI*t/(b[3]*b[3]))*sin(2*M_PI*t/b[3])
                      +b[5]*(2*M_PI*t/(b[3]*b[3]))*cos(2*M_PI*t/b[3]));
      gsl_matrix_set (df, i, 4, -cos(2*M_PI*t/b[3]));
      gsl_matrix_set (df, i, 5, -sin(2*M_PI*t/b[3]));
      gsl_matrix_set (df, i, 6, 
                     -b[7] * (2*M_PI*t/(b[6]*b[6])) * sin(2*M_PI*t/b[6])
                     +b[8] * (2*M_PI*t/(b[6]*b[6])) * cos(2*M_PI*t/b[6]));
      gsl_matrix_set (df, i, 7, -cos(2*M_PI*t/b[6]));
      gsl_matrix_set (df, i, 8, -sin(2*M_PI*t/b[6]));
    }

  return GSL_SUCCESS;
}

int
 enso_fdf(const gsl_vector * x, void *params,
           gsl_vector * f, gsl_matrix * df)
{
  enso_f (x, params, f);
  enso_df (x, params, df);

  return GSL_SUCCESS;
}





