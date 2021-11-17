#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

const size_t brown_N = 20;
const size_t brown_P = 4;

MpIeee brown_X[20][4] =  {
  {24.3485677, 4.71448798, -2.19486633, 2.69405755},
  {22.4116222, 3.93075538, -1.42344852, 2.5233557},
  {17.88886, 2.9290853, .125174936, -3.96823353},
  {17.3237176, 2.99606803, 2.03285653, 2.28992327},
  {17.0906508, 3.02485425, .296995153, .0876226126},
  {16.578006, 3.1036312, -.18617941, .103262914},
  {15.692993, 3.33088442, .0706406887, 1.05923955},
  {14.3232177, 3.85604218, -2.3762839, -3.09486813},
  {14.1279266, 3.97896121, .446109351, 1.40023753},
  {13.6081961, 4.16435075, -1.51250057, -1.52510626},
  {13.4295245, 4.22697223, -.196985195, .532009293},
  {13.0176117, 4.3579261, -.353131208, .301377627},
  {12.2713535, 4.62398535, -.00183585584, .894170703},
  {11.0316144, 5.13967727, -2.38978772, -2.89510064},
  {10.8807981, 5.24558004, .230495952, 1.27315117},
  {10.4029264, 5.41141257, -1.5116632, -1.47615921},
  {10.2574435, 5.46211045, -.299855732, .451893162},
  {9.87863876, 5.57914292, -.368885288, .358086545},
  {9.1894983, 5.82082741, -.230157969, .621476534},
  {8.00589008, 6.27788753, -1.46022815, -1.33468082}
};

MpIeee brown_F[20] =  {
  MpIeee( "2474.05541" ),
  MpIeee( "1924.69004" ),
  MpIeee( "1280.63194" ),
  MpIeee( "1244.81867" ),
  MpIeee( "1190.53739" ),
  MpIeee( "1159.34935" ),
  MpIeee( "1108.44426" ),
  MpIeee( "1090.11073" ),
  MpIeee( "1015.92942" ),
  MpIeee( "1002.43533" ),
  MpIeee( "971.221084" ),
  MpIeee( "949.589435" ),
  MpIeee( "911.359899" ),
  MpIeee( "906.522994" ),
  MpIeee( "840.525729" ),
  MpIeee( "833.950164" ),
  MpIeee( "807.557511" ),
  MpIeee( "791.00924" ),
  MpIeee( "761.09598" ),
  MpIeee( "726.787783" ),
};

MpIeee brown_cov[4][4] =  {
  { 1.8893186910e-01, -4.7099989571e-02,  5.2154168404e-01,  1.6608168209e-02},
  {-4.7099989571e-02,  1.1761534388e-02, -1.2987843074e-01, -4.1615942391e-03},
  { 5.2154168404e-01, -1.2987843074e-01,  1.4653936514e+00,  1.5738321686e-02},
  { 1.6608168209e-02, -4.1615942391e-03,  1.5738321686e-02,  4.2348042340e-02},
};

MpIeee brown_x0[4] =  { MpIeee( "25" ), MpIeee( "5" ), -MpIeee( "5" ), -MpIeee( "1" ) };

int
 brown_f(const gsl_vector * x, void *params, gsl_vector * f)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));
  MpIeee x2=  gsl_vector_get (x, MpIeee( "2" ));
  MpIeee x3=  gsl_vector_get (x, MpIeee( "3" ));
  size_t i;

  for (i = 0; i < 20; i++)
    {
      MpIeee ti=  MpIeee( "0.2" ) * (i + MpIeee( "1" ));
      MpIeee ui=  x0 + x1 * ti - exp (ti);
      MpIeee vi=  x2 + x3 * sin (ti) - cos (ti);

      gsl_vector_set (f, i, ui * ui + vi * vi);
    }

  return GSL_SUCCESS;
}

int
 brown_df(const gsl_vector * x, void *params, gsl_matrix * df)
{
  MpIeee x0=  gsl_vector_get (x, MpIeee( "0" ));
  MpIeee x1=  gsl_vector_get (x, MpIeee( "1" ));
  MpIeee x2=  gsl_vector_get (x, MpIeee( "2" ));
  MpIeee x3=  gsl_vector_get (x, MpIeee( "3" ));
  size_t i;

  for (i = 0; i < 20; i++)
    {
      MpIeee ti=  MpIeee( "0.2" ) * (i + MpIeee( "1" ));
      MpIeee ui=  x0 + x1 * ti - exp (ti);
      MpIeee vi=  x2 + x3 * sin (ti) - cos (ti);

      gsl_matrix_set (df, i, 0, 2 * ui);
      gsl_matrix_set (df, i, 1, 2 * ui * ti);
      gsl_matrix_set (df, i, 2, 2 * vi);
      gsl_matrix_set (df, i, 3, 2 * vi * sin (ti));

    }
  return GSL_SUCCESS;
}

int
 brown_fdf(const gsl_vector * x, void *params,
           gsl_vector * f, gsl_matrix * df)
{
  brown_f (x, params, f);
  brown_df (x, params, df);

  return GSL_SUCCESS;
}

