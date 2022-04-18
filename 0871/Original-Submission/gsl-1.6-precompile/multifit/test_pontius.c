#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

size_t pontius_n = 40;
size_t pontius_p = 3;

MpIeee pontius_x[] =  { MpIeee( "150000" ), MpIeee( "300000" ), MpIeee( "450000" ), MpIeee( "600000" ), MpIeee( "750000" ), MpIeee( "900000" ),
MpIeee( "1050000" ), MpIeee( "1200000" ), MpIeee( "1350000" ), MpIeee( "1500000" ), MpIeee( "1650000" ), MpIeee( "1800000" ), MpIeee( "1950000" ), MpIeee( "2100000" ),
MpIeee( "2250000" ), MpIeee( "2400000" ), MpIeee( "2550000" ), MpIeee( "2700000" ), MpIeee( "2850000" ), MpIeee( "3000000" ), MpIeee( "150000" ), MpIeee( "300000" ),
MpIeee( "450000" ), MpIeee( "600000" ), MpIeee( "750000" ), MpIeee( "900000" ), MpIeee( "1050000" ), MpIeee( "1200000" ), MpIeee( "1350000" ), MpIeee( "1500000" ),
MpIeee( "1650000" ), MpIeee( "1800000" ), MpIeee( "1950000" ), MpIeee( "2100000" ), MpIeee( "2250000" ), MpIeee( "2400000" ), MpIeee( "2550000" ), MpIeee( "2700000" ),
MpIeee( "2850000" ), MpIeee( "3000000" ) };

MpIeee pontius_y[] =  { MpIeee( ".11019" ), MpIeee( ".21956" ), MpIeee( ".32949" ), MpIeee( ".43899" ), MpIeee( ".54803" ), MpIeee( ".65694" ),
MpIeee( ".76562" ), MpIeee( ".87487" ), MpIeee( ".98292" ), MpIeee( "1.09146" ), MpIeee( "1.20001" ), MpIeee( "1.30822" ), MpIeee( "1.41599" ), MpIeee( "1.52399" ),
MpIeee( "1.63194" ), MpIeee( "1.73947" ), MpIeee( "1.84646" ), MpIeee( "1.95392" ), MpIeee( "2.06128" ), MpIeee( "2.16844" ), MpIeee( ".11052" ), MpIeee( ".22018" ),
MpIeee( ".32939" ), MpIeee( ".43886" ), MpIeee( ".54798" ), MpIeee( ".65739" ), MpIeee( ".76596" ), MpIeee( ".87474" ), MpIeee( ".98300" ), MpIeee( "1.09150" ),
MpIeee( "1.20004" ), MpIeee( "1.30818" ), MpIeee( "1.41613" ), MpIeee( "1.52408" ), MpIeee( "1.63159" ), MpIeee( "1.73965" ), MpIeee( "1.84696" ),
MpIeee( "1.95445" ), MpIeee( "2.06177" ), MpIeee( "2.16829" ) };

void
test_pontius ()
{
  size_t i, j;
  {
    gsl_multifit_linear_workspace * work = 
      gsl_multifit_linear_alloc (pontius_n, pontius_p);

    gsl_matrix * X = gsl_matrix_alloc (pontius_n, pontius_p);
    gsl_vector_view y = gsl_vector_view_array (pontius_y, pontius_n);
    gsl_vector * c = gsl_vector_alloc (pontius_p);
    gsl_matrix * cov = gsl_matrix_alloc (pontius_p, pontius_p);
    gsl_vector_view diag;

    MpIeee chisq;

    MpIeee expected_c[3] =  { MpIeee( "0.673565789473684E-03" ),
                             MpIeee( "0.732059160401003E-06" ),
                            -MpIeee( "0.316081871345029E-14" )};

    MpIeee expected_sd[3] =  { MpIeee( "0.107938612033077E-03" ),
                              MpIeee( "0.157817399981659E-09" ),
                              MpIeee( "0.486652849992036E-16" ) };

    MpIeee expected_chisq=  MpIeee( "0.155761768796992E-05" );

    for (i = 0 ; i < pontius_n; i++) 
      {
        for (j = 0; j < pontius_p; j++) 
          {
            gsl_matrix_set(X, i, j, pow(pontius_x[i], j));
          }
      }

    gsl_multifit_linear (X, &y.vector, c, cov, &chisq, work);

    gsl_test_rel (gsl_vector_get(c,0), expected_c[0], 1e-10, "pontius gsl_fit_multilinear c0") ;
    gsl_test_rel (gsl_vector_get(c,1), expected_c[1], 1e-10, "pontius gsl_fit_multilinear c1") ;
    gsl_test_rel (gsl_vector_get(c,2), expected_c[2], 1e-10, "pontius gsl_fit_multilinear c2") ;

    diag = gsl_matrix_diagonal (cov);

    gsl_test_rel (gsl_vector_get(&diag.vector,0), pow(expected_sd[0],2.0), 1e-10, "pontius gsl_fit_multilinear cov00") ;
    gsl_test_rel (gsl_vector_get(&diag.vector,1), pow(expected_sd[1],2.0), 1e-10, "pontius gsl_fit_multilinear cov11") ;
    gsl_test_rel (gsl_vector_get(&diag.vector,2), pow(expected_sd[2],2.0), 1e-10, "pontius gsl_fit_multilinear cov22") ;

    gsl_test_rel (chisq, expected_chisq, 1e-10, "pontius gsl_fit_multilinear chisq") ;

    gsl_vector_free(c);
    gsl_matrix_free(cov);
    gsl_matrix_free(X);
    gsl_multifit_linear_free (work);
  }


  {
    gsl_multifit_linear_workspace * work = 
      gsl_multifit_linear_alloc (pontius_n, pontius_p);

    gsl_matrix * X = gsl_matrix_alloc (pontius_n, pontius_p);
    gsl_vector_view y = gsl_vector_view_array (pontius_y, pontius_n);
    gsl_vector * w = gsl_vector_alloc (pontius_n);
    gsl_vector * c = gsl_vector_alloc (pontius_p);
    gsl_matrix * cov = gsl_matrix_alloc (pontius_p, pontius_p);

    MpIeee chisq;

    MpIeee expected_c[3] =  {  MpIeee( "0.673565789473684E-03" ),
                               MpIeee( "0.732059160401003E-06" ),
                               -MpIeee( "0.316081871345029E-14" )};

    MpIeee expected_chisq=  MpIeee( "0.155761768796992E-05" );

    MpIeee expected_cov[3][3] = { 
      {2.76754385964916e-01 , -3.59649122807024e-07,   9.74658869395731e-14},
      {-3.59649122807024e-07,   5.91630591630603e-13,  -1.77210703526497e-19},
      {9.74658869395731e-14,  -1.77210703526497e-19,   5.62573661988878e-26} };


    for (i = 0 ; i < pontius_n; i++) 
      {
        for (j = 0; j < pontius_p; j++) 
          {
            gsl_matrix_set(X, i, j, pow(pontius_x[i], j));
          }
      }

    gsl_vector_set_all (w, 1.0);

    gsl_multifit_wlinear (X, w, &y.vector, c, cov, &chisq, work);

    gsl_test_rel (gsl_vector_get(c,0), expected_c[0], 1e-10, "pontius gsl_fit_multilinear c0") ;
    gsl_test_rel (gsl_vector_get(c,1), expected_c[1], 1e-10, "pontius gsl_fit_multilinear c1") ;
    gsl_test_rel (gsl_vector_get(c,2), expected_c[2], 1e-10, "pontius gsl_fit_multilinear c2") ;


    for (i = 0; i < pontius_p; i++) 
      {
        for (j = 0; j < pontius_p; j++)
          {
            gsl_test_rel (gsl_matrix_get(cov,i,j), expected_cov[i][j], 1e-10, 
                          "pontius gsl_fit_wmultilinear cov(%d,%d)", i, j) ;
          }
      }

    gsl_test_rel (chisq, expected_chisq, 1e-10, "pontius gsl_fit_multilinear chisq") ;

    gsl_vector_free(w);
    gsl_vector_free(c);
    gsl_matrix_free(cov);
    gsl_matrix_free(X);
    gsl_multifit_linear_free (work);
  }
}
