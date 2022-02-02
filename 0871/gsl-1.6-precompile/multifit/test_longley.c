#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

size_t longley_n = 16;
size_t longley_p = 7;

MpIeee longley_x[] =  {
  MpIeee( "1" ),  MpIeee( "83.0" ),   MpIeee( "234289" ),   MpIeee( "2356" ),     MpIeee( "1590" ),    MpIeee( "107608" ),  MpIeee( "1947" ),
  MpIeee( "1" ),  MpIeee( "88.5" ),   MpIeee( "259426" ),   MpIeee( "2325" ),     MpIeee( "1456" ),    MpIeee( "108632" ),  MpIeee( "1948" ),
  MpIeee( "1" ),  MpIeee( "88.2" ),   MpIeee( "258054" ),   MpIeee( "3682" ),     MpIeee( "1616" ),    MpIeee( "109773" ),  MpIeee( "1949" ),
  MpIeee( "1" ),  MpIeee( "89.5" ),   MpIeee( "284599" ),   MpIeee( "3351" ),     MpIeee( "1650" ),    MpIeee( "110929" ),  MpIeee( "1950" ),
  MpIeee( "1" ),  MpIeee( "96.2" ),   MpIeee( "328975" ),   MpIeee( "2099" ),     MpIeee( "3099" ),    MpIeee( "112075" ),  MpIeee( "1951" ),
  MpIeee( "1" ),  MpIeee( "98.1" ),   MpIeee( "346999" ),   MpIeee( "1932" ),     MpIeee( "3594" ),    MpIeee( "113270" ),  MpIeee( "1952" ),
  MpIeee( "1" ),  MpIeee( "99.0" ),   MpIeee( "365385" ),   MpIeee( "1870" ),     MpIeee( "3547" ),    MpIeee( "115094" ),  MpIeee( "1953" ),
  MpIeee( "1" ), MpIeee( "100.0" ),   MpIeee( "363112" ),   MpIeee( "3578" ),     MpIeee( "3350" ),    MpIeee( "116219" ),  MpIeee( "1954" ),
  MpIeee( "1" ), MpIeee( "101.2" ),   MpIeee( "397469" ),   MpIeee( "2904" ),     MpIeee( "3048" ),    MpIeee( "117388" ),  MpIeee( "1955" ),
  MpIeee( "1" ), MpIeee( "104.6" ),   MpIeee( "419180" ),   MpIeee( "2822" ),     MpIeee( "2857" ),    MpIeee( "118734" ),  MpIeee( "1956" ),
  MpIeee( "1" ), MpIeee( "108.4" ),   MpIeee( "442769" ),   MpIeee( "2936" ),     MpIeee( "2798" ),    MpIeee( "120445" ),  MpIeee( "1957" ),
  MpIeee( "1" ), MpIeee( "110.8" ),   MpIeee( "444546" ),   MpIeee( "4681" ),     MpIeee( "2637" ),    MpIeee( "121950" ),  MpIeee( "1958" ),
  MpIeee( "1" ), MpIeee( "112.6" ),   MpIeee( "482704" ),   MpIeee( "3813" ),     MpIeee( "2552" ),    MpIeee( "123366" ),  MpIeee( "1959" ),
  MpIeee( "1" ), MpIeee( "114.2" ),   MpIeee( "502601" ),   MpIeee( "3931" ),     MpIeee( "2514" ),    MpIeee( "125368" ),  MpIeee( "1960" ),
  MpIeee( "1" ), MpIeee( "115.7" ),   MpIeee( "518173" ),   MpIeee( "4806" ),     MpIeee( "2572" ),    MpIeee( "127852" ),  MpIeee( "1961" ),
  MpIeee( "1" ), MpIeee( "116.9" ),   MpIeee( "554894" ),   MpIeee( "4007" ),     MpIeee( "2827" ),    MpIeee( "130081" ),  MpIeee( "1962" ) } ;

MpIeee longley_y[] =  {MpIeee( "60323" ), MpIeee( "61122" ), MpIeee( "60171" ), MpIeee( "61187" ), MpIeee( "63221" ), MpIeee( "63639" ), MpIeee( "64989" ), MpIeee( "63761" ),
                       MpIeee( "66019" ), MpIeee( "67857" ), MpIeee( "68169" ), MpIeee( "66513" ), MpIeee( "68655" ), MpIeee( "69564" ), MpIeee( "69331" ), MpIeee( "70551" )};


void 
test_longley ()
{     
  size_t i, j;
  {
    gsl_multifit_linear_workspace * work = 
      gsl_multifit_linear_alloc (longley_n, longley_p);

    gsl_matrix_view X = gsl_matrix_view_array (longley_x, longley_n, longley_p);
    gsl_vector_view y = gsl_vector_view_array (longley_y, longley_n);
    gsl_vector * c = gsl_vector_alloc (longley_p);
    gsl_matrix * cov = gsl_matrix_alloc (longley_p, longley_p);
    gsl_vector_view diag;

    MpIeee chisq;

    MpIeee expected_c[7] =  {  -MpIeee( "3482258.63459582" ),
                              MpIeee( "15.0618722713733" ),
                              -MpIeee( "0.358191792925910E-01" ),
                              -MpIeee( "2.02022980381683" ),
                              -MpIeee( "1.03322686717359" ),
                              -MpIeee( "0.511041056535807E-01" ),
                              MpIeee( "1829.15146461355" ) };

    MpIeee expected_sd[7]  =  {  MpIeee( "890420.383607373" ),      
                                MpIeee( "84.9149257747669" ),      
                                MpIeee( "0.334910077722432E-01" ), 
                                MpIeee( "0.488399681651699" ),     
                                MpIeee( "0.214274163161675" ),     
                                MpIeee( "0.226073200069370" ),     
                                MpIeee( "455.478499142212" ) } ;  

    MpIeee expected_chisq=  MpIeee( "836424.055505915" );

    gsl_multifit_linear (&X.matrix, &y.vector, c, cov, &chisq, work);

    gsl_test_rel (gsl_vector_get(c,0), expected_c[0], 1e-10, "longley gsl_fit_multilinear c0") ;
    gsl_test_rel (gsl_vector_get(c,1), expected_c[1], 1e-10, "longley gsl_fit_multilinear c1") ;
    gsl_test_rel (gsl_vector_get(c,2), expected_c[2], 1e-10, "longley gsl_fit_multilinear c2") ;
    gsl_test_rel (gsl_vector_get(c,3), expected_c[3], 1e-10, "longley gsl_fit_multilinear c3") ;
    gsl_test_rel (gsl_vector_get(c,4), expected_c[4], 1e-10, "longley gsl_fit_multilinear c4") ;
    gsl_test_rel (gsl_vector_get(c,5), expected_c[5], 1e-10, "longley gsl_fit_multilinear c5") ;
    gsl_test_rel (gsl_vector_get(c,6), expected_c[6], 1e-10, "longley gsl_fit_multilinear c6") ;

    diag = gsl_matrix_diagonal (cov);

    gsl_test_rel (gsl_vector_get(&diag.vector,0), pow(expected_sd[0],2.0), 1e-10, "longley gsl_fit_multilinear cov00") ;
    gsl_test_rel (gsl_vector_get(&diag.vector,1), pow(expected_sd[1],2.0), 1e-10, "longley gsl_fit_multilinear cov11") ;
    gsl_test_rel (gsl_vector_get(&diag.vector,2), pow(expected_sd[2],2.0), 1e-10, "longley gsl_fit_multilinear cov22") ;
    gsl_test_rel (gsl_vector_get(&diag.vector,3), pow(expected_sd[3],2.0), 1e-10, "longley gsl_fit_multilinear cov33") ;
    gsl_test_rel (gsl_vector_get(&diag.vector,4), pow(expected_sd[4],2.0), 1e-10, "longley gsl_fit_multilinear cov44") ;
    gsl_test_rel (gsl_vector_get(&diag.vector,5), pow(expected_sd[5],2.0), 1e-10, "longley gsl_fit_multilinear cov55") ;
    gsl_test_rel (gsl_vector_get(&diag.vector,6), pow(expected_sd[6],2.0), 1e-10, "longley gsl_fit_multilinear cov66") ;

    gsl_test_rel (chisq, expected_chisq, 1e-10, "longley gsl_fit_multilinear chisq") ;

    gsl_vector_free(c);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free (work);
  }


  {
    gsl_multifit_linear_workspace * work = 
      gsl_multifit_linear_alloc (longley_n, longley_p);

    gsl_matrix_view X = gsl_matrix_view_array (longley_x, longley_n, longley_p);
    gsl_vector_view y = gsl_vector_view_array (longley_y, longley_n);
    gsl_vector * w = gsl_vector_alloc (longley_n);
    gsl_vector * c = gsl_vector_alloc (longley_p);
    gsl_matrix * cov = gsl_matrix_alloc (longley_p, longley_p);

    MpIeee chisq;

    MpIeee expected_c[7] =  {  -MpIeee( "3482258.63459582" ),
                              MpIeee( "15.0618722713733" ),
                              -MpIeee( "0.358191792925910E-01" ),
                              -MpIeee( "2.02022980381683" ),
                              -MpIeee( "1.03322686717359" ),
                              -MpIeee( "0.511041056535807E-01" ),
                              MpIeee( "1829.15146461355" ) };

    MpIeee expected_cov[7][7] =  { { 8531122.56783558,
-166.727799925578, 0.261873708176346, 3.91188317230983,
1.1285582054705, -0.889550869422687, -4362.58709870581},

{-166.727799925578, 0.0775861253030891, -1.98725210399982e-05,
-0.000247667096727256, -6.82911920718824e-05, 0.000136160797527761,
0.0775255245956248},

{0.261873708176346, -1.98725210399982e-05, 1.20690316701888e-08,
1.66429546772984e-07, 3.61843600487847e-08, -6.78805814483582e-08,
-0.00013158719037715},

{3.91188317230983, -0.000247667096727256, 1.66429546772984e-07,
2.56665052544717e-06, 6.96541409215597e-07, -9.00858307771567e-07,
-0.00197260370663974},

{1.1285582054705, -6.82911920718824e-05, 3.61843600487847e-08,
6.96541409215597e-07, 4.94032602583969e-07, -9.8469143760973e-08,
-0.000576921112208274},

{-0.889550869422687, 0.000136160797527761, -6.78805814483582e-08,
-9.00858307771567e-07, -9.8469143760973e-08, 5.49938542664952e-07,
0.000430074434198215},

{-4362.58709870581, 0.0775255245956248, -0.00013158719037715,
-0.00197260370663974, -0.000576921112208274, 0.000430074434198215,
2.23229587481535 }} ;

    MpIeee expected_chisq=  MpIeee( "836424.055505915" );

    gsl_vector_set_all (w, 1.0);

    gsl_multifit_wlinear (&X.matrix, w, &y.vector, c, cov, &chisq, work);

    gsl_test_rel (gsl_vector_get(c,0), expected_c[0], 1e-10, "longley gsl_fit_wmultilinear c0") ;
    gsl_test_rel (gsl_vector_get(c,1), expected_c[1], 1e-10, "longley gsl_fit_wmultilinear c1") ;
    gsl_test_rel (gsl_vector_get(c,2), expected_c[2], 1e-10, "longley gsl_fit_wmultilinear c2") ;
    gsl_test_rel (gsl_vector_get(c,3), expected_c[3], 1e-10, "longley gsl_fit_wmultilinear c3") ;
    gsl_test_rel (gsl_vector_get(c,4), expected_c[4], 1e-10, "longley gsl_fit_wmultilinear c4") ;
    gsl_test_rel (gsl_vector_get(c,5), expected_c[5], 1e-10, "longley gsl_fit_wmultilinear c5") ;
    gsl_test_rel (gsl_vector_get(c,6), expected_c[6], 1e-10, "longley gsl_fit_wmultilinear c6") ;

    for (i = 0; i < longley_p; i++) 
      {
        for (j = 0; j < longley_p; j++)
          {
            gsl_test_rel (gsl_matrix_get(cov,i,j), expected_cov[i][j], 1e-7, 
                          "longley gsl_fit_wmultilinear cov(%d,%d)", i, j) ;
          }
      }

    gsl_test_rel (chisq, expected_chisq, 1e-10, "longley gsl_fit_wmultilinear chisq") ;

    gsl_vector_free(w);
    gsl_vector_free(c);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free (work);
  }
}
