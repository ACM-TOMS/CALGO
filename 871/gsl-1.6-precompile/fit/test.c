#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* These tests are based on the NIST Statistical Reference Datasets
   See http://www.nist.gov/itl/div898/strd/index.html for more
   information. */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_fit.h>

#include <gsl/gsl_ieee_utils.h>

size_t norris_n = 36;

MpIeee norris_x[] =  { MpIeee( "0.2" ), MpIeee( "337.4" ), MpIeee( "118.2" ), MpIeee( "884.6" ), MpIeee( "10.1" ), MpIeee( "226.5" ), MpIeee( "666.3" ), MpIeee( "996.3" ),
                      MpIeee( "448.6" ), MpIeee( "777.0" ), MpIeee( "558.2" ), MpIeee( "0.4" ), MpIeee( "0.6" ), MpIeee( "775.5" ), MpIeee( "666.9" ), MpIeee( "338.0" ), 
                      MpIeee( "447.5" ), MpIeee( "11.6" ), MpIeee( "556.0" ), MpIeee( "228.1" ), MpIeee( "995.8" ), MpIeee( "887.6" ), MpIeee( "120.2" ), MpIeee( "0.3" ), 
                      MpIeee( "0.3" ), MpIeee( "556.8" ), MpIeee( "339.1" ), MpIeee( "887.2" ), MpIeee( "999.0" ), MpIeee( "779.0" ), MpIeee( "11.1" ), MpIeee( "118.3" ),
                      MpIeee( "229.2" ), MpIeee( "669.1" ), MpIeee( "448.9" ), MpIeee( "0.5" ) } ;

MpIeee norris_y[] =  { MpIeee( "0.1" ), MpIeee( "338.8" ), MpIeee( "118.1" ), MpIeee( "888.0" ), MpIeee( "9.2" ), MpIeee( "228.1" ), MpIeee( "668.5" ), MpIeee( "998.5" ),
                      MpIeee( "449.1" ), MpIeee( "778.9" ), MpIeee( "559.2" ), MpIeee( "0.3" ), MpIeee( "0.1" ), MpIeee( "778.1" ), MpIeee( "668.8" ), MpIeee( "339.3" ), 
                      MpIeee( "448.9" ), MpIeee( "10.8" ), MpIeee( "557.7" ), MpIeee( "228.3" ), MpIeee( "998.0" ), MpIeee( "888.8" ), MpIeee( "119.6" ), MpIeee( "0.3" ), 
                      MpIeee( "0.6" ), MpIeee( "557.6" ), MpIeee( "339.3" ), MpIeee( "888.0" ), MpIeee( "998.5" ), MpIeee( "778.9" ), MpIeee( "10.2" ), MpIeee( "117.6" ),
                      MpIeee( "228.9" ), MpIeee( "668.4" ), MpIeee( "449.2" ), MpIeee( "0.2" )};

size_t noint1_n = 11;
MpIeee noint1_x[] =  { MpIeee( "60" ), MpIeee( "61" ), MpIeee( "62" ), MpIeee( "63" ), MpIeee( "64" ), MpIeee( "65" ), MpIeee( "66" ), MpIeee( "67" ), MpIeee( "68" ), MpIeee( "69" ), MpIeee( "70" ) };
MpIeee noint1_y[] =  { MpIeee( "130" ), MpIeee( "131" ), MpIeee( "132" ), MpIeee( "133" ), MpIeee( "134" ), MpIeee( "135" ), MpIeee( "136" ), MpIeee( "137" ), MpIeee( "138" ), MpIeee( "139" ), MpIeee( "140" )};

size_t noint2_n = 3;
MpIeee noint2_x[] =  { MpIeee( "4" ), MpIeee( "5" ), MpIeee( "6" ) } ;
MpIeee noint2_y[] =  { MpIeee( "3" ), MpIeee( "4" ), MpIeee( "4" ) } ;

int
main (void)
{


  MpIeee x[1000];MpIeee  y[1000];MpIeee  w[1000];

  size_t xstride = 2, wstride = 3, ystride = 5;
  size_t i;

  for (i = 0; i < norris_n; i++) 
    {
      x[i*xstride] = norris_x[i];
      w[i*wstride] = MpIeee( "1.0" );
      y[i*ystride] = norris_y[i];
    }

  gsl_ieee_env_setup();

  {
    MpIeee c0;MpIeee  c1;MpIeee  cov00;MpIeee  cov01;MpIeee  cov11;MpIeee  sumsq;
       
    MpIeee expected_c0=  -MpIeee( "0.262323073774029" );
    MpIeee expected_c1=   MpIeee( "1.00211681802045" ); 
    MpIeee expected_cov00=  pow(MpIeee( "0.232818234301152" ), MpIeee( "2.0" ));
    MpIeee expected_cov01=  -MpIeee( "7.74327536339570e-05" );  /* computed from octave */
    MpIeee expected_cov11=  pow(MpIeee( "0.429796848199937E-03" ), MpIeee( "2.0" ));
    MpIeee expected_sumsq=  MpIeee( "26.6173985294224" );
    
    gsl_fit_linear (x, xstride, y, ystride, norris_n, 
                    &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
    
    /* gsl_fit_wlinear (x, xstride, w, wstride, y, ystride, norris_n, 
                     &c0, &c1, &cov00, &cov01, &cov11, &sumsq); */
  
    gsl_test_rel (c0, expected_c0, 1e-10, "norris gsl_fit_linear c0") ;
    gsl_test_rel (c1, expected_c1, 1e-10, "norris gsl_fit_linear c1") ;
    gsl_test_rel (cov00, expected_cov00, 1e-10, "norris gsl_fit_linear cov00") ;
    gsl_test_rel (cov01, expected_cov01, 1e-10, "norris gsl_fit_linear cov01") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "norris gsl_fit_linear cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "norris gsl_fit_linear sumsq") ;
  }

  {
    MpIeee c0;MpIeee  c1;MpIeee  cov00;MpIeee  cov01;MpIeee  cov11;MpIeee  sumsq;
       
    MpIeee expected_c0=  -MpIeee( "0.262323073774029" );
    MpIeee expected_c1=   MpIeee( "1.00211681802045" ); 
    MpIeee expected_cov00=  MpIeee( "6.92384428759429e-02" );  /* computed from octave */
    MpIeee expected_cov01=  -MpIeee( "9.89095016390515e-05" ); /* computed from octave */
    MpIeee expected_cov11=  MpIeee( "2.35960747164148e-07" );  /* computed from octave */
    MpIeee expected_sumsq=  MpIeee( "26.6173985294224" );
    
    gsl_fit_wlinear (x, xstride, w, wstride, y, ystride, norris_n, 
                     &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  
    gsl_test_rel (c0, expected_c0, 1e-10, "norris gsl_fit_wlinear c0") ;
    gsl_test_rel (c1, expected_c1, 1e-10, "norris gsl_fit_wlinear c1") ;
    gsl_test_rel (cov00, expected_cov00, 1e-10, "norris gsl_fit_wlinear cov00") ;
    gsl_test_rel (cov01, expected_cov01, 1e-10, "norris gsl_fit_wlinear cov01") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "norris gsl_fit_wlinear cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "norris gsl_fit_wlinear sumsq") ;
  }

  for (i = 0; i < noint1_n; i++) 
    {
      x[i*xstride] = noint1_x[i];
      w[i*wstride] = MpIeee( "1.0" );
      y[i*ystride] = noint1_y[i];
    }

  {
    MpIeee c1;MpIeee  cov11;MpIeee  sumsq;
       
    MpIeee expected_c1=  MpIeee( "2.07438016528926" ); 
    MpIeee expected_cov11=  pow(MpIeee( "0.165289256198347E-01" ), MpIeee( "2.0" ));  
    MpIeee expected_sumsq=  MpIeee( "127.272727272727" );
    
    gsl_fit_mul (x, xstride, y, ystride, noint1_n, &c1, &cov11, &sumsq);
  
    gsl_test_rel (c1, expected_c1, 1e-10, "noint1 gsl_fit_mul c1") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "noint1 gsl_fit_mul cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "noint1 gsl_fit_mul sumsq") ;
  }

  {
    MpIeee c1;MpIeee  cov11;MpIeee  sumsq;
       
    MpIeee expected_c1=  MpIeee( "2.07438016528926" ); 
    MpIeee expected_cov11=  MpIeee( "2.14661371686165e-05" ); /* computed from octave */
    MpIeee expected_sumsq=  MpIeee( "127.272727272727" );
    
    gsl_fit_wmul (x, xstride, w, wstride, y, ystride, noint1_n, &c1, &cov11, &sumsq);

    gsl_test_rel (c1, expected_c1, 1e-10, "noint1 gsl_fit_wmul c1") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "noint1 gsl_fit_wmul cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "noint1 gsl_fit_wmul sumsq") ;
  }


  for (i = 0; i < noint2_n; i++) 
    {
      x[i*xstride] = noint2_x[i];
      w[i*wstride] = MpIeee( "1.0" );
      y[i*ystride] = noint2_y[i];
    }

  {
    MpIeee c1;MpIeee  cov11;MpIeee  sumsq;
       
    MpIeee expected_c1=  MpIeee( "0.727272727272727" ); 
    MpIeee expected_cov11=  pow(MpIeee( "0.420827318078432E-01" ), MpIeee( "2.0" ));  
    MpIeee expected_sumsq=  MpIeee( "0.272727272727273" );
    
    gsl_fit_mul (x, xstride, y, ystride, noint2_n, &c1, &cov11, &sumsq);
  
    gsl_test_rel (c1, expected_c1, 1e-10, "noint2 gsl_fit_mul c1") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "noint2 gsl_fit_mul cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "noint2 gsl_fit_mul sumsq") ;
  }

  {
    MpIeee c1;MpIeee  cov11;MpIeee  sumsq;
       
    MpIeee expected_c1=  MpIeee( "0.727272727272727" ); 
    MpIeee expected_cov11=  MpIeee( "1.29870129870130e-02" ) ; /* computed from octave */
    MpIeee expected_sumsq=  MpIeee( "0.272727272727273" );
    
    gsl_fit_wmul (x, xstride, w, wstride, y, ystride, noint2_n, &c1, &cov11, &sumsq);

    gsl_test_rel (c1, expected_c1, 1e-10, "noint2 gsl_fit_wmul c1") ;
    gsl_test_rel (cov11, expected_cov11, 1e-10, "noint2 gsl_fit_wmul cov11") ;
    gsl_test_rel (sumsq, expected_sumsq, 1e-10, "noint2 gsl_fit_wmul sumsq") ;
  }

  /* now summarize the results */

  exit (gsl_test_summary ());
}
