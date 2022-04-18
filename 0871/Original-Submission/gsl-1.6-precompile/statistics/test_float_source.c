#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/test_float_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Jim Davies, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

void FUNCTION (test, func) (const size_t stridea, const size_t strideb);

void
FUNCTION (test, func) (const size_t stridea, const size_t strideb)
{
  /* sample sets of doubles */
  size_t i;
  const size_t na = 14, nb = 14;

  const MpIeee rawa[] = 
  {.0421, .0941, .1064, .0242, .1331,
   .0773, .0243, .0815, .1186, .0356,
   .0728, .0999, .0614, .0479};

  const MpIeee rawb[] = 
  {.1081, .0986, .1566, .1961, .1125,
   .1942, .1079, .1021, .1583, .1673,
   .1675, .1856, .1688, .1512};

  const MpIeee raww[] =  
  {.0000, .0000, .0000, 3.000, .0000,
   1.000, 1.000, 1.000, 0.000, .5000,
   7.000, 5.000, 4.000, 0.123};

  MpIeee * sorted;

  MpIeee * groupa=  (MpIeee *) malloc (stridea * na * sizeof(MpIeee));
  MpIeee * groupb=  (MpIeee *) malloc (strideb * nb * sizeof(MpIeee));
  MpIeee * w=  (MpIeee *) malloc (strideb * na * sizeof(MpIeee));

#ifdef BASE_FLOAT
  MpIeee rel=  MpIeee( "1" )e-MpIeee( "6" );
#else
  MpIeee rel=  MpIeee( "1" )e-MpIeee( "10" );
#endif

  for (i = 0 ; i < na ; i++)
    groupa[i * stridea] = (MpIeee) rawa[i] ;

  for (i = 0 ; i < na ; i++)
    w[i * strideb] = (MpIeee) raww[i] ;

  for (i = 0 ; i < nb ; i++)
    groupb[i * strideb] = (MpIeee) rawb[i] ;


  {
    MpIeee mean=  FUNCTION(gsl_stats,mean) (groupa, stridea, na);
    MpIeee expected=  MpIeee( "0.0728" );
    gsl_test_rel (mean, expected, rel, NAME(gsl_stats) "_mean");
  }

  {
    MpIeee mean=  FUNCTION(gsl_stats,mean) (groupa, stridea, na);
    MpIeee var=  FUNCTION(gsl_stats,variance_with_fixed_mean) (groupa, stridea, na, mean);
    MpIeee expected=  MpIeee( "0.00113837428571429" );
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_variance_with_fixed_mean");
  }


  {
    MpIeee mean=  FUNCTION(gsl_stats,mean) (groupa, stridea, na);
    MpIeee var=  FUNCTION(gsl_stats,sd_with_fixed_mean) (groupa, stridea, na, mean);
    MpIeee expected=  MpIeee( "0.0337398026922845" );
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_sd_with_fixed_mean");
  }


  {
    MpIeee var=  FUNCTION(gsl_stats,variance) (groupb, strideb, nb);
    MpIeee expected=  MpIeee( "0.00124956615384615" );
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_variance");
  }

  {
    MpIeee sd=  FUNCTION(gsl_stats,sd) (groupa, stridea, na);
    MpIeee expected=  MpIeee( "0.0350134479659107" );
    gsl_test_rel (sd, expected, rel, NAME(gsl_stats) "_sd");
  }

  {
    MpIeee absdev=  FUNCTION(gsl_stats,absdev) (groupa, stridea, na);
    MpIeee expected=  MpIeee( "0.0287571428571429" );
    gsl_test_rel (absdev, expected, rel, NAME(gsl_stats) "_absdev");
  }

  {
    MpIeee skew=  FUNCTION(gsl_stats,skew) (groupa, stridea, na);
    MpIeee expected=  MpIeee( "0.0954642051479004" );
    gsl_test_rel (skew, expected, rel, NAME(gsl_stats) "_skew");
  }

  {
    MpIeee kurt=  FUNCTION(gsl_stats,kurtosis) (groupa, stridea, na);
    MpIeee expected=  -MpIeee( "1.38583851548909" ) ;
    gsl_test_rel (kurt, expected, rel, NAME(gsl_stats) "_kurtosis");
  }

  {
    MpIeee wmean=  FUNCTION(gsl_stats,wmean) (w, strideb, groupa, stridea, na);
    MpIeee expected=  MpIeee( "0.0678111523670601" );
    gsl_test_rel (wmean, expected, rel, NAME(gsl_stats) "_wmean");
  }

  {
    MpIeee wmean=  FUNCTION(gsl_stats,wmean) (w, strideb, groupa, stridea, na);
    MpIeee wvar=  FUNCTION(gsl_stats,wvariance_with_fixed_mean) (w, strideb, groupa, stridea, na, wmean);
    MpIeee expected=  MpIeee( "0.000615793060878654" );
    gsl_test_rel (wvar, expected, rel, NAME(gsl_stats) "_wvariance_with_fixed_mean");
  }

  {
    MpIeee est_wvar=  FUNCTION(gsl_stats,wvariance) (w, strideb, groupa, stridea, na);
    MpIeee expected=  MpIeee( "0.000769562962860317" );
    gsl_test_rel (est_wvar, expected, rel, NAME(gsl_stats) "_wvariance");
  }

  {
    MpIeee wsd=  FUNCTION(gsl_stats,wsd) (w, strideb, groupa, stridea, na);
    MpIeee expected=  MpIeee( "0.0277409978706664" );
    gsl_test_rel (wsd, expected, rel, NAME(gsl_stats) "_wsd");
  }

  {
    MpIeee wabsdev=  FUNCTION(gsl_stats,wabsdev) (w, strideb, groupa, stridea, na);
    MpIeee expected=  MpIeee( "0.0193205027504008" );
    gsl_test_rel (wabsdev, expected, rel, NAME(gsl_stats) "_wabsdev");
  }

  {
    MpIeee wskew=  FUNCTION(gsl_stats,wskew) (w, strideb, groupa, stridea, na);
    MpIeee expected=  -MpIeee( "0.373631000307076" );
    gsl_test_rel (wskew, expected, rel, NAME(gsl_stats) "_wskew");
  }

  {
    MpIeee wkurt=  FUNCTION(gsl_stats,wkurtosis) (w, strideb, groupa, stridea, na);
    MpIeee expected=  -MpIeee( "1.48114233353963" );
    gsl_test_rel (wkurt, expected, rel, NAME(gsl_stats) "_wkurtosis");
  }

  {
    MpIeee c=  FUNCTION(gsl_stats,covariance) (groupa, stridea, groupb, strideb, nb);
    MpIeee expected=  -MpIeee( "0.000139021538461539" );
    gsl_test_rel (c, expected, rel, NAME(gsl_stats) "_covariance");
  }


  {
    MpIeee pv=  FUNCTION(gsl_stats,pvariance) (groupa, stridea, na, groupb, strideb, nb);
    MpIeee expected=  MpIeee( "0.00123775384615385" );
    gsl_test_rel (pv, expected, rel, NAME(gsl_stats) "_pvariance");
  }

  {
    MpIeee t=  FUNCTION(gsl_stats,ttest) (groupa, stridea, na, groupb, strideb, nb);
    MpIeee expected=  -MpIeee( "5.67026326985851" );
    gsl_test_rel (t, expected, rel, NAME(gsl_stats) "_ttest");
  }

  {
    MpIeee expected=  (MpIeee)MpIeee( "0.1331" );
    gsl_test  (FUNCTION(gsl_stats,max) (groupa, stridea, na) != expected,
               NAME(gsl_stats) "_max (" OUT_FORMAT " observed vs " OUT_FORMAT " expected)", 
               FUNCTION(gsl_stats,max) (groupa, stridea, na), expected);
  }

  {
    MpIeee min=  FUNCTION(gsl_stats,min) (groupa, stridea, na);
    MpIeee expected=  (MpIeee)MpIeee( "0.0242" );
    gsl_test (min != expected,
              NAME(gsl_stats) "_min (" OUT_FORMAT " observed vs " OUT_FORMAT " expected)", 
              min, expected);
  }

  {
    MpIeee min;MpIeee  max;
    MpIeee expected_max=  (MpIeee)MpIeee( "0.1331" );
    MpIeee expected_min=  (MpIeee)MpIeee( "0.0242" );
    
    FUNCTION(gsl_stats,minmax) (&min, &max, groupa, stridea, na);
 
    gsl_test  (max != expected_max,
               NAME(gsl_stats) "_minmax max (" OUT_FORMAT " observed vs " OUT_FORMAT " expected)", 
               max, expected_max);
    gsl_test  (min != expected_min,
               NAME(gsl_stats) "_minmax min (" OUT_FORMAT " observed vs " OUT_FORMAT " expected)", 
               min, expected_min);
  }

  {
    int  max_index=  FUNCTION(gsl_stats,max_index) (groupa, stridea, na);
    int  expected=  4;
    gsl_test (max_index != expected,
              NAME(gsl_stats) "_max_index (%d observed vs %d expected)",
              max_index, expected);
  }

  {
    int  min_index=  FUNCTION(gsl_stats,min_index) (groupa, stridea, na);
    int  expected=  3;
    gsl_test (min_index != expected,
              NAME(gsl_stats) "_min_index (%d observed vs %d expected)",
              min_index, expected);
  }

  {
    size_t min_index, max_index;
    size_t expected_max_index = 4;
    size_t expected_min_index = 3;

    FUNCTION(gsl_stats,minmax_index) (&min_index, &max_index, groupa, stridea, na);

    gsl_test  (max_index != expected_max_index,
               NAME(gsl_stats) "_minmax_index max (%u observed vs %u expected)", 
               max_index, expected_max_index);
    gsl_test  (min_index != expected_min_index,
               NAME(gsl_stats) "_minmax_index min (%u observed vs %u expected)", 
               min_index, expected_min_index);
  }


  sorted = (MpIeee *) malloc(stridea * na * sizeof(MpIeee)) ;
  
  for (i = 0 ; i < na ; i++)
    sorted[stridea * i] = groupa[stridea * i] ;
  
  TYPE(gsl_sort)(sorted, stridea, na) ;

  {
    MpIeee median=  FUNCTION(gsl_stats,median_from_sorted_data)(sorted, stridea, na) ;
    MpIeee expected=  MpIeee( "0.07505" );
    gsl_test_rel  (median,expected, rel,
                   NAME(gsl_stats) "_median_from_sorted_data (even)");
  }

  {
    MpIeee median=  FUNCTION(gsl_stats,median_from_sorted_data)(sorted, stridea, na - MpIeee( "1" )) ;
    MpIeee expected=  MpIeee( "0.0728" );
    gsl_test_rel  (median,expected, rel,
                   NAME(gsl_stats) "_median_from_sorted_data");
  }


  {
    MpIeee zeroth=  FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, na, MpIeee( "0.0" )) ;
    MpIeee expected=  MpIeee( "0.0242" );
    gsl_test_rel  (zeroth,expected, rel,
                   NAME(gsl_stats) "_quantile_from_sorted_data (0)");
  }

  {
    MpIeee top=  FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, na, MpIeee( "1.0" )) ;
    MpIeee expected=  MpIeee( "0.1331" );
    gsl_test_rel  (top,expected, rel,
                   NAME(gsl_stats) "_quantile_from_sorted_data (100)");
  }

  {
    MpIeee median=  FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, na, MpIeee( "0.5" )) ;
    MpIeee expected=  MpIeee( "0.07505" );
    gsl_test_rel  (median,expected, rel,
                   NAME(gsl_stats) "_quantile_from_sorted_data (50even)");
  }

  {
    MpIeee median=  FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, na - MpIeee( "1" ), MpIeee( "0.5" ));
    MpIeee expected=  MpIeee( "0.0728" );
    gsl_test_rel  (median,expected, rel,
                   NAME(gsl_stats) "_quantile_from_sorted_data (50odd)");

  }

  free (sorted);
  free (groupa);
  free (groupb);
  free (w);
}
