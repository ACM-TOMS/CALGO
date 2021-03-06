#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* statistics/test_int_source.c
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
  /* sample sets of integers */
  size_t i;
  const size_t ina = 20, inb = 20;

  const BASE raw1[] = {1, 2, 3, 4, 5, 6} ;
  
  const BASE irawa[] =
  {17, 18, 16, 18, 12,
   20, 18, 20, 20, 22,
   20, 10, 8, 12, 16,
   16, 18, 20, 18, 21};

  const BASE irawb[] =
  {19, 20, 22, 24, 10,
   25, 20, 22, 21, 23,
   20, 10, 12, 14, 12,
   20, 22, 24, 23, 17};

  MpIeee * sorted;

  MpIeee * test1=  (MpIeee *) malloc (stridea * MpIeee( "6" ) * sizeof(MpIeee));
  MpIeee * igroupa=  (MpIeee *) malloc (stridea * ina * sizeof(MpIeee));
  MpIeee * igroupb=  (MpIeee *) malloc (strideb * inb * sizeof(MpIeee));

  MpIeee rel=  MpIeee( "1" )e-MpIeee( "10" ) ;

  for (i = 0 ; i < ina ; i++)
    igroupa[i * stridea] = irawa[i] ;

  for (i = 0 ; i < inb ; i++)
    igroupb[i * strideb] = irawb[i] ;

  for (i = 0 ; i < 6 ; i++)
    test1[i * stridea] = raw1[i] ;



  {
    MpIeee mean=  FUNCTION(gsl_stats,mean) (igroupa, stridea, ina);
    MpIeee expected=  MpIeee( "17.0" );
    gsl_test_rel (mean,expected, rel, NAME(gsl_stats) "_mean (integer)");
  }

  {
    MpIeee mean=  FUNCTION(gsl_stats,mean) (test1, stridea, MpIeee( "6" ));
    MpIeee expected=  MpIeee( "3.5" );
    gsl_test_rel (mean,expected, rel, NAME(gsl_stats) "_mean (fractional)");
  }

  {
    MpIeee mean=  FUNCTION(gsl_stats,mean) (igroupa, stridea, ina);
    MpIeee var=  FUNCTION(gsl_stats,variance_with_fixed_mean) (igroupa, stridea, ina, mean);
    MpIeee expected=  MpIeee( "13.7" );
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_variance_with_fixed_mean");
  }

  {
    MpIeee mean=  FUNCTION(gsl_stats,mean) (igroupa, stridea, ina);
    MpIeee sd=  FUNCTION(gsl_stats,sd_with_fixed_mean) (igroupa, stridea, ina, mean);
    MpIeee expected=  MpIeee( "3.70135110466435" );
    gsl_test_rel (sd, expected, rel, NAME(gsl_stats) "_sd_with_fixed_mean");
  }

  {
    MpIeee var=  FUNCTION(gsl_stats,variance) (igroupa, stridea, ina);
    MpIeee expected=  MpIeee( "14.4210526315789" );
    gsl_test_rel (var, expected, rel, NAME(gsl_stats) "_variance");
  }

  {
    MpIeee sd_est=  FUNCTION(gsl_stats,sd) (igroupa, stridea, ina);
    MpIeee expected=  MpIeee( "3.79750610685209" );
    gsl_test_rel (sd_est, expected, rel, NAME(gsl_stats) "_sd");
  }

  {
    MpIeee absdev=  FUNCTION(gsl_stats,absdev) (igroupa, stridea, ina);
    MpIeee expected=  MpIeee( "2.9" );
    gsl_test_rel (absdev, expected, rel, NAME(gsl_stats) "_absdev");
  }

  {
    MpIeee skew=  FUNCTION(gsl_stats,skew) (igroupa, stridea, ina);
    MpIeee expected=  -MpIeee( "0.909355923168064" );
    gsl_test_rel (skew, expected, rel, NAME(gsl_stats) "_skew");
  }

  {
    MpIeee kurt=  FUNCTION(gsl_stats,kurtosis) (igroupa, stridea, ina);
    MpIeee expected=  -MpIeee( "0.233692524908094" ) ;
    gsl_test_rel (kurt, expected, rel, NAME(gsl_stats) "_kurtosis");
  }

  {
    MpIeee c=  FUNCTION(gsl_stats,covariance) (igroupa, stridea, igroupb, strideb, inb);
    MpIeee expected=  MpIeee( "14.5263157894737" );
    gsl_test_rel (c, expected, rel, NAME(gsl_stats) "_covariance");
  }


  {
    MpIeee pv=  FUNCTION(gsl_stats,pvariance) (igroupa, stridea, ina, igroupb, strideb, inb);
    MpIeee expected=  MpIeee( "18.8421052631579" );
    gsl_test_rel (pv, expected, rel, NAME(gsl_stats) "_pvariance");
  }

  {
    MpIeee t=  FUNCTION(gsl_stats,ttest) (igroupa, stridea, ina, igroupb, strideb, inb);
    MpIeee expected=  -MpIeee( "1.45701922702927" );
    gsl_test_rel (t, expected, rel, NAME(gsl_stats) "_ttest");
  }

  {
    int  max=  FUNCTION(gsl_stats,max) (igroupa, stridea, ina);
    int  expected=  22;
    gsl_test (max != expected,
              NAME(gsl_stats) "_max (%d observed vs %d expected)", max, expected);
  }

  {
    int  min=  FUNCTION(gsl_stats,min) (igroupa, stridea, ina);
    int  expected=  8;
    gsl_test (min != expected,
              NAME(gsl_stats) "_min (%d observed vs %d expected)", min, expected);
  }

  {
    MpIeee min;MpIeee  max;
    MpIeee expected_max=  MpIeee( "22" );
    MpIeee expected_min=  MpIeee( "8" );
    
    FUNCTION(gsl_stats,minmax) (&min, &max, igroupa, stridea, ina);
 
    gsl_test  (max != expected_max,
               NAME(gsl_stats) "_minmax max (" OUT_FORMAT " observed vs " OUT_FORMAT " expected)", 
               max, expected_max);
    gsl_test  (min != expected_min,
               NAME(gsl_stats) "_minmax min (" OUT_FORMAT " observed vs " OUT_FORMAT " expected)", 
               min, expected_min);
  }

  {
    int  max_index=  FUNCTION(gsl_stats,max_index) (igroupa, stridea, ina);
    int  expected=  9 ;
    gsl_test (max_index != expected,
              NAME(gsl_stats) "_max_index (%d observed vs %d expected)",
              max_index, expected);
  }

  {
    int  min_index=  FUNCTION(gsl_stats,min_index) (igroupa, stridea, ina);
    int  expected=  12 ;
    gsl_test (min_index != expected,
              NAME(gsl_stats) "_min_index (%d observed vs %d expected)",
              min_index, expected);
  }

  {
    size_t min_index, max_index;
    size_t expected_max_index = 9;
    size_t expected_min_index = 12;

    FUNCTION(gsl_stats,minmax_index) (&min_index, &max_index, igroupa, stridea, ina);

    gsl_test  (max_index != expected_max_index,
               NAME(gsl_stats) "_minmax_index max (%u observed vs %u expected)", 
               max_index, expected_max_index);
    gsl_test  (min_index != expected_min_index,
               NAME(gsl_stats) "_minmax_index min (%u observed vs %u expected)", 
               min_index, expected_min_index);
  }


  sorted = (MpIeee *) malloc(stridea * ina * sizeof(MpIeee)) ;

  for (i = 0 ; i < ina ; i++)
    sorted[stridea * i] = igroupa[stridea * i] ;


  TYPE(gsl_sort)(sorted, stridea, ina) ;

  {
    MpIeee median=  FUNCTION(gsl_stats,median_from_sorted_data)(sorted, stridea, ina) ;
    MpIeee expected=  MpIeee( "18" );
    gsl_test_rel (median,expected, rel,
                  NAME(gsl_stats) "_median_from_sorted_data (even)");
  }

  {
    MpIeee median=  FUNCTION(gsl_stats,median_from_sorted_data)(sorted, stridea, ina - MpIeee( "1" )) ;
    MpIeee expected=  MpIeee( "18" );
    gsl_test_rel (median,expected, rel,
                  NAME(gsl_stats) "_median_from_sorted_data (odd)");
  }


  {
    MpIeee zeroth=  FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, ina, MpIeee( "0.0" )) ;
    MpIeee expected=  MpIeee( "8" );
    gsl_test_rel (zeroth,expected, rel,
                  NAME(gsl_stats) "_quantile_from_sorted_data (0)");
  }

  {
    MpIeee top=  FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, ina, MpIeee( "1.0" )) ;
    MpIeee expected=  MpIeee( "22" );
    gsl_test_rel (top,expected, rel,
                  NAME(gsl_stats) "_quantile_from_sorted_data (100)");
  }

  {
    MpIeee median=  FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, ina, MpIeee( "0.5" )) ;
    MpIeee expected=  MpIeee( "18" );
    gsl_test_rel (median,expected, rel,
                  NAME(gsl_stats) "_quantile_from_sorted_data (50, even)");
  }

  {
    MpIeee median=  FUNCTION(gsl_stats,quantile_from_sorted_data)(sorted, stridea, ina - MpIeee( "1" ), MpIeee( "0.5" ));
    MpIeee expected=  MpIeee( "18" );
    gsl_test_rel (median,expected, rel,
                  NAME(gsl_stats) "_quantile_from_sorted_data (50, odd)");
  }

  free (sorted);
  free (igroupa);
  free (igroupb);
  free (test1);
}
