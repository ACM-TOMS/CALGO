#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* specfunc/bessel_amp_phase.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include "bessel_amp_phase.h"

/* chebyshev expansions for amplitude and phase
   functions used in bessel evaluations

   These are the same for J0,Y0 and for J1,Y1, so
   they sit outside those functions.
*/
        
static MpIeee bm0_data[21] =  {
   MpIeee( "0.09284961637381644" ),
  -MpIeee( "0.00142987707403484" ),
   MpIeee( "0.00002830579271257" ),
  -MpIeee( "0.00000143300611424" ),
   MpIeee( "0.00000012028628046" ),
  -MpIeee( "0.00000001397113013" ),
   MpIeee( "0.00000000204076188" ),
  -MpIeee( "0.00000000035399669" ),
   MpIeee( "0.00000000007024759" ),
  -MpIeee( "0.00000000001554107" ),
   MpIeee( "0.00000000000376226" ),
  -MpIeee( "0.00000000000098282" ),
   MpIeee( "0.00000000000027408" ),
  -MpIeee( "0.00000000000008091" ),
   MpIeee( "0.00000000000002511" ),
  -MpIeee( "0.00000000000000814" ),
   MpIeee( "0.00000000000000275" ),
  -MpIeee( "0.00000000000000096" ),
   MpIeee( "0.00000000000000034" ),
  -MpIeee( "0.00000000000000012" ),
   MpIeee( "0.00000000000000004" )
}; 
const cheb_series _gsl_sf_bessel_amp_phase_bm0_cs = {
  bm0_data,
  20,
  -1, 1,
  10
};
      
static MpIeee bth0_data[24] =  {
  -MpIeee( "0.24639163774300119" ),
   MpIeee( "0.001737098307508963" ),
  -MpIeee( "0.000062183633402968" ),
   MpIeee( "0.000004368050165742" ),
  -MpIeee( "0.000000456093019869" ),
   MpIeee( "0.000000062197400101" ),
  -MpIeee( "0.000000010300442889" ),
   MpIeee( "0.000000001979526776" ),
  -MpIeee( "0.000000000428198396" ),
   MpIeee( "0.000000000102035840" ),
  -MpIeee( "0.000000000026363898" ),
   MpIeee( "0.000000000007297935" ),
  -MpIeee( "0.000000000002144188" ),
   MpIeee( "0.000000000000663693" ),
  -MpIeee( "0.000000000000215126" ),
   MpIeee( "0.000000000000072659" ),
  -MpIeee( "0.000000000000025465" ),
   MpIeee( "0.000000000000009229" ),
  -MpIeee( "0.000000000000003448" ),
   MpIeee( "0.000000000000001325" ),
  -MpIeee( "0.000000000000000522" ),
   MpIeee( "0.000000000000000210" ),
  -MpIeee( "0.000000000000000087" ),
   MpIeee( "0.000000000000000036" )
};
const cheb_series _gsl_sf_bessel_amp_phase_bth0_cs = {
  bth0_data,
  23,
  -1, 1,
  12
};


static MpIeee bm1_data[21] =  {
   MpIeee( "0.1047362510931285" ), 
   MpIeee( "0.00442443893702345" ),
  -MpIeee( "0.00005661639504035" ),
   MpIeee( "0.00000231349417339" ),
  -MpIeee( "0.00000017377182007" ),
   MpIeee( "0.00000001893209930" ),
  -MpIeee( "0.00000000265416023" ),
   MpIeee( "0.00000000044740209" ),
  -MpIeee( "0.00000000008691795" ),
   MpIeee( "0.00000000001891492" ),
  -MpIeee( "0.00000000000451884" ),
   MpIeee( "0.00000000000116765" ),
  -MpIeee( "0.00000000000032265" ),
   MpIeee( "0.00000000000009450" ),
  -MpIeee( "0.00000000000002913" ),
   MpIeee( "0.00000000000000939" ),
  -MpIeee( "0.00000000000000315" ),
   MpIeee( "0.00000000000000109" ),
  -MpIeee( "0.00000000000000039" ),
   MpIeee( "0.00000000000000014" ),
  -MpIeee( "0.00000000000000005" ),
}; 
const cheb_series _gsl_sf_bessel_amp_phase_bm1_cs = {
  bm1_data,
  20,
  -1, 1,
  10
};

static MpIeee bth1_data[24] =  {
   MpIeee( "0.74060141026313850" ), 
  -MpIeee( "0.004571755659637690" ),
   MpIeee( "0.000119818510964326" ),
  -MpIeee( "0.000006964561891648" ),
   MpIeee( "0.000000655495621447" ),
  -MpIeee( "0.000000084066228945" ),
   MpIeee( "0.000000013376886564" ),
  -MpIeee( "0.000000002499565654" ),
   MpIeee( "0.000000000529495100" ),
  -MpIeee( "0.000000000124135944" ),
   MpIeee( "0.000000000031656485" ),
  -MpIeee( "0.000000000008668640" ),
   MpIeee( "0.000000000002523758" ),
  -MpIeee( "0.000000000000775085" ),
   MpIeee( "0.000000000000249527" ),
  -MpIeee( "0.000000000000083773" ),
   MpIeee( "0.000000000000029205" ),
  -MpIeee( "0.000000000000010534" ),
   MpIeee( "0.000000000000003919" ),
  -MpIeee( "0.000000000000001500" ),
   MpIeee( "0.000000000000000589" ),
  -MpIeee( "0.000000000000000237" ),
   MpIeee( "0.000000000000000097" ),
  -MpIeee( "0.000000000000000040" ),
};
const cheb_series _gsl_sf_bessel_amp_phase_bth1_cs = {
  bth1_data,
  23,
  -1, 1,
  12
};


int
 gsl_sf_bessel_asymp_Mnu_e(const MpIeee nu, const MpIeee x, MpIeee * result)
{
  const MpIeee r=  2.0*nu/x;
  const MpIeee r2=  r*r;
  const MpIeee x2=  x*x;
  const MpIeee term1=  (r2-1.0/x2)/8.0;
  const MpIeee term2=  (r2-1.0/x2)*(r2-9.0/x2)*3.0/128.0;
  const MpIeee Mnu2_c=  2.0/(M_PI) * (1.0 + term1 + term2);
  *result = sqrt(Mnu2_c)/sqrt(x); /* will never underflow this way */
  return GSL_SUCCESS;
}


int
 gsl_sf_bessel_asymp_thetanu_corr_e(const MpIeee nu, const MpIeee x, MpIeee * result)
{
  const MpIeee r=  2.0*nu/x;
  const MpIeee r2=  r*r;
  const MpIeee x2=  x*x;
  const MpIeee term1=  x*(r2 - 1.0/x2)/8.0;
  const MpIeee term2=  x*(r2 - 1.0/x2)*(r2 - 25.0/x2)/384.0;
  *result = (-MpIeee( "0.25" )*M_PI + term1 + term2);
  return GSL_SUCCESS;
}
