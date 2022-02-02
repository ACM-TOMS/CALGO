/* const/test.c
 * 
 * Copyright (C) 2003 Brian Gough
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

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_const.h>
#include <gsl/gsl_test.h>

#include <gsl/gsl_ieee_utils.h>

int
main (void)
{
  gsl_ieee_env_setup ();

  /* Basic check to make sure the header files are functioning */

  {
    double c = GSL_CONST_MKS_SPEED_OF_LIGHT;
    double eps = GSL_CONST_MKS_VACUUM_PERMITTIVITY;
    double mu = GSL_CONST_MKS_VACUUM_PERMEABILITY;

    gsl_test_rel (c, 1.0/sqrt(eps*mu), 1e-6, "speed of light (mks)");
  }

  {
    double ly = GSL_CONST_CGS_LIGHT_YEAR;
    double c = GSL_CONST_CGS_SPEED_OF_LIGHT;
    double y = 365.2425 * GSL_CONST_CGS_DAY;
    
    gsl_test_rel (ly, c * y, 1e-6, "light year (cgs)");
  }

  {
    double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
    double eps = GSL_CONST_MKSA_VACUUM_PERMITTIVITY;
    double mu = GSL_CONST_MKSA_VACUUM_PERMEABILITY;

    gsl_test_rel (c, 1.0/sqrt(eps*mu), 1e-6, "speed of light (mksa)");
  }

  {
    double ly = GSL_CONST_CGSM_LIGHT_YEAR;
    double c = GSL_CONST_CGSM_SPEED_OF_LIGHT;
    double y = 365.2425 * GSL_CONST_CGSM_DAY;
    
    gsl_test_rel (ly, c * y, 1e-6, "light year (cgsm)");
  }

  {
    double micro = GSL_CONST_NUM_MICRO;
    double mega = GSL_CONST_NUM_MEGA;
    double kilo = GSL_CONST_NUM_KILO;

    gsl_test_rel (mega/kilo, 1/(micro*kilo), 1e-10, "kilo (mega/kilo, 1/(micro*kilo))");
  }

  exit (gsl_test_summary ());
}

