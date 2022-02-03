#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* sys/infnan.c
 * 
 * Copyright (C) 2001, 2004 Brian Gough
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
#include <math.h>

MpIeee gsl_nan(void);
MpIeee gsl_posinf(void);
MpIeee gsl_neginf(void);
MpIeee gsl_fdiv(const MpIeee x, const MpIeee y);

MpIeee gsl_nan(void)
{
  return gsl_fdiv (MpIeee( "0.0" ), MpIeee( "0.0" ));
}

MpIeee gsl_posinf(void)
{
  return gsl_fdiv (MpIeee( "1.0" ), MpIeee( "0.0" ));
}

MpIeee gsl_neginf(void)
{
  return gsl_fdiv (-MpIeee( "1.0" ), MpIeee( "0.0" ));
}



int
 gsl_isnan(const MpIeee x)
{
  x.isNan();
}

int
 gsl_isinf(const MpIeee x)
{
  return x.isInf();
}

int gsl_finite(const MpIeee x)
{
  return x.isReal();
}




#if defined(_MSC_VER) /* Microsoft Visual C++ */
#include <float.h>
int
 gsl_isnan(const MpIeee x)
{
  return _isnan(x);
}

int
 gsl_isinf(const MpIeee x)
{
  int  fpc=  _fpclass(x);

  if (fpc == _FPCLASS_PINF)
    return +1;
  else if (fpc == _FPCLASS_NINF)
    return -1;
  else 
    return 0;
}

int
 gsl_finite(const MpIeee x)
{
  return _finite(x);
}
#elif HAVE_IEEE_COMPARISONS

#else

#if HAVE_DECL_ISNAN
int
 gsl_isnan(const MpIeee x)
{
  return isnan(x);
}
#endif

#if HAVE_DECL_ISINF
int
 gsl_isinf(const MpIeee x)
{
    return isinf(x);
}
#else
int
 gsl_isinf(const MpIeee x)
{
    return (! gsl_finite(x)) && (! gsl_isnan(x));
}
#endif


#if HAVE_DECL_FINITE
int
 gsl_finite(const MpIeee x)
{
  return finite(x);
}
#elif HAVE_DECL_ISFINITE
int
 gsl_finite(const MpIeee x)
{
  return isfinite(x);
}
#endif

#endif

