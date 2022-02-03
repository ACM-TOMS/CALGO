#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* sys/ldfrexp.c
 * 
 * Copyright (C) 2002, Gert Van den Eynde
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
#include <gsl/gsl_math.h>

MpIeee gsl_ldexp(const MpIeee x, const int e)
{
  MpIeee p2=  pow (MpIeee( "2.0" ), (MpIeee)e);
  return x * p2;
}

MpIeee gsl_frexp(const MpIeee x, int  *e)
{
  if (x == 0.0)
    {
      *e = 0;
      return MpIeee( "0.0" );
    }
  else
    {
      MpIeee ex=  ceil (log (fabs (x)) / M_LN2);
      int  ei=  (int) ex.toInt();
      MpIeee f=  gsl_ldexp (x, -ei);

      while (fabs (f) >= 1.0)
        {
          ei++;
          f /= MpIeee( "2.0" );
        }
      
      while (fabs (f) < 0.5)
        {
          ei--;
          f *= MpIeee( "2.0" );
        }

      *e = ei;
      return f;
    }
}
