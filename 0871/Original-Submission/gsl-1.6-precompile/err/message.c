#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* err/message.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_message.h>

unsigned int  gsl_message_mask=  GSL_MESSAGE_MASK;

void
gsl_message (const char * reason, const char * file, int  line, 
             unsigned int  mask)
{
  if (mask & gsl_message_mask)
    {
      gsl_stream_printf ("MESSAGE", file, line, reason);
    }
}
