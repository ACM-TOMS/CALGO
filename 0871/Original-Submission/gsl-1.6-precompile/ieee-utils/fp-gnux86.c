#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* ieee-utils/fp-gnux86.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

#include <stdio.h>
#include <fpu_control.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_ieee_utils.h>

  /* Handle libc5, where _FPU_SETCW is not available, suggested by
     OKUJI Yoshinori <okuji@gnu.org> and Evgeny Stambulchik
     <fnevgeny@plasma-gate.weizmann.ac.il> */

#ifndef _FPU_SETCW
#include <i386/fpu_control.h>
#define _FPU_SETCW(cw) __setfpucw(cw)
#endif

int
 gsl_ieee_set_mode(int  precision, int  rounding, int  exception_mask)
{
  unsigned short mode = 0 ;

  switch (precision)
    {
    case GSL_IEEE_SINGLE_PRECISION:
      mode |= _FPU_SINGLE ;
      break ;
    case GSL_IEEE_DOUBLE_PRECISION:
      mode |= _FPU_DOUBLE ;
      break ;
    case GSL_IEEE_EXTENDED_PRECISION:
      mode |= _FPU_EXTENDED ;
      break ;
    default:
      mode |= _FPU_EXTENDED ;
    }

  switch (rounding)
    {
    case GSL_IEEE_ROUND_TO_NEAREST:
      mode |= _FPU_RC_NEAREST ;
      break ;
    case GSL_IEEE_ROUND_DOWN:
      mode |= _FPU_RC_DOWN ;
      break ;
    case GSL_IEEE_ROUND_UP:
      mode |= _FPU_RC_UP ;
      break ;
    case GSL_IEEE_ROUND_TO_ZERO:
      mode |= _FPU_RC_ZERO ;
      break ;
    default:
      mode |= _FPU_RC_NEAREST ;
    }

  if (exception_mask & GSL_IEEE_MASK_INVALID)
    mode |= _FPU_MASK_IM ;

  if (exception_mask & GSL_IEEE_MASK_DENORMALIZED)
    mode |= _FPU_MASK_DM ;

  if (exception_mask & GSL_IEEE_MASK_DIVISION_BY_ZERO)
    mode |= _FPU_MASK_ZM ;

  if (exception_mask & GSL_IEEE_MASK_OVERFLOW)
    mode |= _FPU_MASK_OM ;

  if (exception_mask & GSL_IEEE_MASK_UNDERFLOW)
    mode |= _FPU_MASK_UM ;

  if (exception_mask & GSL_IEEE_TRAP_INEXACT)
    {
      mode &= ~ _FPU_MASK_PM ;
    }
  else
    {
      mode |= _FPU_MASK_PM ;
    }

  _FPU_SETCW(mode) ;

  return GSL_SUCCESS ;
}
