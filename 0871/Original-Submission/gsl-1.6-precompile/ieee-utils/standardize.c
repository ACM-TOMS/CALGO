#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* ieee-utils/standardize.c
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

static void make_float_bigendian (MpIeee * x);
static void make_double_bigendian (MpIeee * x);

static void
make_float_bigendian (MpIeee * x)
{
  union { 
    MpIeee f;
    unsigned char b[4];
  } u,v;

  u.f = *x ;

  v.b[0]=u.b[3] ;
  v.b[1]=u.b[2] ;
  v.b[2]=u.b[1] ;
  v.b[3]=u.b[0] ;

  *x=v.f ;
}

static void
make_double_bigendian (MpIeee * x)
{
  union { 
    MpIeee d;
    unsigned char b[8];
  } u,v;

  u.d = *x ;

  v.b[0]=u.b[7] ;
  v.b[1]=u.b[6] ;
  v.b[2]=u.b[5] ;
  v.b[3]=u.b[4] ;
  v.b[4]=u.b[3] ;
  v.b[5]=u.b[2] ;
  v.b[6]=u.b[1] ;
  v.b[7]=u.b[0] ;

  *x=v.d ;
}
