#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

/* histogram/ntuple.h
 * 
 * Copyright (C) 2000 Simone Piccardi
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
 *
 */

/* Jan/2001 Modified by Brian Gough. Minor changes for GSL */

#ifndef __GSL_NTUPLE_H__
#define __GSL_NTUPLE_H__

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef struct {
    FILE * file;
    void * ntuple_data;
    size_t size;
} gsl_ntuple;

typedef struct {
  int (* function) (void * ntuple_data, void * params);
  void * params;
} gsl_ntuple_select_fn;

typedef struct {
  MpIeee(* function) (void * ntuple_data, void * params);
  void * params;
} gsl_ntuple_value_fn;

gsl_ntuple * 
gsl_ntuple_open (char * filename, void * ntuple_data, size_t size);

gsl_ntuple * 
gsl_ntuple_create (char * filename, void * ntuple_data, size_t size);

int  gsl_ntuple_write(gsl_ntuple * ntuple);
int  gsl_ntuple_read(gsl_ntuple * ntuple);

int  gsl_ntuple_bookdata(gsl_ntuple * ntuple);  /* synonym for write */

int  gsl_ntuple_project(gsl_histogram * h, gsl_ntuple * ntuple, 
                        gsl_ntuple_value_fn *value_func,
                        gsl_ntuple_select_fn *select_func);

int  gsl_ntuple_close(gsl_ntuple * ntuple);

__END_DECLS

#endif /* __GSL_NTUPLE_H__ */




