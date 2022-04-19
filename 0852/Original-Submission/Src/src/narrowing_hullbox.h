/****************************************************************************
 * RealPaver v. 0.4                                                         *
 *--------------------------------------------------------------------------*
 * Author: Laurent Granvilliers                                             *
 * Copyright (c) 1999-2003 Institut de Recherche en Informatique de Nantes  *
 * Copyright (c) 2004      Laboratoire d'Informatique de Nantes Atlantique  *
 *--------------------------------------------------------------------------*
 * RealPaver is distributed WITHOUT ANY WARRANTY. Read the associated       *
 * COPYRIGHT file for more details.                                         *
 *--------------------------------------------------------------------------*
 * narrowing_hullbox.h                                                      *
 ****************************************************************************/

#ifndef __narrowing_hullbox_h
#define __narrowing_hullbox_h


#include "constant.h"

#define IBNarrowFailure  0
#define IBNarrowSuccess  1
#define IBNarrowFlounder 2

#define IBAllocDichotomicSearch 20

/*-- Search for leftmost quasi-zeros */
int IBNarrowShrinkLeft  (IBConstraint *c, int locvar, int globvar, IBDomains d,
                         IBInterval *in, IBInterval *out, double precision);

/*-- Search for rightmost quasi-zeros */
int IBNarrowShrinkRight (IBConstraint *c, int locvar, int globvar, IBDomains d,
                         IBInterval *in, IBInterval *out, double precision);

/*-- Narrowing operator for box consistency over projection (c,globvar) */
int IBNarrowBC3revise   (IBConstraint *c, int locvar, int globvar, IBDomains d,
                         IBInterval *out, double precision);

/*-- Narrowing operators for hull consistency over c */
int IBNarrowHC3revise (IBConstraint *c, IBDomains d);
int IBNarrowHC4revise (IBConstraint *c, IBDomains d);

#endif
