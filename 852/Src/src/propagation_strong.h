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
 * propagation_strong.h                                                     *
 ****************************************************************************/

#ifndef __propagation_strong_h
#define __propagation_strong_h

#include "propagation_local.h"


/* Generic type for propagation algorithms based on local (2B) consistency */
typedef int (* IBLocalPropagation)(IBDomains, IBDmodified *);

#define IBF2Bhc3        IBFilteringHC3decomp
#define IBF2Bhc4        IBFilteringHC4
#define IBF2Bhc4I       IBFilteringHC3
#define IBF2Bhc4Newton  IBFilteringHC4Newton
#define IBF2Bbc3        IBFilteringBC3
#define IBF2Bbc3Newton  IBFilteringBC3Newton
#define IBF2Bbc4        IBFilteringBC4
#define IBF2Bbc5        IBFilteringBC5


/* Generic type for propagation algorithms */
typedef int (* IBPropagation)(IBDomains, IBDmodified *, IBLocalPropagation);


/* 2B consistency-based propagation algorithms -> propagation algorithms */

static inline int IBFhc3       (IBDomains d, IBDmodified *dm, IBLocalPropagation useless) {
  return IBFilteringHC3decomp(d,dm);
}

static inline int IBFhc4I      (IBDomains d, IBDmodified *dm, IBLocalPropagation useless) {
  return IBFilteringHC3(d,dm);
}

static inline int IBFhc4       (IBDomains d, IBDmodified *dm, IBLocalPropagation useless) {
  return IBFilteringHC4(d,dm);
}

static inline int IBFhc4Newton (IBDomains d, IBDmodified *dm, IBLocalPropagation useless) {
  return IBFilteringHC4Newton(d,dm);
}

static inline int IBFbc3       (IBDomains d, IBDmodified *dm, IBLocalPropagation useless) {
  return IBFilteringBC3(d,dm);
}

static inline int IBFbc3Newton (IBDomains d, IBDmodified *dm, IBLocalPropagation useless) {
  return IBFilteringBC3Newton(d,dm);
}

static inline int IBFbc4       (IBDomains d, IBDmodified *dm, IBLocalPropagation useless) {
  return IBFilteringBC4(d,dm);
}

static inline int IBFbc5       (IBDomains d, IBDmodified *dm, IBLocalPropagation useless) {
  return IBFilteringBC5(d,dm);
}

/*--- FILTERING ALGORITHM for << weak 3B(w) CONSISTENCY >> */
int IBFilteringWeak3B(IBDomains d, IBDmodified *dmodified, IBLocalPropagation f2b);


/*--- FILTERING ALGORITHM for 3B(w) CONSISTENCY */
int IBFiltering3B(IBDomains d, IBDmodified *dmodified, IBLocalPropagation f2B);

#define IBF3B IBFiltering3B
#define IBF3Bweak IBFilteringWeak3B

#endif
