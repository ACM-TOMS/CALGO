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
 * search.h                                                                 *
 ****************************************************************************/

#ifndef __search_h
#define __search_h

#include "propagation_strong.h"
#include "list_domains.h"


/* Can the domain d on variable i be bisected, given the precision p
   and the array of variables av ? */
#define IBIsDomainBisectable(av,i,d,p) ( (IBIsBranchVar(av,i)) && \
                                         (!IBCanonicalI(IBDomV(d,i))) && \
                                         (IBWidthI(IBDomV(d,i))>p) )


/*------ Bisection functions */
typedef void (* IBBisectArity)(IBDList *, int, long*);

void IBDListBisect2 (IBDList *dlist, int var, long *nbdom);  /* 2 parts */
void IBDListBisect3 (IBDList *dlist, int var, long *nbdom);  /* 3 parts */


#define IBBsplit2 IBDListBisect2 
#define IBBsplit3 IBDListBisect3


/*------ Bisection Strategies */
#define IBBisectNone           0
#define IBBisectRoundRobin     1
#define IBBisectLargestFirst   2
#define IBBisectMaxNarrow      3


/* generic type for a bisection strategy function */
typedef int (* IBBisectVar)(IBDomains, IBDomains, int);

int IBBisectVariableLF (IBDomains d, IBDomains dold, int var);  /* largest first strategy */
int IBBisectVariableRR (IBDomains d, IBDomains dold, int var);  /* round-robin strategy */
int IBBisectVariableMN (IBDomains d, IBDomains dold, int var);  /* max reduction strategy */

#define IBBlf IBBisectVariableLF
#define IBBrr IBBisectVariableRR
#define IBBmn IBBisectVariableMN


#define IBBisectIterate(subpaving,nbsol,nbdlist,parts,number) \
   ( (nbdlist!=0) && \
     (subpaving ? \
      (nbsol+nbdlist+parts-1 <= number) : \
      (nbsol < number)) )


/*------ BISECTION ALGORITHM: returns the number of output boxes */
int IBBisection (IBDomains d, int Nobisect, int* completeProcess);

#endif
