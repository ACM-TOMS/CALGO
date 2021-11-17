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
 * list_domains.h                                                           *
 ****************************************************************************/

#ifndef __list_domains_h
#define __list_domains_h

#include "domain.h"
#include "variable.h"


/* Doubly linked lists of domains used by bisection algorithms */

typedef struct
{
  int var;         /* previous bisected variable in the round-robin strategy */
  IBDomains d;     /* domains */
} IBDListInfo;

#define IBDListInfoCopy(copy,source) \
  (copy).var = (source).var; \
  IBCopyD((copy).d,(source).d,IBVnb(variables))

#define IBDListInfoWrite(x) \
   IBWriteVdom(stdout,x.d,variables,IBPragmaIntervalDigits,IBPragmaStyleInterval)

#define IBDListInfoAlloc(x) (x).d = IBNewD(IBVnb(variables))

#define IBDListInfoDesalloc(x) free((x).d)

typedef struct
{
  IBDListInfo info;
} IBDListElement;

#define IBDListAllocUnit 10

typedef struct
{
  IBDListElement elem[IBDListAllocUnit];
  int first;
  int last;
  int nb;
} IBDListArrayElements;

#define IBDListArrayNextIndex(x) x = (x+1) % IBDListAllocUnit
#define IBDListArrayPrevIndex(x) x = (((x)==0) ? (IBDListAllocUnit-1) : ((x)-1))
#define IBDListArrayFull(a) a->elem.nb==IBDListAllocUnit


/* Doubly linked lists */
typedef struct IBDListCell
{
  IBDListArrayElements elem;
  struct IBDListCell* prev;
  struct IBDListCell* next;
} IBDListCell;

typedef struct
{
  IBDListCell* first;
  IBDListCell* last;
  int ncell;
} IBDList;


#define IBDListEmpty(l)    ((l->first==NULL) || (l->first->elem.nb==0))
#define IBDListNotEmpty(l) (!IBDListEmpty(l))

#define IBDListNbElements(l) \
   ( (l->ncell==0) ? \
        0 : \
        ( (l->ncell==1) ? \
             (l->first->elem.nb) : \
             ((l->ncell-2)*IBDListAllocUnit + \
                 l->first->elem.nb + \
                 l->last->elem.nb) ) )


/* Functions that return the next domain to consider */
typedef IBDomains   (* IBDListGetDomain)(IBDList *);

static inline IBDomains IBDListGetFirstDomain(IBDList* l) {
  return l->first->elem.elem[l->first->elem.first].info.d;
}

static inline IBDomains IBDListGetLastDomain(IBDList* l) {
  return l->last->elem.elem[l->last->elem.last].info.d;
}


/* Functions that return the variable associated to the next domain to consider */
typedef int         (* IBDListGetVar)(IBDList *);
#define IBDListGetFirstVar(l) l->first->elem.elem[l->first->elem.first].info.var
#define IBDListGetLastVar(l)  l->last->elem.elem[l->last->elem.last].info.var


#define IBDListSetLastVar(l,v) l->last->elem.elem[l->last->elem.last].info.var = v


IBDListCell* IBDListNewCell          ();
IBDList*     IBDListNew              ();
void         IBDListFree             (IBDList* l);
IBDomains    IBDListAddLast          (IBDList* l, IBDListInfo x);
IBDomains    IBDListAddLastSugar     (IBDList* l, int var, IBDomains d);


/* Functions that remove the domain which has been considered */
typedef IBDListInfo (* IBDListRemove)(IBDList *);
typedef IBDomains   (* IBDListRemoveSugar)(IBDList *);
IBDListInfo  IBDListRemoveLast       (IBDList* l);
IBDomains    IBDListRemoveLastSugar  (IBDList* l);
IBDListInfo  IBDListRemoveFirst      (IBDList* l);
IBDomains    IBDListRemoveFirstSugar (IBDList* l);


/* Used for debugging */
void         IBDListWrite            (IBDList* l);

#endif
