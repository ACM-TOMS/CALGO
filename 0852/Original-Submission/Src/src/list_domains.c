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
 * list_domains.c                                                           *
 ****************************************************************************/


#include "list_domains.h"
#include <stdlib.h>


extern IBVariables   variables;           /* array of constrained variables */
extern int    IBPragmaStyleInterval;      /* used for interval printing */
extern int    IBPragmaIntervalDigits;     /* number of digits for interval printing */



IBDListCell* IBDListNewCell()
/***************************************************************************
*  Allocation of a new list cell
*/
{
  int i;
  IBDListCell* lc;
  lc = (IBDListCell*)malloc(sizeof(IBDListCell));
  lc->prev = lc->next = NULL;
  lc->elem.nb = 0;
  for( i=0; i<IBDListAllocUnit; ++i )
  {
    IBDListInfoAlloc(lc->elem.elem[i].info);
  }
  return lc;
}

IBDList* IBDListNew()
/***************************************************************************
*  Allocation of a new list of domains
*/
{
  IBDList* l;
  l = (IBDList*)malloc(sizeof(IBDList));
  l->first = l->last = IBDListNewCell();
  l->ncell = 1;
  return( l );
}

void IBDListFree(IBDList* l)
/***************************************************************************
*  Desallocation of a new list of domains
*/
{
  IBDListCell* cl = l->first, *pl;
  int i;
  while( cl!=NULL )
  {
    pl = cl;
    cl = cl->next;
    for( i=0; i<IBDListAllocUnit; ++i )
    {
      IBDListInfoDesalloc(pl->elem.elem[i].info);
    }
    free(pl);
  }
  free(l);
}

IBDomains IBDListAddLast(IBDList* l, IBDListInfo x)
/***************************************************************************
*  Insertion of an element in the last position
*/
{
  if( IBDListArrayFull(l->last) )
  {
    if( l->last->next==NULL )
    {
      l->last->next = IBDListNewCell();
      l->last->next->prev = l->last;
      l->last = l->last->next;
    }
    else
    {
      l->last = l->last->next;
    }
    l->ncell ++;
  }
  if( l->last->elem.nb==0 )
  {
    l->last->elem.first = 0;
    l->last->elem.last = -1;
  }
  IBDListArrayNextIndex(l->last->elem.last);
  IBDListInfoCopy(l->last->elem.elem[l->last->elem.last].info,x);
  l->last->elem.nb++;

  return l->last->elem.elem[l->last->elem.last].info.d;
}


inline IBDomains IBDListAddLastSugar(IBDList* l, int var, IBDomains d)
/***************************************************************************
*  Insertion of an element in the last position
*/
{
  IBDListInfo info;
  info.var = var;
  info.d = d;
  return IBDListAddLast(l,info);
}


IBDListInfo IBDListRemoveLast(IBDList* l)
/***************************************************************************
*  Returns and removes the element in the last position
*/
{
  IBDListCell* cl;
  IBDListInfo x = l->last->elem.elem[l->last->elem.last].info;
  l->last->elem.nb--;

  if( l->last->elem.nb==0 )
  {
    if( l->last == l->first )
    {
      /* nothing to do */
    }
    else
    {
      l->last = l->last->prev;
      l->ncell --;
    }
  }
  else
  {
    IBDListArrayPrevIndex(l->last->elem.last);
  }
  return( x );
}


inline IBDomains IBDListRemoveLastSugar(IBDList* l)
/***************************************************************************
*  Returns and removes the domain in the last position
*/
{
  IBDListInfo info = IBDListRemoveLast(l);
  return info.d;
}


IBDListInfo IBDListRemoveFirst(IBDList* l)
/***************************************************************************
*  Returns and removes the element in the last position
*/
{
  IBDListInfo x = l->first->elem.elem[l->first->elem.first].info;
  l->first->elem.nb--;

  if( l->first->elem.nb==0 )
  {
    if( l->first == l->last )
    {
      /* nothing to do */
    }
    else
    {
      l->first = l->first->next;
      free(l->first->prev);
      l->first->prev = NULL;
      l->ncell --;
    }
  }
  else
  {
    IBDListArrayNextIndex(l->first->elem.first);
  }
  return( x );
}


inline IBDomains IBDListRemoveFirstSugar(IBDList* l)
/***************************************************************************
*  Returns and removes the domain in the first position
*/
{
  IBDListInfo info = IBDListRemoveFirst(l);
  return info.d;
}


void IBDListWrite(IBDList* l)
/***************************************************************************
*  Display of l, used for debugging
*/
{
  IBDListCell* cl = l->first;
  int i, n;

  while( cl!=NULL )
  {
    for( i=cl->elem.first, n=0; n<cl->elem.nb; IBDListArrayNextIndex(i), ++n)
    {
      /* IBDListInfoWrite(cl->elem.elem[i].info); */
      printf("#");
      printf(" ");
    }
    printf(" @@ ");
    cl = cl->next;
  }
}
