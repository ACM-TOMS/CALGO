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
 * variable.c                                                               *
 ****************************************************************************/

#include <memory.h>
#include <string.h>
#include <stdlib.h>
#include "variable.h"


IBVariables IBNewV()
/***************************************************************************
*  Allocation of an array of variables
*/
{
  IBVariables a;
  a = (struct IBV *)malloc(sizeof(struct IBV));

  a->a     = (IBVar *)malloc(IBVarAllocUnit*sizeof(IBVar));
  a->d     = IBNewD(IBVarAllocUnit);
  a->N     = 0;
  a->Nfree = IBVarAllocUnit;

  return( a );
}


void IBFreeV(IBVariables a)
/***************************************************************************
*  Desallocation of an array of variables
*/
{
  int i;

  for( i=0; i<a->N; i++ )
  {
    free(a->a[i].name);
    if( a->a[i].dep!=NULL )
    {
      free(a->a[i].dep->a);
      free(a->a[i].dep);
    }
  }
  free(a->a);
  IBFreeD(a->d);
  free(a);
}


int IBIsPresentInV(IBVariables a, char *name)
/***************************************************************************
*  Returns -1 if variable name i not present in a, its index otherwise
*/
{
  int i;
  for( i=0; i<a->N; i++ )
  {
    if( strcmp(name,a->a[i].name)==0 ) return( i );
  }
  return( -1 );
}


int IBAddV(IBVariables a, char *name, int status)
/***************************************************************************
*  To add variable name in a
*  Returns its index in a
*/
{
  int i = IBIsPresentInV(a,name);
  if( i>=0 ) return( i );   /* name is already in a */

  if( a->Nfree==0 )         /* a is full */
  {
    a->a = realloc(a->a,(a->N + IBVarAllocUnit)*sizeof(IBVar));
    a->d = realloc(a->d,(a->N + IBVarAllocUnit)*sizeof(IBDom));
    a->Nfree = IBVarAllocUnit;
  }

  a->a[a->N].name = (char *)malloc((1+strlen(name))*sizeof(char));
  strcpy(a->a[a->N].name,name);

  a->a[a->N].dep = NULL;
  a->a[a->N].IsInt = 0;
  a->a[a->N].weight = 1;
  a->a[a->N].status = status;

  IBToLargestI(IBDomV(a->d,a->N));  /* set the domain to [-oo,+oo] */

  a->N++;
  a->Nfree--;
  return( a->N - 1 );
}


inline void IBSetIntV(IBVariables a, int var)
/***************************************************************************
*  var is an integer variable
*/
{
  a->a[var].IsInt = 1;
}


long IBNbFreshVar(IBVariables a)
/***************************************************************************
*  Returns the number of fresh variables in a
*/
{
  int i;
  long n = 0;
  for( i=0; i<a->N; i++ )
  {
    if( IBIsFreshVar(a,i) ) n++;
  }
  return( n );
}

void IBWriteV(FILE *out, IBVariables a)
/***************************************************************************
*  To write a
*/
{
  int i;
  for( i=0; i<a->N; i++ )
  {
    fprintf(out,"%s(%d) ",IBNameV(a,i),i);
  }
  fprintf(out,"\n");
}


void IBWriteVdom(FILE *out, IBDomains d, IBVariables a, int digits, int mode)
/***************************************************************************
*  To write the domains in d
*/
{
  int i, j;
  for( i=0; i<a->N; i++ )
  {
    if( IBIsUserVar(a,i) )  /* user variable, not a fresh variable used in HC3 */
    {
      fprintf(out,"  %s",IBNameV(a,i));
      for( j=0; j<a->maxname-strlen(a->a[i].name); j++ ) fprintf(out," ");

      if( IBIsDoubleI(IBDomV(d,i)) ) fprintf(out," = ");
      else if( mode==IBPrintIntervalBounds ) fprintf(out," in ");
      else fprintf(out," = ");

      IBWriteIverb(out,IBDomV(d,i),digits,mode);
      fprintf(out,"\n");
    }
  }
}


void IBReallocV(IBVariables a)
/***************************************************************************
*  Reallocation of structures; needed after the parsing
*/
{
  int max = 0, i;

  a->a = realloc(a->a,(a->N)*sizeof(IBVar));
  a->d = realloc(a->d,(a->N)*sizeof(IBDom));

  for( i=0; i<a->N; i++ )
  {
    max = IBMax(max,strlen(a->a[i].name));
  }
  a->maxname = max;
}


void IBInitDependencyV(IBVariables a, int n)
/***************************************************************************
*  Initialisation of the list of dependencies of all variables
*  n is the number of constraints
*  Using an unsigned long where one bit = one constraint
*  the number of long which are required is size
*/
{
  int i;
  int size = (int)ceil((double)n / (double)sizeof(unsigned long));

  for(i=0; i<a->N; i++ )
  {
    a->a[i].dep    = (IBDependencyV *)malloc(sizeof(IBDependencyV));
    a->a[i].dep->a = (IBDependencyWord *)malloc(size*sizeof(IBDependencyWord));
    a->a[i].dep->N = 0;
  }
}


void IBReallocDependencyV(IBVariables a)
/***************************************************************************
*  Destruction of unused memory
*/
{
  int i;
  for(i=0; i<a->N; i++ )
  {
    a->a[i].dep->a = realloc(a->a[i].dep->a,
                             (a->a[i].dep->N)*sizeof(IBDependencyWord));
  }
}


void IBAddDependencyV(IBVariables a, int globvar, int cstr)
/***************************************************************************
*  To add cstr in the dependency list of variable globvar
*  cstr must be greater than all the constraints already in the dependency
*  list of globvar
*/
{
  int i = (int)floor((double)cstr / (double)sizeof(unsigned long));
  int p = cstr - sizeof(unsigned long)*i;
  IBDependencyV *dep = a->a[globvar].dep;

  if( dep->N==0 )  /* no constraint in dep(globvar) */
  {
    dep->a[0].index = i;
    dep->a[0].val = (1<<p);
    dep->N = 1;
  }
  else     /* index no i already created */
  {
    if( dep->a[dep->N-1].index==i )
    {
      dep->a[dep->N-1].val |= (1<<p);
    }
    else  /* index no i not created yet -> creation in the next place free in dep */
    {
      dep->a[dep->N].index = i;
      dep->a[dep->N].val = (1<<p);
      dep->N ++;
    }
  }
}


void IBWriteDependencyV(IBVariables a)
/***************************************************************************
* 
*/
{
  int i,j, bitpos;
  IBDependencyV *dep;
  unsigned long val;

  for( i=0; i<a->N; i++ )
  {
    printf("Dep(%s) : ",IBNameV(a,i));
    dep = a->a[i].dep;

    for( j=0; j<dep->N; j++ )
    {
      val = dep->a[j].val;
      while( val!=0 )
      {
        bitpos = (int)floor(log((double)val)/0.6931471805599453);
        val &= ~(1<<bitpos);
        printf("C%d ",bitpos+sizeof(unsigned long)*dep->a[j].index);
      }
    }
    printf("\n");
  }
}
