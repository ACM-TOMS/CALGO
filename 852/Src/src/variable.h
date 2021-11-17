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
 * variable.h                                                               *
 ****************************************************************************/

#ifndef __variable_h
#define __variable_h


#include "domain.h"


/*------ definition of variables dependencies

   Every variable depends on a list of constraints it occurs in:
   this list is an array of words, each bit in a word corresponding
   to exactly one constraint:

   constraint number C corresponds to bit number `C - size(one word)*W'
   in word number `W=C div size(one word)'
----------*/

typedef struct
{
    int index;          /* index of word, has to be memorized since
                           the array is sparse */
    unsigned long val;  /* one word */
}
IBDependencyWord;


typedef struct
{
    IBDependencyWord *a;  /* sparse array of all words */
    int N;                /* number of words */
}
IBDependencyV;


#define IBDVnb(dep)      dep->N
#define IBDVdep(dep)     dep->a

#define IBDVindex(dep,i) dep->a[i].index
#define IBDVval(dep,i)   dep->a[i].val




/*------ definition of variable ------*/

typedef struct
{
    char *name;          /* name of the variable */
    IBDependencyV *dep;  /* dependencies */
    int IsInt;           /* 1 if Integer variable, 0 otherwise */
    int status;
    int IsHidden;        /* 1 if it is a fresh variable, 0 otherwise */
    int weight;          /* used for the nested form */
}
IBVar;

/*------ definition of a structure for variables ------*/
struct IBV
{
    IBVar *a;       /* array of variables */
    IBDomains d;    /* variable domains for the parsing */
    int N;          /* number of variables */
    int Nfree;      /* number of empty places in the array */
    int maxname;    /* maximum length of variable name */
};
typedef struct IBV *IBVariables;

#define IBVarAllocUnit    10
#define IBNameV(av,var)   av->a[var].name
#define IBWeightV(av,var) av->a[var].weight
#define IBStatusV(av,var) av->a[var].status
#define IBDepV(av,var)    av->a[var].dep
#define IBDomVars(av)     av->d
#define IBVnb(av)         av->N

#define IBVstatusFresh  1   /* fresh variable */
#define IBVstatusHidden 2   /* hidden variable */
#define IBVstatusUser   3   /* user's variable */

#define IBIsFreshVar(av,var)  (av->a[var].status==IBVstatusFresh)
#define IBIsUserVar(av,var)   (av->a[var].status==IBVstatusUser)
#define IBIsBranchVar(av,var) ((av->a[var].status==IBVstatusUser) || (av->a[var].status==IBVstatusHidden))

#define IBIsIntegerVar(av,var) (av->a[var].IsInt)


/*------ functions ------*/

/* memory allocation for a new structure for variables */
IBVariables IBNewV        ();

/* memory desallocation */
void        IBFreeV       (IBVariables a);

/* Adds in 'a' a new variable which name is 'name'
   Returns the index of the variable in a
   If a variable which name is 'name' is already present in 'a', no variable is created */
int         IBAddV        (IBVariables a, char *name, int status);

/* var becomes an integer variable */
void        IBSetIntV     (IBVariables a, int var);

/* output of 'a' on 'out', used for debugging */
void        IBWriteV      (FILE *out, IBVariables a);

/* output of variables and domains on 'out */
void        IBWriteVdom   (FILE *out, IBDomains d, IBVariables a, int digits, int mode);

/* destruction of unused memory, needed after parsing */
void        IBReallocV    (IBVariables a);

/* returns the number of fresh variables in 'a' */
long        IBNbFreshVar  (IBVariables a);

/* initialisation of the list of dependencies of all variables */
void IBInitDependencyV    (IBVariables a, int n);

/* destruction of unused memory in the lists of dependencies */
void IBReallocDependencyV (IBVariables a);

/* globvar is declared to depend on cstr */
void IBAddDependencyV     (IBVariables a, int globvar, int cstr);

/* to write the lists of dependencies, used for debugging */
void IBWriteDependencyV   (IBVariables a);

#endif
