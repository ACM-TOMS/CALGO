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
 * constraint.c                                                             *
 ****************************************************************************/

#include "constraint.h"
#include "evaluator.h"
#include <string.h>


extern IBVariables variables;     /* global array of constrained variables */


IBTree *IBTNewExp(int exp)
/***************************************************************************
*  To create a tree's node for the interval [exp,exp]
*/
{
  IBTree *f = (IBTree *)malloc(sizeof(IBTree));

  IBSetI(IBTitv(f),exp,exp);
  IBCopyI(IBTfwd(f),IBTitv(f));
  IBThc4U(f) = IBNewU(2);

  IBTleft(f) = IBTright(f) = NULL;
  IBTtype(f) = IBTNodeItv;

  return( f );
}


IBTree *IBTNewItv(IBItv i)
/***************************************************************************
*  To create a tree's node for the interval i
*/
{
  IBTree *f = (IBTree *)malloc(sizeof(IBTree));

  IBCopyI(IBTitv(f),i);
  IBCopyI(IBTfwd(f),i);
  IBThc4U(f) = IBNewU(2);

  IBTleft(f) = IBTright(f) = NULL;
  IBTtype(f) = IBTNodeItv;

  return( f );
}


IBTree *IBTNewUseless()
/***************************************************************************
* 
*/
{
  IBTree *f  = (IBTree *)malloc(sizeof(IBTree));
  IBToLargestI(IBTfwd(f));      /* must be non-empty for interval evaluation */
  IBTtype(f) = IBTNodeUseless;
  IBTleft(f) = IBTright(f) = NULL;
  IBThc4U(f) = IBNewU(2);
  return( f );
}


IBTree *IBTNewVar(int globvar, int locvar)
/***************************************************************************
*  To create a tree's node for the variable globvar/locvar
*/
{
  IBTree *f = (IBTree *)malloc(sizeof(IBTree));

  IBTleft(f)    = IBTright(f) = NULL;
  IBTtype(f)    = IBTNodeVar;
  IBTglobvar(f) = globvar;
  IBTlocvar(f)  = locvar;
  IBThc4U(f)    = IBNewU(2);

  return( f );
}


IBTree *IBTNewOp(int op, IBTree  *l, IBTree *r)
/***************************************************************************
*  To create a tree's node for the interval i
*/
{
  IBTree *f  = (IBTree *)malloc(sizeof(IBTree));

  IBTop(f)    = op;
  IBTleft(f)  = l;
  IBTright(f) = r;
  IBTtype(f)  = IBTNodeOp;
  IBThc4U(f)    = IBNewU(2);

  return( f );
}


void IBTFree(IBTree *f)
/***************************************************************************
*  To desallocate a tree
*/
{
  if( f==NULL ) return;
  IBFreeU(IBThc4U(f));
  IBTFree(IBTleft(f));
  IBTFree(IBTright(f));
  free(f);
}


IBTree *IBTCopy(IBTree *f)
/***************************************************************************
*  To copy a tree and return the copy
*/
{
  IBTree *l, *r;
  if( f==NULL ) return( NULL );

  switch( IBTtype(f) )
  {
    case IBTNodeVar:
         return( IBTNewVar(IBTglobvar(f),IBTlocvar(f)) );
         break;
    case IBTNodeItv:
         return( IBTNewItv(IBTitv(f)) );
         break;
    case IBTNodeUseless:
         return( IBTNewUseless() );
         break;
    case IBTNodeOp:
         l = IBTCopy(IBTleft(f));
         r = IBTCopy(IBTright(f));
         return( IBTNewOp(IBTop(f),l,r) );
         break;
  }
}


int IBTiszero(IBTree *f)
/***************************************************************************
*  Returns 1 if f=0, 0 otherwise
*/
{
  if( f==NULL )
  {
    return 1;
  }
  else if( IBTtype(f)==IBTNodeItv )
  {
    if( IBIsZeroI(IBTitv(f)) ) return( 1 );
    else return( 0 );
  }
  else return( 0 );
}


int IBTDerivable(IBTree *f)
/***************************************************************************
*  Returns 1 if f is derivable
*/
{
  if( f==NULL ) return( 0 );

  switch( IBTtype(f) )
  {
    case IBTNodeVar:
         return( 1 );
         break;
    case IBTNodeItv:
         return( 1 );
         break;
    case IBTNodeUseless:
         return( 1 );
         break;
    case IBTNodeOp:
      if( IBOpIsDerivable(IBTop(f)) )
      {
	return( IBTDerivable(IBTleft(f)) && IBTDerivable(IBTright(f)) );
      }
      else
      {
        return( 0 );
      }
      break;
  }
}


int IBTWriteNeedsNoBracket(IBTree *f)
/***************************************************************************
*  Auxiliary function used to write a function's tree in order to write
*  the minimum of brackets for the infix operators considering their usual
*  priorities
*/
{
  if( IBTtype(f) != IBTNodeOp )  return( 1 );
  else if( IBTop(f)==IBOpSqrtI ) return( 1 );
  else if( IBTop(f)==IBOpExpI )  return( 1 );
  else if( IBTop(f)==IBOpLogI )  return( 1 );
  else if( ((IBTop(f)==IBOpSqrI) ||
           (IBTop(f)==IBOpPowI)) &&
           (IBTWriteNeedsNoBracket(IBTleft(f))) ) return( 1 );
  else return( 0 );
}


void IBWriteT(FILE *out, IBTree *f)
/***************************************************************************
*  To write f on out
*/
{
  if( f==NULL ) return;

  switch( IBTtype(f) )
  {
    case IBTNodeVar:
         fprintf(out,"%s",IBNameV(variables,IBTglobvar(f)));
         break;
    case IBTNodeItv:
         IBWriteI(out,IBTitv(f),9,IBPrintIntervalBounds);
         break;
    case IBTNodeOp:
         if( IBTop(f)==IBOpAddII )
         { IBWriteT(out,IBTleft(f)); fprintf(out,"+");
           IBWriteT(out,IBTright(f));
         }
         else if( IBTop(f)==IBOpAddRI )
         {
           IBWriteT(out,IBTleft(f)); fprintf(out,"+");
           IBWriteT(out,IBTright(f));
         }
         else if( IBTop(f)==IBOpSubII )
         {
           IBWriteT(out,IBTleft(f)); fprintf(out,"-");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpSubRI )
         {
           IBWriteT(out,IBTleft(f));
           fprintf(out,"-");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpSubIR )
         {
           IBWriteT(out,IBTleft(f));
           fprintf(out,"-");
           IBWriteT(out,IBTright(f));
         }
         else if( IBTop(f)==IBOpNegI )
         {
           fprintf(out,"-");
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTleft(f));
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpMulII )
         {
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTleft(f));
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,")");
           fprintf(out,"*");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpMulRI )
         {
           IBWriteT(out,IBTleft(f));
           fprintf(out,"*");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpMulRnegI )
         {
           IBWriteT(out,IBTleft(f));
           fprintf(out,"*");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpMulRposI )
         {
           IBWriteT(out,IBTleft(f));
           fprintf(out,"*");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpDivII )
         {
          if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTleft(f));
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,")");
           fprintf(out,"/");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpDivIR )
         {
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTleft(f));
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,")");
           fprintf(out,"/");
           IBWriteT(out,IBTright(f));
         }
         else if( IBTop(f)==IBOpDivRI )
         {
           IBWriteT(out,IBTleft(f));
           fprintf(out,"/");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpDivIRneg )
         {
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTleft(f));
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,")");
           fprintf(out,"/");
           IBWriteT(out,IBTright(f));
         }
         else if( IBTop(f)==IBOpDivIRpos )
         {
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTleft(f));
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,")");
           fprintf(out,"/");
           IBWriteT(out,IBTright(f));
         }
         else if( IBTop(f)==IBOpDivRnegI )
         {
           IBWriteT(out,IBTleft(f));
           fprintf(out,"/");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpDivRposI )
         {
           IBWriteT(out,IBTleft(f));
           fprintf(out,"/");
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTright(f));
           if( !IBTWriteNeedsNoBracket(IBTright(f)) ) fprintf(out,")");
         }
         else if( IBTop(f)==IBOpSqrI )
         {
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTleft(f));
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,")");
           fprintf(out,"^2");
         }
         else if( IBTop(f)==IBOpSqrtI )
         {
           fprintf(out,"sqrt(");
           IBWriteT(out,IBTleft(f));
           fprintf(out,")");
         }
         else if( IBTop(f)==IBOpExpI )
         {
           fprintf(out,"exp(");
           IBWriteT(out,IBTleft(f));
           fprintf(out,")");
         }
         else if( IBTop(f)==IBOpLogI )
         {
           fprintf(out,"log(");
           IBWriteT(out,IBTleft(f));
           fprintf(out,")");
         }
         else if( IBTop(f)==IBOpPowI )
         {
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,"(");
           IBWriteT(out,IBTleft(f));
           if( !IBTWriteNeedsNoBracket(IBTleft(f)) ) fprintf(out,")");
           fprintf(out,"^%d",(int)IBMinI(IBTitv(IBTright(f))));
         }
       }
}


int IBLocvarInTree(IBTree *f, int locvar)
/***************************************************************************
*  Returns 1 if locvar is in f, 0 otherwise
*/
{
  if( f==NULL ) return( 0 );
  if( IBTtype(f)==IBTNodeVar )
  {
    if( IBTlocvar(f)==locvar ) return( 1 );
    else return( 0 );
  }
  else if( IBTtype(f)==IBTNodeOp )
  {
    if( IBLocvarInTree(IBTleft(f),locvar) ) return( 1 );
    else return( IBLocvarInTree(IBTright(f),locvar) );
  }
  else return( 0 );
}



struct IBListDepNodes *IBCVCreateLDN(IBTree *f, struct IBListDepNodes *l,
                                     int locvar)
/***************************************************************************
*  To create the list of nodes of f depending on variable locvar
*/
{
  struct IBListDepNodes *p;

  if( f==NULL ) return( l );
  if( IBLocvarInTree(f,locvar) )  /* root(f) is an op. or a var. */
  {
    if( IBTtype(f)==IBTNodeVar )  /* necessary locvar */
    {
      p = (struct IBListDepNodes *)malloc(sizeof(struct IBListDepNodes));
      p->t = f;
      p->next = l;
      return( p );
    }
    else if( IBTtype(f)==IBTNodeOp )
    {
      p = (struct IBListDepNodes *)malloc(sizeof(struct IBListDepNodes));
      p->t = f;
      p->next = l;
      l = IBCVCreateLDN(IBTleft(f),p,locvar);
      l = IBCVCreateLDN(IBTright(f),l,locvar);
      return( l );
    }
  }
  else return( l );
}


int IBCGlobvarToLocvar(IBConstraint *c, int globvar)
/***************************************************************************
*  Returns -1 if global variable globvar is not in c
*  the local index otherwise
*/
{
  int j;
  for( j=0; j<IBCNbVar(c); j++ )
  {
    if( IBCVglobvar(c,j)==globvar ) return( j );
  }
  return( -1 );
}


void IBCCreateLocVars(IBConstraint *c, IBTree *f, int *nbfree, int size)
/***************************************************************************
*  To create the informations about the variables in c
*/
{
  int i, j;
  if( IBTtype(f)==IBTNodeVar )
  {
    /* test if variable is already created in c */
    j = -1;
    for( i=0; i<IBCNbVar(c); i++ )
    {
      if( IBTglobvar(f)==IBCVglobvar(c,i) ) j = i;
    }

    if( j>=0 ) /*variable already created */
    {
      IBCVnbocc(c,j)++;
      IBTlocvar(f) = j;  /* variable in f has local index j */
    }
    else
    {
      if( *nbfree==0 ) /* no place in IBCvars(c) */
      {
        IBCvars(c) = realloc(IBCvars(c),
                             (IBCNbVar(c)+size)*sizeof(IBClocvar));
        *nbfree = size;
      }

      IBCVglobvar(c,IBCNbVar(c)) = IBTglobvar(f);
      IBCVnbocc(c,IBCNbVar(c)) = 1;  /* first time variable is found */
      IBTlocvar(f) = IBCNbVar(c); /* variable in f has local index IBCNbVar(c) */
      IBCNbVar(c)++;
      (*nbfree)--;
    }
  }
  else if( IBTtype(f)==IBTNodeOp )
  {
    IBCCreateLocVars(c,IBTleft(f),nbfree,size);
    IBCCreateLocVars(c,IBTright(f),nbfree,size);
  }
}



void IBCVcomputeNbOcc(IBTree *f, int locvar, int *occ)
/***************************************************************************
*  Compute in *occ the number of occurrences of variable locvar in f
*/
{
  if( f==NULL ) return;
  else if( IBTtype(f)==IBTNodeVar )
  {
    if( IBTlocvar(f)==locvar )
    {
      (*occ) ++;
    }
  }
  else if( IBTtype(f)==IBTNodeOp )
  {
    IBCVcomputeNbOcc(IBTleft(f),locvar,occ);
    IBCVcomputeNbOcc(IBTright(f),locvar,occ);
  }
}


IBConstraint *IBNewConstraint(IBTree *l, int rel, IBTree *r, char *s,
                              int indexctr, int isPartOfModel)
/***************************************************************************
*  To allocate a new constraint l rel r s.t. !(l=0 and r=0) 
*/
{
  IBConstraint *c;
  int nf, i;
  IBTree *f;

  c = (IBConstraint *)malloc(sizeof(IBConstraint));

  IBCctrNumdomGS(c) = -1;
  IBCleft(c)  = l;
  IBCright(c) = r;
  IBCrel(c)   = rel;
  IBCname(c)  = (char *)malloc((1+strlen(s))*sizeof(char));
  strcpy(IBCname(c),s);

  IBCrelfunc(c) = rel;

  IBCPartOfModel(c) = isPartOfModel;

  if( IBTiszero(r) ) IBCfunc(c) = l;
  else if( IBTiszero(l) )  /* relation symbol != IBRelationSET */
  {
    IBCfunc(c) = r;
    /* Inversion of relation symbol */
    if( rel==IBRelationSUP )      IBCrelfunc(c) = IBRelationINF;
    else if( rel==IBRelationINF ) IBCrelfunc(c) = IBRelationSUP;
  }
  else IBCfunc(c) = IBTNewOp(IBOpSubII,l,r);


  /* Creation of local structures for the variables in c */
  IBCvars(c) = (IBClocvar *)malloc(sizeof(IBClocvar));
  IBCNbVar(c) = 0;
  nf = 1;
  IBCCreateLocVars(c,IBCfunc(c),&nf,10);
  IBCvars(c) = realloc(IBCvars(c),IBCNbVar(c)*sizeof(IBClocvar));


  for( i=0; i<IBCNbVar(c); i++ )
    IBCVfunc(c,i) = IBCfunc(c);        /* the natural expression */



  /* Backward interval of the root node of the functional part of c
     is initialized to [1,1]
     This is a constant information which will never change */
  IBSetI(IBTbwd(IBCfunc(c)),1.0,1.0);


  IBCProjOne(c) = (IBProjection **)malloc(IBCNbVar(c)*sizeof(IBProjection *));
  IBCNbProjOne(c) = 0;
  IBCProjMul(c) = (IBProjection **)malloc(IBCNbVar(c)*sizeof(IBProjection *));
  IBCNbProjMul(c) = 0;
  for( i=0; i<IBCNbVar(c); i++ )
  {
    /* Creation of projections of c: partition of projections in two parts
       with respect to variables having one or multiple occurrences in c */
    if( IBCVnbocc(c,i)==1 )  /* one occurrence of local variable i */
    {
      IBConeprj(c,IBCNbProjOne(c)) = (IBProjection *)malloc(sizeof(IBProjection));
      IBCPctr(IBConeprj(c,IBCNbProjOne(c)))    = indexctr;
      IBCPvar(IBConeprj(c,IBCNbProjOne(c)))    = i;
      IBCPasleep(IBConeprj(c,IBCNbProjOne(c))) = 1;
      IBCNbProjOne(c)++;
     }
    else                     /* multiple occurrences of local variable i */
    {
      IBCmulprj(c,IBCNbProjMul(c)) = (IBProjection *)malloc(sizeof(IBProjection));
      IBCPctr(IBCmulprj(c,IBCNbProjMul(c)))    = indexctr;
      IBCPvar(IBCmulprj(c,IBCNbProjMul(c)))    = i;
      IBCPasleep(IBCmulprj(c,IBCNbProjMul(c))) = 1;
      IBCNbProjMul(c)++;
    }

    /* Creation of the list of nodes of IBCVfunc(c,i) depending on variable i */
    IBCVdepnodes(c,i) = IBCVCreateLDN(IBCVfunc(c,i),NULL,i);
  }

  IBCProjMul(c) = realloc(IBCProjMul(c),IBCNbProjMul(c)*sizeof(IBProjection *));
  IBCProjOne(c) = realloc(IBCProjOne(c),IBCNbProjOne(c)*sizeof(IBProjection *));

  return( c );
}


void IBFreeConstraint(IBConstraint *c)
/***************************************************************************
*  To desallocate a constraint
*/
{
  int i;
  struct IBListDepNodes *p, *q;

  /* to desallocate the nested forms
  for( i=0; i<IBCNbVar(c); i++ )
  {
    if( IBCVfunc(c,i)!=IBCfunc(c) )  IBTFree(IBCVfunc(c,i));
  }
  */

  if( (IBCfunc(c)==IBCleft(c)) ||(IBCfunc(c)==IBCright(c))  )
  {
    IBTFree(IBCleft(c));
    IBTFree(IBCright(c));
  }
  else
  {
    IBTFree(IBCleft(c));
    IBTFree(IBCright(c));
    free(IBCfunc(c));   /* to desallocate root node with minus operation */
  }
  for( i=0; i<IBCNbProjOne(c); i++ )
  {
    free(IBCProjOne(c)[i]);
  }
  for( i=0; i<IBCNbProjMul(c); i++ )
  {
    free(IBCProjMul(c)[i]);
  }
  for( i=0; i<IBCNbVar(c); i++ )
  {
    p = IBCVdepnodes(c,i);
    while( p!=NULL )
    {
      q = p;
      p = p->next;
      free(q);
    }
  }
  free(IBCProjOne(c));
  free(IBCProjMul(c));
  free(IBCname(c));
  free(IBCvars(c));
  free(c);
}


void IBWriteC(FILE *out, IBConstraint *c)
/***************************************************************************
*  To write c on out
*/
{
  if( strcmp(IBCname(c),"")!=0 )
        fprintf(out,"%s: ",IBCname(c));
  IBWriteT(out,IBCleft(c));
  if( IBCrel(c)==IBRelationEQU ) fprintf(out," = ");
  else if( IBCrel(c)==IBRelationINF ) fprintf(out," <= ");
  else if( IBCrel(c)==IBRelationSUP ) fprintf(out," >= ");
  else fprintf(out," in ");
  IBWriteT(out,IBCright(c));
}



IBConstraints IBNewConstraints()
/***************************************************************************
*  Allocation of an array of constraints
*/
{
  IBConstraints a;
  a = (struct IBCstr *)malloc(sizeof(struct IBCstr));

  a->a     = (IBConstraint **)malloc(IBCstrAllocUnit*sizeof(IBConstraint *));
  a->Nproj = 0;
  a->N     = 0;
  a->Nfree = IBCstrAllocUnit;

  return( a );
}


void IBFreeConstraints(IBConstraints a)
/***************************************************************************
*  Desallocation of an array of constraints
*/
{
  int i;

  for( i=0; i<a->N; i++ )
  {
    IBFreeConstraint(IBCCtr(a,i));
  }
  free(a->a);
  free(a);
}


void IBAddConstraint(IBConstraints a, IBTree *l, int rel, IBTree *r,
                     char *s, int isPartOfModel)
/***************************************************************************
*  To add constraint `l rel r' with name s in a
*/
{
  if( a->Nfree==0 )         /* a is full */
  {
    a->a = realloc(a->a,(a->N + IBCstrAllocUnit)*sizeof(IBConstraint *));
    a->Nfree = IBCstrAllocUnit;
  }

  a->a[a->N] = IBNewConstraint(l,rel,r,s,a->N,isPartOfModel);
  a->Nproj += a->a[a->N]->Nvar;
  a->N++;
  a->Nfree--;
}


void IBAddIntegerTypeConstraints(IBConstraints a, IBVariables av)
/***************************************************************************
*  To add constraints integer(x) for all the integer variables in av
*/
{
  int i;
  for( i=0; i<IBVnb(av); ++i )
  {
    if( IBIsIntegerVar(av,i) )
    {
      IBAddConstraint(a,IBTNewVar(i,-1),IBRelationINT,IBTNewUseless(),"",1);
    }
  }
}


int IBNbOccGlobVar(IBTree *f, int globvar)
/***************************************************************************
*  Returns the number of occurrences of variable globvar in f
*/
{
  if( f==NULL ) return( 0 );
  if( IBTtype(f)==IBTNodeVar )
  {
    if( IBTglobvar(f)==globvar ) return( 1 );
    else return( 0 );
  }
  else if( IBTtype(f)==IBTNodeOp )
  {
    return( IBNbOccGlobVar(IBTleft(f),globvar) +
            IBNbOccGlobVar(IBTright(f),globvar) );
  }
  else return( 0 );
}


int IBNbOperations(IBTree *f)
/***************************************************************************
*  Returns the number of operations in f
*/
{
  if( f==NULL ) return( 0 );
  if( IBTtype(f)==IBTNodeOp )
  {
    return( 1 + IBNbOperations(IBTleft(f)) + IBNbOperations(IBTright(f)) );
  }
  else return( 0 );
}


int IBAddNewFreshVar(IBVariables av)
/***************************************************************************
*  Create a fresh variable in av generated by the decomposition process
*  It is named _Vn where n is the value of NumFreshVar
*/
{
  char name[20];
  int n;
  static long NumFreshVar = 0;

  sprintf(name,"_V%ld",NumFreshVar);
  NumFreshVar++;

  n = IBAddV(av,name,IBVstatusFresh);

  /*
  IBMinI(IBDomV(IBDomVars(av),n)) = IBMinDouble;
  IBMaxI(IBDomV(IBDomVars(av),n)) = IBMaxDouble;
  */

  IBSetToRealDomain(IBDomV(IBDomVars(av),n));

  return( n );
}



void IBDecompRemoveMultipleVar(IBVariables av, IBTree *f,
                               IBTree *one, IBTree *two)
/***************************************************************************
*  To remove the multiple occurrences of variables in f
*  A multiple occurrence is an occurrence of variable x such that
*  x appears more than once in (one,two)
*  A multiple occurrence is replaced by a new fresh variable added in av
*/
{
  int n, globvar;

  if( f==NULL ) return;
  if( IBTtype(f)==IBTNodeVar )
  {
    globvar = IBTglobvar(f);
    n = IBNbOccGlobVar(one,globvar) +  /* number of occurrences */
        IBNbOccGlobVar(two,globvar);

    if( n > 1 )  /* then remove this occurrence */
    {
   /* addition of a new fresh variable in av; IBTglobvar(f) is modified */
      IBTglobvar(f) = IBAddNewFreshVar(av);

   /* this new variable takes the same domain than globvar */
      IBCopyI(IBDomV(IBDomVars(variables),IBTglobvar(f)),
              IBDomV(IBDomVars(variables),globvar));
    }
  }
  else if( IBTtype(f)==IBTNodeOp )
  {
    IBDecompRemoveMultipleVar(av,IBTleft(f),one,two);
    IBDecompRemoveMultipleVar(av,IBTright(f),one,two);
  }
}


void IBDecompTerm(IBConstraints a, IBVariables av, IBTree *f, int globvar,
                  int isPartOfModel)
/***************************************************************************
*   To decompose constraint `globvar = f'
*/
{
  IBTree *new1, *new2, *new3, *new4;
  char namedecomp[] = "Decomp";
  int nop1, nop2, newvar, newvar1, newvar2;

  new1 = IBTNewVar(globvar,-1);

  if( IBNbOperations(f)<=1 )
  {
    new2 = IBTCopy(f);
    IBAddConstraint(a,new1,IBRelationEQU,new2,namedecomp,isPartOfModel); /* `globvar = f' */
  }
  else  /* decomposition */
  {
    nop1 = IBNbOperations(IBTleft(f));
    nop2 = IBNbOperations(IBTright(f));

    if( nop1==0 )  /* nop2>=2 then IBTright(f) to be decomposed */
    {
      new2 = IBTCopy(IBTleft(f));
      newvar = IBAddNewFreshVar(av);
      new3 = IBTNewVar(newvar,-1);
      new4 = IBTNewOp(IBTop(f),new2,new3);
      IBAddConstraint(a,new1,IBRelationEQU,new4,namedecomp,isPartOfModel);
      IBDecompTerm(a,av,IBTright(f),newvar,isPartOfModel);
    }
    else if( nop2==0 )  /* nop1>=2 then IBTleft(f) to be decomposed */
    {
      new2 = IBTCopy(IBTright(f));
      newvar = IBAddNewFreshVar(av);
      new3 = IBTNewVar(newvar,-1);
      new4 = IBTNewOp(IBTop(f),new3,new2);
      IBAddConstraint(a,new1,IBRelationEQU,new4,namedecomp,isPartOfModel);
      IBDecompTerm(a,av,IBTleft(f),newvar,isPartOfModel);
    }
    else  /* IBTleft(f) and IBTright(f) to be decomposed */
    {
      newvar1 = IBAddNewFreshVar(av);
      new2 = IBTNewVar(newvar1,-1);
      newvar2 = IBAddNewFreshVar(av);
      new3 = IBTNewVar(newvar2,-1);
      new4 = IBTNewOp(IBTop(f),new2,new3);
      IBAddConstraint(a,new1,IBRelationEQU,new4,namedecomp,isPartOfModel);
      IBDecompTerm(a,av,IBTleft(f),newvar1,isPartOfModel);
      IBDecompTerm(a,av,IBTright(f),newvar2,isPartOfModel);
    }
  }
}



void IBDecompConstraint(IBConstraints a, IBVariables av, IBTree *l, int rel,
                                                         IBTree *r, char *s,
                                                         int isPartOfModel)
/***************************************************************************
*  To decompose constraint `l rel r' with name s
*  The resulting primitive constraints are added in a
*/
{
  long nop1, nop2;
  IBTree *new1, *new2;
  char namedecomp[] = "Decomp";
  int newvar1, newvar2;


  /* First procedure: replacing all multiple occurrences of variables
     by new fresh variables */  
  /*  IBDecompRemoveMultipleVar(av,l,l,r);
      IBDecompRemoveMultipleVar(av,r,l,r);
  */


  /* Second procedure: decomposition of constraint `l rel r' */
  nop1 = IBNbOperations(l);
  nop2 = IBNbOperations(r);

  if( nop1+nop2<=1 )  /* constraint `l rel r' is already primitive */
  {
    IBAddConstraint(a,l,rel,r,s,isPartOfModel);
  }
  else if( nop1==0 )
  {
    newvar2 = IBAddNewFreshVar(av);
    new2 = IBTNewVar(newvar2,-1);
    IBAddConstraint(a,l,rel,new2,namedecomp,isPartOfModel);             /* `l rel new2' */
    IBDecompTerm(a,av,r,newvar2,isPartOfModel);   /* Decomposition of r with `new2 = r' */
    IBTFree(r);
  }
  else if( nop2==0 )
  {
    newvar1 = IBAddNewFreshVar(av);
    new1 = IBTNewVar(newvar1,-1);
    IBAddConstraint(a,new1,rel,r,namedecomp,isPartOfModel);   /* `l rel new1' */
    IBDecompTerm(a,av,l,newvar1,isPartOfModel);   /* Decomposition of r with `l = new1' */
    IBTFree(l);
   }
  else  /* general case */
  {
    newvar1 = IBAddNewFreshVar(av);
    newvar2 = IBAddNewFreshVar(av);
    new1 = IBTNewVar(newvar1,-1);
    new2 = IBTNewVar(newvar2,-1);
    IBAddConstraint(a,new1,rel,new2,namedecomp,isPartOfModel);   /* `new2 rel new1' */
    IBDecompTerm(a,av,l,newvar1,isPartOfModel);   /* Decomposition of r with `l = new1' */
    IBDecompTerm(a,av,r,newvar2,isPartOfModel);   /* Decomposition of r with `new2 = r' */
    IBTFree(l);
    IBTFree(r);
  }
}


void IBReallocConstraints(IBConstraints a)
/***************************************************************************
*  Reallocation of structures; needed after the parsing
*/
{
  a->a = realloc(a->a,(a->N)*sizeof(IBConstraint *));
}


void IBWriteConstraints(FILE *out, IBConstraints a)
/***************************************************************************
*  To write on out all the constraints in a
*/
{
  int i;

  for( i=0; i<IBCNbCtr(a); i++ )
  {
    IBWriteC(out,IBCCtr(a,i));
    printf("\n");
  }
}


void IBCreateDependencies(IBConstraints ac, IBVariables av)
/***************************************************************************
*  To create the dependencies between constraints and variables
*/
{
  int ctr, var;

  IBInitDependencyV(av,IBCNbCtr(ac));

  for( ctr=0; ctr<IBCNbCtr(ac); ctr++ )
  {
    for( var=0; var<IBCNbVar(IBCCtr(ac,ctr)); var++ )
    {
      IBAddDependencyV(av,IBCVglobvar(IBCCtr(ac,ctr),var),ctr);
    }
  }

  IBReallocDependencyV(av);
}
