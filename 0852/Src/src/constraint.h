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
 * constraint.h                                                             *
 ****************************************************************************/

#ifndef __constraint_h
#define __constraint_h


#include "variable.h"


/*--- Identifiers of relation symbols ---*/
#define IBRelationEQU 1      /* f == g */
#define IBRelationSUP 2      /* f >= g */
#define IBRelationINF 3      /* f <= g */
#define IBRelationINT 4      /* integer(x) */
#define IBRelationSET 5      /* f in g */



/*--- Identifiers of interval operations symbols ---*/
#define IBOpAddII      0
#define IBOpAddRI      1
#define IBOpSubII      2
#define IBOpSubRI      3
#define IBOpSubIR      4
#define IBOpNegI       5
#define IBOpMulII      6
#define IBOpMulRI      7
#define IBOpMulRnegI   8
#define IBOpMulRposI   9
#define IBOpDivII     10
#define IBOpDivIR     11
#define IBOpDivRI     12
#define IBOpDivIRneg  13
#define IBOpDivIRpos  14
#define IBOpDivRnegI  15
#define IBOpDivRposI  16
#define IBOpSqrI      17
#define IBOpSqrtI     18
#define IBOpPowI      19
#define IBOpExpI      20
#define IBOpLogI      21
#define IBOpMinII     22
#define IBOpMaxII     23
#define IBOpCosI      24
#define IBOpSinI      25
#define IBOpTanI      26
#define IBOpCoshI     27
#define IBOpSinhI     28
#define IBOpTanhI     29
#define IBOpAcosI     30
#define IBOpAsinI     31
#define IBOpAtanI     32
#define IBOpAcoshI    33
#define IBOpAsinhI    34
#define IBOpAtanhI    35

#define IBNbOp        36


static int IBOpDerivable[] = { 1,   /* IBOpAddII */
                               1,   /* IBOpAddRI */
                               1,   /* IBOpSubII */
                               1,   /* IBOpSubRI */
                               1,   /* IBOpSubIR */
                               1,   /* IBOpNegI */
                               1,   /* IBOpMulII */
                               1,   /* IBOpMulRI */
                               1,   /* IBOpMulRnegI */
                               1,   /* IBOpMulRposI */
                               1,   /* IBOpDivII */
                               1,   /* IBOpDivIR */
                               1,   /* IBOpDivRI */
                               1,   /* IBOpDivIRneg */
                               1,   /* IBOpDivIRpos */
                               1,   /* IBOpDivRnegI */
                               1,   /* IBOpDivRposI */
                               1,   /* IBOpSqrI */
                               1,   /* IBOpSqrtI */
                               1,   /* IBOpPowI */
                               1,   /* IBOpExpI */
                               1,   /* IBOpLogI */
                               0,   /* IBOpMinII */
                               0,   /* IBOpMaxII */
                               1,   /* IBOpCosI */
                               1,   /* IBOpSinI */
                               0,   /* IBOpTanI */
                               1,   /* IBOpCoshI */
                               1,   /* IBOpSinhI */
                               1,   /* IBOpTanhI */
                               1,   /* IBOpAcosI */
                               1,   /* IBOpAsinI */
                               1,   /* IBOpAtanI */
                               0,   /* IBOpAcoshI */
                               1,   /* IBOpAsinhI */
                               0 }; /* IBOpAtanhI */


#define IBOpIsDerivable(i) IBOpDerivable[i]


/*--- Trees ---*/

#define IBTNodeVar     0   /* variable */
#define IBTNodeOp      1   /* operation */
#define IBTNodeItv     2   /* interval = constant */
#define IBTNodeUseless 3   /* empty node for unary operations */


union IBValTree
{
    int var[2];  /* var[0]: index of variable in the global array of variables
                    var[1]: index of variable in the local array of variables */
    int op;      /* index of operation symbol */
    IBItv i;     /* interval constant or exponent */
};


struct IBTreeNode
{
    int type;                 /* type of information: var, operation, constant */
    union IBValTree val;      /* value at this node */
    IBItv forward;            /* forward interval evaluation at this node for f */
    IBItv backward;           /* backward interval evaluation at this node for dF */
    IBUnion *u;               /* used for HC4 */
    struct IBTreeNode *left;  /* left subtree */
    struct IBTreeNode *right; /* right subtree */
};
typedef struct IBTreeNode IBTree;


#define IBTleft(f)    ((f)->left)
#define IBTright(f)   ((f)->right)
#define IBTtype(f)    ((f)->type)
#define IBTfwd(f)     ((f)->forward)
#define IBTbwd(f)     ((f)->backward)
#define IBThc4U(f)    ((f)->u)
#define IBTitv(f)     ((f)->val.i)
#define IBTglobvar(f) ((f)->val.var[0])
#define IBTlocvar(f)  ((f)->val.var[1])
#define IBTop(f)      ((f)->val.op)


void    IBTFree       (IBTree *f);
IBTree *IBTCopy       (IBTree *f);
IBTree *IBTNewExp     (int exp);
IBTree *IBTNewItv     (IBItv i);
IBTree *IBTNewUseless ();
IBTree *IBTNewVar     (int globvar, int locvar);
IBTree *IBTNewOp      (int op, IBTree  *l, IBTree *r);
int     IBTiszero     (IBTree *f);
void    IBWriteTop    (FILE *out, int op);
void    IBWriteT      (FILE *out, IBTree *f);
int     IBTDerivable  (IBTree *f);


/*--- Constraints ---*/

typedef struct
{
  int ctr;
  int locvar;
  int asleep;     /* used in filtering algorithms to know if projection has to
                     be used or not */
} IBProjection;   /* projection of ctr over locvar */

#define IBCPctr(p)    p->ctr
#define IBCPvar(p)    p->locvar
#define IBCPasleep(p) p->asleep

struct IBListDepNodes
{
  IBTree *t;
  struct IBListDepNodes *next;
};


/* information on one variable of a given constraint */
typedef struct
{
  int globvar;
  int nbocc;
  IBItv deriv;
  IBTree *function;             /* for example, nested form for box consistency */
  struct IBListDepNodes *list;  /* nodes of the tree depending on this variable */
} IBClocvar;

typedef struct
{
  IBTree *left;
  IBTree *right;
  int rel;             /* constraint is: left rel right */

  IBTree *function;    /* left if right=0, right if left=0, left-right otherwise */
  int relfunc;         /* relation symbol for constraint 'function relfunc 0' */
                       /* relfunc different from rel since it can be inverted */

  IBClocvar *vars;     /* informations concerning variables in the constraint */
  int Nvar;            /* number of variables in the constraint */

  char *name;          /* constraint name */

  IBProjection **pone; /* projections over variables with one occurrence */
  int Npone;           /* the number of such projections */
  IBProjection **pmul; /* projections over variables with multiple occurrences */
  int Npmul;           /* the number of such projections */

  int asleep;          /* flag for propagation in HC4 */

  long numdom;         /* flag to know for which domains the constraint has been
                          evaluated in HC3revise or HC4revise -> then it has not
                          to be evaluated in Gauss-Seidel */

  int IsPartOfModel;   /* 1 if the constraint is part of the model,
                          0 if the constraint is redundant
                          Useful for the test of INNER box */
} IBConstraint;


#define IBCleft(c)          (c)->left
#define IBCright(c)         (c)->right
#define IBCfunc(c)          (c)->function
#define IBCrel(c)           (c)->rel
#define IBCrelfunc(c)       (c)->relfunc
#define IBCname(c)          (c)->name
#define IBCctrAsleep(c)     (c)->asleep
#define IBCctrNumdomGS(c)   (c)->numdom
#define IBCPartOfModel(c)   (c)->IsPartOfModel

#define IBCProjOne(c)       c->pone
#define IBCNbProjOne(c)     c->Npone
#define IBCProjMul(c)       c->pmul
#define IBCNbProjMul(c)     c->Npmul
#define IBCmulprj(c,i)      c->pmul[i]
#define IBConeprj(c,i)      c->pone[i]

#define IBCNbVar(c)         c->Nvar
#define IBCvars(c)          c->vars
#define IBCVglobvar(c,i)    c->vars[i].globvar
#define IBCVnbocc(c,i)      c->vars[i].nbocc
#define IBCVfunc(c,i)       c->vars[i].function
#define IBCVderiv(c,i)      c->vars[i].deriv
#define IBCVdepnodes(c,i)   c->vars[i].list


IBConstraint *IBNewConstraint  (IBTree *l, int rel, IBTree *r, char *s,
                                int indexctr, int isPartOfModel);
void          IBFreeConstraint (IBConstraint *c);
int           IBCGlobvarToLocvar (IBConstraint *c, int globvar);
void          IBCCreateLocVars (IBConstraint *c, IBTree *f, int *nbfree, int size);
void          IBWriteC         (FILE *out, IBConstraint *c);


/*--- array of constraints ---*/
struct IBCstr
{
    IBConstraint **a;  /* array of constraints */
    long Nproj;        /* number of projections in the constraint system */
    int N;             /* number of constraints */
    int Nfree;         /* number of empty places in the array */
};
typedef struct IBCstr *IBConstraints;

#define IBCstrAllocUnit 10

#define IBCNbCtr(c)  c->N
#define IBCNbProj(c) c->Nproj
#define IBCCtr(c,i)  c->a[i]


IBConstraints IBNewConstraints ();
void          IBAddConstraint  (IBConstraints a, IBTree *l, int rel,
                                                 IBTree *r, char *s,
                                                 int isPartOfModel);

void          IBAddIntegerTypeConstraints (IBConstraints a, IBVariables av);

void          IBDecompConstraint (IBConstraints a, IBVariables av, IBTree *l,
                                                   int rel, IBTree *r, char *s,
                                                   int isPartOfModel);

void          IBReallocConstraints (IBConstraints a);
void          IBWriteConstraints   (FILE *out, IBConstraints a);



/* Creation of the dependencies between variables and constraints */
void IBCreateDependencies (IBConstraints ac, IBVariables av);


#endif
