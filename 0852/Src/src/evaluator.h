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
 * evaluator.h                                                              *
 ****************************************************************************/

#ifndef __evaluator_h
#define __evaluator_h


#include "constraint.h"


/*-- interval evaluation of derivatives in backward mode --*/

/* generic type for an operation of derivation */
typedef void (* IBEvalBwdI)(IBTree *);

void IBBwdAddII      (IBTree *f);
void IBBwdAddRI      (IBTree *f);
void IBBwdSubII      (IBTree *f);
void IBBwdSubIR      (IBTree *f);
void IBBwdSubRI      (IBTree *f);
void IBBwdNegI       (IBTree *f);
void IBBwdMulII      (IBTree *f);
void IBBwdMulRI      (IBTree *f);
void IBBwdDivII      (IBTree *f);
void IBBwdDivIR      (IBTree *f);
void IBBwdDivRI      (IBTree *f);
void IBBwdSqrI       (IBTree *f);
void IBBwdPowI       (IBTree *f);
void IBBwdSqrtI      (IBTree *f);
void IBBwdLogI       (IBTree *f);
void IBBwdExpI       (IBTree *f);
void IBBwdMinII      (IBTree *f);
void IBBwdMaxII      (IBTree *f);
void IBBwdCosI       (IBTree *f);
void IBBwdSinI       (IBTree *f);
void IBBwdTanI       (IBTree *f);
void IBBwdCoshI      (IBTree *f);
void IBBwdSinhI      (IBTree *f);
void IBBwdTanhI      (IBTree *f);
void IBBwdAcosI      (IBTree *f);
void IBBwdAsinI      (IBTree *f);
void IBBwdAtanI      (IBTree *f);
void IBBwdAcoshI     (IBTree *f);
void IBBwdAsinhI     (IBTree *f);
void IBBwdAtanhI     (IBTree *f);


/*-- inversion operations for hull consistency in the union of intervals mode */

/* generic type for an operation of inversion */
typedef int (* IBEvalHC4I)(IBTree *);

int IBHC4AddII    (IBTree *f);
int IBHC4AddRI    (IBTree *f);
int IBHC4SubII    (IBTree *f);
int IBHC4SubIR    (IBTree *f);
int IBHC4SubRI    (IBTree *f);
int IBHC4NegI     (IBTree *f);
int IBHC4MulII    (IBTree *f);
int IBHC4MulRI    (IBTree *f);
int IBHC4DivII    (IBTree *f);
int IBHC4DivIR    (IBTree *f);
int IBHC4DivRI    (IBTree *f);
int IBHC4PowI     (IBTree *f);
int IBHC4SqrtI    (IBTree *f);
int IBHC4ExpI     (IBTree *f);
int IBHC4LogI     (IBTree *f);
int IBHC4MinII    (IBTree *f);
int IBHC4MaxII    (IBTree *f);
int IBHC4CosI     (IBTree *f);
int IBHC4SinI     (IBTree *f);
int IBHC4TanI     (IBTree *f);
int IBHC4CoshI    (IBTree *f);
int IBHC4SinhI    (IBTree *f);
int IBHC4TanhI    (IBTree *f);
int IBHC4AcosI    (IBTree *f);
int IBHC4AsinI    (IBTree *f);
int IBHC4AtanI    (IBTree *f);
int IBHC4AcoshI   (IBTree *f);
int IBHC4AsinhI   (IBTree *f);
int IBHC4AtanhI   (IBTree *f);


/*-- inversion operations for hull consistency in the interval mode */

/* generic type for an operation of inversion for hull consistency */
typedef int (* IBEvalHC3I)(IBTree *);

int IBHC3AddII    (IBTree *f);
int IBHC3AddRI    (IBTree *f);
int IBHC3SubII    (IBTree *f);
int IBHC3SubIR    (IBTree *f);
int IBHC3SubRI    (IBTree *f);
int IBHC3NegI     (IBTree *f);
int IBHC3MulII    (IBTree *f);
int IBHC3MulRI    (IBTree *f);
int IBHC3DivII    (IBTree *f);
int IBHC3DivIR    (IBTree *f);
int IBHC3DivRI    (IBTree *f);
int IBHC3PowI     (IBTree *f);
int IBHC3SqrtI    (IBTree *f);
int IBHC3ExpI     (IBTree *f);
int IBHC3LogI     (IBTree *f);
int IBHC3MinII    (IBTree *f);
int IBHC3MaxII    (IBTree *f);
int IBHC3CosI     (IBTree *f);
int IBHC3SinI     (IBTree *f);
int IBHC3TanI     (IBTree *f);
int IBHC3CoshI    (IBTree *f);
int IBHC3SinhI    (IBTree *f);
int IBHC3TanhI    (IBTree *f);
int IBHC3AcosI    (IBTree *f);
int IBHC3AsinI    (IBTree *f);
int IBHC3AtanI    (IBTree *f);
int IBHC3AcoshI   (IBTree *f);
int IBHC3AsinhI   (IBTree *f);
int IBHC3AtanhI   (IBTree *f);


/*-- to each operation symbol is associated a set of evaluation functions */
struct IBO
{
  IBEvalOpI  eval;    /* interval evaluation operation */
  IBEvalBwdI deriv;   /* interval derivation operation */
  IBEvalHC4I inv;     /* interval inversion operation (HC4) */
  IBEvalHC3I invhc3;  /* interval inversion operation (HC3) */
};
typedef struct IBO *IBOperations;

#define IBTevalOp(a,op)    a[op].eval
#define IBTevalBwd(a,op)   a[op].deriv
#define IBTevalHC4(a,op)   a[op].inv
#define IBTevalHC3(a,op)   a[op].invhc3



/* memory management */
IBOperations IBOperationsInit ();
void         IBOperationsFree (IBOperations a);


/* Interval evaluation of expression f with domains d */
void IBTevalAll (IBTree *f, IBDomains d);


/* Interval evaluation of expression f with domains d
   and domain domvar for variable globvar */
void IBTeval (IBTree *f, IBItv domvar, int globvar, IBDomains d);


/* Interval evaluation of expression f with domains d
   and domain domvar for variable globvar ; only the nodes depending on
   this variable are considered */
void IBTevalOnevar(IBTree *f, IBItv domvar, int globvar,
                   struct IBListDepNodes *list, IBDomains d);


/* returns true if f is variable free */
int  IBTIsConstant   (IBTree *f);

/* Interval evaluation of constant (free variable) expression f */
void IBTevalConstant(IBTree *f);


/* Evaluates with IBTevalConstant and replaces all free variables
   subtrees in f and returns the new tree */
IBTree *IBRemoveConstantSubtrees (IBTree *f);


/* Interval evaluation of partial derivatives of the functional part of c
   the interval evaluation is supposed to be already computed */
void IBCderiv(IBConstraint *c);

#endif
