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
 * parser.y                                                                 *
 *     -> needs GNU bison parser generator                                  *
 ****************************************************************************/

%{

#include "search.h"

#define YYINITDEPTH 180000  /* SIZE OF THE PARSER STACK */
#define YYMAXDEPTH  180000
extern char *yytext;
extern FILE *yyin;


/*--- global variables ---*/
extern IBVariables variables;      /* global array of constrained variables */
extern IBConstants constants;      /* global array of constants */
extern IBConstraints constraints;  /* global array of constraints */

extern double             IBPragmaPrecision;
extern int                IBPragmaBisection;
extern int                IBPragmaNumberBisection; 
extern long               IBPragmaMaxTime;
extern IBBisectArity      IBsplit;
extern int                IBargBisectNo;
extern int                IBPragmaHullMode; 
extern int                IBPragmaIntervalDigits;
extern int                IBPragmaStyleInterval;
extern unsigned long      IBPragmaMaxSolution;
extern IBLocalPropagation IBfilter2B;
extern int                IBcompute3B;
extern int                IBcomputeWeak3B;
extern double             IBPragmaPrecisionShrink;
extern double             IBPragmaPrecision3B; 
extern double             IBPragmaImprovement; 
extern int                IBPragmaSubpaving;


/*--- variables used for parsing ---*/
IBItv itv;
IBInterval *itv2;
double x1, x2;
char s[40];
IBTree *f1, *f2;
int i,
    newVar,
    the_index,
    IBIndexArrayLeft,
    IBIndexArrayRight;
long lt;
%}


%union
{
   int  nvar;
   char str_num[40];
   char str_unum[40];
   char str_sign[2];
   char str_integ[40];
   char str_float[40];
   char str_var[20];
   int  n, vtype, bbound;
   double lbound, rbound;
   void *fun;
}

%type <str_var>   ConstraintName
%type <str_var>   ConstraintIdent
%type <n>         Relation

%type <fun>       Expr
%type <fun>       ExprMul
%type <fun>       ExprExp
%type <fun>       ExprUnit
%type <fun>       Ident
%type <n>         Exposant

%type <str_var>   Variable
%type <n>         VariableArray
%type <n>         IsNewData
%type <vtype>     VarType
%type <n>         PrevSuccNumber
%type <lbound>    ExprLeftBound
%type <rbound>    ExprRightBound
%type <bbound>    BracketBound

%type <str_var>   ConstName
%type <str_var>   IdentName
%type <str_var>   IdentArray

%type <str_num>   Number
%type <str_integ> Integer
%type <str_float> Float


%token PRAGMACONSTANTS
%token PRAGMADOMAINS
%token PRAGMACONSTRAINTS

%token PRAGMABISECTION
%token BISECTIONCHOICE
%token BISECTIONPARTS
%token NOBISECTION
%token BISECTIONCHOICERR
%token BISECTIONCHOICELF
%token BISECTIONCHOICEMN
%token BISECTIONNUMBER
%token BISECTIONSUBPAVING
%token BISECTIONPOINTS

%token PRAGMAMAXTIME

%token MODE

%token PRAGMAOUTPUT
%token OUTPUTHULLMODE
%token OUTPUTUNIONMODE
%token OUTPUTDIGITS
%token OUTPUTSTYLE
%token OUTPUTBOUNDSTYLE
%token OUTPUTMIDPOINTSTYLE
%token OUTPUTSOLUTION
%token OUTPUTALLSOLUTION

%token PRAGMACONSISTENCY
%token LOCALCONSISTENCY
%token STRONGCONSISTENCY

%token CONSISTENCYBC3
%token CONSISTENCYBC3Newton
%token CONSISTENCYBC4
%token CONSISTENCYBC5
%token CONSISTENCYHC3
%token CONSISTENCYHC4
%token CONSISTENCYHC4I
%token CONSISTENCYHC4Newton
%token PRECISION2B
%token IMPROVEMENT2B

%token CONSISTENCY3B
%token CONSISTENCYWEAK3B
%token PRECISION3B

%token PRECISION

%token INDOM
%token COMMA
%token TWOPOINTS
%token SCOLON
%token COLON
%token UNDERSCORE
%token LSBR
%token RSBR
%token LBR
%token RBR
%token IDENT
%token INTEGER
%token FLOAT
%token REALPOSINFINITY
%token REALNEGINFINITY
%token ADD
%token SUB
%token MUL
%token DIV
%token POW
%token SQRT
%token LOG
%token EXP
%token MINIMUM
%token MAXIMUM
%token COS
%token SIN
%token TAN
%token COSH
%token SINH
%token TANH
%token ACOS
%token ASIN
%token ATAN
%token ACOSH
%token ASINH
%token ATANH
%token IDENT
%token PREV
%token SUCC
%token INTEGERTYPE
%token REALTYPE

%token NEWDATA

%token EQU
%token INF
%token SUP


%expect 6

%start First

%%
/****************************************************************************
 *                        PROGRAM = LIST OF PRAGMAS                         *
 ****************************************************************************/
 First: Pragma SCOLON NextPragma
      ;
 NextPragma: 
           | First 
           ;
 Pragma: PRAGMADOMAINS Domains
       | PRAGMACONSTRAINTS Constraints
       | PRAGMACONSTANTS Constants
       | PRAGMABISECTION Bisection
       | PRAGMAMAXTIME MaxTime
       | PRAGMAOUTPUT Output
       | PRAGMACONSISTENCY Consistency
       ;


/****************************************************************************
 *                        MAXIMUM COMPUTATION TIME                          *
 ****************************************************************************/
 MaxTime: EQU Number
          {
            lt = strtoul($2,NULL,10);
            if( lt>0 )
            {
              IBPragmaMaxTime = lt;
	    }
            else
	    {
	      yyerror("invalid time.");
	      YYABORT;                    
	    }
	  }
        ;


/****************************************************************************
 *                                 OUTPUT                                   *
 ****************************************************************************/
 Output: OutputArgument NextOutput
       ;

 NextOutput:
           | COMMA Output
           ;

 OutputArgument:
               | MODE EQU OutputTypeMode
               | OUTPUTDIGITS EQU Integer
                 {
                    i = IBStrToInt($3);
		    if( i>=0 )
		    {
                      IBPragmaIntervalDigits = i;
                    }
                    else
		    {
	              yyerror("invalid number of digits.");
	              YYABORT;
		    }
		 }
               | OUTPUTSTYLE EQU OutputTypeStyle
	       ;

 OutputTypeMode: OUTPUTHULLMODE
                 {
                   IBPragmaHullMode = 1;
                 }
               | OUTPUTUNIONMODE
                 {
                   IBPragmaHullMode = 0;
                 }
               ;

 OutputTypeStyle: OUTPUTBOUNDSTYLE     { IBPragmaStyleInterval = IBPrintIntervalBounds; }
                | OUTPUTMIDPOINTSTYLE  { IBPragmaStyleInterval = IBPrintIntervalMidError; }
                ;


/****************************************************************************
 *                               CONSISTENCY                                *
 ****************************************************************************/
 Consistency: ConsistencyArgument NextConsistency
            ;

 NextConsistency:
                | COMMA Consistency
                ;

 ConsistencyArgument:
                   | LOCALCONSISTENCY EQU ConsistencyLocalType
                   | STRONGCONSISTENCY EQU ConsistencyStrongType
                   | PRECISION2B EQU Number
                     {
                       x2 = IBStrToDouble($3);
                       if( x2>=0 )
		       {
		         IBPragmaPrecisionShrink = x2;  /* precision of local (box) consistency */
		       }
                       else
		       {
	                 yyerror("invalid precision.");
	                 YYABORT;     
		       }
		     }
                   | PRECISION3B EQU Number
                     {
                       x2 = IBStrToDouble($3);
                       if( x2>=0 )
		       {
		         IBPragmaPrecision3B = x2;  /* precision of (weak) 3B consistency */
		       }
                       else
		       {
	                 yyerror("invalid precision.");
	                 YYABORT;     
		       }
		     }
                   | IMPROVEMENT2B EQU Number
                     {
                       x2 = IBStrToDouble($3);
                       if( (x2>=0) && (x2<=100.0) )
		       {
                         if( x2==0.0 )
			 {
		           IBPragmaImprovement = 1.0;
			 }
                         else
                         {
		           IBPragmaImprovement = 1 - (x2/100); /* improvement factor of 2B consistency */
			 }
		       }
                       else
		       {
	                 yyerror("invalid improvement factor.");
	                 YYABORT;     
		       }
		     }
		   ;

 ConsistencyLocalType: CONSISTENCYBC3
                       {
			 IBfilter2B = IBF2Bbc3;
		       }
                     | CONSISTENCYBC3Newton
                       {
			 IBfilter2B = IBF2Bbc3Newton;
		       }
                     | CONSISTENCYBC4
                       {
			 IBfilter2B = IBF2Bbc4;
		       }
                     | CONSISTENCYBC5
                       {
			 IBfilter2B = IBF2Bbc5;
		       }
                     | CONSISTENCYHC3
                       {
			 IBfilter2B = IBF2Bhc3;
		       }
                     | CONSISTENCYHC4
                       {
			 IBfilter2B = IBF2Bhc4;
		       }
                     | CONSISTENCYHC4I
                       {
			 IBfilter2B = IBF2Bhc4I;
		       }
                     | CONSISTENCYHC4Newton
                       {
			 IBfilter2B = IBF2Bhc4Newton;
		       }
                     ;

ConsistencyStrongType: CONSISTENCY3B
                       {
                         IBcompute3B = 1;
		       }
                     | CONSISTENCYWEAK3B
                       {
                         IBcomputeWeak3B = 1;
		       }
                     ;
 
/****************************************************************************
 *                                BISECTION                                 *
 ****************************************************************************/
 Bisection: BisectionArgument NextBisection
          ;

 NextBisection:
              | COMMA Bisection;

 BisectionArgument:
                  | NOBISECTION
                    {
                      IBPragmaNumberBisection = 0;
                      IBargBisectNo = 1;
		    }

                  | PRECISION EQU Number
                    {
                      x2 = IBStrToDouble($3);
                      if( x2>=0 )
		      {
                        IBPragmaPrecision = x2;
		      }
                      else
		      {
	                yyerror("invalid precision.");
	                YYABORT;     
		      }
		    }
                  | BISECTIONCHOICE EQU BisectionTypeChoice
                  | BISECTIONPARTS EQU Integer
                    {
                      i = IBStrToInt($3);
		      if( i==2 )
		      {
                        IBsplit = IBBsplit2;
                        IBPragmaNumberBisection = 2;
		      }
		      else if( i==3 )
		      {
                        IBsplit = IBBsplit3;
                        IBPragmaNumberBisection = 3;
		      }
		    }
                  | BISECTIONNUMBER EQU BisectionNumber
                  | MODE EQU BisectionMode
		  ;

 BisectionTypeChoice: BISECTIONCHOICERR { IBPragmaBisection = IBBisectRoundRobin; }
                    | BISECTIONCHOICELF { IBPragmaBisection = IBBisectLargestFirst; }
                    | BISECTIONCHOICEMN { IBPragmaBisection = IBBisectMaxNarrow; }
                    ;

 BisectionNumber: Integer
                  {
                    IBPragmaMaxSolution = strtoul($1,NULL,10);
 
	 	    if( IBPragmaMaxSolution<1 )
		    {
	              yyerror("invalid number of solutions.");
	              YYABORT;
		    }
		  }
                | REALPOSINFINITY
                  {
		    IBPragmaMaxSolution = ~0;
		  }
                ;

 BisectionMode: BISECTIONPOINTS    { IBPragmaSubpaving = 0; }
              | BISECTIONSUBPAVING { IBPragmaSubpaving = 1; }
              ;


/****************************************************************************
 *                         DEFINITION OF CONSTANTS                          *
 ****************************************************************************/
 Constants: Constant NextConstants
          ;

 NextConstants:
              | COMMA Constants
              ;

 Constant: ConstName EQU Expr
           {
             if( !IBAddConstant(constants,$1,$3) )
	     {
	       yyerror("a constant expression evaluates in the empty interval.");
	       YYABORT;                    
	     }
	   }
         | ConstName EQU LSBR Expr DomainSep Expr RSBR
           {
             f1 = (IBTree *)$4;
             IBTevalConstant(f1);
	     if( IBEmptyI(IBTfwd(f1)) )
	     {
	       yyerror("a constant expression evaluates in the empty interval.");
	       YYABORT;
	     }

             f2 = (IBTree *)$6;
             IBTevalConstant(f2);
	     if( IBEmptyI(IBTfwd(f2)) )
	     {
	       yyerror("a constant expression evaluates in the empty interval.");
	       YYABORT;
	     }

             IBSetI(itv,IBMinI(IBTfwd(f1)),IBMaxI(IBTfwd(f2)));
             IBTFree(f1);
             IBTFree(f2);
             IBAddConstant(constants,$1,IBTNewItv(itv));
	   }
         ;

 ConstName: IDENT { strcpy($$,yytext); }
          ;


/****************************************************************************
 *                        DEFINITION OF CONSTRAINTS                         *
 ****************************************************************************/
 Constraints: Constraint NextConstraints
            ;

 NextConstraints:
                | COMMA Constraints
                ;

 Constraint: IsNewData ConstraintName Expr Relation Expr
             {
               if( !((IBTiszero($3)) && (IBTiszero($5))) )
	       {
		  if( (f1=IBRemoveConstantSubtrees($3))==NULL )
		  {
	            yyerror("a constant expression evaluates in the empty interval.");
	            YYABORT;
		  }
                  if( (f2=IBRemoveConstantSubtrees($5))==NULL )
		  {
	            yyerror("a constant expression evaluates in the empty interval.");
	            YYABORT;
		  }

                  if( IBfilter2B==IBF2Bhc3 )
		  {
                    IBDecompConstraint(constraints,variables,f1,$4,f2,$2,!($1));
		  }
                  else
		  {
                    IBAddConstraint(constraints,f1,$4,f2,$2,!($1));
		  }
	       }
             }
           ;

 ConstraintName: { strcpy($$,""); }
               | ConstraintIdent COLON { strcpy($$,$1); }
               ;

 ConstraintIdent: UNDERSCORE IDENT { strcpy($$,yytext); }
                ;

 Relation: EQU    { $$ = IBRelationEQU; }
         | INF    { $$ = IBRelationINF; }
         | SUP    { $$ = IBRelationSUP; }
         | INDOM  { $$ = IBRelationSET; }
         ;


/****************************************************************************
 *                 DEFINITION OF VARIABLES AND DOMAINS                      *
 ****************************************************************************/
 Domains: Domain NextDomains
        ;

 NextDomains:
            | COMMA Domains
            ;

 Domain: VarType IsNewData Variable VariableArray
         INDOM
         BracketBound ExprLeftBound DomainSep ExprRightBound BracketBound
         {
	   if( ($7>$9) || (($7==$9) && (($6==2) || ($10==1))) ) {
	     yyerror("the interval is empty.");
	     YYABORT;
	   }

           if( $2 )
  	   {
	     newVar = IBVstatusHidden;
	   }
           else
	   {
	     newVar = IBVstatusUser;
	   }

           if( $4 )  /* array of variables */
	   {
             if( (IBIndexArrayLeft>=0) && (IBIndexArrayRight>=IBIndexArrayLeft) )
	     {
               for( i=IBIndexArrayLeft; i<=IBIndexArrayRight; i++ )
	       {
                 sprintf(s,"%s[%d]",$3,i);
                 the_index = IBAddV(variables,s,newVar);
                 IBSetD(IBDomVars(variables),the_index,$7,$9);
                 if( $1==3 )  /* integer variable */
	         {
                   IBSetIntV(variables,the_index);
	         }
	       }
	     }
             else
	     {
               yyerror("wrong index of variables array.");
	       YYABORT;
	     }
	   }
           else      /* one variable */
	   {
             the_index = IBAddV(variables,$3,newVar);
             IBSetD(IBDomVars(variables),the_index,$7,$9);
             if( $1==3 )  /* integer variable */
	     {
               IBSetIntV(variables,the_index);
	     }
	   }
         }
       ;

 VariableArray: { $$ = 0; }
              | LSBR Integer DomainSep Integer RSBR
                {
                  IBIndexArrayLeft  = IBStrToInt($2);
                  IBIndexArrayRight = IBStrToInt($4);
                  $$ = 1;
                }
              ;

 DomainSep: COMMA
          | TWOPOINTS
          ;

 VarType:             { $$ = 1; }
        | REALTYPE    { $$ = 2; }
        | INTEGERTYPE { $$ = 3; }
	;

 BracketBound: LSBR { $$ = 1; }   /* [ */
             | RSBR { $$ = 2; }   /* ] */
             ;

 ExprLeftBound: PrevSuccNumber Expr
                {
                  f1 = (IBTree *)$2;
                  IBTevalConstant(f1);
		  if( IBEmptyI(IBTfwd(f1)) )
		  {
	            yyerror("the left bound is not defined by a constant expression.");
	            YYABORT;
		  }
                  else
		  {
                    x1 = IBMinI(IBTfwd(f1));

	            if ($1==1) {
                      $$ = IBPrevDouble(x1);
	            }
	            else if ($1==2) {
                      $$ = IBNextDouble(x1);
	            }
                    else
		    {
		      $$ = x1;
		    }
                    IBTFree(f1);
		  }
                }
              ;

 ExprRightBound: PrevSuccNumber Expr
                 {
                   f1 = (IBTree *)$2;
                   IBTevalConstant(f1);
		   if( IBEmptyI(IBTfwd(f1)) )
		   {
	             yyerror("the right bound is not defined by a constant expression.");
	             YYABORT;
	 	   }
                   else
		   {
                     x1 = IBMaxI(IBTfwd(f1));
	             if ($1==1) {
                       $$ = IBPrevDouble(x1);
	             }
	             else if ($1==2) {
                       $$ = IBNextDouble(x1);
	             }
                     else
		     {
		       $$ = x1;
		     }
                     IBTFree(f1);
		   }
                 }
              ;


/****************************************************************************
 *                                EXPRESSIONS                               *
 ****************************************************************************/
 Expr: Expr ADD ExprMul
       {
         if( IBTtype((IBTree *)$1)==IBTNodeItv )
	 {
           if( IBIsDoubleI(IBTitv((IBTree *)$1)) )
             $$ = IBTNewOp(IBOpAddRI,(IBTree *)$1,(IBTree *)$3);
           else
             $$ = IBTNewOp(IBOpAddII,(IBTree *)$1,(IBTree *)$3);
	 }
         else if( IBTtype((IBTree *)$3)==IBTNodeItv )
	 {
           if( IBIsDoubleI(IBTitv((IBTree *)$3)) )
             $$ = IBTNewOp(IBOpAddRI,(IBTree *)$3,(IBTree *)$1);
           else
             $$ = IBTNewOp(IBOpAddII,(IBTree *)$3,(IBTree *)$1);
	 }
         else
           $$ = IBTNewOp(IBOpAddII,(IBTree *)$1,(IBTree *)$3);
       }

     | Expr SUB ExprMul
       {
         if( IBTtype((IBTree *)$1)==IBTNodeItv )
	 {
           if( IBIsDoubleI(IBTitv((IBTree *)$1)) )
             $$ = IBTNewOp(IBOpSubRI,(IBTree *)$1,(IBTree *)$3);
           else
             $$ = IBTNewOp(IBOpSubII,(IBTree *)$1,(IBTree *)$3);
	 }
         else if( IBTtype((IBTree *)$3)==IBTNodeItv )
	 {
           if( IBIsDoubleI(IBTitv((IBTree *)$3)) )
             $$ = IBTNewOp(IBOpSubIR,(IBTree *)$1,(IBTree *)$3);
           else
             $$ = IBTNewOp(IBOpSubII,(IBTree *)$1,(IBTree *)$3);
	 }
         else
           $$ = IBTNewOp(IBOpSubII,(IBTree *)$1,(IBTree *)$3);
       }

     | ExprMul
       { 
         $$ = (IBTree*)$1; 
       }
     ;

 ExprMul: ExprMul MUL ExprExp
          {
            if( IBTtype((IBTree *)$1)==IBTNodeItv )
            {
              if( IBIsDoubleI(IBTitv((IBTree *)$1)) )
              {
                if( IBMinI(IBTitv((IBTree *)$1))>=0 )
                     $$ = IBTNewOp(IBOpMulRposI,(IBTree *)$1,(IBTree *)$3);
                else
                     $$ = IBTNewOp(IBOpMulRnegI,(IBTree *)$1,(IBTree *)$3);
              }
              else
                $$ = IBTNewOp(IBOpMulII,(IBTree *)$1,(IBTree *)$3);
            }
            else if( IBTtype((IBTree *)$3)==IBTNodeItv )
            {
              if( IBIsDoubleI(IBTitv((IBTree *)$3)) )
              {
                if( IBMinI(IBTitv((IBTree *)$3))>=0 )
                     $$ = IBTNewOp(IBOpMulRposI,(IBTree *)$3,(IBTree *)$1);
                else
                     $$ = IBTNewOp(IBOpMulRnegI,(IBTree *)$3,(IBTree *)$1);
              }
              else
                $$ = IBTNewOp(IBOpMulII,(IBTree *)$1,(IBTree *)$3);
            }
            else
              $$ = IBTNewOp(IBOpMulII,(IBTree *)$1,(IBTree *)$3);
          }

        | ExprMul DIV ExprExp
          {
            if( IBTtype((IBTree *)$1)==IBTNodeItv )
            {
              if( IBIsDoubleI(IBTitv((IBTree *)$1)) )
              {
                if( IBMinI(IBTitv((IBTree *)$1))>=0 )
                     $$ = IBTNewOp(IBOpDivRposI,(IBTree *)$1,(IBTree *)$3);
                else
                     $$ = IBTNewOp(IBOpDivRnegI,(IBTree *)$1,(IBTree *)$3);
              }
              else
                $$ = IBTNewOp(IBOpDivII,(IBTree *)$1,(IBTree *)$3);
            }
            else if( IBTtype((IBTree *)$3)==IBTNodeItv )
            {
              if( IBIsDoubleI(IBTitv((IBTree *)$3)) )
              {
                if( IBMinI(IBTitv((IBTree *)$3))>=0 )
                     $$ = IBTNewOp(IBOpDivIRpos,(IBTree *)$1,(IBTree *)$3);
                else
                     $$ = IBTNewOp(IBOpDivIRneg,(IBTree *)$1,(IBTree *)$3);
              }
              else
                $$ = IBTNewOp(IBOpDivII,(IBTree *)$1,(IBTree *)$3);
            }
            else
            {
              $$ = IBTNewOp(IBOpDivII,(IBTree *)$1,(IBTree *)$3);

            }
          }

        | ExprExp
          { 
            $$ = (IBTree *)$1; 
          }
        ;

 ExprExp: ExprUnit Exposant
          {
            if( $2>1 )
            {
              if( $2==2 ) $$ = IBTNewOp(IBOpSqrI,(IBTree *)$1,IBTNewExp(2));
              else        $$ = IBTNewOp(IBOpPowI,(IBTree *)$1,IBTNewExp($2));

            }
          }
        ;

 Exposant:             { $$ = 1; }
         | POW Integer { $$ = IBStrToInt($2); }
         | POW Ident   { f1 = (IBTree *)$2;
	                 if( (IBTtype(f1)==IBTNodeItv) && (IBIsIntegerI(IBTitv(f1))) )
	                 {
                           $$ = (int)(IBMinI(IBTitv(f1)));
	                 }
                         else
                         {
	                   yyerror("invalid exponent.");
	                   YYABORT;
                         }
	               }
         ;

 ExprUnit: LBR Expr RBR
           { $$ = $2; }

         | Ident
           { $$ = $1; }

         | Number
           {
             if( strcmp($1,"+oo")==0 )
	     {
               IBSetI(itv,IBPosInfinity,IBPosInfinity);
	     }
	     else if( strcmp($1,"-oo")==0 )
	     {
               IBSetI(itv,IBNegInfinity,IBNegInfinity);
	     }
             else
	     {
	       IBStrToI($1,itv);
	     }
             $$ = IBTNewItv(itv);
           }

         | SUB ExprMul
           { $$ = IBTNewOp(IBOpNegI,(IBTree *)$2,IBTNewUseless()); }

         | ADD ExprMul
           { $$ = $2; }

         | SQRT LBR Expr RBR
           { $$ =IBTNewOp(IBOpSqrtI,(IBTree *)$3,IBTNewUseless()); }

         | EXP LBR Expr RBR
           { $$ =IBTNewOp(IBOpExpI,(IBTree *)$3,IBTNewUseless()); }

         | LOG LBR Expr RBR
           { $$ =IBTNewOp(IBOpLogI,(IBTree *)$3,IBTNewUseless()); }

         | MINIMUM LBR Expr COMMA Expr RBR
           { $$ =IBTNewOp(IBOpMinII,(IBTree *)$3,(IBTree *)$5); }

         | MAXIMUM LBR Expr COMMA Expr RBR
           { $$ =IBTNewOp(IBOpMaxII,(IBTree *)$3,(IBTree *)$5); }

         | COS LBR Expr RBR
           { $$ =IBTNewOp(IBOpCosI,(IBTree *)$3,IBTNewUseless()); }

         | SIN LBR Expr RBR
           { $$ =IBTNewOp(IBOpSinI,(IBTree *)$3,IBTNewUseless()); }

         | TAN LBR Expr RBR
           { $$ =IBTNewOp(IBOpTanI,(IBTree *)$3,IBTNewUseless()); }

         | COSH LBR Expr RBR
           { $$ =IBTNewOp(IBOpCoshI,(IBTree *)$3,IBTNewUseless()); }

         | SINH LBR Expr RBR
           { $$ =IBTNewOp(IBOpSinhI,(IBTree *)$3,IBTNewUseless()); }

         | TANH LBR Expr RBR
           { $$ =IBTNewOp(IBOpTanhI,(IBTree *)$3,IBTNewUseless()); }

         | ACOS LBR Expr RBR
           { $$ =IBTNewOp(IBOpAcosI,(IBTree *)$3,IBTNewUseless()); }

         | ASIN LBR Expr RBR
           { $$ =IBTNewOp(IBOpAsinI,(IBTree *)$3,IBTNewUseless()); }

         | ATAN LBR Expr RBR
           { $$ =IBTNewOp(IBOpAtanI,(IBTree *)$3,IBTNewUseless()); }

         | ACOSH LBR Expr RBR
           { $$ =IBTNewOp(IBOpAcoshI,(IBTree *)$3,IBTNewUseless()); }

         | ASINH LBR Expr RBR
           { $$ =IBTNewOp(IBOpAsinhI,(IBTree *)$3,IBTNewUseless()); }

         | ATANH LBR Expr RBR
           { $$ =IBTNewOp(IBOpAtanhI,(IBTree *)$3,IBTNewUseless()); }
         ;


/****************************************************************************
 *                          IDENTIFIERS AND NUMBERS                         *
 ****************************************************************************/
 Ident: IdentName IdentArray
        {
          sprintf(s,"%s%s",$1,$2);

	  if( s[0]=='@' )  /* predefined constant */
	  {
            itv2 = IBGetConstant(constants,s);
            if( itv2==NULL )
	    {
	      strcat(s,": predefined constant unknown.");
              yyerror(s);
	      YYABORT;                      
	    }
            else             /* s is the name of a constant */
	    {
              $$ = IBTNewItv(IBNewCopyI(itv2));
	    }
	  }
          else
	  {
            itv2 = IBGetConstant(constants,s);
            if( itv2==NULL )    /* s is the name of a variable */
	    {
              if ((i = IBIsPresentInV(variables,s))>=0)
              {
                $$ = IBTNewVar(i,-1);
              }
	      else {
		strcat(s,": identifier unknown.");
		yyerror(s);
                YYABORT;
	      }
	    }
            else             /* s is the name of a constant */
	    {
              $$ = IBTNewItv(IBNewCopyI(itv2));
	    }
	  }
        }
      ;

 IdentName: IDENT { strcpy($$,yytext); }
          ;

 IdentArray: { strcpy($$,""); }
           | LSBR Integer RSBR
             {
               sprintf(s,"[%s]",$2);
               strcpy($$,s);
             }
           ;

 Variable: IDENT { strcpy($$,yytext); }
         ;

 PrevSuccNumber:       { $$ = 0; }
               | PREV { $$ = 1; }
               | SUCC { $$ = 2; }
               ;

 Number: Integer         { strcpy($$,$1); }
       | Float           { strcpy($$,$1); }
       | REALPOSINFINITY { strcpy($$,"+oo"); }
       | REALNEGINFINITY { strcpy($$,"-oo"); }
       ;

 Float: FLOAT { strcpy($$,yytext); }
      ;

 Integer: INTEGER { strcpy($$,yytext); }
        ;

 IsNewData:         { $$ = 0; }
          | NEWDATA { $$ = 1; }
          ;
%%
