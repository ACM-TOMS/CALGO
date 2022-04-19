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
 * main.c                                                                   *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "narrowing_hullbox.h"
#include "narrowing_newton.h"
#include "search.h"
#include "main.h"


/* Structures used to represent a constraint system */
IBVariables   variables;           /* array of constrained variables */
IBConstants   constants;           /* array of constants */
IBConstraints constraints;         /* array of constraints */
IBOperations  operations;          /* array of operations */


/* Propagation and bisection algorithms */
IBPropagation      IBfilter;       /* propagation algorithm */
IBLocalPropagation IBfilter2B;     /* 2B-like propagation algorithm used by IBfilter */
IBBisectVar        IBbisect;       /* choice function for the bisected variable */
IBBisectArity      IBsplit;        /* bisection step: generation of sub-domains */


/* Lists used in propagation algorithms, only one allocation in the main function */
IBPropagationList    *IBpropaglist,
                     *IBpropaglistsave;
IBPropagationListCtr *IBpropaglistctr;
IBPropagationGlobal  *IBpropaglistglobal;


/* Matrices and domains used in the Interval Newton method
   => only one allocation in the main function */
IBMInterval *IBMIzero,               /* zero interval matrix */
            *IBMIjacobian,           /* interval Jacobian matrix */
            *IBMIfinalsystem;        /* matrix after preconditionning */
IBMDouble   *IBMDmidjacobian,        /* midpoint of Jacobian */
            *IBMDzero,               /* zero real matrix */
            *IBMDidentity,           /* identity real matrix */
            *IBMDinverse;            /* used in the inversion process for preconditionners */
IBDomains   IBdnwt1,
            IBdnwt2,
            IBdnwt3,
            IBdnwt4;
int         IBComputableIntervalNewton; /* true if Interval Newton can be applied */


/* Array used in the 3B consistency algorithm */
double *IBwidth3B;


/* Global pragmas */
long   IBPragmaNbGeneratedDomains; /* number of generated boxes */
unsigned long  IBPragmaMaxSolution; /* maximum number of output boxes */
double IBPragmaPrecision;          /* maximum width of output boxes */
int    IBPragmaBisection;          /* bisection strategy */
int    IBPragmaNumberBisection;    /* arity of each bisection operation */
double IBPragmaImprovement;        /* parameter for propagation algorithms */
double IBPragmaPrecisionShrink;    /* precision of box consistency */
double IBPragmaPrecision3B;        /* precision of 3B consistency */
int    IBPragmaStyleInterval;      /* interval printing style */
int    IBPragmaIntervalDigits;     /* number of digits for interval printing */
int    IBPragmaHullMode;           /* true the result is the hull of all the output boxes */
long   IBPragmaMaxTime;            /* stop after IBPragmaMaxTime milliseconds */
int    IBargBisectNo;              /* 1 if no biscetion */
int    IBcompute3B;                /* 1 if 3B consistency is enforced */
int    IBcomputeWeak3B;            /* 1 if weak 3B consistency is enforced */
int    IBPragmaSubpaving;          /* 1 if a subpaving is computed */

/* Variables and functions used for parsing */
extern FILE* yyin;                 /* input file */
int    IBNbParsingError = 0;       /* number of parsing errors */
int    IByyline;                   /* current line in the input file */
char   IBParsingError[100] = "";   /* error message */



int yyerror(char* s)
/***************************************************************************
*  Called after a parsing error
*     => copy of s in IBParsinError
*/
{
  IBNbParsingError++;
  strcpy(IBParsingError,s);
}


int yywarning(char* s)
/***************************************************************************
*  Called after a warning in the parsing process
*/
{
  printf("  !! l.%d: warning: %s\n",IByyline,s);
}


int IBparser(char *namefile)
/***************************************************************************
*  Parser of the Software
*/
{
  yyin = fopen(namefile, "r");  /* the existence of the file has already been checked */
  yyparse();
  fclose(yyin);

  if( IBNbParsingError==0 ) return( 1 );
  else return( 0 );
}


void IBhelp()
/***************************************************************************
*   Help of the Software
*/
{
  printf("SYNOPSIS\n");
  printf("   realpaver [options...] filename\n\n");
  printf("   Where options include:\n");
  printf("      -union            return the list of computed boxes\n");
  printf("      -hull             return the hull of computed boxes\n");
  printf("      -bound            display of intervals: [a,b]\n");
  printf("      -midpoint         display of intervals: m + [-e,+e]\n");
  printf("\n");
  printf("      -precision value  minimum width of divided boxes\n"); 
  printf("      -nosplit          disable the splitting process\n");
  printf("      -split2           domains split in two parts\n");
  printf("      -split3           domains split in three parts\n");
  printf("      -rr               round-robin splitting strategy\n");
  printf("      -lf               largest-first splitting strategy\n");
  printf("      -mn               maximum narrowing splitting strategy\n");
  printf("      -paving           compute a paving of the solution space\n");
  printf("                        => Breadth-First-Search and largest_first choice\n");
  printf("      -points           minimum width of divided boxes\n"); 
  printf("                        => Depth-First-Search\n");
  printf("      -number +oo       compute all possible boxes\n");
  printf("      -number value     maximum number of computed boxes\n");
  printf("\n");
  printf("      -hc3              hull consistency over decomposed constraints\n");
  printf("      -hc4              hull consistency over user' constraints\n");
  printf("      -hc4_newton       hc4 combined with the interval Newton method\n");
  printf("      -bc3              box consistency\n");
  printf("      -bc3_newton       bc3 combined with the interval Newton method\n");
  printf("      -bc4              hc4 combined with bc3\n");
  printf("      -bc5              bc4 combined with the interval Newton method\n");
  printf("      -3B value         3B(value) consistency\n");
  printf("      -weak3B value     weak 3B(value) consistency\n");

  printf("\n   All these options can be parameterized in filename. See\n");
  printf("   the documentation for that.\n");

  printf("\n\nDESCRIPTION\n");

  printf("    RealPaver is a real constraint paver.\n\n");
  printf("    Given a system of nonlinear constraints over the real numbers,\n");
  printf("    it computes a paving of the solution set. A paving is a union\n");
  printf("    of Cartesian product of interval domains (boxes) enclosing\n");
  printf("    the solution set.\n\n");

  printf("    The solving engine of RealPaver implements a branch-and-prune\n");
  printf("    algorithm that alternates reduction steps of variable domains\n");
  printf("    by means of constraint consistency techniques and splitting steps.\n");

  printf("\n\nINPUT FILE\n");
  printf("    The input file is composed of a nonlinear constraint system\n");
  printf("    (constants, variables and domains, constraints) and flags used\n");
  printf("    to parameterize the resolution.\n\n");
  printf("    See the files benchmarks/abc or benchmarks/circle\n");

  printf("\n\nREFERENCES\n");
  printf("\n    [1] F. Benhamou, F. Goualard, L. Granvilliers and J.-F. Puget.\n");
  printf("        Revising Hull and Box Consistency. Procs. of ICLP'99,\n");
  printf("        Las Cruces, USA, 1999. The MIT Press.\n");

  printf("\n    [2] E. Davis. Constraint Propagation with Interval Labels.\n");
  printf("        Artificial Intelligence, 32:281-331, 1987\n");

  printf("\n    [3] L. Granvilliers. On the Combination of Interval Constraint Solvers.\n");
  printf("        Reliable Computing, 7(6):467:483, 2001.\n");

  printf("\n    [4] R. B. Kearfott. Some Tests of Generalized Bisection.\n");
  printf("        ACM TOMS, 13(3):197-220, 1987.\n");

  printf("\n    [5] O. Lhomme. Consistency Techniques for Numeric CSPs.\n");
  printf("        Procs. of IJCAI'93, Chambery, France, 1993. Morgan Kaufman\n");

  printf("\n    [6] R. E. Moore. Interval Analysis. Prentice-Hall,\n");
  printf("        Englewood Cliffs, NJ, 1966.\n");

  printf("\n    [7] P. Van Hentenryck, L. Michel and Y. Deville. Numerica:\n");
  printf("        a Modeling Language for Global Optimization. MIT Press, 1997.\n");

  printf("\n");
}


main(int argc, char **argv)
/***************************************************************************
*  The main function of the Software
*/
{
  char version[] = VERSION,          /* version of the Software from config.h */
       software_run[] = PACKAGE,     /* name of the Software from config.h */
       software_name[] = SOFTWARE_NAME, /* name of the Software from config.h */
       IBarguments[100] = "",        /* used for printing wrong arguments */
       slong[30],                    /* used for printing long integers */
       *IBinfile,                    /* name of input file from argv */
       IBInfoConsistency[50];        /* name of consistency enforced */

  FILE *faux;                        /* used to open the input file */

  int  argBisect,                    /* type of bisection strategy */
       argHelp,                      /* 1 if the help of the Software is printed */
       argVerbose,                   /* 1 if verbose mode */
       argPhi,                       /* 1 if phi of box_phi consistency is parameterized */
       argFile,                      /* 1 if an input file exists */
       slongesp,                     /* used for printing long integers */
       sysVar,                       /* information on the variables occurring in constraints */
       sysCtr,                       /* information on the constraints */
       nEq,                          /* number of equations */
       nIneq,                        /* number of inequations */
       nOccOne,                      /* number of variables occurring once */
       nOccMul,                      /* number of variables occurring more than once */
       completeProcess,              /* 1 if the used algorithm is complete */
       nbsol,                        /* number of output boxes */
       i, n;

  long IBnblocvar,                   /* sum_(c in constraints) cardinal({x occurring in c}) */
       time;                         /* used to implement the online mode */

  double x;


/* Initialisation of structures */
  variables        = IBNewV();
  constants        = IBNewConstants();
  constraints      = IBNewConstraints();
  operations       = IBOperationsInit();


/* Initialization of pragmas */
  IBPragmaNbGeneratedDomains = 1;
  IBPragmaMaxSolution        = 1024;
  IBPragmaPrecision          = 1.0e-8;
  IBPragmaBisection          = IBBisectRoundRobin;
  IBPragmaNumberBisection    = 3;
  IBPragmaImprovement        = 0.9;    /* 10% */
  IBPragmaPrecisionShrink    = 0.0;
  IBPragmaPrecision3B        = 0.001;
  IBPragmaStyleInterval      = IBPrintIntervalBounds;
  IBPragmaIntervalDigits     = 16;
  IBPragmaHullMode           = 0; /* default: union mode, returns the list of output boxes */
  IBPragmaMaxTime            = 1000000000;   /* default: 1 million seconds (11.57 days) */
  IBPragmaSubpaving          = 0; /* default: no subpaving */

/* Initialization of other variables */
  IBComputableIntervalNewton = 0;
  IBfilter2B                 = IBF2Bbc5;
  IBsplit                    = IBBsplit3;
  IBcompute3B                = 0;
  IBcomputeWeak3B            = 0;
  IBwidth3B                  = NULL; 

  IBpropaglist               = NULL;
  IBpropaglistsave           = NULL;
  IBpropaglistctr            = NULL;
  IBpropaglistglobal         = NULL;

  completeProcess = 1;

/* Initialization of modules */
  IBClockInit();  /* clock management */
  IBInitIA();     /* interval arithmetic */


/* Initialization of arguments */
  argHelp = argFile = 0;
  argBisect = IBBisectNone;
  IBargBisectNo = 0;
  argVerbose = 0;
  argPhi = 0;

  printf("\n%s v. %s (c) LINA 2004\n",software_name,version);
  /* printf("[report bugs to Laurent.Granvilliers@lina.univ-nantes.fr]\n"); */
  printf("\n");


  /* PARSING OF THE ARGUMENTS ON THE COMMAND LINE */
  i=1;
  while( i<argc )
  {
    if( (faux=fopen(argv[i],"r"))!=NULL )
    {
      fclose(faux);
      IBinfile = argv[i];
      argFile = 1;
    }
    else
    {
      /* CONSISTENCY */
      if( strcmp(argv[i],"-bc3")==0 )
      {
        IBfilter2B = IBF2Bbc3;
      }
      else if( strcmp(argv[i],"-bc3_newton")==0 )
      {
        IBfilter2B = IBF2Bbc3Newton;
      }
      else if( strcmp(argv[i],"-bc4")==0 )
      {
        IBfilter2B = IBF2Bbc4;
      }
      else if( strcmp(argv[i],"-bc5")==0 )
      {
        IBfilter2B = IBF2Bbc5;
      }
      else if( strcmp(argv[i],"-hc3")==0 )
      {
        IBfilter2B = IBF2Bhc3;
      }
      else if( strcmp(argv[i],"-hc4")==0 )
      {
        IBfilter2B = IBF2Bhc4;
      }
      else if( strcmp(argv[i],"-hc4I")==0 )
      {
        IBfilter2B = IBF2Bhc4I;
      }
      else if( strcmp(argv[i],"-hc4_newton")==0 )
      {
        IBfilter2B = IBF2Bhc4Newton;
      }
      else if( strcmp(argv[i],"-3B")==0 )
      {
        IBcompute3B = 1;              /* 3B consistency is enforced */
        if( i<argc-1 )
	{
          if( (x=strtod(argv[i+1],NULL))>0.0 )
          {
            IBPragmaPrecision3B = x;  /* precision of 3B */
            i++;
	  }
	}
      }
      else if( strcmp(argv[i],"-weak3B")==0 )
      {
        IBcomputeWeak3B = 1;          /* weak 3B consistency is enforced */
        if( i<argc-1 )
	{
          if( (x=strtod(argv[i+1],NULL))>0.0 )
          {
            IBPragmaPrecision3B = x;  /* precision of 3B */
            i++;
	  }
	}
      }

      /* PRECISION */
      else if( strcmp(argv[i],"-precision")==0 )
      {
        if( i<argc-1 )
	{
          if( (x=strtod(argv[i+1],NULL))>0.0 )
          {
            IBPragmaPrecision = x;
            i++;
	  }
	}
      }

      /* BISECTION */
      else if( strcmp(argv[i],"-paving")==0 )
      {
        IBPragmaSubpaving = 1;
      }
      else if( strcmp(argv[i],"-points")==0 )
      {
        IBPragmaSubpaving = 0;
      }
      else if( strcmp(argv[i],"-rr")==0 )
      {
        argBisect = IBBisectRoundRobin;
      }
      else if( strcmp(argv[i],"-lf")==0 )
      {
        argBisect = IBBisectLargestFirst;
      }
      else if( strcmp(argv[i],"-mn")==0 )
      {
        argBisect = IBBisectMaxNarrow;
      }
      else if( strcmp(argv[i],"-bisect2")==0 )
      {
        IBsplit = IBBsplit2;
        IBPragmaNumberBisection = 2;
      }
      else if( strcmp(argv[i],"-split2")==0 )
      {
        IBsplit = IBBsplit2;
        IBPragmaNumberBisection = 2;
      }
      else if( strcmp(argv[i],"-bisect3")==0 )
      {
        IBPragmaNumberBisection = 3;
        IBsplit = IBBsplit3;
      }
      else if( strcmp(argv[i],"-split3")==0 )
      {
        IBPragmaNumberBisection = 3;
        IBsplit = IBBsplit3;
      }
      else if( strcmp(argv[i],"-nobisect")==0 )
      {
        IBPragmaNumberBisection = 0;
        IBargBisectNo = 1;
      }
      else if( strcmp(argv[i],"-nosplit")==0 )
      {
        IBPragmaNumberBisection = 0;
        IBargBisectNo = 1;
      }
      else if( strcmp(argv[i],"-number")==0 )
      {
        if( i<argc-1 )
	{
          if( strcmp(argv[i+1],"+oo")==0 )
	  {
            IBPragmaMaxSolution = ~0;    /* stop after "all" output boxes */
            i++;
	  }
          else if( (n=atoi(argv[i+1]))>0 )
          {
            IBPragmaMaxSolution = n;   /* maximum number of output boxes */
            i++;
	  }
          else
	  {
	    if (strcmp(IBarguments,"")!=0) strcat(IBarguments,", ");
            strcat(IBarguments,argv[i]);
	  }
	}
        else
        {
	  if (strcmp(IBarguments,"")!=0) strcat(IBarguments,", ");
          strcat(IBarguments,argv[i]);
        }
      }

      else if( strcmp(argv[i],"-hull")==0 )
      {
        IBPragmaHullMode = 1;
      }
      else if( strcmp(argv[i],"-union")==0 )
      {
        IBPragmaHullMode = 0;
      }

      /* PARAMETERS OF PROPAGATION ALGORITHMS */
      else if( strcmp(argv[i],"-improve")==0 )
      {
        if( i<argc-1 )
	{
          x = 0.0;
          x = strtod(argv[i+1],NULL);
          if( x>0.0 )
          {
            IBPragmaImprovement = 1 - (x/100);
            i++;
	  }
          else if( x==0.0 )
	  {
            IBPragmaImprovement = 1.0;
            i++;
	  }
	}
      }
      else if( strcmp(argv[i],"-phi")==0 )
      {
        if( i<argc-1 )
	{
          x = 0.0;
          x = atof(argv[i+1]);
          if( x>=0.0 )
          {
            IBPragmaPrecisionShrink = x;
            argPhi = 1;
            i++;
	  }
	}
        else
	{
	  if (strcmp(IBarguments,"")!=0) strcat(IBarguments,", ");
          strcat(IBarguments,argv[i]);
	}
      }

      /* INTERVAL PRINTING STYLE */
      else if( strcmp(argv[i],"-bound")==0 )
      {
        IBPragmaStyleInterval = IBPrintIntervalBounds;
      }
      else if( strcmp(argv[i],"-midpoint")==0 )
      {
        IBPragmaStyleInterval = IBPrintIntervalMidError; 
      }

      /* VERBOSE MODE AND HELP */
      else if( (strcmp(argv[i],"-verb")==0) || (strcmp(argv[i],"-verbose")==0) )
      {
        argVerbose = 1;
      }
      else if( (strcmp(argv[i],"-h")==0) || (strcmp(argv[i],"-help")==0) )
      {
        argHelp = 1;
	i = argc;
      }

      /* ONLINE MODE: MAX TIME */
      else if( strcmp(argv[i],"-time")==0 )
      {
        if( i<argc-1 )
	{
          time = 0;
          time = strtoul(argv[i+1],NULL,10);
          if( time>0 )
          {
            IBPragmaMaxTime = time;
            i++;
	  }
	}
        else
	{
	  if (strcmp(IBarguments,"")!=0) strcat(IBarguments,", ");
          strcat(IBarguments,argv[i]);
	}
      }

      /* WRONG ARGUMENT */
      else
      {
	  if (strcmp(IBarguments,"")!=0) {
            strcat(IBarguments,", ");
	  }
          strcat(IBarguments,argv[i]);
      }
    }
    i++;
  }


  if( argHelp )
  {
    IBhelp();

    IBFreeConstraints(constraints);
    IBFreeV(variables);
    IBOperationsFree(operations);
    IBClockFree();
    return;
  }

  if( argVerbose )
  {
    printf("  !! sorry, no verbose mode...\n\n");
  }

  if( !argFile )
  {
    printf("  !! error: input file not found [use '%s -h' for help]\n\n",software_run);
    return;
  }

  if (strcmp(IBarguments,"")!=0) {
    printf("  !! wrong arguments on the command line: %s\n\n",IBarguments);
  }


  /* PARSING
   *------------------------------------------------------------------------*/
  IBClockBegin(IBClockParse);
  IByyline = 1;


#if SOFTWARE_PROFILE
  printf("PARSING\n");
  printf("  File: %s\n",IBinfile);
#endif

  if( IBparser(IBinfile) )
  {
    IBClockEnd(IBClockParse);


#if SOFTWARE_PROFILE
    _IBprintlong(slong,IBClockGet(IBClockParse),1);
    printf("  Elapsed time: %s ms\n",slong);
#endif


    /* end of parsing: creation of constraints integer(x) */
    IBAddIntegerTypeConstraints(constraints,variables);

    /* after parsing: reallocation of structures */
    IBReallocV(variables);
    IBReallocConstraints(constraints);
    IBFreeConstants(constants);  /* constants are no more used */

    /* creation of the dependencies between variables and constraints */
    IBCreateDependencies(constraints,variables);

    /* creation of propagation lists */
    IBpropaglist       = IBPLnew(IBCNbProj(constraints));
    IBpropaglistsave   = IBPLnew(IBVnb(variables));
    IBpropaglistctr    = IBPLCnew(IBCNbCtr(constraints));
    IBpropaglistglobal = IBPGlobalNew(IBCNbCtr(constraints));


    /* Type of the constraint system
     *------------------------------*/
    IBnblocvar = nOccOne = nOccMul = nEq = nIneq = 0;
    for( i=0; i<IBCNbCtr(constraints); i++ )
    {
      if( IBCNbProjOne(IBCCtr(constraints,i)) ) nOccOne++;
      if( IBCNbProjMul(IBCCtr(constraints,i)) ) nOccMul++;
      if( IBCrel(IBCCtr(constraints,i))==IBRelationEQU ) nEq++;
      else nIneq++;
      IBnblocvar += IBCNbVar(IBCCtr(constraints,i));
    }
    if( nOccMul==0 )      sysVar = IBSystemOccOne;
    else if( nOccOne==0 ) sysVar = IBSystemOccMul;
    else                  sysVar = IBSystemOccAll;

    if( nEq==0 )          sysCtr = IBSystemCtrIneq;
    else if( nIneq==0 )   sysCtr = IBSystemCtrEq;
    else                  sysCtr = IBSystemCtrAll;


    /* Is the Interval Newton method computable ?
     *-------------------------------------------*/
    if( (IBCNbCtr(constraints)>=IBVnb(variables)) && (IBCNbCtr(constraints)>1) )
    {
      IBComputableIntervalNewton=1;
      i=0;
      while( IBComputableIntervalNewton && (i<IBVnb(variables)) )
      {
        if( (IBCrelfunc(IBCCtr(constraints,i))!=IBRelationEQU) ||
            (!IBTDerivable(IBCfunc(IBCCtr(constraints,i)))) )
	{
          IBComputableIntervalNewton = 0;
	}
        else
        {
          i++;
	}
      }
      if( IBComputableIntervalNewton )
      {
        IBMIjacobian    = IBMIntervalNew(IBVnb(variables));
        IBMIfinalsystem = IBMIntervalNew(IBVnb(variables));
        IBMIzero        = IBMIntervalNewZero(IBVnb(variables));
        IBMDmidjacobian = IBMDoubleNew(IBVnb(variables));
        IBMDinverse     = IBMDoubleNew(IBVnb(variables));
        IBMDzero        = IBMDoubleNewZero(IBVnb(variables));
        IBMDidentity    = IBMDoubleNewIdentity(IBVnb(variables));
        IBdnwt1         = IBNewD(IBVnb(variables));
        IBdnwt2         = IBNewD(IBVnb(variables));
        IBdnwt3         = IBNewD(IBVnb(variables));
        IBdnwt4         = IBNewD(IBVnb(variables));
      }
    }
    else IBComputableIntervalNewton = 0;

    /* Reasoning about the best 2B consistency-based propagation algorithm
     *--------------------------------------------------------------------*/
    if( IBfilter2B==IBF2Bbc5 )
    {
      if( sysVar==IBSystemOccOne )  /* bc5 equivalent to hc4 */
      {
        IBfilter2B = IBF2Bhc4Newton;
      }
      else
      if( sysVar==IBSystemOccMul )  /* bc5 equivalent to bc3 */
      {
        IBfilter2B = IBF2Bbc3Newton;
      }
    }


    /* If Interval Newton cannot be computed, BC5, BC3+Newton and HC4+Newton
     * must be avoided
     *--------------------------------------------------------------------*/
    if (!IBComputableIntervalNewton) {
      if( IBfilter2B==IBF2Bbc5 ) {
        IBfilter2B = IBF2Bbc4;            /* BC4 instead of BC5 */
      }
      else if ( IBfilter2B==IBF2Bhc4Newton) {
        IBfilter2B = IBF2Bhc4;            /* HC4 instead of HC4+Newton */
      }
      else if ( IBfilter2B==IBF2Bbc3Newton ) {
        IBfilter2B = IBF2Bbc3;            /* BC3 instead of BC3+Newton */
      }
    }


    /* Exact box consistency is computed if Interval Newton is not used
     * and if no phi argument on the command line or in the input file is specified
     * Otherwise, box(1.0e-8) consistency is computed
     *-----------------------------------------------------------------*/
    if( (argPhi==0) && ((IBfilter2B!=IBF2Bbc3) && (IBfilter2B!=IBF2Bbc4)) )
    {
      IBPragmaPrecisionShrink = 1.0e-8;
    }


    /* Reasoning about the choice algorithm for the bisected variable's domain
     *------------------------------------------------------------------------*/

    /* Force the strategy largest_first if a subpaving is computed */
    if( IBPragmaSubpaving )
    {
      IBPragmaBisection = IBBisectLargestFirst;
      IBbisect = IBBlf;
    }
    else
    {
      if( argBisect!=IBBisectNone )       /* bisection specified on comand line */
           IBPragmaBisection = argBisect;

      if( IBPragmaBisection==IBBisectRoundRobin )
      {
        IBbisect = IBBrr;
      }
      else if( IBPragmaBisection==IBBisectLargestFirst )
      {
        IBbisect = IBBlf;
      }
      else   /* MAX REDUCTION strategy, used only if Interval Newton is applied */
      {
        if( IBComputableIntervalNewton &&
            ((IBfilter2B==IBF2Bbc5) ||
             (IBfilter2B==IBF2Bhc4Newton) ||
             (IBfilter2B==IBF2Bbc3Newton)) )
        {
          IBbisect = IBBmn;
          IBPragmaBisection = IBBisectMaxNarrow;
        }
        else
        {
          IBbisect = IBBrr;
          IBPragmaBisection = IBBisectRoundRobin;
        }
      }
    }


    /* Association 2B consistency-based propagation algorithm (IBfilter2B) /
       propagation algorithm (IBfilter)
    *----------------------------------------------------------------------*/

    if( IBcompute3B )
    {
      IBfilter=IBF3B;
      IBwidth3B = (double *)malloc(IBVnb(variables)*sizeof(double));
    }
    else if( IBcomputeWeak3B )
    {
      IBfilter=IBF3Bweak;
      IBwidth3B = (double *)malloc(IBVnb(variables)*sizeof(double));
    }
    else
    {
      if( IBfilter2B==IBF2Bhc3 )
      {
        IBfilter=IBFhc3;
      }
      else if (IBfilter2B==IBF2Bhc4 )
      {
        IBfilter=IBFhc4;
      }
      else if (IBfilter2B==IBF2Bhc4I )
      {
        IBfilter=IBFhc4I;
      }
      else if (IBfilter2B==IBF2Bhc4Newton )
      {
        IBfilter=IBFhc4Newton;
      }
      else if (IBfilter2B==IBF2Bbc3 )
      {
        IBfilter=IBFbc3;
      }
      else if (IBfilter2B==IBF2Bbc3Newton )
      {
        IBfilter=IBFbc3Newton;
      }
      else if (IBfilter2B==IBF2Bbc4 )
      {
        IBfilter=IBFbc4;
      }
      else if (IBfilter2B==IBF2Bbc5 )
      {
        IBfilter=IBFbc5;
      }
    }
  

    /* Propagation algorithm used => message
     *--------------------------------------*/
    if( IBfilter==IBFbc5 )
    {
      sprintf(IBInfoConsistency,"BC5(%.4g)",IBPragmaPrecisionShrink);
    }
    else if( IBfilter==IBFbc4 )
    {
      sprintf(IBInfoConsistency,"BC4(%.4g)",IBPragmaPrecisionShrink);
    }
    else if( IBfilter==IBFbc3 )
    {
      sprintf(IBInfoConsistency,"BC3(%.4g)",IBPragmaPrecisionShrink);
    }
    else if( IBfilter==IBFbc3Newton )
    {
      sprintf(IBInfoConsistency,"BC3(%.4g) + Newton",IBPragmaPrecisionShrink);
    }
    else if( IBfilter==IBFhc4Newton )
    {
      sprintf(IBInfoConsistency,"HC4 + Newton");
    }
    else if( IBfilter==IBFhc3 )
    {
      sprintf(IBInfoConsistency,"HC3");
    }
    else if( IBfilter==IBFhc4I )
    {
      sprintf(IBInfoConsistency,"HC4/Intervals");
    }
    else if( IBfilter==IBFhc4 )
    {
      sprintf(IBInfoConsistency,"HC4/Union");
    }
    else if( (IBfilter==IBF3B) || (IBfilter==IBF3Bweak) )
    {
      if( IBfilter==IBF3B )
      {
        sprintf(IBInfoConsistency,"3B(%.4g) / ",IBPragmaPrecision3B);
      }
      else
      {
        sprintf(IBInfoConsistency,"weak 3B(%.4g) / ",IBPragmaPrecision3B);
      }

      if( IBfilter2B==IBF2Bhc3 )
      {
        strcat(IBInfoConsistency,"HC3");
      }
      else if (IBfilter2B==IBF2Bhc4 )
      {
        strcat(IBInfoConsistency,"HC4");
      }
      else if (IBfilter2B==IBF2Bhc4I )
      {
        strcat(IBInfoConsistency,"HC4I");
      }
      else if (IBfilter2B==IBF2Bhc4Newton )
      {
        strcat(IBInfoConsistency,"HC4 + Newton");
      }
      else if (IBfilter2B==IBF2Bbc3 )
      {
        strcat(IBInfoConsistency,"BC3");
      }
      else if (IBfilter2B==IBF2Bbc3Newton )
      {
        strcat(IBInfoConsistency,"BC3 + Newton");
      }
      else if (IBfilter2B==IBF2Bbc4 )
      {
        strcat(IBInfoConsistency,"BC4");
      }
      else if (IBfilter2B==IBF2Bbc5 )
      {
        strcat(IBInfoConsistency,"BC5");
      }
    }

    /* Information on solving algorithms
     *----------------------------------*/
#if SOFTWARE_PROFILE
    printf("\nSOLVING\n");

    printf("  System: %d x %d",IBCNbCtr(constraints),
                                                 IBVnb(variables));
    if( sysCtr==IBSystemCtrEq )        printf(" of equations");
    else if( sysCtr==IBSystemCtrIneq ) printf(" of inequations");
    else                               printf(" of equations/inequations");

    printf(", density: %.3g\n",
             (double)IBnblocvar/
             ((double)IBCNbCtr(constraints)*(double)IBVnb(variables)));

    printf("  Reduction: %s\n",IBInfoConsistency);

    if( !IBargBisectNo ) {
      printf("  Split: ");
      if( IBPragmaSubpaving ) printf("paving");
      else                    printf("points");

      if( IBPragmaBisection==IBBisectRoundRobin )         printf(", round-robin");
      else if( IBPragmaBisection==IBBisectLargestFirst )  printf(", largest-first");
      else if( IBPragmaBisection==IBBisectMaxNarrow )     printf(", max-narrow");

      printf(", %d parts",IBPragmaNumberBisection);
      printf(", precision: %g\n\n",IBPragmaPrecision);
    }
#endif
    

    /* SOLVING */
    IBClockBegin(IBClockSolve);
    if( nbsol=IBBisection(IBDomVars(variables),IBargBisectNo,&completeProcess) )
    {
      IBClockEnd(IBClockSolve);
      printf("\nEND OF SOLVING\n");
    }
    else
    {
      printf("\nEND OF SOLVING\n");
      IBClockEnd(IBClockSolve);
      if( completeProcess )
      {
        printf("  Property:     no solution in the initial box\n");
      }
      else
      {
        _IBprintlong(slong,IBClockGet(IBClockSolve),1);
        printf("\n  Property:     no consistent box at %g found in %s ms\n",IBPragmaPrecision,slong);
      }
    }


    /* Information on resolution
    *--------------------------*/
    if( nbsol>0 )
    {
      if( completeProcess )
      {
        printf("  Property:     reliable process (no solution is lost)\n");
      }
      else
      {
        printf("  Property:     non reliable process (some solutions may be lost)\n");
      }
    }


    _IBprintlong(slong,IBClockGet(IBClockSolve),1);
    printf("  Elapsed time: %s ms\n",slong);

 
#if SOFTWARE_PROFILE
    printf("\nPROFILING\n");
   if( IBPragmaNbGeneratedDomains<=1 ) {
      printf("  Split: no, 1 box examined\n");
    }
    else {
      _IBprintlong(slong,(IBPragmaNbGeneratedDomains-1)/IBPragmaNumberBisection,1);
      printf("  Split: %s",slong);

      _IBprintlong(slong,IBPragmaNbGeneratedDomains,1);
      printf(", %s boxes examined\n",slong);
    }


    printf("  Solving time in consistency algorithms:\n");


    /* The solving times for HC3 and HC4 may be wrong due to rounding errors...
       if there are a lot of iterations such that each one is very short */

    if( IBClockGet(IBClockHC3)>0 )  /* HC3 is used */
    {
      if( IBClockGet(IBClockHC3) +
          IBClockGet(IBClockBC3) +
          IBClockGet(IBClockINwt) > IBClockGet(IBClockSolve) )
      {
        IBClockSet(IBClockHC3,IBClockGet(IBClockSolve) -
                              IBClockGet(IBClockINwt)  -
                              IBClockGet(IBClockBC3));
      }
    }

    /* Note that HC3 and HC4 cannot be used together */
    if( IBClockGet(IBClockHC4)>0 )  /* HC4 is used */
    {
      if( IBClockGet(IBClockHC4) +
          IBClockGet(IBClockBC3) +
          IBClockGet(IBClockINwt) > IBClockGet(IBClockSolve) )
      {
        IBClockSet(IBClockHC4,IBClockGet(IBClockSolve) -
                              IBClockGet(IBClockINwt)  -
                              IBClockGet(IBClockBC3));
      }
    }

    slongesp=0;
    if (slongesp<(i=_IBNbDigits(IBClockGet(IBClockBC3))))  slongesp=i;
    if (slongesp<(i=_IBNbDigits(IBClockGet(IBClockHC3))))  slongesp=i;
    if (slongesp<(i=_IBNbDigits(IBClockGet(IBClockHC4))))  slongesp=i;
    if (slongesp<(i=_IBNbDigits(IBClockGet(IBClockINwt)))) slongesp=i;
    slongesp+=(slongesp-1)/3;

    _IBprintlong(slong,IBClockGet(IBClockBC3),slongesp);
    printf("     in BC3:    %s ms\n",slong);

    _IBprintlong(slong,IBClockGet(IBClockHC3),slongesp);
    printf("     in HC3:    %s ms\n",slong);

    _IBprintlong(slong,IBClockGet(IBClockHC4),slongesp);
    printf("     in HC4:    %s ms\n",slong);

    _IBprintlong(slong,IBClockGet(IBClockINwt),slongesp);
    printf("     in Newton: %s ms\n",slong);

    if( IBfilter2B==IBF2Bhc3 )
    {
      _IBprintlong(slong,IBNbFreshVar(variables),1);
      printf("\n  Number of fresh variables in HC3: %s\n",slong);
    }

    IBProfileIA();
#endif


    printf("\n");
  }
  else
  {
    printf("  !! l.%d: error: %s\n\n",IByyline,IBParsingError);
  }


  /* the end: desallocation of global structures --*/
  if( IBComputableIntervalNewton )
  {
    IBMIntervalFree(IBMIjacobian);
    IBMIntervalFree(IBMIfinalsystem);
    IBMIntervalFree(IBMIzero);
    IBMDoubleFree(IBMDmidjacobian);
    IBMDoubleFree(IBMDinverse);
    IBMDoubleFree(IBMDzero);
    IBMDoubleFree(IBMDidentity);
    IBFreeD(IBdnwt1);
    IBFreeD(IBdnwt2);
    IBFreeD(IBdnwt3);
    IBFreeD(IBdnwt4);
  }

  if( IBpropaglist!=NULL )       IBPLfree(IBpropaglist);
  if( IBpropaglistsave!=NULL )   IBPLfree(IBpropaglistsave);
  if( IBpropaglistctr!=NULL )    IBPLCfree(IBpropaglistctr);
  if( IBpropaglistglobal!=NULL ) IBPGlobalFree(IBpropaglistglobal);

  if( IBwidth3B!=NULL ) free(IBwidth3B);

  IBFreeConstraints(constraints);
  IBFreeV(variables);
  IBOperationsFree(operations);
  IBClockFree();
}
