
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse         ppParse_bisonparse
#define yylex           ppParse_bisonlex
#define yyerror         ppParse_bisonerror
#define yylval          ppParse_bisonlval
#define yychar          ppParse_bisonchar
#define yydebug         ppParse_bisondebug
#define yynerrs         ppParse_bisonnerrs


/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 3 "ppParse_bison.y"

#include "ppParse.h"

int declarationType = -1, isDefined = 0, lastExpWasNumber = 0;
int definingType = -1, typeFound = -1, locFound = -1;
extern int ppUsingParamHomotopy, ppAddAllNums;
extern preProcessArray *ppArray;
int *endOfFile;
int *lineNum;



/* Line 189 of yacc.c  */
#line 94 "src/ppParse_bison.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ppEND = 258,
     ppHOMVARGP = 259,
     ppVARGP = 260,
     ppPATHVAR = 261,
     ppVAR = 262,
     ppPARAM = 263,
     ppCONST = 264,
     ppFUNC = 265,
     ppSUBFUNC = 266,
     ppI = 267,
     ppPi = 268,
     ppNUMBER = 269,
     ppSCE = 270,
     ppPOW = 271,
     ppNAME = 272,
     ppNEG = 273
   };
#endif
/* Tokens.  */
#define ppEND 258
#define ppHOMVARGP 259
#define ppVARGP 260
#define ppPATHVAR 261
#define ppVAR 262
#define ppPARAM 263
#define ppCONST 264
#define ppFUNC 265
#define ppSUBFUNC 266
#define ppI 267
#define ppPi 268
#define ppNUMBER 269
#define ppSCE 270
#define ppPOW 271
#define ppNAME 272
#define ppNEG 273




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 16 "ppParse_bison.y"

  char *name;



/* Line 214 of yacc.c  */
#line 172 "src/ppParse_bison.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 184 "src/ppParse_bison.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  18
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   91

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  28
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  9
/* YYNRULES -- Number of rules.  */
#define YYNRULES  39
/* YYNRULES -- Number of states.  */
#define YYNSTATES  64

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   273

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      25,    26,    20,    18,    27,    19,     2,    21,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    24,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    23,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    22
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     4,     6,     9,    13,    15,    17,    19,
      21,    23,    27,    31,    35,    39,    43,    46,    49,    53,
      55,    57,    59,    61,    63,    65,    70,    77,    82,    86,
      88,    90,    92,    94,    96,    98,   100,   102,   104,   108
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      29,     0,    -1,    -1,     3,    -1,    35,    36,    -1,    30,
      24,    31,    -1,    17,    -1,    14,    -1,    15,    -1,    12,
      -1,    13,    -1,    31,    18,    31,    -1,    31,    19,    31,
      -1,    31,    20,    31,    -1,    31,    21,    31,    -1,    31,
      23,    31,    -1,    18,    31,    -1,    19,    31,    -1,    25,
      31,    26,    -1,    32,    -1,    33,    -1,    14,    -1,    17,
      -1,    12,    -1,    13,    -1,    15,    25,    31,    26,    -1,
      16,    25,    31,    27,    31,    26,    -1,    17,    25,    34,
      26,    -1,    34,    27,    31,    -1,    31,    -1,     4,    -1,
       5,    -1,     6,    -1,     7,    -1,     8,    -1,     9,    -1,
      10,    -1,    11,    -1,    36,    27,    17,    -1,    17,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    32,    32,    36,    41,    45,    52,    87,    93,    99,
     105,   114,   118,   122,   126,   130,   134,   138,   146,   150,
     151,   158,   169,   199,   203,   210,   218,   226,   244,   248,
     255,   261,   267,   273,   279,   285,   291,   297,   306,   323
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ppEND", "ppHOMVARGP", "ppVARGP",
  "ppPATHVAR", "ppVAR", "ppPARAM", "ppCONST", "ppFUNC", "ppSUBFUNC", "ppI",
  "ppPi", "ppNUMBER", "ppSCE", "ppPOW", "ppNAME", "'+'", "'-'", "'*'",
  "'/'", "ppNEG", "'^'", "'='", "'('", "')'", "','", "$accept", "line",
  "lval", "expression", "elem_expression", "function_call", "arguments",
  "declaration", "declist", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,    43,    45,
      42,    47,   273,    94,    61,    40,    41,    44
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    28,    29,    29,    29,    29,    30,    30,    30,    30,
      30,    31,    31,    31,    31,    31,    31,    31,    31,    31,
      31,    32,    32,    32,    32,    33,    33,    33,    34,    34,
      35,    35,    35,    35,    35,    35,    35,    35,    36,    36
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     1,     2,     3,     1,     1,     1,     1,
       1,     3,     3,     3,     3,     3,     2,     2,     3,     1,
       1,     1,     1,     1,     1,     4,     6,     4,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     3,    30,    31,    32,    33,    34,    35,    36,    37,
       9,    10,     7,     8,     6,     0,     0,     0,     1,     0,
      39,     4,    23,    24,    21,     0,     0,    22,     0,     0,
       0,     5,    19,    20,     0,     0,     0,     0,    16,    17,
       0,     0,     0,     0,     0,     0,    38,     0,     0,    29,
       0,    18,    11,    12,    13,    14,    15,    25,     0,    27,
       0,     0,    28,    26
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,    15,    16,    31,    32,    33,    50,    17,    21
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -29
static const yytype_int8 yypact[] =
{
      30,   -29,   -29,   -29,   -29,   -29,   -29,   -29,   -29,   -29,
     -29,   -29,   -29,   -29,   -29,     5,   -12,     9,   -29,     6,
     -29,     0,   -29,   -29,   -29,     3,     4,    21,     6,     6,
       6,    68,   -29,   -29,    36,     6,     6,     6,    25,    25,
      41,     6,     6,     6,     6,     6,   -29,    50,    31,    68,
     -16,   -29,   -17,   -17,    25,    25,    25,   -29,     6,   -29,
       6,    59,    68,   -29
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -29,   -29,   -29,   -28,   -29,   -29,   -29,   -29,   -29
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      38,    39,    40,    43,    44,    18,    45,    47,    48,    49,
      59,    60,    19,    52,    53,    54,    55,    56,    22,    23,
      24,    25,    26,    27,    28,    29,    20,    34,    35,    36,
      61,    30,    62,     1,     2,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    37,    14,    45,    41,
      42,    43,    44,    46,    45,     0,     0,     0,    58,    41,
      42,    43,    44,     0,    45,     0,     0,    51,    41,    42,
      43,    44,     0,    45,     0,     0,    57,    41,    42,    43,
      44,     0,    45,     0,     0,    63,    41,    42,    43,    44,
       0,    45
};

static const yytype_int8 yycheck[] =
{
      28,    29,    30,    20,    21,     0,    23,    35,    36,    37,
      26,    27,    24,    41,    42,    43,    44,    45,    12,    13,
      14,    15,    16,    17,    18,    19,    17,    27,    25,    25,
      58,    25,    60,     3,     4,     5,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    25,    17,    23,    18,
      19,    20,    21,    17,    23,    -1,    -1,    -1,    27,    18,
      19,    20,    21,    -1,    23,    -1,    -1,    26,    18,    19,
      20,    21,    -1,    23,    -1,    -1,    26,    18,    19,    20,
      21,    -1,    23,    -1,    -1,    26,    18,    19,    20,    21,
      -1,    23
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    17,    29,    30,    35,     0,    24,
      17,    36,    12,    13,    14,    15,    16,    17,    18,    19,
      25,    31,    32,    33,    27,    25,    25,    25,    31,    31,
      31,    18,    19,    20,    21,    23,    17,    31,    31,    31,
      34,    26,    31,    31,    31,    31,    31,    26,    27,    26,
      27,    31,    31,    26
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

/* Line 1455 of yacc.c  */
#line 32 "ppParse_bison.y"
    { // we are at the end of the file
          *endOfFile = 1;
          lastExpWasNumber = 0;
        }
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 37 "ppParse_bison.y"
    { // we are at the end of the file
          *endOfFile = 1;
          lastExpWasNumber = 0;
        }
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 42 "ppParse_bison.y"
    { // declares things
          lastExpWasNumber = 0;
        }
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 46 "ppParse_bison.y"
    { // defines things
          lastExpWasNumber = 0;
        }
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 53 "ppParse_bison.y"
    { // determine if it is in ppArray
     	  if (lookup_preProcessArray(&typeFound, &locFound, ppArray, (yyvsp[(1) - (1)].name)))
          { // found the item in the array - it must be a function that is not defined!
            if (!ppArray->types[typeFound].isDefined[locFound])
            { // declared but undefined
              if (typeFound == FUNCTIONTYPE || typeFound == PARAMETERTYPE || typeFound == CONSTANTTYPE)
  	      { // it is a function, parameter or constant - mark as defined
	        ppArray->types[typeFound].isDefined[locFound] = 1;
              }
              else
              { // not a function
                printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", (yyvsp[(1) - (1)].name), *lineNum);
                errorStatement();
              }
            }
	    else
            { // trying to redefine!
              printf("ERROR: Trying to redefine %s in system statement %d.\n", (yyvsp[(1) - (1)].name), *lineNum);
              errorStatement();
            }
          } 
          else
          { // add as an inline subfunction that is defined here
	    typeFound = INLINESUBFUNCTIONTYPE;
	    addEntry_preProcessArray_type(&ppArray->types[INLINESUBFUNCTIONTYPE], (yyvsp[(1) - (1)].name), 1, *lineNum);
          }
          // setup what we are defining
          definingType = typeFound;
          lastExpWasNumber = 0;

          // free memory
          free((yyvsp[(1) - (1)].name));
          (yyvsp[(1) - (1)].name) = NULL;
	}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 88 "ppParse_bison.y"
    { // invalid lval
	  printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", (yyvsp[(1) - (1)].name), *lineNum);
          errorStatement();
          lastExpWasNumber = 1;
        }
    break;

  case 8:

/* Line 1455 of yacc.c  */
#line 94 "ppParse_bison.y"
    { // invalid lval
          printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", (yyvsp[(1) - (1)].name), *lineNum);
    	  errorStatement();
          lastExpWasNumber = 0;
        }
    break;

  case 9:

/* Line 1455 of yacc.c  */
#line 100 "ppParse_bison.y"
    { // invalid lval
	  printf("ERROR: Invalid left-hand value (I) in system statement %d.\n", *lineNum);
    	  errorStatement();
          lastExpWasNumber = 0;
        }
    break;

  case 10:

/* Line 1455 of yacc.c  */
#line 106 "ppParse_bison.y"
    { // invalid lval
	  printf("ERROR: Invalid left-hand value (Pi) in system statement %d.\n", *lineNum);
    	  errorStatement();
          lastExpWasNumber = 0;
        }
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 115 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 119 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 123 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 127 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 131 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 135 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 139 "ppParse_bison.y"
    {
          if (ppAddAllNums && lastExpWasNumber)
          { // negate last number entered
            negate_lastEntry(&ppArray->types[NUMBERTYPE]);
          }
          lastExpWasNumber = 0;
        }
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 147 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 152 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 159 "ppParse_bison.y"
    { // determine if this number has been used before
          if (ppAddAllNums || !lookup_preProcessArray_type(&locFound, &ppArray->types[NUMBERTYPE], (yyvsp[(1) - (1)].name)))
          { // add it to the list as defined
      	    addEntry_preProcessArray_type(&ppArray->types[NUMBERTYPE], (yyvsp[(1) - (1)].name), 1, *lineNum);
          }
	  // free memory
	  free((yyvsp[(1) - (1)].name));
	  (yyvsp[(1) - (1)].name) = NULL;
          lastExpWasNumber = 1;
	}
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 170 "ppParse_bison.y"
    { // determine what it is
          if (!lookup_preProcessArray(&typeFound, &locFound, ppArray, (yyvsp[(1) - (1)].name)))
          { // not found
            printf("ERROR: Attempting to use the unknown symbol '%s' in system statement %d.\n", (yyvsp[(1) - (1)].name), *lineNum);
    	    errorStatement();
          }

	  // verify that is has been defined
	  if (!ppArray->types[typeFound].isDefined[locFound])
          { // not defined
            printf("ERROR: Attempting to use the symbol '%s' before it has been defined in system statement %d.\n", (yyvsp[(1) - (1)].name), *lineNum);
    	    errorStatement();
          }

          // verify it is okay
          verify_dependent_name(definingType, typeFound, (yyvsp[(1) - (1)].name), *lineNum);

          // verify not defining an expression based on itself
          if (ppArray->types[typeFound].lineNumber[locFound] == *lineNum)
          { // not defined
            printf("ERROR: Attempting to define the symbol '%s' in terms of itself in system statement %d.\n", (yyvsp[(1) - (1)].name), *lineNum);
            errorStatement();
          }

          // free memory
          free((yyvsp[(1) - (1)].name));
	  (yyvsp[(1) - (1)].name) = NULL;
          lastExpWasNumber = 0;
	}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 200 "ppParse_bison.y"
    { // I -- ignore
         lastExpWasNumber = 0;
        }
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 204 "ppParse_bison.y"
    { // Pi -- ignore
          lastExpWasNumber = 0;
        }
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 211 "ppParse_bison.y"
    { // sin, cos, tan, or exp
          lastExpWasNumber = 0;

          // free memory
          free((yyvsp[(1) - (4)].name));
          (yyvsp[(1) - (4)].name) = NULL;
        }
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 219 "ppParse_bison.y"
    { // pow
          lastExpWasNumber = 0;

          // free memory
          free((yyvsp[(1) - (6)].name));
	  (yyvsp[(1) - (6)].name) = NULL;
        }
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 227 "ppParse_bison.y"
    { // determine if this corresponds to a defined subfunction
          if (!lookup_preProcessArray_type(&locFound, &ppArray->types[DEFINEDSUBFUNCTIONTYPE], (yyvsp[(1) - (4)].name)))
          { // not found
            printf("ERROR: Invalid name of a defined subfunction (%s) in system statement %d.\n", (yyvsp[(1) - (4)].name), *lineNum);
    	    errorStatement();
          }

          // verify it is okay
          verify_dependent_name(definingType, DEFINEDSUBFUNCTIONTYPE, (yyvsp[(1) - (4)].name), *lineNum);

          // free memory
          free((yyvsp[(1) - (4)].name));
	  (yyvsp[(1) - (4)].name) = NULL;
          lastExpWasNumber = 0;
        }
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 245 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 249 "ppParse_bison.y"
    {
          lastExpWasNumber = 0;
        }
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 256 "ppParse_bison.y"
    { // setup declarationType
	  declarationType = HOMVARIABLEGROUPTYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
	}
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 262 "ppParse_bison.y"
    { // setup declarationType
          declarationType = VARIABLEGROUPTYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
	}
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 268 "ppParse_bison.y"
    { // setup declarationType
          declarationType = PATHVARIABLETYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
        }
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 274 "ppParse_bison.y"
    { // setup declarationType
          declarationType = VARIABLETYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
        }
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 280 "ppParse_bison.y"
    { // setup declarationType
          declarationType = PARAMETERTYPE;
	  isDefined = ppUsingParamHomotopy;
          lastExpWasNumber = 0;
        }
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 286 "ppParse_bison.y"
    { // setup declarationType
          declarationType = CONSTANTTYPE;
	  isDefined = 0;
          lastExpWasNumber = 0;
        }
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 292 "ppParse_bison.y"
    { // setup declarationType
          declarationType = FUNCTIONTYPE;
	  isDefined = 0;
          lastExpWasNumber = 0;
        }
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 298 "ppParse_bison.y"
    { // setup declarationType
          declarationType = DEFINEDSUBFUNCTIONTYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
        }
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 307 "ppParse_bison.y"
    { // determine if it is already in the list
          if (lookup_preProcessArray(&typeFound, &locFound, ppArray, (yyvsp[(3) - (3)].name)))
          { // redefining!
            printf("ERROR: Trying to redefine %s in system statement %d.\n", (yyvsp[(3) - (3)].name), *lineNum);
    	    errorStatement();
          }
          else
          { // add to list
    	    addEntry_preProcessArray_type(&ppArray->types[declarationType], (yyvsp[(3) - (3)].name), isDefined, *lineNum);
          }

          // free memory
          free((yyvsp[(3) - (3)].name));
	  (yyvsp[(3) - (3)].name) = NULL;
          lastExpWasNumber = 0;
	}
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 324 "ppParse_bison.y"
    { // determine if it is already in the list
          if (lookup_preProcessArray(&typeFound, &locFound, ppArray, (yyvsp[(1) - (1)].name)))
          { // redefining!
            printf("ERROR: Trying to redefine %s in system statement %d.\n", (yyvsp[(1) - (1)].name), *lineNum);
    	    errorStatement();
          }
          else
          { // add to list
            addEntry_preProcessArray_type(&ppArray->types[declarationType], (yyvsp[(1) - (1)].name), isDefined, *lineNum);
          }

          // free memory
          free((yyvsp[(1) - (1)].name));
	  (yyvsp[(1) - (1)].name) = NULL;
          lastExpWasNumber = 0;
	}
    break;



/* Line 1455 of yacc.c  */
#line 1900 "src/ppParse_bison.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 342 "ppParse_bison.y"


void verify_dependent_name(int defineType, int currType, char *name, int line)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify that 'defineType' can depend upon 'currType'    *
\***************************************************************/
{
  if (defineType == CONSTANTTYPE)
  { // must depend on other constants
    if (!(currType == CONSTANTTYPE || currType == NUMBERTYPE))
    { // not a constant
      printf("ERROR: Attempting to use the symbol '%s' when defining a constant in system statement %d.\nConstants must be defined using other constants or numbers.\n", name, line);
      errorStatement();
    }
  }
  else if (defineType == PARAMETERTYPE)
  { // must only depend on constants, pathvariables, parameters, or inline subfunction
    if (!(currType == PATHVARIABLETYPE || currType == CONSTANTTYPE || currType == NUMBERTYPE))
    { // not valid
      printf("ERROR: Attempting to use the symbol '%s' when defining a parameter in system statement %d.\nParameters must be defined using pathvariables, constants, or numbers.\n", name, line);
      errorStatement();
    }
  }
  else if (defineType == FUNCTIONTYPE)
  { // cannot depend on path variables or other functions
    if (currType == PATHVARIABLETYPE || currType == FUNCTIONTYPE)
    { // is a pathvariable
      printf("ERROR: Attempting to use the symbol '%s' when defining a function in system statement %d.\nFunctions cannot directly depend upon path variables or other functions.\n", name, line);
      errorStatement();
    }
  }

  return;
}


