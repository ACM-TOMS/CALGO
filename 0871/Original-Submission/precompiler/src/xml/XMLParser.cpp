#define YY_XMLParser_h_included

/*  A Bison++ parser, made from src/xml/xml.yacc  */

 /* with Bison++ version bison++ Version 1.21-8, adapted from GNU bison by coetmeur@icdc.fr
  */


#line 1 "/usr/local/lib/bison.cc"
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Bob Corbett and Richard Stallman

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 1, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

/* HEADER SECTION */
#if defined( _MSDOS ) || defined(MSDOS) || defined(__MSDOS__) 
#define __MSDOS_AND_ALIKE
#endif
#if defined(_WINDOWS) && defined(_MSC_VER)
#define __HAVE_NO_ALLOCA
#define __MSDOS_AND_ALIKE
#endif

#ifndef alloca
#if defined( __GNUC__)
#define alloca __builtin_alloca

#elif (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc)  || defined (__sgi)
#include <alloca.h>

#elif defined (__MSDOS_AND_ALIKE)
#include <malloc.h>
#ifndef __TURBOC__
/* MS C runtime lib */
#define alloca _alloca
#endif

#elif defined(_AIX)
#include <malloc.h>
#pragma alloca

#elif defined(__hpux)
#ifdef __cplusplus
extern "C" {
void *alloca (unsigned int);
};
#else /* not __cplusplus */
void *alloca ();
#endif /* not __cplusplus */

#endif /* not _AIX  not MSDOS, or __TURBOC__ or _AIX, not sparc.  */
#endif /* alloca not defined.  */
#ifdef c_plusplus
#ifndef __cplusplus
#define __cplusplus
#endif
#endif
#ifdef __cplusplus
#ifndef YY_USE_CLASS
#define YY_USE_CLASS
#endif
#else
#ifndef __STDC__
#define const
#endif
#endif
#include <stdio.h>
#define YYBISON 1  

/* #line 73 "/usr/local/lib/bison.cc" */
#line 85 "src/xml/XMLParser.cpp"
#line 15 "src/xml/xml.yacc"

  #include <stdio.h>
  #include <string> 
  #include <list>
  #include <iostream>
  #include <fstream>
  #include <FlexLexer.h>
  #include "xmlnode.h"
  #include "textnode.h"
  #include "elementnode.h"
  #include "pinode.h"
  #include "commentnode.h"
    

#line 35 "src/xml/xml.yacc"
typedef union {
    string*   str;
  } yy_XMLParser_stype;
#define YY_XMLParser_STYPE yy_XMLParser_stype
#define YY_XMLParser_MEMBERS  \
  virtual ~XMLParser();                       \
  void init();                                \
  int parse(ifstream* in);                    \
  int parse();                                \
  int errorLine();                            \
  int errorColumn();                          \
  int syntaxError(const string&);             \
  int semanticError(const string&);           \
  void warning(const string&);                \
  string charpToStr(char*);                   \
  string stripWhiteAndQuotes(const string&);  \
  string stripWhites(const string&);          \
  void checkVersion(const string&);           \
  void openTag();                             \
  void singleTag();                           \
  void closeTag(const string&);               \
  void showTagStack();                        \
  void newTextNode(const string&);            \
  void newWhiteTextNode(const string&);       \
  void setIncludeWhites(bool);                \
  yyFlexLexer xmlLexer;                       \
  bool bSemError;                             \
  bool bSynError;                             \
  bool bSingleTag;                            \
  bool bIncludeWhites;                        \
  list<string> tagStack;                      \
  bool emptyAllowed;                          \
  string tagName;                             \
  XMLNode *tree,*treePos;                     \
  vector<Attribute> attributes;        
#define YY_XMLParser_LEX_BODY  { return xmlLexer.yylex(); }
#define YY_XMLParser_ERROR_BODY  {  } 
#define YY_XMLParser_CONSTRUCTOR_PARAM 
#define YY_XMLParser_CONSTRUCTOR_CODE  { \
  init(); \
  bIncludeWhites=1; \
}

#line 73 "/usr/local/lib/bison.cc"
/* %{ and %header{ and %union, during decl */
#define YY_XMLParser_BISON 1
#ifndef YY_XMLParser_COMPATIBILITY
#ifndef YY_USE_CLASS
#define  YY_XMLParser_COMPATIBILITY 1
#else
#define  YY_XMLParser_COMPATIBILITY 0
#endif
#endif

#if YY_XMLParser_COMPATIBILITY != 0
/* backward compatibility */
#ifdef YYLTYPE
#ifndef YY_XMLParser_LTYPE
#define YY_XMLParser_LTYPE YYLTYPE
#endif
#endif
#ifdef YYSTYPE
#ifndef YY_XMLParser_STYPE 
#define YY_XMLParser_STYPE YYSTYPE
#endif
#endif
#ifdef YYDEBUG
#ifndef YY_XMLParser_DEBUG
#define  YY_XMLParser_DEBUG YYDEBUG
#endif
#endif
#ifdef YY_XMLParser_STYPE
#ifndef yystype
#define yystype YY_XMLParser_STYPE
#endif
#endif
/* use goto to be compatible */
#ifndef YY_XMLParser_USE_GOTO
#define YY_XMLParser_USE_GOTO 1
#endif
#endif

/* use no goto to be clean in C++ */
#ifndef YY_XMLParser_USE_GOTO
#define YY_XMLParser_USE_GOTO 0
#endif

#ifndef YY_XMLParser_PURE

/* #line 117 "/usr/local/lib/bison.cc" */
#line 192 "src/xml/XMLParser.cpp"

#line 117 "/usr/local/lib/bison.cc"
/*  YY_XMLParser_PURE */
#endif

/* section apres lecture def, avant lecture grammaire S2 */

/* #line 121 "/usr/local/lib/bison.cc" */
#line 201 "src/xml/XMLParser.cpp"

#line 121 "/usr/local/lib/bison.cc"
/* prefix */
#ifndef YY_XMLParser_DEBUG

/* #line 123 "/usr/local/lib/bison.cc" */
#line 208 "src/xml/XMLParser.cpp"

#line 123 "/usr/local/lib/bison.cc"
/* YY_XMLParser_DEBUG */
#endif


#ifndef YY_XMLParser_LSP_NEEDED

/* #line 128 "/usr/local/lib/bison.cc" */
#line 218 "src/xml/XMLParser.cpp"

#line 128 "/usr/local/lib/bison.cc"
 /* YY_XMLParser_LSP_NEEDED*/
#endif



/* DEFAULT LTYPE*/
#ifdef YY_XMLParser_LSP_NEEDED
#ifndef YY_XMLParser_LTYPE
typedef
  struct yyltype
    {
      int timestamp;
      int first_line;
      int first_column;
      int last_line;
      int last_column;
      char *text;
   }
  yyltype;

#define YY_XMLParser_LTYPE yyltype
#endif
#endif
/* DEFAULT STYPE*/
      /* We used to use `unsigned long' as YY_XMLParser_STYPE on MSDOS,
	 but it seems better to be consistent.
	 Most programs should declare their own type anyway.  */

#ifndef YY_XMLParser_STYPE
#define YY_XMLParser_STYPE int
#endif
/* DEFAULT MISCELANEOUS */
#ifndef YY_XMLParser_PARSE
#define YY_XMLParser_PARSE yyparse
#endif
#ifndef YY_XMLParser_LEX
#define YY_XMLParser_LEX yylex
#endif
#ifndef YY_XMLParser_LVAL
#define YY_XMLParser_LVAL yylval
#endif
#ifndef YY_XMLParser_LLOC
#define YY_XMLParser_LLOC yylloc
#endif
#ifndef YY_XMLParser_CHAR
#define YY_XMLParser_CHAR yychar
#endif
#ifndef YY_XMLParser_NERRS
#define YY_XMLParser_NERRS yynerrs
#endif
#ifndef YY_XMLParser_DEBUG_FLAG
#define YY_XMLParser_DEBUG_FLAG yydebug
#endif
#ifndef YY_XMLParser_ERROR
#define YY_XMLParser_ERROR yyerror
#endif
#ifndef YY_XMLParser_PARSE_PARAM
#ifndef __STDC__
#ifndef __cplusplus
#ifndef YY_USE_CLASS
#define YY_XMLParser_PARSE_PARAM
#ifndef YY_XMLParser_PARSE_PARAM_DEF
#define YY_XMLParser_PARSE_PARAM_DEF
#endif
#endif
#endif
#endif
#ifndef YY_XMLParser_PARSE_PARAM
#define YY_XMLParser_PARSE_PARAM void
#endif
#endif
#if YY_XMLParser_COMPATIBILITY != 0
/* backward compatibility */
#ifdef YY_XMLParser_LTYPE
#ifndef YYLTYPE
#define YYLTYPE YY_XMLParser_LTYPE
#else
/* WARNING obsolete !!! user defined YYLTYPE not reported into generated header */
#endif
#endif
#ifndef YYSTYPE
#define YYSTYPE YY_XMLParser_STYPE
#else
/* WARNING obsolete !!! user defined YYSTYPE not reported into generated header */
#endif
#ifdef YY_XMLParser_PURE
#ifndef YYPURE
#define YYPURE YY_XMLParser_PURE
#endif
#endif
#ifdef YY_XMLParser_DEBUG
#ifndef YYDEBUG
#define YYDEBUG YY_XMLParser_DEBUG 
#endif
#endif
#ifndef YY_XMLParser_ERROR_VERBOSE
#ifdef YYERROR_VERBOSE
#define YY_XMLParser_ERROR_VERBOSE YYERROR_VERBOSE
#endif
#endif
#ifndef YY_XMLParser_LSP_NEEDED
#ifdef YYLSP_NEEDED
#define YY_XMLParser_LSP_NEEDED YYLSP_NEEDED
#endif
#endif
#endif
#ifndef YY_USE_CLASS
/* TOKEN C */

/* #line 236 "/usr/local/lib/bison.cc" */
#line 331 "src/xml/XMLParser.cpp"
#define	PROLOG_START	258
#define	VERSION	259
#define	EQ	260
#define	NR	261
#define	PROLOG_END	262
#define	DOCTYPE	263
#define	COMMENT_BEGIN	264
#define	COMMENT_END	265
#define	COMMENT_DATA	266
#define	OPEN_TAG_BEGIN	267
#define	OPEN_TAG_END	268
#define	ATTRIBUTE	269
#define	SINGLE_TAG_END	270
#define	PROCESS_TAG_TARGET	271
#define	PROCESS_TAG_DATA	272
#define	PROCESS_TAG_END	273
#define	CLOSE_TAG	274
#define	CLOSE_TAG_END	275
#define	COMMENT	276
#define	TEXTNODE	277
#define	WHITESPACE	278


#line 236 "/usr/local/lib/bison.cc"
 /* #defines tokens */
#else
/* CLASS */
#ifndef YY_XMLParser_CLASS
#define YY_XMLParser_CLASS XMLParser
#endif
#ifndef YY_XMLParser_INHERIT
#define YY_XMLParser_INHERIT
#endif
#ifndef YY_XMLParser_MEMBERS
#define YY_XMLParser_MEMBERS 
#endif
#ifndef YY_XMLParser_LEX_BODY
#define YY_XMLParser_LEX_BODY  
#endif
#ifndef YY_XMLParser_ERROR_BODY
#define YY_XMLParser_ERROR_BODY  
#endif
#ifndef YY_XMLParser_CONSTRUCTOR_PARAM
#define YY_XMLParser_CONSTRUCTOR_PARAM
#endif
#ifndef YY_XMLParser_CONSTRUCTOR_CODE
#define YY_XMLParser_CONSTRUCTOR_CODE
#endif
#ifndef YY_XMLParser_CONSTRUCTOR_INIT
#define YY_XMLParser_CONSTRUCTOR_INIT
#endif
/* choose between enum and const */
#ifndef YY_XMLParser_USE_CONST_TOKEN
#define YY_XMLParser_USE_CONST_TOKEN 0
/* yes enum is more compatible with flex,  */
/* so by default we use it */ 
#endif
#if YY_XMLParser_USE_CONST_TOKEN != 0
#ifndef YY_XMLParser_ENUM_TOKEN
#define YY_XMLParser_ENUM_TOKEN yy_XMLParser_enum_token
#endif
#endif

class YY_XMLParser_CLASS YY_XMLParser_INHERIT
{
public: 
#if YY_XMLParser_USE_CONST_TOKEN != 0
/* static const int token ... */

/* #line 280 "/usr/local/lib/bison.cc" */
#line 402 "src/xml/XMLParser.cpp"
static const int PROLOG_START;
static const int VERSION;
static const int EQ;
static const int NR;
static const int PROLOG_END;
static const int DOCTYPE;
static const int COMMENT_BEGIN;
static const int COMMENT_END;
static const int COMMENT_DATA;
static const int OPEN_TAG_BEGIN;
static const int OPEN_TAG_END;
static const int ATTRIBUTE;
static const int SINGLE_TAG_END;
static const int PROCESS_TAG_TARGET;
static const int PROCESS_TAG_DATA;
static const int PROCESS_TAG_END;
static const int CLOSE_TAG;
static const int CLOSE_TAG_END;
static const int COMMENT;
static const int TEXTNODE;
static const int WHITESPACE;


#line 280 "/usr/local/lib/bison.cc"
 /* decl const */
#else
enum YY_XMLParser_ENUM_TOKEN { YY_XMLParser_NULL_TOKEN=0

/* #line 283 "/usr/local/lib/bison.cc" */
#line 432 "src/xml/XMLParser.cpp"
	,PROLOG_START=258
	,VERSION=259
	,EQ=260
	,NR=261
	,PROLOG_END=262
	,DOCTYPE=263
	,COMMENT_BEGIN=264
	,COMMENT_END=265
	,COMMENT_DATA=266
	,OPEN_TAG_BEGIN=267
	,OPEN_TAG_END=268
	,ATTRIBUTE=269
	,SINGLE_TAG_END=270
	,PROCESS_TAG_TARGET=271
	,PROCESS_TAG_DATA=272
	,PROCESS_TAG_END=273
	,CLOSE_TAG=274
	,CLOSE_TAG_END=275
	,COMMENT=276
	,TEXTNODE=277
	,WHITESPACE=278


#line 283 "/usr/local/lib/bison.cc"
 /* enum token */
     }; /* end of enum declaration */
#endif
public:
 int YY_XMLParser_PARSE (YY_XMLParser_PARSE_PARAM);
 virtual void YY_XMLParser_ERROR(char *msg) YY_XMLParser_ERROR_BODY;
#ifdef YY_XMLParser_PURE
#ifdef YY_XMLParser_LSP_NEEDED
 virtual int  YY_XMLParser_LEX (YY_XMLParser_STYPE *YY_XMLParser_LVAL,YY_XMLParser_LTYPE *YY_XMLParser_LLOC) YY_XMLParser_LEX_BODY;
#else
 virtual int  YY_XMLParser_LEX (YY_XMLParser_STYPE *YY_XMLParser_LVAL) YY_XMLParser_LEX_BODY;
#endif
#else
 virtual int YY_XMLParser_LEX() YY_XMLParser_LEX_BODY;
 YY_XMLParser_STYPE YY_XMLParser_LVAL;
#ifdef YY_XMLParser_LSP_NEEDED
 YY_XMLParser_LTYPE YY_XMLParser_LLOC;
#endif
 int   YY_XMLParser_NERRS;
 int    YY_XMLParser_CHAR;
#endif
#if YY_XMLParser_DEBUG != 0
 int YY_XMLParser_DEBUG_FLAG;   /*  nonzero means print parse trace     */
#endif
public:
 YY_XMLParser_CLASS(YY_XMLParser_CONSTRUCTOR_PARAM);
public:
 YY_XMLParser_MEMBERS 
};
/* other declare folow */
#if YY_XMLParser_USE_CONST_TOKEN != 0

/* #line 314 "/usr/local/lib/bison.cc" */
#line 490 "src/xml/XMLParser.cpp"
const int YY_XMLParser_CLASS::PROLOG_START=258;
const int YY_XMLParser_CLASS::VERSION=259;
const int YY_XMLParser_CLASS::EQ=260;
const int YY_XMLParser_CLASS::NR=261;
const int YY_XMLParser_CLASS::PROLOG_END=262;
const int YY_XMLParser_CLASS::DOCTYPE=263;
const int YY_XMLParser_CLASS::COMMENT_BEGIN=264;
const int YY_XMLParser_CLASS::COMMENT_END=265;
const int YY_XMLParser_CLASS::COMMENT_DATA=266;
const int YY_XMLParser_CLASS::OPEN_TAG_BEGIN=267;
const int YY_XMLParser_CLASS::OPEN_TAG_END=268;
const int YY_XMLParser_CLASS::ATTRIBUTE=269;
const int YY_XMLParser_CLASS::SINGLE_TAG_END=270;
const int YY_XMLParser_CLASS::PROCESS_TAG_TARGET=271;
const int YY_XMLParser_CLASS::PROCESS_TAG_DATA=272;
const int YY_XMLParser_CLASS::PROCESS_TAG_END=273;
const int YY_XMLParser_CLASS::CLOSE_TAG=274;
const int YY_XMLParser_CLASS::CLOSE_TAG_END=275;
const int YY_XMLParser_CLASS::COMMENT=276;
const int YY_XMLParser_CLASS::TEXTNODE=277;
const int YY_XMLParser_CLASS::WHITESPACE=278;


#line 314 "/usr/local/lib/bison.cc"
 /* const YY_XMLParser_CLASS::token */
#endif
/*apres const  */
YY_XMLParser_CLASS::YY_XMLParser_CLASS(YY_XMLParser_CONSTRUCTOR_PARAM) YY_XMLParser_CONSTRUCTOR_INIT
{
#if YY_XMLParser_DEBUG != 0
YY_XMLParser_DEBUG_FLAG=0;
#endif
YY_XMLParser_CONSTRUCTOR_CODE;
};
#endif

/* #line 325 "/usr/local/lib/bison.cc" */
#line 528 "src/xml/XMLParser.cpp"


#define	YYFINAL		49
#define	YYFLAG		-32768
#define	YYNTBASE	24

#define YYTRANSLATE(x) ((unsigned)(x) <= 278 ? yytranslate[x] : 35)

static const char yytranslate[] = {     0,
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
     2,     2,     2,     2,     2,     1,     2,     3,     4,     5,
     6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
    16,    17,    18,    19,    20,    21,    22,    23
};

#if YY_XMLParser_DEBUG != 0
static const short yyprhs[] = {     0,
     0,     2,     5,     8,    12,    14,    17,    23,    30,    33,
    35,    38,    40,    42,    45,    47,    49,    51,    54,    56,
    58,    61,    65,    67,    70,    74,    76
};

static const short yyrhs[] = {    27,
     0,    26,    27,     0,    25,    27,     0,    26,    25,    27,
     0,     8,     0,    23,    25,     0,     3,     4,     5,     6,
     7,     0,     3,     4,     5,     6,    32,     7,     0,    23,
    26,     0,    28,     0,    28,    27,     0,    29,     0,    33,
     0,    19,    20,     0,    22,     0,    31,     0,    23,     0,
    12,    30,     0,    13,     0,    15,     0,    14,    30,     0,
    16,    32,    18,     0,    17,     0,    32,    17,     0,     9,
    34,    10,     0,    11,     0,    34,    11,     0
};

#endif

#if YY_XMLParser_DEBUG != 0
static const short yyrline[] = { 0,
   141,   142,   143,   144,   150,   151,   154,   159,   165,   170,
   171,   175,   178,   179,   183,   187,   188,   196,   205,   208,
   211,   220,   234,   238,   246,   259,   263
};

static const char * const yytname[] = {   "$","error","$illegal.","PROLOG_START",
"VERSION","EQ","NR","PROLOG_END","DOCTYPE","COMMENT_BEGIN","COMMENT_END","COMMENT_DATA",
"OPEN_TAG_BEGIN","OPEN_TAG_END","ATTRIBUTE","SINGLE_TAG_END","PROCESS_TAG_TARGET",
"PROCESS_TAG_DATA","PROCESS_TAG_END","CLOSE_TAG","CLOSE_TAG_END","COMMENT","TEXTNODE",
"WHITESPACE","xml_doc","doctype","prolog","internal_parts","part","open_tag",
"open_tag_end","process_tag","process_data","comment","comment_text",""
};
#endif

static const short yyr1[] = {     0,
    24,    24,    24,    24,    25,    25,    26,    26,    26,    27,
    27,    28,    28,    28,    28,    28,    28,    29,    30,    30,
    30,    31,    32,    32,    33,    34,    34
};

static const short yyr2[] = {     0,
     1,     2,     2,     3,     1,     2,     5,     6,     2,     1,
     2,     1,     1,     2,     1,     1,     1,     2,     1,     1,
     2,     3,     1,     2,     3,     1,     2
};

static const short yydefact[] = {     0,
     0,     5,     0,     0,     0,     0,    15,    17,     0,     0,
     1,    10,    12,    16,    13,     0,    26,     0,    19,     0,
    20,    18,    23,     0,    14,     0,     6,     9,    17,     3,
    17,     0,     2,    11,     0,    25,    27,    21,    24,    22,
     0,     4,     0,     7,     0,     8,     0,     0,     0
};

static const short yydefgoto[] = {    47,
    27,    28,    11,    12,    13,    22,    14,    24,    15,    18
};

static const short yypact[] = {    -1,
    26,-32768,    21,    32,    -7,    23,-32768,    11,    19,    17,
-32768,    19,-32768,-32768,-32768,    45,-32768,     2,-32768,    32,
-32768,-32768,-32768,    31,-32768,     1,-32768,-32768,-32768,-32768,
    29,    19,-32768,-32768,    46,-32768,-32768,-32768,-32768,-32768,
    -3,-32768,    10,-32768,    27,-32768,    51,    53,-32768
};

static const short yypgoto[] = {-32768,
     6,    54,    -9,-32768,-32768,    35,-32768,    13,-32768,-32768
};


#define	YYLAST		56


static const short yytable[] = {    30,
    33,     1,    34,     1,     2,     9,     2,     3,     2,    23,
     4,    36,    37,     1,     5,    32,    44,     6,     2,    41,
     7,     8,    42,    26,     2,     3,    23,     3,     4,    16,
     4,    17,     5,    46,     5,     6,     2,     6,     7,    31,
     7,    29,    25,    39,    19,    20,    21,    39,    40,    35,
    48,    43,    49,    10,    38,    45
};

static const short yycheck[] = {     9,
    10,     3,    12,     3,     8,     0,     8,     9,     8,    17,
    12,    10,    11,     3,    16,    10,     7,    19,     8,    23,
    22,    23,    32,    23,     8,     9,    17,     9,    12,     4,
    12,    11,    16,     7,    16,    19,     8,    19,    22,    23,
    22,    23,    20,    17,    13,    14,    15,    17,    18,     5,
     0,     6,     0,     0,    20,    43
};

#line 325 "/usr/local/lib/bison.cc"
 /* fattrs + tables */

/* parser code folow  */


/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

/* Note: dollar marks section change
   the next  is replaced by the list of actions, each action
   as one case of the switch.  */ 

#if YY_XMLParser_USE_GOTO != 0
/* 
 SUPRESSION OF GOTO : on some C++ compiler (sun c++)
  the goto is strictly forbidden if any constructor/destructor
  is used in the whole function (very stupid isn't it ?)
 so goto are to be replaced with a 'while/switch/case construct'
 here are the macro to keep some apparent compatibility
*/
#define YYGOTO(lb) {yy_gotostate=lb;continue;}
#define YYBEGINGOTO  enum yy_labels yy_gotostate=yygotostart; \
                     for(;;) switch(yy_gotostate) { case yygotostart: {
#define YYLABEL(lb) } case lb: {
#define YYENDGOTO } } 
#define YYBEGINDECLARELABEL enum yy_labels {yygotostart
#define YYDECLARELABEL(lb) ,lb
#define YYENDDECLARELABEL  };
#else
/* macro to keep goto */
#define YYGOTO(lb) goto lb
#define YYBEGINGOTO 
#define YYLABEL(lb) lb:
#define YYENDGOTO
#define YYBEGINDECLARELABEL 
#define YYDECLARELABEL(lb)
#define YYENDDECLARELABEL 
#endif
/* LABEL DECLARATION */
YYBEGINDECLARELABEL
  YYDECLARELABEL(yynewstate)
  YYDECLARELABEL(yybackup)
/* YYDECLARELABEL(yyresume) */
  YYDECLARELABEL(yydefault)
  YYDECLARELABEL(yyreduce)
  YYDECLARELABEL(yyerrlab)   /* here on detecting error */
  YYDECLARELABEL(yyerrlab1)   /* here on error raised explicitly by an action */
  YYDECLARELABEL(yyerrdefault)  /* current state does not do anything special for the error token. */
  YYDECLARELABEL(yyerrpop)   /* pop the current state because it cannot handle the error token */
  YYDECLARELABEL(yyerrhandle)  
YYENDDECLARELABEL
/* ALLOCA SIMULATION */
/* __HAVE_NO_ALLOCA */
#ifdef __HAVE_NO_ALLOCA
int __alloca_free_ptr(char *ptr,char *ref)
{if(ptr!=ref) free(ptr);
 return 0;}

#define __ALLOCA_alloca(size) malloc(size)
#define __ALLOCA_free(ptr,ref) __alloca_free_ptr((char *)ptr,(char *)ref)

#ifdef YY_XMLParser_LSP_NEEDED
#define __ALLOCA_return(num) \
            return( __ALLOCA_free(yyss,yyssa)+\
		    __ALLOCA_free(yyvs,yyvsa)+\
		    __ALLOCA_free(yyls,yylsa)+\
		   (num))
#else
#define __ALLOCA_return(num) \
            return( __ALLOCA_free(yyss,yyssa)+\
		    __ALLOCA_free(yyvs,yyvsa)+\
		   (num))
#endif
#else
#define __ALLOCA_return(num) return(num)
#define __ALLOCA_alloca(size) alloca(size)
#define __ALLOCA_free(ptr,ref) 
#endif

/* ENDALLOCA SIMULATION */

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (YY_XMLParser_CHAR = YYEMPTY)
#define YYEMPTY         -2
#define YYEOF           0
#define YYACCEPT        __ALLOCA_return(0)
#define YYABORT         __ALLOCA_return(1)
#define YYERROR         YYGOTO(yyerrlab1)
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL          YYGOTO(yyerrlab)
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do                                                              \
  if (YY_XMLParser_CHAR == YYEMPTY && yylen == 1)                               \
    { YY_XMLParser_CHAR = (token), YY_XMLParser_LVAL = (value);                 \
      yychar1 = YYTRANSLATE (YY_XMLParser_CHAR);                                \
      YYPOPSTACK;                                               \
      YYGOTO(yybackup);                                            \
    }                                                           \
  else                                                          \
    { YY_XMLParser_ERROR ("syntax error: cannot back up"); YYERROR; }   \
while (0)

#define YYTERROR        1
#define YYERRCODE       256

#ifndef YY_XMLParser_PURE
/* UNPURE */
#define YYLEX           YY_XMLParser_LEX()
#ifndef YY_USE_CLASS
/* If nonreentrant, and not class , generate the variables here */
int     YY_XMLParser_CHAR;                      /*  the lookahead symbol        */
YY_XMLParser_STYPE      YY_XMLParser_LVAL;              /*  the semantic value of the */
				/*  lookahead symbol    */
int YY_XMLParser_NERRS;                 /*  number of parse errors so far */
#ifdef YY_XMLParser_LSP_NEEDED
YY_XMLParser_LTYPE YY_XMLParser_LLOC;   /*  location data for the lookahead     */
			/*  symbol                              */
#endif
#endif


#else
/* PURE */
#ifdef YY_XMLParser_LSP_NEEDED
#define YYLEX           YY_XMLParser_LEX(&YY_XMLParser_LVAL, &YY_XMLParser_LLOC)
#else
#define YYLEX           YY_XMLParser_LEX(&YY_XMLParser_LVAL)
#endif
#endif
#ifndef YY_USE_CLASS
#if YY_XMLParser_DEBUG != 0
int YY_XMLParser_DEBUG_FLAG;                    /*  nonzero means print parse trace     */
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif
#endif



/*  YYINITDEPTH indicates the initial size of the parser's stacks       */

#ifndef YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif


#if __GNUC__ > 1                /* GNU C and GNU C++ define this.  */
#define __yy_bcopy(FROM,TO,COUNT)       __builtin_memcpy(TO,FROM,COUNT)
#else                           /* not GNU C or C++ */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */

#ifdef __cplusplus
static void __yy_bcopy (char *from, char *to, int count)
#else
#ifdef __STDC__
static void __yy_bcopy (char *from, char *to, int count)
#else
static void __yy_bcopy (from, to, count)
     char *from;
     char *to;
     int count;
#endif
#endif
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}
#endif

int
#ifdef YY_USE_CLASS
 YY_XMLParser_CLASS::
#endif
     YY_XMLParser_PARSE(YY_XMLParser_PARSE_PARAM)
#ifndef __STDC__
#ifndef __cplusplus
#ifndef YY_USE_CLASS
/* parameter definition without protypes */
YY_XMLParser_PARSE_PARAM_DEF
#endif
#endif
#endif
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YY_XMLParser_STYPE *yyvsp;
  int yyerrstatus;      /*  number of tokens to shift before error messages enabled */
  int yychar1=0;          /*  lookahead token as an internal (translated) token number */

  short yyssa[YYINITDEPTH];     /*  the state stack                     */
  YY_XMLParser_STYPE yyvsa[YYINITDEPTH];        /*  the semantic value stack            */

  short *yyss = yyssa;          /*  refer to the stacks thru separate pointers */
  YY_XMLParser_STYPE *yyvs = yyvsa;     /*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YY_XMLParser_LSP_NEEDED
  YY_XMLParser_LTYPE yylsa[YYINITDEPTH];        /*  the location stack                  */
  YY_XMLParser_LTYPE *yyls = yylsa;
  YY_XMLParser_LTYPE *yylsp;

#define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;

#ifdef YY_XMLParser_PURE
  int YY_XMLParser_CHAR;
  YY_XMLParser_STYPE YY_XMLParser_LVAL;
  int YY_XMLParser_NERRS;
#ifdef YY_XMLParser_LSP_NEEDED
  YY_XMLParser_LTYPE YY_XMLParser_LLOC;
#endif
#endif

  YY_XMLParser_STYPE yyval;             /*  the variable used to return         */
				/*  semantic values from the action     */
				/*  routines                            */

  int yylen;
/* start loop, in which YYGOTO may be used. */
YYBEGINGOTO

#if YY_XMLParser_DEBUG != 0
  if (YY_XMLParser_DEBUG_FLAG)
    fprintf(stderr, "Starting parse\n");
#endif
  yystate = 0;
  yyerrstatus = 0;
  YY_XMLParser_NERRS = 0;
  YY_XMLParser_CHAR = YYEMPTY;          /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YY_XMLParser_LSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
YYLABEL(yynewstate)

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YY_XMLParser_STYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YY_XMLParser_LSP_NEEDED
      YY_XMLParser_LTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef YY_XMLParser_LSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if yyoverflow is a macro.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yyls1, size * sizeof (*yylsp),
		 &yystacksize);
#else
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yystacksize);
#endif

      yyss = yyss1; yyvs = yyvs1;
#ifdef YY_XMLParser_LSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  YY_XMLParser_ERROR("parser stack overflow");
	  __ALLOCA_return(2);
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
      yyss = (short *) __ALLOCA_alloca (yystacksize * sizeof (*yyssp));
      __yy_bcopy ((char *)yyss1, (char *)yyss, size * sizeof (*yyssp));
      __ALLOCA_free(yyss1,yyssa);
      yyvs = (YY_XMLParser_STYPE *) __ALLOCA_alloca (yystacksize * sizeof (*yyvsp));
      __yy_bcopy ((char *)yyvs1, (char *)yyvs, size * sizeof (*yyvsp));
      __ALLOCA_free(yyvs1,yyvsa);
#ifdef YY_XMLParser_LSP_NEEDED
      yyls = (YY_XMLParser_LTYPE *) __ALLOCA_alloca (yystacksize * sizeof (*yylsp));
      __yy_bcopy ((char *)yyls1, (char *)yyls, size * sizeof (*yylsp));
      __ALLOCA_free(yyls1,yylsa);
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YY_XMLParser_LSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YY_XMLParser_DEBUG != 0
      if (YY_XMLParser_DEBUG_FLAG)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YY_XMLParser_DEBUG != 0
  if (YY_XMLParser_DEBUG_FLAG)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  YYGOTO(yybackup);
YYLABEL(yybackup)

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* YYLABEL(yyresume) */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    YYGOTO(yydefault);

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (YY_XMLParser_CHAR == YYEMPTY)
    {
#if YY_XMLParser_DEBUG != 0
      if (YY_XMLParser_DEBUG_FLAG)
	fprintf(stderr, "Reading a token: ");
#endif
      YY_XMLParser_CHAR = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (YY_XMLParser_CHAR <= 0)           /* This means end of input. */
    {
      yychar1 = 0;
      YY_XMLParser_CHAR = YYEOF;                /* Don't call YYLEX any more */

#if YY_XMLParser_DEBUG != 0
      if (YY_XMLParser_DEBUG_FLAG)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(YY_XMLParser_CHAR);

#if YY_XMLParser_DEBUG != 0
      if (YY_XMLParser_DEBUG_FLAG)
	{
	  fprintf (stderr, "Next token is %d (%s", YY_XMLParser_CHAR, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef YYPRINT
	  YYPRINT (stderr, YY_XMLParser_CHAR, YY_XMLParser_LVAL);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    YYGOTO(yydefault);

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	YYGOTO(yyerrlab);
      yyn = -yyn;
      YYGOTO(yyreduce);
    }
  else if (yyn == 0)
    YYGOTO(yyerrlab);

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YY_XMLParser_DEBUG != 0
  if (YY_XMLParser_DEBUG_FLAG)
    fprintf(stderr, "Shifting token %d (%s), ", YY_XMLParser_CHAR, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (YY_XMLParser_CHAR != YYEOF)
    YY_XMLParser_CHAR = YYEMPTY;

  *++yyvsp = YY_XMLParser_LVAL;
#ifdef YY_XMLParser_LSP_NEEDED
  *++yylsp = YY_XMLParser_LLOC;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  YYGOTO(yynewstate);

/* Do the default action for the current state.  */
YYLABEL(yydefault)

  yyn = yydefact[yystate];
  if (yyn == 0)
    YYGOTO(yyerrlab);

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
YYLABEL(yyreduce)
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YY_XMLParser_DEBUG != 0
  if (YY_XMLParser_DEBUG_FLAG)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


/* #line 811 "/usr/local/lib/bison.cc" */
#line 1151 "src/xml/XMLParser.cpp"

  switch (yyn) {

case 1:
#line 141 "src/xml/xml.yacc"
{;
    break;}
case 2:
#line 142 "src/xml/xml.yacc"
{;
    break;}
case 3:
#line 143 "src/xml/xml.yacc"
{;
    break;}
case 4:
#line 144 "src/xml/xml.yacc"
{;
    break;}
case 5:
#line 150 "src/xml/xml.yacc"
{;
    break;}
case 6:
#line 151 "src/xml/xml.yacc"
{;
    break;}
case 7:
#line 154 "src/xml/xml.yacc"
{
                                            checkVersion(*yyvsp[-1].str);
                                            delete yyvsp[-1].str; //NR is string*
                                          ;
    break;}
case 8:
#line 159 "src/xml/xml.yacc"
{
                                                  checkVersion(*yyvsp[-2].str);
                                                  //warning("Unused prolog data:"+*$5);
                                                  delete yyvsp[-2].str;
                                                  delete yyvsp[-1].str;
                                                ;
    break;}
case 9:
#line 165 "src/xml/xml.yacc"
{;
    break;}
case 10:
#line 170 "src/xml/xml.yacc"
{;
    break;}
case 11:
#line 171 "src/xml/xml.yacc"
{;
    break;}
case 12:
#line 175 "src/xml/xml.yacc"
{
                                openTag();
                              ;
    break;}
case 13:
#line 178 "src/xml/xml.yacc"
{;
    break;}
case 14:
#line 179 "src/xml/xml.yacc"
{
                                closeTag(*yyvsp[-1].str);
                                delete(yyvsp[-1].str);
                              ;
    break;}
case 15:
#line 183 "src/xml/xml.yacc"
{
                                newTextNode(*yyvsp[0].str);
                                delete yyvsp[0].str;
                              ;
    break;}
case 16:
#line 187 "src/xml/xml.yacc"
{;
    break;}
case 17:
#line 188 "src/xml/xml.yacc"
{
                                if(bIncludeWhites) newWhiteTextNode(*yyvsp[0].str);
                                delete yyvsp[0].str;
                              ;
    break;}
case 18:
#line 196 "src/xml/xml.yacc"
{
                                    tagName=string(*yyvsp[-1].str);
                                    delete yyvsp[-1].str;
                                    
                                  ;
    break;}
case 19:
#line 205 "src/xml/xml.yacc"
{
                                bSingleTag=false;
                              ;
    break;}
case 20:
#line 208 "src/xml/xml.yacc"
{
                                bSingleTag=true;
                              ;
    break;}
case 21:
#line 211 "src/xml/xml.yacc"
{
                                Attribute a(*yyvsp[-1].str);
                                attributes.push_back(a);
                                delete yyvsp[-1].str; 
                              ;
    break;}
case 22:
#line 220 "src/xml/xml.yacc"
{
                                                          PINode* pNode=new PINode(*yyvsp[-2].str,*yyvsp[-1].str);
                                                          if(treePos==0){
                                                            semanticError("You need to define a root tag before using a processing instruction");
                                                          }
                                                          else{
                                                            treePos->appendChild(pNode);
                                                          }
                                                          delete yyvsp[-2].str;
                                                          delete yyvsp[-1].str;
                                                        ;
    break;}
case 23:
#line 234 "src/xml/xml.yacc"
{
                        yyval.str=new string(*yyvsp[0].str);
                        delete yyvsp[0].str;
                     ;
    break;}
case 24:
#line 238 "src/xml/xml.yacc"
{
                                    yyval.str=new string(*yyvsp[-1].str+*yyvsp[0].str);
                                    delete yyvsp[-1].str;
                                    delete yyvsp[0].str;
                                   ;
    break;}
case 25:
#line 246 "src/xml/xml.yacc"
{
                                              CommentNode* cNode=new CommentNode(*yyvsp[-1].str);
                                              if(treePos==0){
                                                semanticError("no comments allowed outside of a root");
                                              }
                                              else{
                                                treePos->appendChild(cNode);
                                              }
                                              delete yyvsp[-1].str;
                                            ;
    break;}
case 26:
#line 259 "src/xml/xml.yacc"
{
                                    yyval.str=new string(*yyvsp[0].str);
                                    delete yyvsp[0].str;
                                  ;
    break;}
case 27:
#line 263 "src/xml/xml.yacc"
{
                                    yyval.str=new string(*yyvsp[-1].str+*yyvsp[0].str);
                                    delete yyvsp[-1].str;
                                    delete yyvsp[0].str;
                                  ;
    break;}
}

#line 811 "/usr/local/lib/bison.cc"
   /* the action file gets copied in in place of this dollarsign  */
  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YY_XMLParser_LSP_NEEDED
  yylsp -= yylen;
#endif

#if YY_XMLParser_DEBUG != 0
  if (YY_XMLParser_DEBUG_FLAG)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YY_XMLParser_LSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = YY_XMLParser_LLOC.first_line;
      yylsp->first_column = YY_XMLParser_LLOC.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  YYGOTO(yynewstate);

YYLABEL(yyerrlab)   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++YY_XMLParser_NERRS;

#ifdef YY_XMLParser_ERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
	  for (x = (yyn < 0 ? -yyn : 0);
	       x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (yyn < 0 ? -yyn : 0);
		       x < (sizeof(yytname) / sizeof(char *)); x++)
		    if (yycheck[x + yyn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, yytname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      YY_XMLParser_ERROR(msg);
	      free(msg);
	    }
	  else
	    YY_XMLParser_ERROR ("parse error; also virtual memory exceeded");
	}
      else
#endif /* YY_XMLParser_ERROR_VERBOSE */
	YY_XMLParser_ERROR("parse error");
    }

  YYGOTO(yyerrlab1);
YYLABEL(yyerrlab1)   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (YY_XMLParser_CHAR == YYEOF)
	YYABORT;

#if YY_XMLParser_DEBUG != 0
      if (YY_XMLParser_DEBUG_FLAG)
	fprintf(stderr, "Discarding token %d (%s).\n", YY_XMLParser_CHAR, yytname[yychar1]);
#endif

      YY_XMLParser_CHAR = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;              /* Each real token shifted decrements this */

  YYGOTO(yyerrhandle);

YYLABEL(yyerrdefault)  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) YYGOTO(yydefault);
#endif

YYLABEL(yyerrpop)   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YY_XMLParser_LSP_NEEDED
  yylsp--;
#endif

#if YY_XMLParser_DEBUG != 0
  if (YY_XMLParser_DEBUG_FLAG)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

YYLABEL(yyerrhandle)

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    YYGOTO(yyerrdefault);

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    YYGOTO(yyerrdefault);

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	YYGOTO(yyerrpop);
      yyn = -yyn;
      YYGOTO(yyreduce);
    }
  else if (yyn == 0)
    YYGOTO(yyerrpop);

  if (yyn == YYFINAL)
    YYACCEPT;

#if YY_XMLParser_DEBUG != 0
  if (YY_XMLParser_DEBUG_FLAG)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = YY_XMLParser_LVAL;
#ifdef YY_XMLParser_LSP_NEEDED
  *++yylsp = YY_XMLParser_LLOC;
#endif

  yystate = yyn;
  YYGOTO(yynewstate);
/* end loop, in which YYGOTO may be used. */
  YYENDGOTO
}

/* END */

/* #line 1010 "/usr/local/lib/bison.cc" */
#line 1531 "src/xml/XMLParser.cpp"
#line 268 "src/xml/xml.yacc"



#include "XMLParseFunctions.h"
