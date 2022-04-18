#ifndef YY_XMLParser_h_included
#define YY_XMLParser_h_included

#line 1 "/usr/local/lib/bison.h"
/* before anything */
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
#endif
#include <stdio.h>

/* #line 14 "/usr/local/lib/bison.h" */
#line 21 "include/xml/XMLParser.h"
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

#line 14 "/usr/local/lib/bison.h"
 /* %{ and %header{ and %union, during decl */
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
/* WARNING obsolete !!! user defined YYLTYPE not reported into generated header */
/* use %define LTYPE */
#endif
#endif
#ifdef YYSTYPE
#ifndef YY_XMLParser_STYPE 
#define YY_XMLParser_STYPE YYSTYPE
/* WARNING obsolete !!! user defined YYSTYPE not reported into generated header */
/* use %define STYPE */
#endif
#endif
#ifdef YYDEBUG
#ifndef YY_XMLParser_DEBUG
#define  YY_XMLParser_DEBUG YYDEBUG
/* WARNING obsolete !!! user defined YYDEBUG not reported into generated header */
/* use %define DEBUG */
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

/* #line 63 "/usr/local/lib/bison.h" */
#line 133 "include/xml/XMLParser.h"

#line 63 "/usr/local/lib/bison.h"
/* YY_XMLParser_PURE */
#endif

/* #line 65 "/usr/local/lib/bison.h" */
#line 140 "include/xml/XMLParser.h"

#line 65 "/usr/local/lib/bison.h"
/* prefix */
#ifndef YY_XMLParser_DEBUG

/* #line 67 "/usr/local/lib/bison.h" */
#line 147 "include/xml/XMLParser.h"

#line 67 "/usr/local/lib/bison.h"
/* YY_XMLParser_DEBUG */
#endif
#ifndef YY_XMLParser_LSP_NEEDED

/* #line 70 "/usr/local/lib/bison.h" */
#line 155 "include/xml/XMLParser.h"

#line 70 "/usr/local/lib/bison.h"
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

/* TOKEN C */
#ifndef YY_USE_CLASS

#ifndef YY_XMLParser_PURE
extern YY_XMLParser_STYPE YY_XMLParser_LVAL;
#endif


/* #line 143 "/usr/local/lib/bison.h" */
#line 233 "include/xml/XMLParser.h"
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


#line 143 "/usr/local/lib/bison.h"
 /* #defines token */
/* after #define tokens, before const tokens S5*/
#else
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

/* #line 182 "/usr/local/lib/bison.h" */
#line 299 "include/xml/XMLParser.h"
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


#line 182 "/usr/local/lib/bison.h"
 /* decl const */
#else
enum YY_XMLParser_ENUM_TOKEN { YY_XMLParser_NULL_TOKEN=0

/* #line 185 "/usr/local/lib/bison.h" */
#line 329 "include/xml/XMLParser.h"
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


#line 185 "/usr/local/lib/bison.h"
 /* enum token */
     }; /* end of enum declaration */
#endif
public:
 int YY_XMLParser_PARSE(YY_XMLParser_PARSE_PARAM);
 virtual void YY_XMLParser_ERROR(char *msg) YY_XMLParser_ERROR_BODY;
#ifdef YY_XMLParser_PURE
#ifdef YY_XMLParser_LSP_NEEDED
 virtual int  YY_XMLParser_LEX(YY_XMLParser_STYPE *YY_XMLParser_LVAL,YY_XMLParser_LTYPE *YY_XMLParser_LLOC) YY_XMLParser_LEX_BODY;
#else
 virtual int  YY_XMLParser_LEX(YY_XMLParser_STYPE *YY_XMLParser_LVAL) YY_XMLParser_LEX_BODY;
#endif
#else
 virtual int YY_XMLParser_LEX() YY_XMLParser_LEX_BODY;
 YY_XMLParser_STYPE YY_XMLParser_LVAL;
#ifdef YY_XMLParser_LSP_NEEDED
 YY_XMLParser_LTYPE YY_XMLParser_LLOC;
#endif
 int YY_XMLParser_NERRS;
 int YY_XMLParser_CHAR;
#endif
#if YY_XMLParser_DEBUG != 0
public:
 int YY_XMLParser_DEBUG_FLAG;	/*  nonzero means print parse trace	*/
#endif
public:
 YY_XMLParser_CLASS(YY_XMLParser_CONSTRUCTOR_PARAM);
public:
 YY_XMLParser_MEMBERS 
};
/* other declare folow */
#endif


#if YY_XMLParser_COMPATIBILITY != 0
/* backward compatibility */
#ifndef YYSTYPE
#define YYSTYPE YY_XMLParser_STYPE
#endif

#ifndef YYLTYPE
#define YYLTYPE YY_XMLParser_LTYPE
#endif
#ifndef YYDEBUG
#ifdef YY_XMLParser_DEBUG 
#define YYDEBUG YY_XMLParser_DEBUG
#endif
#endif

#endif
/* END */

/* #line 236 "/usr/local/lib/bison.h" */
#line 407 "include/xml/XMLParser.h"
#endif
