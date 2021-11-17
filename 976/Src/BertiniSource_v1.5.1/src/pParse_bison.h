
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
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


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     pEND = 258,
     pHOMVARGP = 259,
     pVARGP = 260,
     pPATHVAR = 261,
     pVAR = 262,
     pPARAM = 263,
     pCONST = 264,
     pFUNC = 265,
     pSUBFUNC = 266,
     pI = 267,
     pPi = 268,
     pNUMBER = 269,
     pSCE = 270,
     pPOW = 271,
     pNAME = 272,
     pNEG = 273
   };
#endif
/* Tokens.  */
#define pEND 258
#define pHOMVARGP 259
#define pVARGP 260
#define pPATHVAR 261
#define pVAR 262
#define pPARAM 263
#define pCONST 264
#define pFUNC 265
#define pSUBFUNC 266
#define pI 267
#define pPi 268
#define pNUMBER 269
#define pSCE 270
#define pPOW 271
#define pNAME 272
#define pNEG 273




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 1676 of yacc.c  */
#line 19 "pParse_bison.y"

  char *name;
  int memLoc;



/* Line 1676 of yacc.c  */
#line 95 "src/pParse_bison.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE pParse_bisonlval;


