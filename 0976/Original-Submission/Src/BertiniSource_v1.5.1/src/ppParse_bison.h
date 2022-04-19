
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

/* Line 1676 of yacc.c  */
#line 16 "ppParse_bison.y"

  char *name;



/* Line 1676 of yacc.c  */
#line 94 "src/ppParse_bison.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE ppParse_bisonlval;


