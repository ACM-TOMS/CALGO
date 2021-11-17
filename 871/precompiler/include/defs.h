/*=============================================================================
author        :Walter Schreppers
filename      :defs.h
created       :/
modified      :6/04/2001
version       :3
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef DEFS_H
#define DEFS_H

#define sym_eof         1
#define sym_other       2

#define comment         3
#define directives      4

#define quoted_string   5
#define sym_word        6

#define assignment      7
#define comparison      8
#define stream_op       9
#define mult            10
#define plus            11
#define minus           12
#define div             13
#define double_colon    14
#define pointer_op      15
#define bool_op         16
#define percent         17

#define open_bracket    18
#define close_bracket   19
#define open_par        20
#define close_par       21
#define open_brace      22
#define close_brace     23
#define semi_colon      24
#define colon           25
#define white_space     26

#define integer         27
#define decimal         28
#define sym_and         29
#define comma           30
#define tilde           31

#define function_argument 32

#endif

