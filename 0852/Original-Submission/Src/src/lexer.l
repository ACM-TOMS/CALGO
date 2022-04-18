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
 * lexer.l                                                                  *
 *     -> needs GNU flex                                                    *
 ****************************************************************************/

/*
float4               {sign}{float1}
float5               {sign}{float2}
float6               {sign}{float3}
{float4}          {return(FLOAT);}
{float5}          {return(FLOAT);}
{float6}          {return(FLOAT);}
integer2             {sign}{integer1}
{integer2}         {return(INTEGER);}
*/

%{
#include "parser.h"
#include <stdio.h>
extern int IByyline;
%}
enter                "\n"
space                [ \t\r]+

lbr                  "("
rbr                  ")"
lsbr                 "["
rsbr                 "]"
dollar               "$"
twopoints            ".."
realposinfinity      "+oo"
realneginfinity      "-oo"
prev                 "pred"
succ                 "succ"

num                  [0-9]
add                  "+"
sub                  "-"
sign                 ({add}|{sub})
integer1             {num}+
exp                  (E|e)

mk                   "."
float1               {num}*{mk}{num}+
float2               {num}*{mk}{num}+{exp}{num}+
float3               {num}*{mk}{num}+{exp}{sign}{num}+


div                  "/"
pow                  "**"|"^"
psup                 ">="
pinf                 "<="

comment              "%".*

pragmamaxtime        "Time"|"TIME"|"time"
pragmaconstants      "Constants"|"CONSTANTS"|"constants"|"Constant"|"CONSTANT"|"constant"
pragmadomains        "Domains"|"DOMAINS"|"domains"|"Variables"|"VARIABLES"|"variables"|"Variable"|"variable"|"VARIABLE"
pragmaconstraints    "Constraints"|"CONSTRAINTS"|"constraints"|"Constraint"|"CONSTRAINT"|"constraint"
pragmabisection      "Bisection"|"BISECTION"|"bisection"|"Branch"|"branch"|"BRANCH"
pragmaoutput         "Output"|"OUTPUT"|"output"

bisectionchoice      "Choice"|"CHOICE"|"choice"
bisectionparts       "Parts"|"PARTS"|"parts"
nobisection          "None"|"NONE"|"none"

bisectionchoicerr    "rr"|"RR"|"round_robin"|"ROUND_ROBIN"|"Round_Robin"
bisectionchoicelf    "lf"|"LF"|"largest_first"|"LARGEST_FIRST"|"Largest_First"
bisectionchoicemn    "mn"|"MN"|"max_narrow"|"MAX_NARROW"|"Max_Narrow"

bisectionnumber      "Number"|"number"|"NUMBER"

bisectionsubpaving   "paving"|"Paving"|"PAVING"
bisectionpoints      "points"|"Points"|"POINTS"

precision            "Precision"|"PRECISION"|"precision"|"Prec"|"PREC"|"prec"

mode                 "Mode"|"MODE"|"mode"

outputhullmode       "Hull"|"HULL"|"hull"
outputunionmode      "Union"|"UNION"|"union"
outputdigits         "Digits"|"DIGITS"|"digits"
outputstyle          "Style"|"STYLE"|"style"
outputboundstyle     "Bound"|"BOUND"|"bound"|"Bounds"|"BOUNDS"|"bounds"
outputmidpointstyle  "Midpoint"|"MIDPOINT"|"midpoint"|"Mid"|"MID"|"mid"
outputsolution       "Solutions"|"SOLUTIONS"|"solutions"|"Sol"|"SOL"|"sol"|"Solution"|"SOLUTION"|"solution"
outputallsolution    "All"|"ALL"|"all"

pragmaconsistency    "Consistency"|"CONSISTENCY"|"consistency"
localconsistency     "Local"|"LOCAL"|"local"
strongconsistency    "Strong"|"STRONG"|"strong"

consistencybc3       "BC3"|"bc3"|"box"|"BOX"|"Box"
consistencybc3newton "BC3_NEWTON"|"bc3_newton"|"BOX_NEWTON"|"box_newton"|"Box_Newton"|"NUMERICA"|"Numerica"|"numerica"
consistencybc4       "BC4"|"bc4"|"hull_box"|"HULL_BOX"|"Hull_Box"
consistencybc5       "BC5"|"bc5"|"hull_box_newton"|"HULL_BOX_NEWTON"|"Hull_Box_Newton"|"Best"|"BEST"|"best"
consistencyhc3       "HC3"|"hc3"|"hull_decomp"|"HULL_DECOMP"|"Hull_Decomp"
consistencyhc4       "HC4"|"hc4"|"hull"|"HULL"|"Hull"
consistencyhc4i      "HC4I"|"hc4i"|"hc4I"
consistencyhc4newton "HC4_NEWTON"|"HC4_Newton"|"hc4_Newton"|"hc4_newton"|"hull_newton"|"HULL_NEWTON"|"Hull_Newton"

consistency3b        "3B"|"3b"
consistencyweak3b    "weak3B"|"WEAK3B"|"Weak3B"|"weak3b"|"WEAK3b"|"Weak3b"
precision3b          "width"|"WIDTH"|"Width"
precision2b          "phi"|"Phi"|"PHI"

improvement2b        "improve"|"IMPROVE"|"Improve"

minimum              "min"
maximum              "max"
integertype          "int"
realtype             "real"

ident            [a-zA-Z@][a-zA-Z0-9_]*
%%

{pragmamaxtime}        {return(PRAGMAMAXTIME);}
{pragmaconstants}      {return(PRAGMACONSTANTS);}
{pragmadomains}        {return(PRAGMADOMAINS);}
{pragmaconstraints}    {return(PRAGMACONSTRAINTS);}
{pragmabisection}      {return(PRAGMABISECTION);}
{pragmaoutput}         {return(PRAGMAOUTPUT);}

{bisectionchoice}      {return(BISECTIONCHOICE);}
{bisectionparts}       {return(BISECTIONPARTS);}
{nobisection}          {return(NOBISECTION);}

{bisectionchoicerr}    {return(BISECTIONCHOICERR);}
{bisectionchoicelf}    {return(BISECTIONCHOICELF);}
{bisectionchoicemn}    {return(BISECTIONCHOICEMN);}

{bisectionnumber}      {return(BISECTIONNUMBER);}
{bisectionsubpaving}   {return(BISECTIONSUBPAVING);}
{bisectionpoints}      {return(BISECTIONPOINTS);}

{precision}            {return(PRECISION);}

{mode}                 {return(MODE);}
{outputhullmode}       {return(OUTPUTHULLMODE);}
{outputunionmode}      {return(OUTPUTUNIONMODE);}
{outputdigits}         {return(OUTPUTDIGITS);}

{outputstyle}          {return(OUTPUTSTYLE);}
{outputboundstyle}     {return(OUTPUTBOUNDSTYLE);}
{outputmidpointstyle}  {return(OUTPUTMIDPOINTSTYLE);}

{outputsolution}       {return(OUTPUTSOLUTION);}
{outputallsolution}    {return(OUTPUTALLSOLUTION);}

{pragmaconsistency}    {return(PRAGMACONSISTENCY);}
{localconsistency}     {return(LOCALCONSISTENCY);}
{strongconsistency}    {return(STRONGCONSISTENCY);}


{consistencybc3}       {return(CONSISTENCYBC3);}
{consistencybc3newton} {return(CONSISTENCYBC3Newton);}
{consistencybc4}       {return(CONSISTENCYBC4);}
{consistencybc5}       {return(CONSISTENCYBC5);}
{consistencyhc3}       {return(CONSISTENCYHC3);}
{consistencyhc4}       {return(CONSISTENCYHC4);}
{consistencyhc4i}      {return(CONSISTENCYHC4I);}
{consistencyhc4newton} {return(CONSISTENCYHC4Newton);}
{precision2b}          {return(PRECISION2B);}
{consistency3b}        {return(CONSISTENCY3B);}
{consistencyweak3b}    {return(CONSISTENCYWEAK3B);}
{precision3b}          {return(PRECISION3B);}
{improvement2b}        {return(IMPROVEMENT2B);}

{minimum}         {return(MINIMUM);}
{maximum}         {return(MAXIMUM);}

{integertype}     {return(INTEGERTYPE);}
{realtype}        {return(REALTYPE);}

{twopoints}       {return(TWOPOINTS);}

{realposinfinity} {return(REALPOSINFINITY);}
{realneginfinity} {return(REALNEGINFINITY);}
{prev}            {return(PREV);}
{succ}            {return(SUCC);}

{dollar}          {return(NEWDATA);}

in                {return(INDOM);}
sqrt              {return(SQRT);}
log               {return(LOG);}
exp               {return(EXP);}
cos               {return(COS);}
sin               {return(SIN);}
tan               {return(TAN);}
cosh              {return(COSH);}
sinh              {return(SINH);}
tanh              {return(TANH);}
acos              {return(ACOS);}
asin              {return(ASIN);}
atan              {return(ATAN);}
acosh             {return(ACOSH);}
asinh             {return(ASINH);}
atanh             {return(ATANH);}

{ident}           {return(IDENT);}
{lsbr}            {return(LSBR);}
{rsbr}            {return(RSBR);}
{lbr}             {return(LBR);}
{rbr}             {return(RBR);}
[;]               {return(SCOLON);}
[:]               {return(COLON);}
[,]               {return(COMMA);}
[_]               {return(UNDERSCORE);}
[=]               {return(EQU);}
{psup}            {return(SUP);}
{pinf}            {return(INF);}

{add}             {return(ADD);}
{sub}             {return(SUB);}
[*]               {return(MUL);}
{div}             {return(DIV);}
{pow}             {return(POW);}

{integer1}         {return(INTEGER);}
{float1}          {return(FLOAT);}
{float2}          {return(FLOAT);}
{float3}          {return(FLOAT);}




"/*"  {
        register int c;
        for( ;; )
        {
          while( (c=input()) != '*' && c!=EOF )
          {
            if( c=='\n' ) IByyline ++;
          }
          ; /* eat up text of comment */

          if( c=='*' )
          {
            while( (c=input()) == '*' )
            {
              if( c=='\n' ) IByyline ++;
            }
            ;
            if( c=='/' ) break;    /* end of comment */
            if( c=='\n' ) IByyline ++;
          }
          if( c==EOF ) { break; }  /* illegal EOF */
        }
      }

{enter}           { IByyline ++; }
{space}           {}
{comment}         {}
%%

int yywrap()
{
  return 1;
}
