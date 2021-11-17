  /* =============================================================================*/
  /* author        :Walter Schreppers                                             */
  /* filename      :pre.lex                                                       */
  /* created       :/                                                             */
  /* modified      :13/04/2001                                                    */
  /* version       :3                                                             */
  /* copyright     :Walter Schreppers                                             */
  /* bugreport(log):/                                                             */
  /* =============================================================================*/

%option noyywrap
%{
  #include <iostream>
  #include <string>
  //using namespace std;

  #include <math.h>

  #include "defs.h"
  #include "ParserConf.h"

  #include <stdio.h>

  ParserConf* pConf=new ParserConf();
%}


%x COMMENTA
%x COMMENTB
%x DIRECTIVES

%%

  /* state COMMENTA copies comments of the type // ... */
"//" {
	      BEGIN(COMMENTA);
	      yytext="//";
	      return comment;
	    }
<COMMENTA>\n  {
                    BEGIN(INITIAL);
                    yytext="\n";
    	              return comment;
              }
<COMMENTA>.	{	return comment;}

  
  
  /* state COMMENTB copies multiple line comments */
"/*" 		{
          BEGIN(COMMENTB);
          yytext=("/*");
          return comment;
        }
<COMMENTB>\n	{
                  yytext="\n";
                  return comment;
              }
<COMMENTB>.	{	return comment;}

<COMMENTB>"*/"  {
                    BEGIN(INITIAL);
                    yytext="*/";
                    return comment;
                }


  /* state DIRECTIVES copies directives ex. #include, #define,... */ 
"#" {
      BEGIN(DIRECTIVES);
      return directives;
    }

<DIRECTIVES>[^\n]*	{
                      #ifdef _DEBUG_
                        cerr<<"found directive: "<<yytext<<endl;
                      #endif
                      return directives;
                    }
<DIRECTIVES>\n	{
                    yytext="\n";
                    BEGIN(INITIAL);
                    return directives;
                  }


  /* this finds a quoted string */
"\""([^"]|"\\\"")*"\""        {   //quoted string 
                                  #ifdef _DEBUG_
                                    cerr<<"found following quoted string"<<yytext<<endl;
                                  #endif
                                  return quoted_string;
                              }


[0-9][0-9]*[fFlL]?  {
                      return integer;
                    }

 /* this decimal rule also accepts '2.' as a decimal string */
[0-9][0-9]*(".")[0-9]*([eE][+-]?[0-9][0-9]*)?[fFlL]?  {
                                                        return decimal;
                                                      }
 /* this variation accepts '.2' as a decimal string */
[0-9]*(".")[0-9][0-9]*([eE][+-]?[0-9][0-9]*)?[fFlL]?  {
                                                        return decimal;
                                                      }

[a-zA-Z_][a-zA-Z_0-9]*  {
                          return sym_word;
                        }

[+\-\*/]?[=]      {
                    return assignment;
                  }
[<>][=]?|("!=")|("==")  {
                          return comparison;
                        }
("&&")|("|")|("||")|("!") {
                            return bool_op;
                          }


"," {
      return comma;
    }

"~" {
      return tilde;
    }

"&" {
      return sym_and;
    }

"%" {
      return percent;
    }
                                
("<<")|(">>") {
                return stream_op;
              }

"->"  {
        return pointer_op;
      }
"*"   {
        return mult;
      }
"+"   {
        return plus;
      }
"-"   {
        return minus;
      }
"/"   {
        return div;
      }
":"   {
        return colon;
      }
      
"::"    {
          return double_colon;
        }



"["   {
        return open_bracket;  
      }

"]"   {
        return close_bracket;
      }
"("   {
			  return open_par;
      }
")"   {
			  return close_par;
      }

"{"		{
			  return open_brace; 
		  }

"}"		{
			  return close_brace;
		  }

";"   {
        return semi_colon;
      }                                

[\n \t]*  {
            return white_space;
          }

.       {
          return sym_other;
        }

<<EOF>> {
          return sym_eof;
        }

%%


