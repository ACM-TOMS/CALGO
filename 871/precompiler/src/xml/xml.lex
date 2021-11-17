  /*==================================================================*/
  /* description  : xml lexer                                         */
  /* author       : Walter Schreppers                                 */
  /* created      : 9/1/2002                                          */
  /* modified     : /                                                 */
  /* bugs         : /                                                 */
  /*==================================================================*/  



%option noyywrap

  /* generate c++ class */
%option c++


%{
  /*==================================================================*/
  /* includes                                                         */
  /*==================================================================*/
  #include "XMLParser.h"
  #include <errno.h>



  
  /*==================================================================*/
  /* variables and functions for lexer                                */
  /*==================================================================*/
  int LineNr=1; //keep track of line 
  int column=1; //keep track of the column
  
  void resetXMLLineCol(){
    LineNr=1;
    column=1;
  }

  XMLParser* parser;
  void setXMLParser(XMLParser* p){
    parser=p;
  }

  int getXMLLineNr(){
    return LineNr;
  } 
  
   
  int getXMLColumn(){
    return column;
  }



  //---------------------------------------------
  // keep track of line and column 
  //
  void count(char* lextext){
    int i;
    for (i = 0; lextext[i] != '\0'; i++){
      if (lextext[i] == '\n'){
        column = 1;
        LineNr++;
      }
      else if (lextext[i] == '\t'){
        column += 8 - (column % 8);
      }
      else{
        column++;
      }
    }
  }

 
%}

  /*==================================================================*/
  /* STATES AND SOME LEXER DEFS                                       */
  /*==================================================================*/




D [0-9]

L [a-zA-Z_0-9]|"."|"-"|":"

E [DdEe][+-]?{D}+
S [ \t\v\f\n]*
Q (("\"")|("'"))

V [a-zA-Z0-9_.:] 
VNUM ({V}|'-')+


%x S_COMMENT 
%x S_TAG_OPEN 
%x S_TAG_CLOSE 
%x S_TAG
%x S_PROCESS_TAG
%x S_PROCESS_TAG_DATA
%x S_PROLOG
%x S_PROLOG_END
%x S_DOCTYPE

%%

  /*==================================================================*/
  /* TOKEN RECOGNITION                                                */
  /*==================================================================*/

  /*------------------------- comments -------------------------------*/

  /*last minute change, since comments are apperantly */ 
  /*allowed outside document root, we don't return anything */
  /*to the parser here!!! Which means comment nodes are just ommitted and not put in the tree*/

"<!--"        {
                count(yytext);
                BEGIN(S_COMMENT);
                //return XMLParser::COMMENT_BEGIN;
              }

<S_COMMENT>"-->"      {                                            
                        count(yytext);
                        BEGIN(INITIAL);                                 
                        //return XMLParser::COMMENT_END;
                      }
<S_COMMENT>\n         {
                        LineNr++;
                        //parser->yylval.str=new string("\n");
                        //return XMLParser::COMMENT_DATA;
                      }                                                                 
<S_COMMENT>.          {                                              
                        count(yytext);
                        //parser->yylval.str=new string(yytext);
                        //return XMLParser::COMMENT_DATA;
                      }                                           


  /*-----------------------simplified doctype spec--------------------*/
("<!DOCTYPE")|("<!doctype")  {
                                  count(yytext);
                                  BEGIN(S_DOCTYPE);
                             }
<S_DOCTYPE>">"               {
                                  count(yytext);
                                  BEGIN(INITIAL);
                                  return XMLParser::DOCTYPE;
                                }
<S_DOCTYPE>.                    {
                                  count(yytext);
                                }

  /*------------------------- tags -----------------------------------*/


("</"){S}             {
                        count(yytext);
                        BEGIN(S_TAG_CLOSE);
                      }

<S_TAG_CLOSE>{L}*         {
                            count(yytext);
                            parser->yylval.str=new string(yytext);
                            return XMLParser::CLOSE_TAG;
                          }


<S_TAG_CLOSE>{S}(">")     {
                            count(yytext);
                            BEGIN(INITIAL);
                            return XMLParser::CLOSE_TAG_END;
                          }
<S_TAG_CLOSE>.            {
                            count(yytext);
                            string errorStr="invalid char in closing tag '"+string(yytext)+"'";
                            parser->syntaxError(errorStr);
                            return EOF;
                          }

("<"){S}   {
                BEGIN(S_TAG_OPEN);
              }


<S_TAG_OPEN>{L}*  {
                    count(yytext);
                    parser->yylval.str=new string(yytext);
                    BEGIN(S_TAG);
                    return XMLParser::OPEN_TAG_BEGIN;
                  }


<S_TAG_OPEN>. {
                count(yytext);
                string errorStr="invalid char in open tag '"+string(yytext)+"'";
                parser->syntaxError(errorStr);
                return EOF;
              }

<S_TAG>{L}*{S}("="){S}"\""([^"]|"\\\"")*"\""{S} {
                                                  count(yytext);
                                                  parser->yylval.str=new string(yytext);
                                                  return XMLParser::ATTRIBUTE; 
                                                }
<S_TAG>{L}*{S}("="){S}"'"[^']*"'"{S}            {
                                                  count(yytext);
                                                  parser->yylval.str=new string(yytext);
                                                  return XMLParser::ATTRIBUTE; 
                                                }
<S_TAG>{S}              { count(yytext); }
<S_TAG>{S}(">")         {
                          count(yytext);
                          BEGIN(INITIAL);
                          return XMLParser::OPEN_TAG_END;
                        }
<S_TAG>{S}("/>")        {
                          count(yytext);
                          BEGIN(INITIAL);
                          return XMLParser::SINGLE_TAG_END;
                        }

<S_TAG>.      {
                count(yytext);
                string errorStr="invalid tag argument '"+string(yytext)+"'";
                parser->syntaxError(errorStr);
                return EOF;
              }





  /*------------- simple document prolog/xml version -----------------*/

("<\?")(X|x)(M|m)(L|l){S}       {
                                  count(yytext);
                                  BEGIN(S_PROLOG);
                                  return XMLParser::PROLOG_START;
                                }

<S_PROLOG>{S}"version"{S}       {
                                  count(yytext);
                                  return XMLParser::VERSION;
                                }
<S_PROLOG>{S}"="{S}             {
                                  count(yytext);
                                  return XMLParser::EQ;
                                }                                
<S_PROLOG>{Q}{VNUM}{Q}          {
                                  count(yytext);
                                  parser->yylval.str=new string(yytext);
                                  BEGIN(S_PROLOG_END);
                                  return XMLParser::NR;
                                }
<S_PROLOG_END>("\?>")           {
                                  count(yytext);
                                  BEGIN(INITIAL);
                                  return XMLParser::PROLOG_END;
                                }
<S_PROLOG_END>.                 {
                                   count(yytext);
                                   parser->yylval.str=new string(yytext);
                                   return XMLParser::PROCESS_TAG_DATA;
                                }

  /*------------------------processing instructions-------------------*/


("<\?")               {
                        count(yytext);
                        BEGIN(S_PROCESS_TAG);
                      }
<S_PROCESS_TAG>{L}*   {
                        //target string...
                        count(yytext);
                        parser->yylval.str=new string(yytext);
                        //BEGIN(S_PROCESS_TAG_DATA);
                        return XMLParser::PROCESS_TAG_TARGET;                        
                      }
<S_PROCESS_TAG>{S}    {
                        BEGIN(S_PROCESS_TAG_DATA);
                      }
<S_PROCESS_TAG>.      {
                        count(yytext);
                        string errorStr="invalid tag argument '"+string(yytext)+"'";
                        parser->syntaxError(errorStr);
                        return EOF;
                      }


<S_PROCESS_TAG_DATA>.                    {
                                            count(yytext);
                                            parser->yylval.str=new string(yytext);
                                            return XMLParser::PROCESS_TAG_DATA;
                                          }
<S_PROCESS_TAG_DATA>\n                    {
                                            LineNr++;
                                            parser->yylval.str=new string("\n");
                                            return XMLParser::PROCESS_TAG_DATA;
                                          }
                                          
<S_PROCESS_TAG_DATA>("\?>")            {
                                            count(yytext);
                                            BEGIN(INITIAL);
                                            return XMLParser::PROCESS_TAG_END;
                                          }





  /*--------------- formatting, errors and textnodes ----------------*/

{S}             {
                  count(yytext);
                  parser->yylval.str=new string(yytext);
                  return XMLParser::WHITESPACE;
                }

(([^<])|(\"))*  {
                  count(yytext);
                  parser->yylval.str=new string(yytext);
                  return XMLParser::TEXTNODE;
                }

%%

