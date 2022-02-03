/*=============================================================================
author        :Walter Schreppers
filename      :precompile.cpp
created       :/
modified      :15/05/2001
version       :4
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#include <iostream>
#include <fstream> 

using namespace std;

/*
#ifdef _WIN32
# pragma warning(disable:4786) // We know basic_string generates long names :-)
#include "pre.h"
#endif
*/

#include <stdio.h>
#include <math.h>

#include "ConfigFile.h"
#include "ParseCmdLine.h"
#include "Parser.h"
#include "ParserConf.h"
#include "PreParser.h"
#include "finder.h"


/*=============================================================================
  these are external variables used by the lexer
=============================================================================*/
extern int yylex();
extern FILE *yyin;
extern char *yytext;
extern int yylen;
extern ParserConf* pConf;


/*-----------------------------------------------------------------------------
name        :Precompile
description :convert cpp file into a cpp file which uses bigint, Mpieee
             (or rational) 
parameters  :ofstream& out, const ParseCmdLine& commands
return      :/
exceptions  :/
algorithm   :- set some options for parser using ParseCmdLine& commands
             - write the includes with parser
             - start retreiving tokens from the lexer and execute the
               corresponding parser.parse which will produce output in
               the outputfile out
-----------------------------------------------------------------------------*/
void Precompile(ofstream& out,const ParseCmdLine& commands){
    
    Parser parser(pConf);

    //set some options for the parser
    parser.setRadix(commands.getRadix());
    parser.setPrecision(commands.getPrecision());
    parser.setExprange(commands.getExprangeLow(),commands.getExprangeUp());
    parser.setRounding(commands.getRounding());
    parser.setIoModeStr(commands.getIoModeStr());
    parser.setNoDefault(!commands.getDefault());
    parser.setNoPrintf(commands.getNoPrintf());    
    parser.setInit( commands.getInit(), commands.getInitFileName() );

    parser.writeIncludes(out);
    int token=yylex();
    while(token!=sym_eof){
      parser.parse(out,yytext,token);
			token = yylex();
    }//while end
    
    parser.dumpq(out); //dump any pending stuff in the parseQueue
}



/*-----------------------------------------------------------------------------
name        :Preparse
description :build a list of all constructed variables with their function and
             types
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :
             - start retreiving tokens from the lexer and execute the
               corresponding preparse.parse which will produce output in
               the outputfile out
-----------------------------------------------------------------------------*/
void Preparse(ofstream& out){
  PreParser preparse(pConf);
  int token = yylex();
  while( token != sym_eof ){
    preparse.parse( out, yytext, token );
    token = yylex();
  }//while end
  preparse.dumpq(out);
}



void findConstants( ofstream& out ){
  Finder finder;
  int token = yylex();
  while( token != sym_eof ){
    finder.find( out, yytext, token );
    token = yylex();
  }
}

/*-----------------------------------------------------------------------------
name        :main
description :We parse the input and open the needed input and output file and pass
             them to the Precompile function which does the convertion. 
parameters  :int argc , char*argv[]
return      :int
exceptions  :/
algorithm   :- use ParseCmdLine to parse the command line
             - open the inputfile for reading and 
               assign it to yyin (input for the lexer)
             - open the outputfile for writing
             - set the options for the lexer
             - call Precompile with the outputfile,ParseCmdLine commands
               to set it's options
-----------------------------------------------------------------------------*/
int main(int argc,char *argv[]){
  try{
	  ParseCmdLine commands;
	  if( commands.parse(argc,argv) ){
	    yyin = fopen( commands.getInFileName(),"r" ); //input file
	    ofstream out( commands.getOutFileName() );    //destination precompiled output file 
      
      if( commands.getConstants() ){
        out << "CONSTANTS IN FILE : '" << commands.getInFileName() << "'" << endl;
        findConstants( out );        
        out << endl;
      }
      else if( commands.getPreParse() ){
        //the xml configuration file
        if ( commands.getConvertFile() ) pConf->loadConvertFile(commands.getConvertFileName());
        //use PreParser instead of Parser!!!
        Preparse( out );
      }
      else{
        //set some options for the parser
	      pConf->setSkipFor(commands.getNoFor());
	      
        //file for the variables and functions to be skipped
	      if (commands.getConfigFile()) pConf->loadConfigFile(commands.getConfFileName()); 

        //the xml configuration file
        if (commands.getConvertFile()) pConf->loadConvertFile(commands.getConvertFileName());
	      Precompile(out,commands);
	    }
	  }
	  else{
	    cout<<commands.getUsage();
	  } 
  }
  catch(Error e){
    cerr<<e.fError<<" at line:"<<pConf->getLine()<<endl;
  }
  catch(...){
    cerr<<"Precompile main caught an exception at line:"<<pConf->getLine()<<endl;
  }
  return 0;
} // main end

		
