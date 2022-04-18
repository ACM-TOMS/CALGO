/*=============================================================================
author        :Walter Schreppers
filename      :ParserConf.h
created       :/
modified      :24/06/2000
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef PARSERCONF_H
#define PARSERCONF_H

#include <string>
#include <fstream>
#include <stack>

#ifdef _WIN32
# pragma warning(disable:4786) // We know basic_string generates long names :-)
#endif

#include "defs.h"
#include "Variables.h"
#include "Var.h"
#include "Reserved.h"
#include "ConfigFile.h"
#include "XMLParser.h"


class ParserConf{
    public:
      //constructor/destructor
      //======================
      ParserConf();                 //default constructor here all reserved words will also be created
      ~ParserConf();

      //members
      //=======
      void newScope();
      void newFunction(
                        const string& name,
                        types ftype,
                        const string& matchType,
                        const string& srcType,
                        const string& targetType
                      );
    
      void newForLoop();
      void closeScope();
      void semiColon();
      void decOpenParentheses();
      void incOpenParentheses();
      void decOpenBracket();
      void incOpenBracket();
      
      void updateLocals(const string& s); //depricated!!!!, instead use Parser::updateLocals, this is faster!
      void loadConfigFile(char* filename);
      void loadConvertFile(char* filename);
      void updateLineNr(const string& s);


      //some getters and setters
      string getStrCurFunction();
      string getStrTypCurFunction();
      types getTypCurFunction();
      
      Var getCurrentFunction();
      
      void setSkipFor(bool b);
      bool getSkipFor(){return bSkipFor;}

      void lineInc(){line_nr++;}
      int getLine(){return line_nr;}
      
      Variables* getVars(){return vars;}
      Variables* getFunctions(){return functions;}
      Reserved* getReswords(){return reswords;}
      ConfigFile* getSkipVar(){return skipVar;}

      bool getInsideFuncParam(){return insideFuncParam;}
      void setInsideFuncParam(bool b){insideFuncParam=b;}
      
      bool getInsideForParam(){return insideForParam;}
      //void setInsideForParam(bool b){insideForParam=b;}
      
      int getOpenParentheses(){return openParentheses;}
      int getOpenBracket(){return openBracket;}
      void skipQuoted(const string&, int&);
      XMLNode* convertTree;
        

    private:

      //private locals
      //==============
      int line_nr;
      int openParentheses;   //this is the number of '('
      int openBracket;         //this is the number of '{'
      bool insideFuncParam;  //=TRUE if we're inside the parameter declaration of a function
      bool insideForParam;   //=TRUE in declaration of for loop
      bool insideFunction;

      enum scopetypes {sSemi,sBracket,sFunction};
      stack<scopetypes> scopes;
      
      bool bSkipFor;  
  
      Variables *vars;                  //the table with all the Variables and scope info
      Variables *functions;         //table with declared functions
  
      Reserved *reswords;       //with reswords we can check for c++ reserved words to avoid adding variables of type tOther
  
      ConfigFile* skipVar;      //to be able to skip convertion in certain variables or procedures
                                          //see also the loadConfigFile function below !
      Var currentFunc; //contains information of current function we are in. (for converting return statements)

      //private members
      //===============
      string scopeName( scopetypes );

};

#endif


