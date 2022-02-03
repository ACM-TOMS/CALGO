/*=============================================================================
author        :Walter Schreppers
filename      :Convert.h
created       :17/04/2000
modified      :25/02/2002
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef CONVERT_H
#define CONVERT_H


#include "ConvertConfig.h"
#include "Matcher.h"
#include "Var.h"
#include "ParserConf.h"

#include "error.h"

class Convert{
	public:
	  
	  //constructor/destructor
	  //======================
    Convert(const Var&,Matcher*,ParserConf*); //the var to be converted, the matcher (from xml file), ParserConf*	  
    ~Convert(){}
	 
    //members
    //=======
    string getValue(); //returns the new value (converted according to rules specified in matcher and convertconfig)
    
  private:
    //locals
    //======
    string strValue;
    types  varType;
    int fLine;
    
    Var fVar;
    Matcher* fMatcher;
    ParserConf* fConf;
		
    //private members
    //===============

    string getFunctionParams( int&, const string&);
    string getFullAssignType( const Var&, unsigned int& );    
    
};    

#endif

