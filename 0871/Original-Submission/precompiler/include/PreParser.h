/*=============================================================================
author        :Walter Schreppers
filename      :PreParser.h
created       :15/05/2001
modified      :15/05/2001
version       :3
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef PREPARSER_H
#define PREPARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>


#include "Parser.h"


class PreParser : public Parser {
	public:
	  
	  //constructor/destructor
	  //======================
	  PreParser();
    PreParser(ParserConf*); 
	  ~PreParser();

	  
		//overloaded members
		//==================
    void dumpq( ofstream& );


	private:
	
	  //private members
	  //===============
	  void writeVar(ofstream& out,const string&);
	
		//private overloaded members
		//==========================
    void parseVar         ( ofstream&, string );
		void parseAssignment  ( ofstream& );
    void setFloatPrecision( ofstream& );

    void moveToken        ( ofstream& );  
    void moveTokens       ( unsigned int&, ofstream& );
    void copyToken        ( ofstream& );  
    void writeStrType     ( ofstream&, const string& );
    void writeOriginalType( ofstream& );
        
    void convertPrintf    ( ofstream& );
    
};

#endif

