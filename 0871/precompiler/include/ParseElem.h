/*=============================================================================
author        :Walter Schreppers
filename      :ParseElem.h
created       :14/04/2001
modified      :7/3/2002
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef PARSE_ELEM_H
#define PARSE_ELEM_H

#include <string>
using namespace std;

class ParseElem {
	public:
	  
	  //constructor/destructor
	  //======================
    ParseElem() { token=0;text=""; }
    ParseElem( int, const string& );
    ParseElem( const ParseElem& );
    ~ParseElem();
    
    
    //public members
    //==============
    ParseElem& operator=( const ParseElem& );
 
 	  string text;
	  int token;

    
	private:
	  //locals
	  //======

};


#endif

