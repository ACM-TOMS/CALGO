/*=============================================================================
author        :Walter Schreppers
filename      :Variables.h
created       :/
modified      :18/04/2000
version       :2
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef VARIABLES_H
#define VARIABLES_H

#include <string>
#include <list>
#include <iostream>

#include "error.h"

#include "Var.h"

class Variables:public list<Var>{
	public:
    	  
	  //constructor/destructor
	  //======================
	  Variables(){}
	  ~Variables(){}
	  
		//members
		//=======
	  void insert(const Var&);
		void insert(types,const string&,const string&); // insert var with a type,pointer/ref spec string, name, typestring
		void insert(types,const string&,const string&,const string&); // insert var with a type,pointer/ref spec string, name, match type
	  void insert( types,
                 const string&,
                 const string&,
                 const string&,
                 const string&,
                 const string& 
               ); // insert var with a type,pointer/ref spec string, name, match type, source type, target type
    		
		void openBracket();	    //add a Var with type tOpenBracket to the end of the list
		void closeBracket();    //remove all Variables till the first occurence of tOpenBracket or if the pseudo stack-list is empty
		void expandLastValue(const string&);
    void expandLastValue(const ParseElem&);
		void expandLastName(const string&);
		void expandLastArrayArg(const string&);
		void show();
		Var lookup(const string&); //find first occurence of a Var with this name
};

#endif

