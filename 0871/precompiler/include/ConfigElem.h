/*=============================================================================
author        :Walter Schreppers
filename      :ConfigElem.h
created       :/
modified      :15/05/2000
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef CONFIGELEM_H
#define CONFIGELEM_H

#include <string>
#include <iostream>
using namespace std;

class ConfigElem {
	public:
		//constructor/destructor
		//======================		
		ConfigElem();
		ConfigElem(const string&,const string&);
		
		~ConfigElem(){};
		
		//members
		//=======
		string getVar();
		string getProcedure();
		
		void setVar(const string&);
		void setVar(char*);

		void setProcedure(const string&);
		void setProcedure(char*);
		
	private:
	  //locals
	  //======
		string fVar;        //variable name to skip
		string fProcedure;  //procedure name in whitch all or some variables must be skipped
		
		//private members
		//===============
		string charpToStr(char *c); //converts char* to string
};

#endif

