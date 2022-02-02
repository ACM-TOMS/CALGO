/*=============================================================================
author        :Walter Schreppers
filename      :ConfigFile.h
created       :/
modified      :15/05/2000
version       :2
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef CONFIGFILE_H
#define CONFIGFILE_H


#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "ConfigElem.h"
#include "XMLParser.h"

class ConfigFile:public vector<ConfigElem>{
  public:
    //exception class
    //===============       
    struct Error{
      string fError;
      Error( string s ): fError(s) {};
    };

  	//constructor/destructor
	  //======================
    ConfigFile(char* filename);   // read a file with formatted as procedurename::variablename 
    ConfigFile();	                // default constructor initialize to empty vector
    ~ConfigFile(){this->clear();}
	  
		//members
		//=======
		
    bool locate( const string& proc, const string& var ); //look for procedure and/or variable name 
    void loadFile (char* );                               // read a file with formatted as procedurename::variablename 
  
  private:
    void readTree( XMLNode* );
    void addElement( ElementNode* );
    void addLine( string );
    string stripSpaces( const string& );
};

#endif


