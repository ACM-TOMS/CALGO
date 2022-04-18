/*=============================================================================
author        :Walter Schreppers
filename      :PrintConverter.h
created       :05/05/2001
modified      :/
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#ifndef PRINTCONVERTER_H
#define PRINTCONVERTER_H

#include <iostream>
#include <vector>

#include <ctype.h> //for isdigit checks
#include <stdio.h>


#include "PrintArg.h"


class PrintConverter{
	public:
	  
	  //constructor/destructor
	  //======================
	  PrintConverter();
    PrintConverter(const string&,const vector<PrintArg>&); 
	  ~PrintConverter();


		//members
		//=======
		void init();
		PrintConverter& operator=(const PrintConverter&);
		string getCoutStatement();
		

    // friends
    // =======
      
    friend ostream& operator <<(ostream&,const PrintConverter&);


  protected:
    
    void insertPercent();    
    bool conversionOperation(string&);
    string defaultSettings();
    string customFlags();
    string getDigitString();
    string getWidthStatement(); //an alternative to getWidth !
		string getWidth();
		string getPrecision();
		string convertFlagChars();
    bool isFormatSpec(string&);

    
	  // locals
	  //=============
	  string format;
	  vector<PrintArg> args;
	  vector<string> flags;
	  unsigned int pos,argpos;
	  bool bNumberSignFlag;
};

#endif



