/*=============================================================================
author        :Walter Schreppers
filename      :Reserved.h
created       :/
modified      :18/04/2000
version       :2
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef RESERVED_H
#define RESERVED_H

#ifdef _WIN32
# pragma warning(disable:4786) // We know basic_string generates long names :-)
#endif

#include <string>
#include <algorithm>  //for using find
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

class Reserved : public vector<string> {
	public:

		//constructor/destructor
		//======================
		Reserved();                 //default constructor here all reserved words will also be created
		Reserved(char* fileName);   //read file with reserved words (is depricated now)
		~Reserved(){};              //destructor
		
		
		//members
		//=======
		bool isReserved(char* word); 
		bool isReserved(string word);
		void printAll();


	private:
	  
		//private members
		//===============
    string stripSpaces(const string&);
};

#endif

