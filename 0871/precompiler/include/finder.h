/*=============================================================================
author        : Walter Schreppers
filename      : finder.h
created       : 12/5/2003 at 16:50:57
modified      : 
version       : 
copyright     : Walter Schreppers
bugreport(log): 
=============================================================================*/

#ifndef FINDER_H
#define FINDER_H


#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class Finder {

  public:
    //constructor & destructor
    //==========================
    Finder();
    ~Finder();

    //public members
    //==============
    void updatePos( const string& );
    void find( ofstream&, char*, int );
    
    int getLine(){ return line; }
    int getCol(){  return col;  }
    
  private:
    //private members:
    //================
    void init();

    //private locals:
    //===============
    int line;
    int col;

}; //end of class Finder

#endif

