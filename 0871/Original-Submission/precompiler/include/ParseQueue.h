/*=============================================================================
author        :Walter Schreppers
filename      :ParseQueue.h
created       :15/04/2001
modified      :/
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef PARSE_QUEUE_H
#define PARSE_QUEUE_H

#include <fstream>
#include <vector>

#include "ParseElem.h"

class ParseQueue:public vector<ParseElem> {
	public:
	  //constructor/destructor
	  //======================
    ParseQueue();
    ~ParseQueue();
    ParseQueue::ParseQueue(const ParseQueue& p);
  
    //operators
    //=========
    ParseQueue& ParseQueue::operator=(const ParseQueue& p);
    
    //public members
    //==============
    void add(int,string);
    void dump(ofstream&);
    void popfirst();
    string poptext();
    ParseElem pop();
    
    ParseElem get(unsigned int);
    int getToken(unsigned int);
    string getText(unsigned int);
    
    bool changed;
    
	private:
	  //locals
	  //======

};


#endif

