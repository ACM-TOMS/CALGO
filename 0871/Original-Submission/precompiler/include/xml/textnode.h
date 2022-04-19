#ifndef TEXTNODE_H
#define TEXTNODE_H


#include "xmlnode.h"
#include <string>
#include <iostream>

class TextNode:public XMLNode{
  public:
    //constructors/destructor
    //=======================
    TextNode();
    TextNode(const string&);    //provide data
    TextNode(const TextNode&);  //copy constructor
    
    
    //public members
    //==============
    void setData(const string&); //replaces &nbsp; etc...
    void setDataLiteral(const string& s) {fData=s;}
    string data();

    XMLNode::NodeType nodeType();
    void show(int indent);    

    //operators
    //=========
    TextNode& operator=(const TextNode&);

    //public vars
    //===========    

  
  private:
    void replaceAll(string& , const string& ,const string& );
    string fData;

};


#endif

