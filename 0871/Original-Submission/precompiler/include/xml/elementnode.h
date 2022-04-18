#ifndef ELEMENTNODE_H
#define ELEMENTNODE_H


#include "xmlnode.h"
#include "attribute.h"
#include <vector>
#include <string>
#include <iostream>

class ElementNode:public XMLNode{
  public:
    //constructors/destructor
    //=======================
    ElementNode();
    ElementNode(const string&); //provide name
    ElementNode(const ElementNode&);
    
    //public members
    //==============
    void setTagName(const string&);
    string tagName();
    void show(int indent);
    
    void insertAttribute(const string&); //string containing name="value"
    void insertAttribute(const string& name,const string& value);
    string attribute(const string&,const string&);
    
    XMLNode::NodeType nodeType();
    
    
    //operators
    //=========
    ElementNode& operator=(const ElementNode&);

    //public vars
    //===========    
    vector<Attribute> attributes; //easy access to attribs (no need for getters/setters :)

  private:
    string name;

};


#endif

