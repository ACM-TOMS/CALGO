#ifndef PINODE_H
#define PINODE_H


#include "xmlnode.h"
#include <string>
#include <iostream>

class PINode:public XMLNode{
  public:
    //constructors/destructor
    //=======================
    PINode();
    PINode(const string&,const string&); //provide target,data
    PINode(const PINode&);
    
    //public members
    //==============
    void setTarget(const string&);
    string target();
    void setData(const string&);
    string data();
    
    void show(int indent);
    XMLNode::NodeType nodeType();

    
    //operators
    //=========
    PINode& operator=(const PINode&);

  private:
    string fTarget;
    string fData;
};


#endif

