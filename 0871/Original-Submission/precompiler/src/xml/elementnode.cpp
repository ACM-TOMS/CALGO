#include "elementnode.h"


ElementNode::ElementNode(){
  name="";
  attributes.clear();
}


ElementNode::ElementNode(const string& n){
  name=n;
  attributes.clear();
}


ElementNode::ElementNode(const ElementNode& a):XMLNode(a){
  if( this != &a ){
    name=a.name;
    for(unsigned int i=0;i<a.attributes.size();i++){
      attributes.push_back( a.attributes[i] );
    }
    //copy children
    for(XMLNode::const_iterator i=a.begin();i!=a.end();i++){
      push_back(*i);
    }
  }
}



ElementNode& ElementNode::operator=(const ElementNode& a){
  name=a.name;
  for(unsigned int i=0;i<a.attributes.size();i++){
    attributes.push_back( a.attributes[i] );
  }
  return *this;
}



void ElementNode::setTagName(const string& n){
  name=n;
}

string ElementNode::tagName(){
  return name;
}

void ElementNode::insertAttribute(const string& s){
  Attribute a(s);
  attributes.push_back(a);
}

void ElementNode::insertAttribute(const string& name,const string& value){
  Attribute a(name,value);
  attributes.push_back(a);
}

string ElementNode::attribute(const string& name, const string& defval){
  for(vector<Attribute>::const_iterator i=attributes.begin();i!=attributes.end();++i){
    if(name==i->name) return i->value;
  }
  return defval;
}


XMLNode::NodeType ElementNode::nodeType(){
    return XMLNode::ElementNode;
}

void ElementNode::show(int indent){
  string s="";
  for(int j=0;j<indent;j++) s+="  ";
  s+="<"+name;
  cout<<s<<" : ";
  for(unsigned int i=0;i<attributes.size();i++){
    cout<<attributes[i].name<<"=";
    cout<<attributes[i].value<<",";
  }
  cout<<">"<<endl;  
}



