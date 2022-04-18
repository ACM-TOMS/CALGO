#include "pinode.h"


PINode::PINode(){
  fTarget="";
  fData="";
}


PINode::PINode(const string& t,const string& d){
  fTarget=t;
  fData=d;
}


PINode::PINode(const PINode& a):XMLNode(a){
  if( this != &a ){
    fTarget=a.fTarget;
    fData=a.fData;
  }
}



PINode& PINode::operator=(const PINode& a){
  fTarget=a.fTarget;
  fData=a.fData;
  return *this;
}



void PINode::setTarget(const string& n){
  fTarget=n;
}

string PINode::target(){
  return fTarget;
}


void PINode::setData(const string& n){
  fData=n;
}

string PINode::data(){
  return fData;
}


XMLNode::NodeType PINode::nodeType(){
    return XMLNode::PINode;
}

void PINode::show(int indent){
  string s="";
  for(int j=0;j<indent;j++) s+="  ";
  string out="";
  out=s+"process instruction target: "+fTarget;
  out+="\n"+s+"process instruction data: "+fData;
  cout<<out<<endl;  
}



