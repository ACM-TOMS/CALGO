#include "commentnode.h"


CommentNode::CommentNode(){
  fData="";
}


CommentNode::CommentNode(const string& s){
  fData=s;
}


CommentNode::CommentNode(const CommentNode& a):XMLNode(a){
  if( this != &a ){
    fData=a.fData;
  }
}



CommentNode& CommentNode::operator=(const CommentNode& a){
  fData=a.fData;
  return *this;
}




void CommentNode::setData(const string& s){
  fData=s;
}

string CommentNode::data(){
  return fData;
}


XMLNode::NodeType CommentNode::nodeType(){
    return XMLNode::CommentNode;
}

void CommentNode::show(int indent){
  string s="";
  for(int i=0;i<indent;i++) s+="  ";
  cout<<s+"comment:'"<<fData<<"'"<<endl;
}
