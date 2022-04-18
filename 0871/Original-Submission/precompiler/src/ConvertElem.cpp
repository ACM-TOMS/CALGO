#include "ConvertElem.h"


ConvertElem::ConvertElem(){
  mtype=mEmpty;
  hasOperation=false;
}


ConvertElem::ConvertElem(const ConvertElem& m){
  *this=m;
}






void ConvertElem::setSource(XMLNode::iterator i){
  if( (*i)->nodeType()==XMLNode::ElementNode){
    ElementNode* e=(ElementNode*) *i;
    if( e->tagName() == "source"){
      sourceType=cSource;
      sourceName=e->attribute("name","");
    }
    else if( e->tagName()=="target" ){
      sourceType=cTarget;
      sourceName=e->attribute("name","");
    }
    else if( e->tagName()=="rhs" ){
      sourceType=cAssignment;
      sourceName=e->attribute("name","");
    }
    else{
      Warning( "ConvertElem:setSource: found undefined tag '"+e->tagName()+"' in convert rule." );
    }
  }
  else{
    Warning( "Warning convert element contains invalid node." );
  }
}



void ConvertElem::setOperation(XMLNode::iterator i){
  if( (*i)->nodeType() == XMLNode::ElementNode){
    ElementNode* e=(ElementNode*) *i;
    if( e->tagName() =="operation" ){
      operationArg = e->attribute("argument",""); //additional, optional argument for operations          
      if(e->hasChildren()){
        XMLNode::iterator t=e->begin();
        if( (*t)->nodeType()== XMLNode::TextNode){
          TextNode* opNode=(TextNode*) *t;
          
          operation=opNode->data();
          this->hasOperation=true;
        }
      }
    }
    else{
      cerr<<"Warning:The node '"<<e->tagName()<<"' should be named operation"<<endl;
    }
  }
}





ConvertElem& ConvertElem::operator=(const ConvertElem& m){
  if(this!=&m){
    name=m.name;
    sourceName=m.sourceName;
    sourceType=m.sourceType;
    target=m.target;
    operation=m.operation;
    operationArg=m.operationArg;
    hasOperation=m.hasOperation;
    mtype=m.mtype;
  }
  return *this;
}



