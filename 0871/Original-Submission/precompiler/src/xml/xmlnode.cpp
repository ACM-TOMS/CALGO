#include "xmlnode.h"


XMLNode::XMLNode(){
  clear();
  parent=0;
}

XMLNode::XMLNode(XMLNode* p){
  clear();
  setParent(p);
}


void XMLNode::setParent(XMLNode* p){
  parent=p;
}

XMLNode::XMLNode( const XMLNode& x ):list<XMLNode*>(x){
  if( this != &x ){
    name=x.name;
    fType=x.fType;
  }
}

XMLNode::~XMLNode(){
  destroy(this);
}

//
//recursively walk down tree and delete every node bottom up
//
void XMLNode::destroy(XMLNode* node){
  if( (node!=0) && (node->size()>0) ) {
    for( XMLNode::iterator i=node->begin(); i!= node->end(); ++i ){
      
      destroy( *i );
      //cout<<"deleting:"<<(*i)->getName()<<endl;
      
      (*i)->clear(); //free children
      //delete ( *i ); //free mem
    }
  }
}


void XMLNode::setName(const string& s){
  name=s;
}

string XMLNode::getName() const{
  return name;
}


void XMLNode::show(int indent){
  string s="";
  for(int i=0;i<indent;i++) s+="  ";
  s+=getName();
  cout<<s<<endl;
}

//
//recursively walk through tree and show node names with indentation
//
void XMLNode::showTree(XMLNode* node, int indent) const{
  indent++;
  if( (node!=0) && (node->size()>0) ){
    for( XMLNode::const_iterator i=node->begin(); i!= node->end(); ++i ){
      (*i)->show( indent );
      showTree( *i , indent );
    }
  }
}



void XMLNode::appendChild(XMLNode* node){
  node->setParent(this);
  push_back(node);
}

void XMLNode::appendSibling(XMLNode* node){
  node->parent=parent;
  if(parent!=0){
    parent->push_back(node);
  }
}

XMLNode::iterator XMLNode::lookup(){
  if(parent!=0){
    XMLNode::iterator i=parent->begin();
    while((*i!=this)&&(i!=parent->end())){
      ++i;
    }
    return i;
  }
  return end();
}



//returns the nextSibling
XMLNode* XMLNode::nextSibling(){
  if (parent) {
    XMLNode::iterator i = lookup();
    if ( i!=parent->end() ){
      ++i;
      // must check after i has been incremented
      // to make sure i isn't at the end before
      // returning the contained pointer value
      if (i!=parent->end()){
        return *i;
      }
    }
  }
  return NULL;
}



XMLNode* XMLNode::prevSibling(){
  if (parent){
    XMLNode::iterator i = lookup();
    // Must make sure we aren't at beginning as well
    // or we can crash when decrementing since we shouldn't
    // decrement before the start of the list
    if ( ( i!=parent->end() ) && ( i!=parent->begin() ) ){
      --i;
      return *i;
    }
  }
  return NULL;
}




// returns first child of a node
XMLNode* XMLNode::firstChild(){
        // Must make sure we aren't empty first!!!
        if (!empty())
        {
                XMLNode::const_iterator child = begin();
                return *child;
        }
        return NULL;
}



//returns last child of a node
XMLNode* XMLNode::lastChild()
{
        // Must make sure we aren't empty first and use rbegin() and
        // a reverse iterator...
        if (!empty())
        {
                XMLNode::const_reverse_iterator child= rbegin();
                return *child;
        }
        return NULL;
}


