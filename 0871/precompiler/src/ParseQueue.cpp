/*=============================================================================
author        :Walter Schreppers
filename      :ParseQueue.cpp
created       :15/04/2001
modified      :/
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#include "ParseQueue.h"

/*-----------------------------------------------------------------------------
name        :ParseQueue
description :constructor
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseQueue::ParseQueue(){
  changed=0;
  this->clear();
}



/*-----------------------------------------------------------------------------
name        :ParseQueue
description :copy constructor, use operator=
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseQueue::ParseQueue(const ParseQueue& p):vector<ParseElem>(p){
  if(this!=&p) *this=p;
}

/*-----------------------------------------------------------------------------
name        :operator=
description :operator=
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseQueue& ParseQueue::operator=(const ParseQueue& p){
  this->clear();
  
  for( ParseQueue::const_iterator i=p.begin(); i!=p.end(); ++i ){
    this->push_back( *i );
  }
  
  return *this;
}


/*-----------------------------------------------------------------------------
name        :~ParseQueue
description :destructor
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseQueue::~ParseQueue(){
  this->clear();
}

/*-----------------------------------------------------------------------------
name        :add
description :add element into queue
parameters  :int token,string text
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
void ParseQueue::add(int token,string text){
  changed=1;
  ParseElem p( token, text );
  this->push_back( p );
}


/*-----------------------------------------------------------------------------
name        :dump
description :dump queue onto ofstream
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
void ParseQueue::dump(ofstream& out){
  if(this->size()>0) changed=1;
  while( this->size() > 0 ){
    out << (*this)[0].text;
    this->erase( this->begin() );
  }
}


/*-----------------------------------------------------------------------------
name        :popfirst
description :remove first item from queue
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
void ParseQueue::popfirst(){ 
  if( this->size() > 0 ){
    this->erase( this->begin() );
    changed=1;
  }
  else{
    //cerr<<"ParseQueue::popping empty queue!"<<endl;
  }
}


/*-----------------------------------------------------------------------------
name        :poptext
description :remove first item from queue and return it's text
             return empty string if queue is empty
parameters  :/
return      :string
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
string ParseQueue::poptext(){
  string retStr=""; 
  
  if( this->size() > 0 ){
    retStr=(*this)[0].text;
    this->erase( this->begin() );
    changed=1;
  }
  
  return retStr;
}



/*-----------------------------------------------------------------------------
name        :pop
description :remove first item from queue and return it
             return empty element if queue is empty
parameters  :/
return      :ParseElem
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseElem ParseQueue::pop(){
  ParseElem retElem(-1,""); 
  
  if( this->size() > 0 ){
    retElem=*( this->begin() );
    this->erase( this->begin() );
    changed=1;
  }
  
  return retElem;
}



/*-----------------------------------------------------------------------------
name        :get
description :return ParseElem at pos 
             return empty element with token = -1 if out of range
parameters  :unsigned int pos
return      :ParseElem
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseElem ParseQueue::get(unsigned int pos){
  if( pos < this->size() ){
    return (*this)[pos];
  }
  else{
    return ParseElem( -1, "" );
  }
}



/*-----------------------------------------------------------------------------
name        :getToken
description :return token at pos, and -1 if pos is out of range
parameters  :unsigned int pos
return      :int
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
int ParseQueue::getToken(unsigned int pos){

  if( pos < this->size() ){
    return (*this)[pos].token;
  }
  else{
    return -1; //a non defined token
  }

}




/*-----------------------------------------------------------------------------
name        :getText
description :return string text at pos, and "" if pos is out of range
parameters  :unsigned int pos
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
string ParseQueue::getText(unsigned int pos){
  if( pos < this->size() ){
    return (*this)[pos].text;
  }
  else{
    return "";
  }
}

