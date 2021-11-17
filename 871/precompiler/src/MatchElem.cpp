#include "MatchElem.h"


void MatchElem::init(){
  name="";
  keyword="";
  token="";
  mtype=mEmpty;
}


MatchElem::MatchElem(){
  init();
}

MatchElem::MatchElem(matchtypes t){
  init();
  mtype=t;
}


MatchElem::MatchElem(const string& name){
  init();
  this->name=name;
}


MatchElem::MatchElem(const MatchElem& m){
  *this=m;
}


MatchElem::MatchElem(const string& name, const string& keyword){
  init();
  this->name=name;
  this->keyword=keyword;
  mtype=mKeyword;
}


MatchElem::MatchElem(const string& name, const string& val,matchtypes t){
  init();
  this->name=name;
  if(t==mToken){
    this->token=val;
  }
  else if(t==mKeyword){
    this->keyword=val;
  }
  mtype=t;
}


void MatchElem::setKeyword(const string& s){
  mtype=mKeyword;
  keyword= s;
}

MatchElem& MatchElem::operator=(const MatchElem& m){
  if(this!=&m){
    name=m.name;
    keyword=m.keyword;
    token=m.token;
    mtype=m.mtype;
  }
  return *this;
}



