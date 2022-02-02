
#include "attribute.h"

Attribute::Attribute(const string& str){
  string::size_type pos=str.find("=");
  if(pos!=string::npos){
    name=str.substr(0,pos);
    value=str.substr(pos+1, str.length()-pos -1);
  }
  else{
    name=str;
  }
}

Attribute::Attribute(const string& n, const string& v){
  name=n;
  value=v;
}


Attribute::Attribute(const Attribute& a){
  name=a.name;
  value=a.value;
}



Attribute& Attribute::operator=(const Attribute& a){
  name=a.name;
  value=a.value;
  return *this;
}


