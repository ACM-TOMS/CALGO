/*=============================================================================
author        :Walter Schreppers
filename      :Var.cpp
created       :/
modified      :25/02/2002
version       :4
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/
#include "Var.h"


/*-----------------------------------------------------------------------------
name        :Var
description :default constructor
parameters  :
return      :Var
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
Var::Var(){
  init();  
};

/*-----------------------------------------------------------------------------
name        :init
description :initialize private locals
parameters  :
return      :
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
void Var::init(){
  fName="";
  fValue="";
  fSpec="";
  fArrayArg="";

  fType=tEmpty;
  fSrcType="";
  fTargType="";
  fMatchType="";
  valQueue.clear();
}



/*-----------------------------------------------------------------------------
name        :Var
description :constructor
parameters  :const string& name,const string& value, types t
return      :Var
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var::Var(const string& name,const string& value,types t){
  init();
	fName=name;
	fValue=value;
	fType=t;
}


/*-----------------------------------------------------------------------------
name        :Var
description :constructor
parameters  :const string& name,const string& value, const string& spec,types t
return      :Var
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var::Var(const string& name,const string& value,const string& spec,types t){
  init();
	fName=name;
	fValue=value;
	fSpec=spec;
	fType=t;
}

/*-----------------------------------------------------------------------------
name        :Var
description :constructor
parameters  :const string& name, 
             const string& value, 
             const string& spec,
             const string& arrayarg,
             types t
return      :Var
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var::Var( const string& name,
          const string& value,
          const string& spec,
          const string& arrayarg,
          types t){

  init();
	fName=name;
	fValue=value;
	fSpec=spec;
	fArrayArg=arrayarg;
	fType=t;

}


/*-----------------------------------------------------------------------------
name        :Var
description :constructor
parameters  :const string& name, 
             const string& value, 
             const string& spec,
             const string& arrayarg,
             types t,
             const string& match
             
return      :Var
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var::Var( const string& name,
          const string& value,
          const string& spec,
          const string& arrayarg,
          types t,
          const string& match){

  init(); 
	fName=name;
	fValue=value;
	fSpec=spec;
	fArrayArg=arrayarg;
	fType=t;
  fMatchType=match;

}


/*-----------------------------------------------------------------------------
name        :Var
description :constructor
parameters  :const string& name, 
             const string& value, 
             const string& spec,
             const string& arrayarg,
             types t,
             const string& match,
             const string& src,
             const string& trg
             
return      :Var
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var::Var( const string& name,
          const string& value,
          const string& spec,
          const string& arrayarg,
          types t,
          const string& match,
          const string& src,
          const string& trg ){

  //init(); //not needed, everything is initialised here
	fName=name;
	fValue=value;
	fSpec=spec;
	fArrayArg=arrayarg;
	fType=t;
  fMatchType=match;
  fSrcType=src;
  fTargType=trg;
}


/*-----------------------------------------------------------------------------
name        :Var
description :constructor
parameters  :types vType,
             const string& vSpec,
             const string& vName,
             const string& matchType,
             const string& src,
             const string& trg
return      :Var
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var::Var( types vType, 
          const string& vSpec,
          const string& vName,
          const string& matchType,
          const string& src,
          const string& trg ){
  init();
  fType = vType;
  fName = vName;
  fSpec = vSpec;
  fMatchType = matchType;
  fSrcType = src;
  fTargType = trg;
}





/*-----------------------------------------------------------------------------
name        :Var
description :copy constructor
parameters  :const Var&
return      :Var
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var::Var(const Var& v){
  if(&v != this){
    *this=v; //use =operator
  }
}

/*-----------------------------------------------------------------------------
name        :Var
description :operator=
parameters  :const Var&
return      :Var&
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var& Var::operator=(const Var& v){
  fName=v.fName;
  fValue=v.fValue;
  fSpec=v.fSpec;
  fType=v.fType;
  fArrayArg=v.fArrayArg;
  fMatchType=v.fMatchType;
  fSrcType=v.fSrcType;
  fTargType=v.fTargType;
  valQueue=v.valQueue;
  
  return *this;
}


/*-----------------------------------------------------------------------------
name        : operator<<
description : output operator
parameters  : ostream& os, const Gen& m
return      : ostream&
exceptions  :
algorithm   : trivial
-----------------------------------------------------------------------------*/
ostream& operator<<( ostream& os, const Var& v ){
  os << " name="      << v.fName;
  os << " value="     << v.fValue;
  os << " spec="      << v.fSpec;
  os << " type="      << v.fType;
  os << " arrayArg="  << v.fArrayArg;
  os << " matchType=" << v.fMatchType;
  os << " srcType="   << v.fSrcType;
  os << " targType="  << v.fTargType;
//  os << " valQueue="  << v.valQueue;

  return os;
}


/*-----------------------------------------------------------------------------
name        :Var
description :destructor
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var::~Var(){
}


/*-----------------------------------------------------------------------------
name        :getName
description :get name of the variable
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Var::getName() const{
	return fName;
}


/*-----------------------------------------------------------------------------
name        :getValue
description :get value of the variable as a string
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Var::getValue() const {
	return fValue;
}


/*-----------------------------------------------------------------------------
name        :getSpec
description :get specification string (pointer/references) of the variable
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Var::getSpec() const{
	return fSpec;
}

/*-----------------------------------------------------------------------------
name        :getArrayArg
description :get the array args ex. float b[5]=... then arrayarg= "[5]"
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Var::getArrayArg() const{
	return fArrayArg;
}


/*-----------------------------------------------------------------------------
name        :getMatchType
description :get matcher type string
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Var::getMatchType() const{
	return fMatchType;
}



/*-----------------------------------------------------------------------------
name        :getTargType
description :get target type string
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Var::getTargetType() const{
	return fTargType;
}


/*-----------------------------------------------------------------------------
name        :getSrcType
description :get source type string
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Var::getSrcType() const{
	return fSrcType;
}




/*-----------------------------------------------------------------------------
name        :getType
description :get type of the variable
parameters  :/
return      :types
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
types Var::getType() const{
	return fType;
}

/*-----------------------------------------------------------------------------
name        :setName
description :set name of variable
parameters  :const string&
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setName(const string& s){
	fName=s;
}


/*-----------------------------------------------------------------------------
name        :setName
description :set name of variable
parameters  :char*
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setName(char* cp){
	fName=string(cp);
}


/*-----------------------------------------------------------------------------
name        :setValue
description :set value of variable
parameters  :const string&
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setValue(const string& s){
	fValue=s;
}

/*-----------------------------------------------------------------------------
name        :setValue
description :set value of variable
parameters  :char*
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setValue(char* cp){
	fValue=string(cp);
}


/*-----------------------------------------------------------------------------
name        :setSpec
description :set Specification of variable
parameters  :const string&
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setSpec(const string& s){
	fSpec=s;
}

/*-----------------------------------------------------------------------------
name        :setSpec
description :set Specification of variable
parameters  :char*
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setSpec(char* cp){
	fSpec=string(cp);
}



/*-----------------------------------------------------------------------------
name        :setArrayArg
description :set ArrayArg of variable
parameters  :const string&
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setArrayArg(const string& s){
	fArrayArg=s;
}

/*-----------------------------------------------------------------------------
name        :setArrayArg
description :set ArrayArg of variable
parameters  :char*
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setArrayArg(char* cp){
	fArrayArg=string(cp);
}


/*-----------------------------------------------------------------------------
name        :setStrType
description :set type string of variable
parameters  :const string&
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setMatchType(const string& s){
	fMatchType=s;
}



/*-----------------------------------------------------------------------------
name        :setType
description :set type of the variable
parameters  :types
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Var::setType(types t){
	fType=t;
}

