/*=============================================================================
author        :Walter Schreppers
filename      :Variables.cpp
created       :/
modified      :25/02/2002
version       :2
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/
#include "Variables.h"


/*-----------------------------------------------------------------------------
name        :insert
description :insert a variable into the list
parameters  :const Var&
return      :/
exceptions  :catches any error that might occur
             and throws Error with an explanation string
algorithm   :trivial 
-----------------------------------------------------------------------------*/
void Variables::insert(const Var& v){
	try {
	  push_back(v); 
  }
  catch(...){
    throw Error("Variables::insert : Can't insert variable");
  }
}



/*-----------------------------------------------------------------------------
name        :insert
description :insert a variable into the list
parameters  :types vType,const string& vSpec, const string& vName
return      :/
exceptions  :catches any error that might occur
             and throws Error with an explanation string
algorithm   :trivial 
-----------------------------------------------------------------------------*/
void Variables::insert(types vType,const string& vSpec,const string& vName){
	try {
	  Var v(vName, "", vSpec, vType);
	  push_back(v); 
  }
  catch(...){
    throw Error("Variables::insert : Can't insert variable");
  }
}



/*-----------------------------------------------------------------------------
name        :insert
description :insert a variable into the list
parameters  :types vType,const string& vSpec, const string& vName,
             const string& matchType
return      :/
exceptions  :catches any error that might occur
             and throws Error with an explanation string
algorithm   :trivial 
-----------------------------------------------------------------------------*/
void Variables::insert( types vType,
                        const string& vSpec,
                        const string& vName,
                        const string& matchType ){
	try {
	  Var v( vName, "", vSpec, "", vType, matchType );
	  push_back(v); 
  }
  catch(...){
    throw Error("Variables::insert : Can't insert variable");
  }
}


/*-----------------------------------------------------------------------------
name        :insert
description :insert a variable into the list
parameters  :types vType,const string& vSpec, const string& vName,
             const string& matchType, const string& src, const string& trg
return      :/
exceptions  :catches any error that might occur
             and throws Error with an explanation string
algorithm   :trivial 
-----------------------------------------------------------------------------*/
void Variables::insert( types vType,
                        const string& vSpec,
                        const string& vName,
                        const string& matchType,
                        const string& src,
                        const string& trg ){
	try {
	  Var v( vName, "", vSpec, "", vType, matchType, src, trg );
	  push_back(v); 
  }
  catch(...){
    throw Error("Variables::insert : Can't insert variable");
  }
}





/*-----------------------------------------------------------------------------
name        :openBracket
description :insert tOpenBracket into list
parameters  :/
return      :/
exceptions  :catches insert Error and throws a modified Error
algorithm   :trivial 
-----------------------------------------------------------------------------*/
void Variables::openBracket(){
  try{
	  insert(tOpenBracket, "", "");
	}
	catch(...){
    throw Error("Variables::openBracket : Cant insert tOpenBracket");
  }
}


/*-----------------------------------------------------------------------------
name        :closeBracket
description :remove all variables upto the first occurance of an openbracket or
             untill list is empty
parameters  :/
return      :/
exceptions  :throw exception Error with an explanation string if
             closeBracket is called on an empty list
algorithm   :trivial 
-----------------------------------------------------------------------------*/
void Variables::closeBracket(){
	  Variables::iterator iter = end();
	  if (size()==0){
	    throw Error("Variables::closeBracket : Cant remove Variable list is empty!");
	  }
	  else if(size()==1){
	    clear();
	  }
	  else{
	    while ( (iter!=begin()) && (iter->getType()!=tOpenBracket) ) {
        iter--;
        if (iter!=begin()){
          pop_back();
        }
		  }
	  }
}


/*-----------------------------------------------------------------------------
name        :expandLastValue
description :this function catenates a string to the value string of the
             variable that was last inserted into the list
parameters  :const string&
return      :Var
exceptions  :throw error Variables::expandLastValue if unsuccessfull
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Variables::expandLastValue(const string& t){
  try{
    Variables::iterator p=end();
    p--;
    string tmp=p->getValue();
    tmp=tmp+t;
    p->setValue(tmp);
  }
  catch(...){
    throw("Variables::expandLastValue: can't add variable value into list");
  }
}



/*-----------------------------------------------------------------------------
name        :expandLastValue
description :this copies a ParseElem to the parseque of the 
             variable that was last inserted into the list
parameters  :const ParseElem&
return      :Var
exceptions  :throw error Variables::expandLastValue if unsuccessfull
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Variables::expandLastValue(const ParseElem& t){
  try{
    Variables::iterator p=end();
    p--;
    p->valQueue.push_back(ParseElem(t));
  }
  catch(...){
    throw("Variables::expandLastValue: can't add variable value into list");
  }
}


/*-----------------------------------------------------------------------------
name        :expandLastName
description :this function catenates a string to the name of the
             variable that was last inserted into the list
parameters  :const string&
return      :Var
exceptions  :throw error Variables::expandLastValue if unsuccessfull
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Variables::expandLastName(const string& t){
  try{
    Variables::iterator p=end();
    p--;
    string tmp=p->getName();
    tmp=tmp+t;
    p->setName(tmp);
  }
  catch(...){
    throw("Variables::expandLastValue: can't add variable value into list");
  }
}


	
/*-----------------------------------------------------------------------------
name        :expandLastArrayArg
description :this function catenates a string to the ArrayArg string of the
             variable that was last inserted into the list
parameters  :char *name,char* value,types t
return      :Var
exceptions  :throw error Variables::expandLastValue if unsuccessfull
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Variables::expandLastArrayArg(const string& t){
  try{
    Variables::iterator p=end();
    p--;
    string tmp=p->getArrayArg();
    tmp=tmp+t;
    p->setArrayArg(tmp);
  }
  catch(...){
    throw("Variables::expandLastValue: can't add variable value into list");
  }
}




/*-----------------------------------------------------------------------------
name        :show
description :show contents of the list 
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial 
-----------------------------------------------------------------------------*/
void Variables::show(){
	Variables::iterator iter = end();
	while (iter!=begin()) {
		iter--;
		if(iter->getType()==tOpenBracket) cout<<"Open bracket";
		cout << iter->getName() << "\n";
	}
}







/*-----------------------------------------------------------------------------
name        :lookup
description :find a variable, return a Var with type tError if it's not found
parameters  :const string& name
return      :/
exceptions  :Variables::lookup in scope failed
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var Variables::lookup(const string& name){
  try{
  	Var result("","",tError);
  	bool bStop = 0;
  	Variables::iterator iter = end();
  	while ((bStop==0) && (iter!=begin())) {
  		iter--;
  		if (iter->getName()==name) {
  			result=*iter;
  			bStop = 1;
  		}
  	}
  	return result;
  }
  catch(...){
    throw("Variables::lookup: Lookup of variable in scope failed");
  }
}



