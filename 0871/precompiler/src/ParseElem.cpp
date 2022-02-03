/*=============================================================================
author        :Walter Schreppers
filename      :ParseElem.cpp
created       :14/04/2001
modified      :/
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#include "ParseElem.h"

/*-----------------------------------------------------------------------------
name        :ParseElem
description :constructor
parameters  :int token, string text
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseElem::ParseElem(int token,const string& text){
  this->token=token;
  this->text=text;
}


/*-----------------------------------------------------------------------------
name        :ParseElem
description :copy constructor
parameters  :const ParseElem&
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseElem::ParseElem(const ParseElem& p){
  if(this!=&p){
    *this=p;
  }
}






/*-----------------------------------------------------------------------------
name        :operator=
description :
parameters  :const ParseElem&
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseElem& ParseElem::operator=(const ParseElem& p){
    text=p.text;
    token=p.token;
    return *this;
}


/*-----------------------------------------------------------------------------
name        :~ParseElem
description :destructor
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
ParseElem::~ParseElem(){
  text="";
}




