/*=============================================================================
author        :Walter Schreppers
filename      :PrintArg.cpp
created       :/
modified      :11/05/2001
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#include "PrintArg.h"


/*-----------------------------------------------------------------------------
name        :PrintArg
description :constructor
parameters  :const string& arg, pargtype t
return      :PrintArg
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
PrintArg::PrintArg(const string& arg,pargtype t){
  this->argument=arg;
  this->atype=t;
}



/*-----------------------------------------------------------------------------
name        :PrintArg
description :destructor
parameters  :/
return      :PrintArg
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
PrintArg::~PrintArg(){
}




/*-----------------------------------------------------------------------------
name        :getType()
description :return argument type
parameters  :/
return      :ostream&
exceptions  :/
algorithm   :/
-----------------------------------------------------------------------------*/
pargtype PrintArg::getType() const{
  return this->atype;
}


/*-----------------------------------------------------------------------------
name        :getArg()
description :return argument (this is a string representation)
parameters  :/
return      :string
exceptions  :/
algorithm   :/
-----------------------------------------------------------------------------*/
string PrintArg::getArg() const{
  return this->argument;
}

