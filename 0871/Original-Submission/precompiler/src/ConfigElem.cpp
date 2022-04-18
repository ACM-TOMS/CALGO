/*=============================================================================
author        :Walter Schreppers
filename      :ConfigElem.cpp
created       :/
modified      :15/05/2000
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#include "ConfigElem.h"

/*-----------------------------------------------------------------------------
name        :ConfigElem
description :Constructor
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
ConfigElem::ConfigElem(){
	fVar="";
	fProcedure="";
}

/*-----------------------------------------------------------------------------
name        :ConfigElem
description :Constructor
parameters  :const string& proc, const string& var
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
ConfigElem::ConfigElem(const string& proc, const string& var){
  fVar=var;
  fProcedure=proc;
}


/*-----------------------------------------------------------------------------
name        :setProcedure
description :set the fProcedure string
parameters  :const string& s
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ConfigElem::setProcedure(const string& s){
  fProcedure=s;
}
		
		

/*-----------------------------------------------------------------------------
name        :setProcedure
description :set the fProcedure string
parameters  :char* s
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ConfigElem::setProcedure(char* s){
  fProcedure=string(s);
}







/*-----------------------------------------------------------------------------
name        :setVar
description :set the fVar string
parameters  :const string& s
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ConfigElem::setVar(const string& s){
  fVar=s;
}
		
		

/*-----------------------------------------------------------------------------
name        :setVar
description :set the fVar string
parameters  :char* s
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ConfigElem::setVar(char* s){
  fVar=string(s);
}




/*-----------------------------------------------------------------------------
name        :getVar
description :return fVar string
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string ConfigElem::getVar(){
  return fVar;
}


/*-----------------------------------------------------------------------------
name        :getProcedure
description :return fProcedure string
parameters  :/
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string ConfigElem::getProcedure(){
  return fProcedure;
}



