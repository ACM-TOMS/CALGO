/*=============================================================================
author        :Walter Schreppers
filename      :PreParser.cpp
created       :15/05/2001
modified      :15/05/2001
version       :2
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#include "PreParser.h"


/*-----------------------------------------------------------------------------
name        :PreParser
description :constructor
parameters  :ParserConf*
return      :/
exceptions  :/
algorithm   :use constructor of Parser
-----------------------------------------------------------------------------*/
PreParser::PreParser():Parser(){
}



/*-----------------------------------------------------------------------------
name        :PreParser
description :constructor
parameters  :ParserConf*
return      :/
exceptions  :/
algorithm   :use constructor of Parser
-----------------------------------------------------------------------------*/
PreParser::PreParser(ParserConf *parserConf):Parser(parserConf){
}



/*-----------------------------------------------------------------------------
name        :PreParser
description :destructor
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
PreParser::~PreParser(){
}


/*-----------------------------------------------------------------------------
name        :writeVar
description :
parameters  :ofstream& out, const string& varname
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
void PreParser::writeVar(ofstream& out,const string& varname ){
  if(pConf->getStrCurFunction().size()==0){
    out<<"/\t\t"; // a '/' as function name meaning it is empty (global var)
  }
  else{
    string func=pConf->getStrCurFunction();
    if(func.find("::") != string::npos) func.erase(0,func.find("::")+2); //strip classname::
    out<<func<<"\t\t";
  }
  string var=varname;
  if(var.find("::") != string::npos) var.erase(0,var.find("::")+2);
  out<<var<<"\t\t";
  out<<getOriginalType();
}




/*-----------------------------------------------------------------------------
name        :parseVar
description :overloading from base class Parser
parameters  :ofstream& out, outputType t
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
void PreParser::parseVar(ofstream& out,string ){
  for(Variables::iterator p=constructs.begin();p!=constructs.end();p++){
    writeVar(out,p->getName());
    out<<p->getSpec()<<endl;
  }
}





/*-----------------------------------------------------------------------------
name        :parseAssignment
description :overloading from base class Parser
parameters  :ofstream& out, bool forloop
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::parseAssignment( ofstream& ){
  //nothing to be done  
}





/*-----------------------------------------------------------------------------
name        :setFloatPrecision
description :overloading from base class Parser
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::setFloatPrecision( ofstream& ){
  //nothing to be done  
}



/*-----------------------------------------------------------------------------
name        :moveToken
description :overloading from base class Parser
             here we just pop the token and don't write it to out
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::moveToken( ofstream& ){
  updateLocals();
  pQueue.popfirst();
}

/*-----------------------------------------------------------------------------
name        :moveTokens
description :overloaded from base class Parser here we don't write to out
parameters  :int& pos,ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::moveTokens( unsigned int& pos, ofstream& out ){
  while( ( pos > 0 ) && ( pQueue.size() > 0 ) ){
    moveToken(out);
    pos--;
  }
}




/*-----------------------------------------------------------------------------
name        :copyToken
description :overloading from base class Parser
             we don't write anything to out
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::copyToken( ofstream& ){
  //do nothing
}

/*-----------------------------------------------------------------------------
name        :dumpq
description :dump the queue as in Parser but don't put it onto ofstream out
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::dumpq( ofstream& out ){
  while(pQueue.size()>0) moveToken(out);
}


/*-----------------------------------------------------------------------------
name        :writeStrType
description :overloaded virtual member from base class Parser here it has
             to do nothing
parameters  :ofstream& out,types t
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::writeStrType( ofstream&, const string& ){
  //do nothing
}

/*-----------------------------------------------------------------------------
name        :writeOriginalType
description :overloaded virtual member from base class Parser here it has
             to do nothing
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::writeOriginalType( ofstream& ){
  //do nothing
}






/*-----------------------------------------------------------------------------
name        :convertPrintf
description :overloading from base class Parser
             we don't write anything to out and just remove the statement
             from the pQueue stack
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PreParser::convertPrintf( ofstream& ){
  //we just dump the printf statement
  while( pQueue.size() > 0 ) pQueue.popfirst();
}


