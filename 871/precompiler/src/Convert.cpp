/*=============================================================================
author        :Walter Schreppers
filename      :Convert.h
created       :17/04/2000
modified      :25/02/2002
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#include "Convert.h"


/*-----------------------------------------------------------------------------
name        :Convert
description :constructor, set the private locals for the new Convert class
parameters  :const Var& v, Matcher* varMatcher, int line
return      :Convert
exceptions  :/
algorithm   :trivial 
-----------------------------------------------------------------------------*/
Convert::Convert(const Var& v,Matcher* varMatcher, ParserConf* pConf){
  fVar=v;              //set private local fVar (the var to be converted)
  fMatcher=varMatcher; //set private local the matcher (to do conversion)
  strValue=v.getValue();
  varType=v.getType();
  fLine=pConf->getLine();
  fConf=pConf;
}




/*-----------------------------------------------------------------------------
name        :getFullAssignType
description : give back a string representing the type such as 'long double'
              which can be found in the rhs of an assignment.
parameters  :const fVar&, unsigned int typePos& (returns the pos up to where type is matched)
return      :string
exceptions  :/
algorithm   :inspired by the types Parser::setVarType(ofstream& out) member
which was for the source lookup. This does almost the same but for the RHS of assignments.
this expands the string to as large as possible to find suitable convert rule
so that things as long double are matched...
-----------------------------------------------------------------------------*/
string Convert::getFullAssignType( const Var& fVar, unsigned int& typePos ){
  unsigned int pos = typePos;
  unsigned int prevPos = pos;
  
  Matcher varMatcher();
  string testType = fVar.valQueue[pos].text;
  string vType = testType;
  pos++;

  MatchElem m=fMatcher->assignLookup( testType );
  MatchElem match;

  while(m.mtype!=mError){
    prevPos = pos-1;
    vType=testType;
    typePos = pos;
    match=m;
  
    pos++;
    if(fVar.valQueue[pos].token==white_space){
      pos++;
      testType+=" "+fVar.valQueue[pos].text; //WARNING we add any token here watch out for bad lookups...
    }
    else{
      testType+=fVar.valQueue[pos].text;
    }
    m=fMatcher->assignLookup(testType);
  }

  typePos = prevPos;
  return vType;
}


/*-----------------------------------------------------------------------------
name        :getValue
description :conversion of the strValue string to string which 
             uses BigInt,MpIeee,...
parameters  :/
return      :string
exceptions  :/
algorithm   :Walk through string and copy variable names literally. 
             Copy functions literally including their parameters by
              using getFunctionParams
             
             constant value's have to be changed by using getConstant
             copy operands +,*,-,/ literally
             
             we take brackets and append them to previous newText, 
             so for instance 'array' '[' '3' '+' '2' ']'
             will be appended into "array[3+2]" and passed to convConf              
-----------------------------------------------------------------------------*/
string Convert::getValue(){
  string targetName=fVar.getTargetType();
  string newValue="";
    
  //for each token find a conversion path leading to targetName
  vector<ConvertElem> converts;
  converts.clear();

  ConvertConfig convConf(fMatcher,fConf);
  
  unsigned int pos=0,typePos=0;
  string newText;
  
  while( pos<fVar.valQueue.size() ){
    newText=fVar.valQueue[pos].text;
    
    typePos=pos;
    
    if( fVar.valQueue[pos+1].token == open_bracket ){ //copy everything between bracket's literally
      pos++;
      while( ( fVar.valQueue[pos].token != close_bracket ) && ( pos<fVar.valQueue.size() ) ){
        newText+=fVar.valQueue[pos].text;
        pos++;
      }
      if( pos<fVar.valQueue.size() ) newText+=fVar.valQueue[pos].text;
      newValue+=newText;
      newText="";
    }
    
    if( fVar.valQueue[pos+1].token == function_argument ){ //copy function arguments literally
      pos++;
      while( ( fVar.valQueue[pos].token == function_argument ) && ( pos<fVar.valQueue.size() ) ){
        newText+=fVar.valQueue[pos].text;
        pos++;
      }
      if( pos<fVar.valQueue.size() ) newText+=fVar.valQueue[pos].text;
      newValue+=newText;
      newText="";
    }
    #ifdef _DEBUG_
      cout<<"'"<<newText<<"' ->"<<fVar.valQueue[typePos].token<<endl;
    #endif
    
    if( newText!="" ){
    try{ 
      if( fVar.valQueue[typePos].token == sym_word ){ // added 11/5/2005, expanding value as large as possible
        converts=convConf.path( getFullAssignType( fVar, typePos ), targetName );
        pos = typePos;
        #ifdef _DEBUG_
          cout << "pos="<<pos<<" typePos="<<typePos<<endl;
        #endif
      } 
      else{
        converts=convConf.path( fVar.valQueue[typePos], targetName );
      }

      for(vector<ConvertElem>::iterator i=converts.begin();i!=converts.end();i++){
        newText=convConf.execute( *i, newText );
      }
    }
    catch(Warning w){
      cerr<<"Warning : "<<w.fWarning<<" at line: "<<fLine<<endl;
    }

    newValue+=newText; 
    }
    
    pos++;
    
  }
  return newValue;
}








/*-----------------------------------------------------------------------------
name        :getFunctionParams
description :copies everything until the close of the function paramers literally
parameters  :int& pos,const string& inVal
return      :/
exceptions  :/
algorithm   :keep count of openbrakkets stop copying when it becomes 0
-----------------------------------------------------------------------------*/
string Convert::getFunctionParams(int& pos,const string& inVal){
  string strParams="";
  int openBraces=1;
  while( ( pos < (int) inVal.length() ) && ( openBraces > 0 ) ){
    if(inVal[pos]=='(') openBraces++;
    if(inVal[pos]==')') openBraces--;
    strParams+=inVal[pos];
    pos++;
  }
  return strParams;  

}




