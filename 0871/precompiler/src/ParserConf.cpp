/*=============================================================================
author        :Walter Schreppers
filename      :ParserConf.cpp
created       :/
modified      :25/02/2002
version       :5
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/
#include "ParserConf.h"


/*-----------------------------------------------------------------------------
name        :ParserConf
description :constructor
parameters  :/
return      :/
exceptions  :/
algorithm   :initialize the variables and pointers
-----------------------------------------------------------------------------*/
ParserConf::ParserConf(){
    line_nr=1;
    openParentheses=0;        //this is the number of '('
    openBracket=0;       //this is the number of '{'
    insideFuncParam=false;   //=TRUE if we're inside the parameter declaration of a function
    insideForParam=false;
    insideFunction=false;

    bSkipFor=false;  //false=replace int in for with BigInt, true=don't replace. Also look at setSkipFor function	
    while(scopes.size()>0) scopes.pop();  
    vars=new Variables;
    
    reswords=new Reserved; //with reswords we can check for c++ reserved words to avoid adding variables of type tOther
  
    skipVar=new ConfigFile;   //to be able to skip convertion in certain variables or procedures
                              //see also the loadConfigFile function below !
  
    functions=new Variables;
    currentFunc=Var("","",tOther);
    convertTree = 0;
}



/*-----------------------------------------------------------------------------
name        :ParserConf
description :destructor
parameters  :/
return      :/
exceptions  :/
algorithm   :release memory, delete used pointers
-----------------------------------------------------------------------------*/
ParserConf::~ParserConf(){
  vars->clear();      
  delete vars;
  
  reswords->clear();
  delete reswords;
  
  skipVar->clear();   
  delete skipVar;
  
  functions->clear(); 
  delete functions;
  
  if( convertTree!= 0 ) delete convertTree;
}





		
/*-----------------------------------------------------------------------------
name        :newScope
description :Start a new scope by inserting a openBracket in vars
             unless this has already been done indicated by the boolean
             insideFuncParam.
             This function should be executed whenever a { is found
parameters  :/
return      :/
exceptions  :vars->openBracket might throw an exception Error but this will
             be caught in the precompile.cpp main() so that lexing stops
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::newScope(){
  #ifdef _DEBUG_
    cerr<<"found open '{' at "<<getLine()<<endl;
  #endif
  
  if( ( scopes.size()!=0 ) && ( scopes.top() == sSemi ) ){ //function or forloop where scope is already opened
    scopes.pop();
    scopes.push( sBracket );
  }
  else{
    scopes.push( sBracket );
    vars->openBracket();
    openBracket++;
  }
  
  insideFuncParam=false;
  insideForParam=false;
}
  
  
  
  
  
/*-----------------------------------------------------------------------------
name        :newFunction
description :Also start a new scope by inserting a openBracket in vars
               set the insideFuncParam=1 to notify newScope
parameters  :/
return      :/
exceptions  :vars->openBracket might throw an exception Error but this will
               be caught in precompile.cpp main() so that lexing stops
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::newFunction(const string& name,
                             types ftype,
                             const string& matchType,
                             const string& srcType,
                             const string& targetType){
  scopes.push(sFunction);

  //name,value,spec,arrayarg,types, match, src, target
  Var vFunc(name, "", "", "", ftype,matchType,srcType,targetType );
  functions->push_back(vFunc);

  if( !insideFunction ){
    currentFunc=vFunc;
    insideFunction = true;
  }
  
#ifdef _DEBUG_
  cerr<<"found function: "<<name <<" at line: "<<line_nr<<endl;
#endif

  vars->openBracket();   /* the parameters of a function belong to the scope of that function*/
  openBracket++;
  
  insideFuncParam=true;
  scopes.push(sSemi);
}

/*-----------------------------------------------------------------------------
name        :newForLoop
description :We use the same principle as for functions SEE: newFunction()
             except that the function name is not set here for obvious reasons.
             and the bool insideForLoop is set to true
parameters  :/
return      :/
exceptions  :vars->openBracket might throw an exception Error but this will
               be caught in precompile.cpp main() so that lexing stops
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::newForLoop(){
  vars->openBracket();  //we implement the new ANSI c++ standard in which variables of FOR loops
                        //belong to the scope of the FOR loop itself like functions
  openBracket++;
  insideForParam=true;
  
  scopes.push(sSemi);
  
#ifdef _DEBUG_
  cerr<<"found for loop: at line: "<<line_nr<<endl;
#endif

}


string ParserConf::scopeName( scopetypes t ){
  switch( t ){
    case sSemi     : return "';'"; break;
    case sBracket  : return "'{'"; break;
    case sFunction : return "'function'"; break;
    default        : return "UNKNOWN"; break;
  }
}

/*-----------------------------------------------------------------------------
name        :closeScope
description :Remove variables of last scope by calling vars->closeBracket()
             also close pending scopes given in closeScopesOnSemi
parameters  :/
return      :/
exceptions  :vars->closeBracket might throw an exception Error but this will
             be caught in precompile.cpp main() so that lexing stops
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::closeScope(){
  
  decOpenBracket(); //close scope of bracket
  vars->closeBracket();
  
  if(scopes.size()>0){
    if(scopes.top() != sBracket){
      cerr<<"Warning: ParserConf::closeScope: there might be a scoping problem at line "<<getLine()<< " scope top = " << scopeName( scopes.top() ) << endl; 
    }
    scopes.pop(); //erase the bracket  
    while( ( scopes.size() > 0 ) &&
           ( scopes.top() == sSemi )){ //close pending semi's scopes
      decOpenBracket();
      vars->closeBracket();  
      scopes.pop();
    }
  }
    
  openParentheses=0; //we reset openParentheses also 
  insideFuncParam=false; //we are not in a function parameter declaration anymore
  
#ifdef _DEBUG_
  cerr<<"closed brakket and removed local vars\n";
#endif

  if(openBracket==0){ //we are not in a function anymore, so currentFunc set to empty
    currentFunc=Var("","",tOther);
    insideFuncParam=false;
    insideFunction=false;
  }
  
  if( ( scopes.size()>0 ) && (scopes.top()==sFunction) ){
    scopes.pop();
    insideFunction = false;
  } 
}


/*-----------------------------------------------------------------------------
name        :semiColon
description :When the closeScopesOnSemi counter is larger than 0 then count
             times a scope is closed. insideFuncParam is set to false.
parameters  :/
return      :/
exceptions  :vars->closeBracket might throw an exception Error but this will
             be caught in precompile.cpp main() so that lexing stops
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::semiColon(){
  if( !insideForParam ) {
    while( ( scopes.size() > 0 ) &&
           ( scopes.top() == sSemi )){ //close pending semi scopes 
      decOpenBracket();
      vars->closeBracket();
      scopes.pop();
    }
  }
  
  insideFuncParam=false;

  while( (scopes.size() > 0) && (scopes.top()==sFunction)){
    scopes.pop();
  }  
}



/*-----------------------------------------------------------------------------
name        :getStrCurFunction
description :return name of current function
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string ParserConf::getStrCurFunction(){
  return currentFunc.getName();
}



/*-----------------------------------------------------------------------------
name        :getTypCurFunction
description :return type of current function
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
types ParserConf::getTypCurFunction(){
  return currentFunc.getType();
}

/*-----------------------------------------------------------------------------
name        :getStrTypCurFunction
description :return matcher typestring of current function
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string ParserConf::getStrTypCurFunction(){
  return currentFunc.getMatchType();
}


/*-----------------------------------------------------------------------------
name        :getStrTypCurFunction
description :return current function as a Var containing all info
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
Var ParserConf::getCurrentFunction(){
  return currentFunc;
}


  
/*-----------------------------------------------------------------------------
name        :setSkipFor
description :set the global boolean bSkipFor. 
             0=convert int in for loops to BigInt
             1=don't convert int in for loops
parameters  :bool b
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::setSkipFor(bool b){
  bSkipFor=b;
}



/*-----------------------------------------------------------------------------
name        :loadConfigFile
description :load a config file which contains var's and functions that should
             be skipped.
parameters  :char* filename
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::loadConfigFile(char* filename){
  skipVar->loadFile(filename);
}



/*-----------------------------------------------------------------------------
name        :loadConvertFile
description :load a convert file with xml parser...
parameters  :char* filename
return      :/
exceptions  :Error (if file does not exist)
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::loadConvertFile(char* filename){
  convertTree=0;  
  ifstream from(filename);
  if(from){
    XMLParser xml;
    xml.setIncludeWhites(false); //skip white space textnodes 
    xml.parse(&from);
    convertTree=xml.tree;
  }
  else{
    throw Error("Error: ParserConf::loadConvertFile: Conversion file '"+string(filename)+"' does not exist.");
  }
}




  
/*-----------------------------------------------------------------------------
name        :updateLineNr
description :in some states a catch all symbols '.' or patterns which may
             include CR's. We use this function here to increase the
             global int line_nr whenever a CR is found in yytext
parameters  :const string&
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void ParserConf::updateLineNr(const string& s){
  for( unsigned int i=0; i<s.length(); i++ ){
    if(s[i]=='\n') line_nr++;
  }
}
  
  


/*-----------------------------------------------------------------------------
name        :skipQuoted
description :skip quoted strings
parameters  :const string&
return      :/
exceptions  :/
algorithm   :skip anything between " or ' also check for escaped \" 
-----------------------------------------------------------------------------*/
void ParserConf::skipQuoted(const string& s, int& i){
  if(s[i]=='"'){
    i++;
    bool escaped=0;
    while( ( i < (int)s.length() ) && ( (s[i]!='"') || escaped ) ){
      if(s[i]=='\\'){
        escaped=1;
      }
      else{
        escaped = 0;
      }
      i++;
    }
    i++;
  }
  else if(s[i]=='\''){
    i++;
    bool escaped=0;
    while( ( i < (int)s.length() ) && ( ( s[i]!='\'')||escaped ) ){
      if(s[i]=='\\'){
        escaped = 1;
      }
      else{
        escaped = 0;
      }
      i++;
    }
    i++;
  }
}

/*-----------------------------------------------------------------------------
name        :updateLocals
description :update the counts for brackets,braces 
parameters  :const string&
return      :/
exceptions  :/
algorithm   :walk through given string and skip quoted strings and something
             between '
             look for Parentheses in string and call incParentheses or decParentheses 
             accordingly
-----------------------------------------------------------------------------*/
void ParserConf::updateLocals(const string& s){
  int i=0;
  while( i < (int)s.length() ){
    //skipQuoted(s,i); //Parser already checks this
    if( i < (int)s.length() ){
      if(s[i]=='(') incOpenParentheses();
      if(s[i]==')') decOpenParentheses();
      if(s[i]=='{') newScope();
      if(s[i]=='}') closeScope();
      if(s[i]==';') insideFuncParam=0;
    }
    i++;
  }
  
  
}
  


/*-----------------------------------------------------------------------------
name        :decOpenParentheses
description :decrease the openParentheses var
             if it is zero or less we know we can set insideForLoop to false
parameters  :const string&
return      :/
exceptions  :/
algorithm   :/
-----------------------------------------------------------------------------*/
void ParserConf::decOpenParentheses(){
  openParentheses--;
  if (openParentheses<=0){
    openParentheses=0;
    insideForParam=false;
  }
  #ifdef _DEBUG_ 
    cerr<<"openbraces="<<openParentheses<<" at: "<<line_nr<<endl;
  #endif
}


/*-----------------------------------------------------------------------------
name        :incOpenParentheses
description :increase the openParentheses var
parameters  :const string&
return      :/
exceptions  :/
algorithm   :/
-----------------------------------------------------------------------------*/
void ParserConf::incOpenParentheses(){
  openParentheses++;
  #ifdef _DEBUG_ 
    cerr<<"openbraces="<<openParentheses<<" at: "<<line_nr<<endl;
  #endif
}


/*-----------------------------------------------------------------------------
name        :decOpenBracket
description :decrease the openBracket var
parameters  :const string&
return      :/
exceptions  :/
algorithm   :/
-----------------------------------------------------------------------------*/
void ParserConf::decOpenBracket(){
  openBracket--;
  if (openBracket<0){
    openBracket=0;
  }
  #ifdef _DEBUG_ 
    cerr<<"openbracket="<<openBracket<<" at: "<<line_nr<<endl;
  #endif
}



/*-----------------------------------------------------------------------------
name        :incOpenBracket
description :increase the openBracket var
parameters  :const string&
return      :/
exceptions  :/
algorithm   :/
-----------------------------------------------------------------------------*/
void ParserConf::incOpenBracket(){
  openBracket++;
  #ifdef _DEBUG_ 
    cerr<<"openbracket="<<openBracket<<" at: "<<line_nr<<endl;
  #endif
}


