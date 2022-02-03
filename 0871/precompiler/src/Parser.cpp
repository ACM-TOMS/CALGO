/*=============================================================================
author        :Walter Schreppers
filename      :Parser.cpp
created       :22/05/2000
modified      :25/02/2002
version       :3
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#include "Parser.h"



/*-----------------------------------------------------------------------------
name        :Parser
description :a default constructor this is only here so that
             derived class PreParser can use it
parameters  :ParserConf*
return      :Parser
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
Parser::Parser(){
  //nothing to be done...
}




/*-----------------------------------------------------------------------------
name        :Parser
description :constructor
parameters  :ParserConf*
return      :Parser
exceptions  :/
algorithm   :initialize the locals with init function
-----------------------------------------------------------------------------*/
Parser::Parser(ParserConf *parserConf){
  init(parserConf);
}



/*-----------------------------------------------------------------------------
name        :init
description :constructor
parameters  :ParserConf*
return      :/
exceptions  :/
algorithm   :initialize the local vars needed by parser 
             from ParserConf pointer, also set needed
             precision and booleans needed for parsing
-----------------------------------------------------------------------------*/
void Parser::init(ParserConf *parserConf){
  if(parserConf!=NULL){
    skipVar=parserConf->getSkipVar();
    pConf=parserConf;
    if((pConf->convertTree!=0)&&(pConf->convertTree!=NULL)){
      varMatcher=new Matcher(pConf->convertTree);
    }
    else{
      varMatcher=new Matcher();
    }
  }
  
  bNoDefault=0;  //if true no precision settings will be set after the main
  bInit=0;       //if true init file is opened and contents copied after main
  bNoPrintf=0;   //if true printf statements are not converted
  bPrintUsed=0;  //set to true when printf instruction is converted
  bMain=0;       //true when main function is being declared
  fRadix=2;
  fPrecision=24;
  fExprangeLow=-126;
  fExprangeUp=127;
  fRounding=tRoundNearest;
  fIoModeStr="";

  initFileName ="";   //holds name of initFile
  strOriginalType=""; //holds the original type when setVarType was called
  
  bArrayConstruct=false; //used for Array constructors

  ArrayType=tOther;
  pQueue.clear(); //the parser queue
  
  constructs.clear(); //clear the constructs variable vector
}




/*-----------------------------------------------------------------------------
name        :Parser
description :destructor
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
Parser::~Parser(){
  if( varMatcher!=NULL ) delete varMatcher;
  if( pConf!=NULL ) delete pConf;
}


/*-----------------------------------------------------------------------------
name        :writeIncludes
description :write some #include directives
             write iostream and iomanip if -noprintf was NOT used
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :traverse xml tree and write contents of include tags
-----------------------------------------------------------------------------*/
void Parser::writeIncludes(ofstream& out){
  if(!bNoPrintf){        //    && (bPrintUsed) ){ no good, it is set after this function is called
    out<<"#include <iostream>"<<endl;
    out<<"#include <iomanip>"<<endl;
    out<<"using namespace std;"<<endl;
  }
  
  XMLNode* tree=pConf->convertTree;
  if( (tree==0) || (tree==NULL) ) return; //empty tree
  
  //write given includes in xml conversion file
  bool bFirstInclude = true;

  for(XMLNode::iterator i=tree->begin();i!=tree->end(); ++i){
    if( (*i)->nodeType() == XMLNode::ElementNode ){
      ElementNode* e=(ElementNode*) *i; //it was an element so cast it
      if( e->tagName() == "include" ) {
        if(e->size()>0){
          for(XMLNode::iterator j=e->begin();j!=e->end(); ++j){
            
            if( (*j)->nodeType() == XMLNode::TextNode){
            
              if( bFirstInclude ){
                out << endl;
                bFirstInclude = false;
              }
              
              TextNode* t=(TextNode*) *j;
              if( e->attribute("global","") == "true" ){
                out<<"#include <"<<t->data()<<">"<<endl;
              }
              else{
                out<<"#include \""<<t->data()<<"\""<<endl;
              }
              
            }
            
          }
        }//has children
        else{
          cerr<<"Warning::include tag has no data"<<endl;
        }
      }//is include element
    }//is element
  }//for
  
  if( !bFirstInclude ) out << endl;
}


/*-----------------------------------------------------------------------------
name        :writeInit
description :copy code in initfile onto ofstream
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
void Parser::writeInit(ofstream& out){
  if( bInit ){
    ifstream in( initFileName );
    if(in){
      char c='\n';
      while( !in.eof() ){
        out << c;
        c=in.get();
      }
    }
    else{
      throw Error("Error: Parser::writeInit: init file '"+string(initFileName)+"' does not exist.");
    }
  }
}



/*-----------------------------------------------------------------------------
name        :writePrecision
description :copy code in <init> tag
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :traverse xml tree and copy pcdata in precision tag to ofstream& out
-----------------------------------------------------------------------------*/
void Parser::writePrecision(ofstream& out){
  XMLNode* tree=pConf->convertTree;
  if( (tree==0) || (tree==NULL) ) return; //empty tree
  
  for(XMLNode::iterator i=tree->begin();i!=tree->end(); ++i){
    if( (*i)->nodeType() == XMLNode::ElementNode ){
      ElementNode* e=(ElementNode*) *i; //it was an element so cast it
      if( e->tagName() == "init" ) {
        if(e->size()>0){
        
          out<<endl<<endl;
          out<<"    "; //indentation
          for(XMLNode::iterator j=e->begin();j!=e->end(); ++j){
            if( (*j)->nodeType() == XMLNode::TextNode){
              TextNode* t=(TextNode*) *j;
              out<<t->data()<<endl;
            }
            else{
              cerr<<"Warning: Parser::writePrecision : <init> tag contains invalid PCDATA"<<endl;
            }
          }
          
        }//has children
        else{
          cerr<<"Warning::<init> tag has no data"<<endl;
        }
      }//is include element
    }//is element
  }//for
}

/*-----------------------------------------------------------------------------
name        :setFloatPrecision
description :insert radix,precision, exponent range and rounding settings
             into ofstream& out.
             ADDED 1/7/2000 don't insert settings if bNoDefault=1
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :use private local vars to set right values
-----------------------------------------------------------------------------*/
void Parser::setFloatPrecision(ofstream& out){
  if(!bNoDefault){
    #ifdef _DEBUG_
  	  cerr<<"PARSER::Inserting precision at: "<<pConf->getLine()<<endl;
  	#endif
    out<<endl<<endl;
    out<<"  MpIeee::fpEnv.setRadix("<< fRadix <<");"<<endl;
    out<<"  MpIeee::fpEnv.setPrecision("<< fPrecision <<");"<<endl;
    out<<"  MpIeee::fpEnv.setExponentRange("<< fExprangeLow <<","<< fExprangeUp <<");"<<endl;
    switch(fRounding){
      case tRoundNearest:
        out<<"  MpIeee::fpEnv.setRound(FP_RN);"<<endl;
        break;
      case tRoundPlus:
        out<<"  MpIeee::fpEnv.setRound(FP_RP);"<<endl;
        break;
      case tRoundMin:
        out<<"  MpIeee::fpEnv.setRound(FP_RM);"<<endl;
        break;
      case tRoundZero:
        out<<"  MpIeee::fpEnv.setRound(FP_RZ);"<<endl;
        break;
    }
  
    if(fIoModeStr.size()>0){
      out<< "  ArithmosIO::setIoMode("<<fIoModeStr<<");"<<endl;
    }
    else{
      out<< "  ArithmosIO::setIoMode(ARITHMOS_IO_DEFAULT);"<<endl;
    }
  }
}



/*-----------------------------------------------------------------------------
name        :getStrType
description :return a string with the with the correct type
             meaning Rational,MpIeee or BigInt
parameters  :types t
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Parser::getStrType(const string& sourceName){
  ConvertElem c=varMatcher->convertLookup(sourceName);
  
  if( c.isError() ){
    //no type found
    return "";
  }
  
  if( c.target.isKeyword() ){
    if(c.hasOperation){
      //cout<<"we need to call ConvertConfig here!"<<endl;
      return strOriginalType;
    }
    else{
      return c.target.keyword;
    }
  }
  else{
    cerr<<"Warning: Conversion rule "<<c.name<<" has a target which is not a keyword!"<<endl;
    return "";
  }
}





/*-----------------------------------------------------------------------------
name        :getStrType
description :return a string with the with the correct type
             meaning Rational,MpIeee or BigInt
parameters  :Variables::iterator
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Parser::getStrType(Variables::iterator p){
  return getStrType( p->getMatchType() );
}



/*-----------------------------------------------------------------------------
name        :getOriginalType
description :when conversion of a variable must be skipped we use this function to get
             the original type of the variable
parameters  :/
return      :string
exceptions  :/
algorithm   :return the private local strOriginalType this is set whenever
             setVarType was called
-----------------------------------------------------------------------------*/
string Parser::getOriginalType(){
  return strOriginalType;
}



/*-----------------------------------------------------------------------------
name        :setOriginalType
description :move pos elements from pQueue into strOriginalType 
             this is used by setVarType
parameters  :unsigned int pos
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setOriginalType(unsigned int pos){
  strOriginalType="";
  for(unsigned int i=0;i<pos;i++) strOriginalType+=pQueue.poptext();
}


/*-----------------------------------------------------------------------------
name        :parseVar
description :output the BigInt ,MpIeee or Rational in proper format 
parameters  :ofstream& out, outputType t
return      :/
exceptions  :/
algorithm   :- get the type string
             - iterate through constructs = var list
             - if skipVar->locate returns true this var is not converted
             - use Convert class to convert assignments
             - if skipVar->locate returns true this var is not converted
             
             little note: the extra space inserted into out between the type and spec
             is necessaray ex: float a,b -> MpIeee a; MpIeeeb if space is not 
             inserted !
-----------------------------------------------------------------------------*/
void Parser::parseVar(ofstream& out,string ending){
  for(Variables::iterator p=constructs.begin();p!=constructs.end();++p){
    if(
        ( !skipVar->locate( pConf->getStrCurFunction(), p->getName() ) )
      ){
      
      out << getStrType(p)<<" "
          << p->getSpec()
          << p->getName()
          << p->getArrayArg();

      if(p->getValue()!=""){
        Convert varval(*p,varMatcher,pConf);
        out << "= " << varval.getValue();
      }
    }
    else{ //output original value (if it has one)
      out << getOriginalType() << " " 
          << p->getSpec()
          << p->getName() 
          << p->getArrayArg();
      
      if(p->getValue()!=""){
        out << "= " <<p->getValue();
      }
    }
    out << ending;
  }
}



/*-----------------------------------------------------------------------------
name        :parseAssignment
description :change format of an assignment string to conform 
             to MpIeee or BigInt
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :- get the first (and only) element in constructs
             - use Convert class to create proper output
-----------------------------------------------------------------------------*/
void Parser::parseAssignment(ofstream& out){
  Variables::iterator p=constructs.begin(); // there is only one assignment per call of this function
  if( !skipVar->locate( pConf->getStrCurFunction(), p->getName() ) ) {
    Convert val(*p,varMatcher,pConf);

    #ifdef _DEBUG_
      string newval = val.getValue();
      out  << newval;
      cout << "NAME ="<<p->getName() << ", TYPE =" << getStrType(p) << ", OLD value=" << p->getValue() << endl;
      cout << "NAME ="<<p->getName() << ", TYPE =" << getStrType(p) << ", NEW value=" << newval << endl;
    #else
      out << val.getValue();
    #endif

  }
  else{
    out<<p->getValue();
  }
}





/*-----------------------------------------------------------------------------
name        :setNoDefault
description :set the bool bNoDefault
             1=don't use setFloatPrecision
             0=do use it to set the settings after main(){}
             ADDED 1/7/2000
parameters  :bool b
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setNoDefault(bool b){
  bNoDefault=b;
}


/*-----------------------------------------------------------------------------
name        :setInit
description :1= use initfile and copy contents to outstream after main funciton
             0= skip copying of initfile
parameters  :bool b, char* filename
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setInit(bool b, char* filename){
  bInit=b;
  initFileName=filename;
}


/*-----------------------------------------------------------------------------
name        :setNoPrintf
description :set the bool bNoPrintf
             if true the parser will not convert printf statements
             ADDED 24/6/2001
parameters  :bool b
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setNoPrintf(bool b){
  bNoPrintf=b;
}




/*-----------------------------------------------------------------------------
name        :setRadix
description :set the local variable fRadix
parameters  :int
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setRadix(int i){
  fRadix=i;
}


/*-----------------------------------------------------------------------------
name        :setPrecision
description :set the local variable fPrecision
parameters  :int
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setPrecision(int i){
  fPrecision=i;
}

/*-----------------------------------------------------------------------------
name        :setExprange
description :set the local variable fExprangeLow and fExprangeUp
parameters  :int low, int up
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setExprange(long low,long up){
  fExprangeLow=low;
  fExprangeUp=up;
}


/*-----------------------------------------------------------------------------
name        :setRounding
description :set the local variable fRounding
parameters  :roundType
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setRounding(roundType r){
  fRounding=r;
}


/*-----------------------------------------------------------------------------
name        :setIoModeStr
description :set the fIoModeStr
parameters  :string s
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::setIoModeStr(string s){
  fIoModeStr=s;
}







/*-----------------------------------------------------------------------------
name        :updateLocals
description :Keep track of open_pars, etc 
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::updateLocals(){
  if(pQueue.size()>0){
    switch(pQueue.getToken(0)){
      case open_par:    pConf->incOpenParentheses();   break;
      case close_par:   pConf->decOpenParentheses();   break;
      case open_brace:  pConf->newScope();        break;
      case close_brace: pConf->closeScope();      break;
      case semi_colon:  pConf->semiColon();       break;
      default:          break;
    }
  }
}




/*-----------------------------------------------------------------------------
name        :moveToken
description :Keep track of open_pars, etc by calling updateLocals
             remove a token from the pQueue and write it's string onto the
             ofstream out
             
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::moveToken(ofstream& out){
  updateLocals();
  if(pQueue.size()>0) out<<pQueue.poptext()<<flush;
}



/*-----------------------------------------------------------------------------
name        :moveTokens
description :write pos number of tokens onto ofstream out and remove them from
             pQueue
parameters  :unsigned int& pos, ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::moveTokens(unsigned int& pos,ofstream& out){
  while((pos>0)&&(pQueue.size()>0)){
    moveToken(out);
    pos--;
  }
}



/*-----------------------------------------------------------------------------
name        :copyToken
description :write string of first token in pQueue onto ofstream out
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::copyToken(ofstream& out){
  out<<pQueue.getText(0);
}


/*-----------------------------------------------------------------------------
name        :dumpq
description :dump the pQueue literally to out
parameters  :ofstream& out 
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::dumpq(ofstream& out){
  while(pQueue.size()>0) moveToken(out);
  pQueue.clear();
}

/*-----------------------------------------------------------------------------
name        :writeStrType
description :write MpIeee,BigInt,Rational to out corresponding to
             t
parameters  :ofstream& out ,types t
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::writeStrType(ofstream& out,const string& name){
  out<<getStrType(name);
}

/*-----------------------------------------------------------------------------
name        :writeOriginalType
description :write the original variable type
parameters  :ofstream& out 
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::writeOriginalType(ofstream& out){
  out<<getOriginalType();
}



/*-----------------------------------------------------------------------------
name        :checkStatementEnd
description :return true if the pQueue contains a statement ending
             which is open_brace, semi_colon or close_brace
             Added check for bArrayConstruct to keep adding tokens if
             we found the closing scope, for pending new variable declarations.
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool Parser::checkStatementEnd(){
  int tok=pQueue[pQueue.size()-1].token;
  if( (tok == open_brace) || 
      (tok == semi_colon) ||
      ( (tok == close_brace) && !bArrayConstruct )      
    ) return true;
  else return false;
}





/*-----------------------------------------------------------------------------
name        :setVarType
description :given the pQueue decide what type of variable is being declared
             The return types can be tInt,tFloat,tDouble meaning 
             one which has to be converted or 
             return tOther if it's a type not to be converted 
             or something else which has to be handled elsewhere
parameters  :ofstream&
return      :types
exceptions  :/
algorithm   :
             we use the varMatcher to find a matching source tag.
             When we have found the matching source tag we look if
             there is a conversion rule for this type, if so we return
             tKnown in the other cases we return tOther.
-----------------------------------------------------------------------------*/
types Parser::setVarType(ofstream& out){
  #ifdef _DEBUG_
    cerr<<"we are in setVarType at "<<pConf->getLine()<<" ";
    showQueue();
  #endif
  skipToWord(out);
  unsigned int pos=0;
  
  string testType=pQueue.getText(pos);
  string vType="";
  MatchElem m=varMatcher->sourceLookup( testType );
  MatchElem match;

  while(m.mtype!=mError){
    vType=testType;
    match=m;
    
    pos++;
    if(pQueue.getToken(pos)==white_space) pos++;
    testType+=" "+pQueue.getText(pos);
    m=varMatcher->sourceLookup(testType);
  }

  if(vType!="") { //we found a good type
    ConvertElem cElem=varMatcher->convertLookup(match.name);
    if( cElem.mtype!=mError ){ //we found a suitable conversion rule
      typeMatch=match;
      setOriginalType(pos);
      return tKnown;
    }
  }
  //no suitable conversion rule or not a defined type.
  for( unsigned int i=0; i<pos; i++ ) moveToken( out ); 
  return tOther;
}


/*-----------------------------------------------------------------------------
name        :skipWhite
description :if token on pQueue at position pos is white_space increase it
parameters  :unsigned int& pos
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::skipWhite(unsigned int& pos){
  if(pQueue.getToken(pos)==white_space) pos++;
}

/*-----------------------------------------------------------------------------
name        :getVarSpec
description :return the specification string of a variable declaration 
             this can be a reference,pointer or simply a white space
             9/2/2005 added support for const pointer declaration : ex. double * const a;
parameters  :unsigned int& pos
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Parser::getVarSpec(unsigned int& pos){
  string varspec="";
  while( ( pos < pQueue.size() ) && 
         ( 
            ( pQueue.getToken(pos) == mult        )||
            ( pQueue.getToken(pos) == sym_and     )||
            ( pQueue.getToken(pos) == white_space )||
            ( ( pQueue.getToken(pos) == sym_word ) && ( pQueue.getText(pos) == "const" ) )
         )
       ){
    varspec+=pQueue.getText(pos);
    pos++;
  }
  return varspec;
}


/*-----------------------------------------------------------------------------
name        :getVarName
description :return the name of a variable or function,
             the name may be preceded by a namespace:: or a classname::
parameters  :unsigned int& pos
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Parser::getVarName( unsigned int& pos ){
  string varname="";  
  if(pQueue.getToken(pos)==sym_word){
    varname=pQueue.getText(pos);
    pos++;
    skipWhite(pos);
    //get possible classname::name
    if(pQueue.getToken(pos)==double_colon){
      int start=pos;
      pos++;
      skipWhite(pos);
      if(pQueue.getToken(pos)==tilde) pos++;
      if(pQueue.getToken(pos)==sym_word){
        pos++;
        for( unsigned int i=start; i<pos; i++ ) varname+=pQueue.getText(i);
      }
      else{
        //this is not right
        cerr<<"Parser::possible wrong convertion at: "<<pConf->getLine()<<endl;
        varname="";
      }
    }
  }
  return varname;
}




/*-----------------------------------------------------------------------------
name        :handleOperatorOverloading
description :special case of function declaration
parameters  :unsigned int& pos,const string& funcname
return      :string
exceptions  :/
algorithm   :not used now still have to discuss this!!! leads to problems
             such as int operator=(3+5); 
             we can regard this as operator overloading if it's in a class definition
             but also just as an int with name operator ...
             
             see also handleFunction
-----------------------------------------------------------------------------*/
void Parser::handleOperatorOverloading( unsigned int& pos,const string& funcname ){
  if(funcname=="operator"){
    unsigned int temp=pos;
    skipWhite( temp );
    if((pQueue.getToken(temp)==assignment)||
        (pQueue.getToken(temp)==comparison)||
        (pQueue.getToken(temp)==bool_op)){
      temp++;
      skipWhite( temp );
      if( ( pQueue.getToken( temp ) == open_par ) &&
          ( pConf->getOpenBracket()==0) )
      {
        for(unsigned int i=pos;i<temp;i++) constructs.expandLastName( pQueue.getText(i) );
        pos=temp; //set pos to this brace
      }
    }
  }
}



/*-----------------------------------------------------------------------------
name        :handleFunction
description :handle new functions
parameters  :unsigned int& pos,const string& funcname,types functype
return      :string
exceptions  :/
algorithm   :we just read a variable name if the next token is an open brace '('
             and we are not inside a scope (meaning scope depth is zero or
             getOpenBracket()==0)
             we know it's a new function being declared.
             we add it, by calling pConf->newFunction
             we also check if the name is main in this case we set the
             local bMain to true
-----------------------------------------------------------------------------*/
bool Parser::handleFunction(unsigned int& pos,const string& funcname,types functype){
  //handleOperatorOverloading(pos,funcname);
  
  //for function pointers inside or outside structs this is needed!
  if(!pConf->getInsideFuncParam()){
    skipWhite(pos);
    int braces=pConf->getOpenParentheses();
    while(( pQueue.getToken(pos) == close_par )&&( braces!=0 )){
      pos++;
      skipWhite(pos);
      braces--;
    }
  }
  
  if(
     (pQueue.getToken(pos)==open_par) &&
     //(pConf->getOpenBracket()==0)&& //only good check when not in class or struct
     (!pConf->getInsideFuncParam())
    ){

    if(funcname=="main"){
      bMain=1;
      functype=tOther;
    }
    
    string origType=getOriginalType();
    string targType=getStrType( typeMatch.name );
    pConf->newFunction(funcname,functype,typeMatch.name,origType,targType);
    return 1;
  }
  else{
    return 0;
  }
}




/*-----------------------------------------------------------------------------
name        :getEvenParentheses
description :skip all tokens in pQueue until the braces are even again
             the pos after the this procedure is on the last brace or
             if not found it will be on last token
parameters  :unsigned int& pos
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
string Parser::getEvenParentheses(unsigned int& pos){
  string retStr=pQueue.getText( pos );
  pos++;
  int braces=1;
  while( ( pos < pQueue.size() ) && ( braces != 0 ) ){
    if(pQueue.getToken(pos)==open_par){
      braces++;
    }
    if(pQueue.getToken(pos)==close_par){
      braces--;
    }
    retStr+=pQueue.getText(pos);
    pos++;
  }
  pos--;
  return retStr;
}




/*-----------------------------------------------------------------------------
name        :handleAssignment
description :if there is an assignment to a variable add this to the
             value string in the constructs.
parameters  :unsigned int& pos,const types vartype
return      :void
exceptions  :/
algorithm   :expand value string of last variable in constructs until
             we find a semi-colon,comma or close brace (skipping even braces)
             when finding a '{' we set bArrayConstruct to true and ArrayType
             to the current type so that we may get the rest of this assignment
             in the next fill of the pQueue (since the pQueue stopped filling
             at the open scope )
-----------------------------------------------------------------------------*/
void Parser::handleAssignment(unsigned int& pos,const types vartype){
  #ifdef _DEBUG_
    cerr<<"handleAssignment at pos "<<pos<<" : "; showQueue(); cerr<<endl;
  #endif
  int prevToken = pQueue.getToken( pos );
  if(pQueue.getToken(pos)==assignment){
    //here we could find something like sin( ... which could be converted ...
    pos++;
    while(
          ( pos < pQueue.size() )&&
          ( pQueue.getToken(pos) != comma )&&
          ( pQueue.getToken(pos) != semi_colon  )&&
          ( pQueue.getToken(pos) != close_par )
         ){
      if( pQueue.getToken(pos) == open_par ){
        /*cout<<"open brace in assignent"<<endl;
        //special case to not convert function arguments in assignment to vars of known type
        if( ( vartype != tOther ) && ( prevToken == sym_word ) ){
          cout<<"found function_args at line "<<pConf->getLine()<<endl;
          int startpos = pos;
          constructs.expandLastValue( getEvenParentheses(pos) ); //value string
          for( unsigned int i=startpos; i<=pos; i++ ){
            ParseElem p = pQueue.get( i );
            p.token = function_argument;
            constructs.expandLastValue( p );
          }
        }
        */
        if( vartype != tOther ){
          int startpos=pos;
          constructs.expandLastValue( getEvenParentheses(pos) ); //value string
          for(unsigned int i=startpos; i<=pos; i++) { //value queue
            constructs.expandLastValue( pQueue.get(i) );
          }
        }
        else{ //no need to add for tOther vars
          string dummy=getEvenParentheses(pos);
        }
      }
      //added this for array constructors! { 1,2,3,.... }
      else if( pQueue.getToken(pos)==open_brace ){
        bArrayConstruct=1;
        ArrayType=vartype;
        if(vartype!=tOther){
          constructs.expandLastValue(pQueue.getText(pos)); //value string
          constructs.expandLastValue(pQueue.get(pos));      //value queue
        }
      }
      else{
        if(vartype!=tOther) { 
          constructs.expandLastValue( pQueue.getText(pos) ); //value string
          constructs.expandLastValue( pQueue.get(pos) ); //value queue
        }
      }
      prevToken = pQueue.getToken( pos );
      pos++;
    }//while
  }//if
  skipWhite(pos);
  #ifdef _DEBUG_
    cerr<<"after handleAssignment at pos "<<pos<<endl;
  #endif
}


/*-----------------------------------------------------------------------------
name        :dumpHandled
description :dump the part of pQueue that has been handled by addVar
             for tOther types we copy it to out literally
parameters  :unsigned int& pos,ofstream& out, const types vartype
return      :string
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::dumpHandled(unsigned int& pos,ofstream& out,const types vartype){
  if(vartype==tOther){
    moveTokens(pos,out);
  }
  else{
    //pop all but the last one, this will probably be a ; or (
    while((pos>0)&&(pQueue.size()>0)){
      updateLocals();
      pQueue.popfirst();
      pos--;
    }
  }
}

/*-----------------------------------------------------------------------------
name        :handleArrayArg
description :copy array args literally (so everything between [ .. ] )
             put then in arrayarg of constructs if it's for a type to be
             converted
parameters  :unsigned int& pos, const types vartype
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::handleArrayArg(unsigned int& pos,const types vartype){
  int start=pos;
  while( ( pQueue.getToken(pos) == open_bracket ) &&
         ( pos < pQueue.size() )
       ){ //skip the array indexing!
    pos++;
    int opened=1;
    while( (opened>0) && (pos < pQueue.size()) ){
      if(pQueue.getToken(pos)==open_bracket) opened++;
      if(pQueue.getToken(pos)==close_bracket) opened--;
      pos++;
    }
  }
  skipWhite(pos);
  string arrayarg="";
  
  for( unsigned int i=start; i<pos; i++ ){
    arrayarg+=pQueue.getText(i);
  }
  if(vartype!=tOther){
    constructs.expandLastArrayArg(arrayarg);
  }
  
}


/*-----------------------------------------------------------------------------
name        :handleArrayAssign
description :this will copy the pQueue entirely when bArrayConstruct was
             set to true in handleAssignment
parameters  :types vartype
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::handleArrayAssign(const types vartype){
  if(bArrayConstruct){
    if(vartype!=tOther){
      while((pQueue.size()>0)&&(pQueue.getToken(0)!=close_brace)){
        constructs.expandLastValue( pQueue.getText(0) ); //value string
        constructs.expandLastValue( pQueue.pop() );      //value queue
      }
      if( pQueue.getToken(0)==close_brace ){
        constructs.expandLastValue( pQueue.getText(0) ); //value string '}'
        constructs.expandLastValue( pQueue.pop() );      //value queue
        if( pQueue.getToken(0)==comma ){
          pQueue.popfirst(); //lose the comma, we will split up the constructs with ';'
        }
      }
    }
    bArrayConstruct=false;
  }
}


/*-----------------------------------------------------------------------------
name        :addVar
description :Add newly constructed variables to constructs
             and to pConf->getVars()
             return 0 if more variables are still on pQueue
parameters  :ofstream& out, currentVarType
return      :bool
exceptions  :/
algorithm   :

-----------------------------------------------------------------------------*/
bool Parser::addVar(ofstream& out,const types currentVarType){
  handleArrayAssign(currentVarType);
  
  if(bMain) return 1;
  
  unsigned int pos=0;

  string currentVarSpec=getVarSpec(pos);
  string varname=getVarName(pos);
  if(varname==""){
    if(currentVarType==tOther) moveTokens(pos,out);
    return 1; //no var to be added
  }

  string origType=getOriginalType();
  string targType=getStrType( typeMatch.name );
  
  if( currentVarType != tOther ){
    Var v( currentVarType, 
           currentVarSpec, 
           varname, 
           typeMatch.name,
           origType,
           targType 
         );
#ifdef _DEBUG_
  cout << "adding var to constructs and parser conf: " << endl << v << endl;
#endif
    constructs.insert( v );
    pConf->getVars()->insert( v );
  }
  else{ //currentVarType == tOther
    pConf->getVars()->insert( currentVarType, currentVarSpec, varname, "", "", "" );  
  }
  
  skipWhite(pos);
  
  handleArrayArg(pos,currentVarType);
  
  //check for a function
  if(handleFunction(pos,varname,currentVarType)){
    if(bMain) return 1; //queue will be literally dumped to out
    dumpHandled(pos,out,currentVarType);
    return 1;//done adding
  }
  
  handleAssignment(pos,currentVarType);

  if(pQueue.getToken(pos)==comma){
    if(pConf->getInsideFuncParam()){
      dumpHandled(pos,out,currentVarType); //will remove everything but the comma for known types from pQueue
      return 1; //done with this type
    }
    else{
      pos++;
      dumpHandled(pos,out,currentVarType); //for known types it will also dump the comma (it will become ';' for known types)
      return 0; //more var's may follow
    }
  }
  
  dumpHandled(pos,out,currentVarType);
  return 1; //we are done with adding variable's of this type
}




/*-----------------------------------------------------------------------------
name        :skipToWord
description :copy everything literally to out until we find a sym_word
             and update the locals of ParserConf
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::skipToWord(ofstream& out){
  while((pQueue.size()>0)&&(pQueue.getToken(0)!=sym_word)){
    moveToken(out);
  }
}





/*-----------------------------------------------------------------------------
name        :convertForLoop
description :convertion of types used by for loop
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :-check if skipfor option was used if so just copy the for loop
              statement literally;
             -else handle for loop same as function meaning the variables
              declared belong to the scope of the for loop.
              
              
              CHANGED  27/06/2001 !!!!!!!
              we do not check getSkipFor here but do it in the 
              convertStatement function !!!!! 
              This is better since the variable will then also be added as a
              tOther type ...
-----------------------------------------------------------------------------*/
void Parser::convertForLoop(ofstream& out){
  moveToken(out); //copy "for" to out

    //handle for loops by setting the right params etc...
    if(pQueue.getToken(0)==white_space) moveToken(out);
    if(pQueue.getToken(0)==open_par){
      //okay we definately have a forloop
      moveToken(out);
      pConf->newForLoop();
      #ifdef _DEBUG_
        showQueue();
      #endif
    }
    else{
      cerr<<"Parser::possible wrong conversion of FOR loop at: "<<pConf->getLine();
      pConf->newForLoop();
    }
}



/*-----------------------------------------------------------------------------
name        :convertReturn
description :convertion of return statements
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :convert the expression after return
             if it's not a function of tOther type
             and skipVar returns false for this function
             and if we are not in the main function!
-----------------------------------------------------------------------------*/
void Parser::convertReturn(ofstream& out){
    if( (pConf->getTypCurFunction()!=tOther) &&
        (!skipVar->locate(pConf->getStrCurFunction(),""))&&
        (!skipVar->locate(pConf->getStrCurFunction(),pConf->getStrCurFunction()))
       ){

      moveToken(out); //write the 'return' to out
      constructs.clear();
      Var retVar=pConf->getCurrentFunction(); //same name and type as function
      retVar.setName("return");
      constructs.insert( retVar );

      while((pQueue.size()>0)&&(pQueue.getToken(0)!=semi_colon)){
        constructs.expandLastValue( pQueue.getText(0) );
        constructs.expandLastValue( pQueue.pop() );
      }
      parseAssignment(out);
    }
    else{//don't convert
      while((pQueue.size()>0)&&(pQueue.getToken(0)!=semi_colon)){
        moveToken(out);
      }
    }
}


/*-----------------------------------------------------------------------------
name        :convertReserved
description :convertion of reserved words
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :if the reserved word is for then call convertForLoop else copy
             the word literally
-----------------------------------------------------------------------------*/
void Parser::convertReserved(ofstream& out){
  #ifdef _DEBUG_
    cerr<<"inside convertReserved at "<<pConf->getLine();
  #endif
  if(pQueue.getText(0)=="for"){
    convertForLoop(out);
  }
  else if(pQueue.getText(0)=="return"){
    convertReturn(out);
  }
  else if(pQueue.getText(0)=="class"){
    if(pQueue.getToken(pQueue.size()-1)==open_brace){
      pConf->decOpenBracket(); //so that functions will be found inside class definitions
    }
  }
  else{
    //just copy literally
    #ifdef _DEBUG_
      cerr<<" copying a word literally from queue:";
      showQueue();
    #endif
    moveToken(out);
  }
}




/*-----------------------------------------------------------------------------
name        :convertPrintf
description :convert printf statements if bNoPrintf is false
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :get the format string and the arguments of the printf statement
             from the pQueue. Use the the class PrintConverter to write out
             an equivalent cout statement.
-----------------------------------------------------------------------------*/
void Parser::convertPrintf(ofstream& out){
  if(bNoPrintf){
    dumpq(out); //copy the printf statement literally
  }
  else{
    bPrintUsed=1;
    //remove the 'printf' etc upto the quoted string
    //or if there is no string upto the first argument
    pQueue.popfirst();
    while((pQueue.size()>0)&&
          (pQueue.getToken(0)!=quoted_string)&&
          (pQueue.getToken(0)!=sym_word)){
      pQueue.popfirst();
    }
  
    //store the format string
    string formatStr="";
    if(pQueue.getToken(0)==quoted_string){
      formatStr=pQueue.poptext();
    }
  

    //store the arguments 
    vector<PrintArg> vArgs;
    string text="";
    pargtype ptype=unknown;

    bool bExpressionArg = false;

    while(pQueue.size()>0){
      #ifdef _DEBUG_
        cout << "PrintArg text = " << text << endl;
      #endif
      if( pQueue.getToken(0) == open_par ){
        text+=pQueue.poptext();
        unsigned int pos = 0;
        text+=getEvenParentheses( pos );
        while( pos-- > 1 ) pQueue.popfirst();
        bExpressionArg = true;
      }
      if((pQueue.getToken(0) == comma)||(pQueue.getToken(0)==close_par)){
        if((ptype!=unknown) && (text!="") ){
        #ifdef _DEBUG_
          cout << "PrintArg added = " << text << " with type " << ptype << endl;
        #endif
          vArgs.push_back( PrintArg( text, ptype ) );
          text="";
          ptype=unknown;
        }
        pQueue.popfirst();
      }
      
      if( bExpressionArg ){
        bExpressionArg = false;
        ptype = pvar;
      }
      else{
        text+=pQueue.getText(0);
      }
      
      switch(pQueue.getToken(0)){
        case sym_word:        ptype=pvar;       break;
        case quoted_string:   ptype=pstring;    break;
        case integer:         ptype=pnumber;    break;
        case decimal:         ptype=pdecnumber; break;
        default: break;
      }
      pQueue.popfirst();
      
    }//while


  #ifdef _DEBUG_
    cerr << "Calling print converter class with formatStr='"<< formatStr <<"'"<<endl;
  #endif

    //write the cout statement to out  
    PrintConverter p(formatStr,vArgs);
    out << "{" << p << "}"; //add {,} because it may contain multiple cout calls (to simulate 1 printf)
  }
}




/*-----------------------------------------------------------------------------
name        :convertKnownVar
description :this is a variable we found that was declared earlier
             if there is an assignment or comparison after it we will convert
             it so it comply's with the new types BigInt,MpIeee or Rational
parameters  :ofstream& out 
return      :/
exceptions  :/
algorithm   : next token must be a comparison or assignment...
              then it has to be converted otherwize left alone
              so name = 45 will become ... somethin like name = BigInt("45")
              
              for now we will use the old style of convertion ...
              
              also added checking for printf, since a devious programmer
              could have done something like int printf=5;
              then later on printf would be recognized as a knownvar but
              it could also be a printf statement! this needs to be checked
              here!
              
              A POSSIBLE TODO:
              Here we could actually completely parse the expression and
              try changing sin and cos etc functions into the MpIeee equivalents.
              
              But this requires a lot of work and is error prone
              since someone might have defined a custom sin function (possibly even in
              another file which was included). This sin function could
              do something totally different than compute the sine of a
              variable. When we then convert it to someting like MpIeee.sin(...
              we would get completely wrong results. To truly get this right
              we would also have to check all used #include statements.
              
              Conclusion we would have to rewrite a whole C/C++ compiler and 
              this is definately beyond the scope of this project!
-----------------------------------------------------------------------------*/
void Parser::convertKnownVar( ofstream& out ){
#ifdef _DEBUG_  
  cerr<<"found known var "<<pQueue.getText(0)<<" at line: "<<pConf->getLine()<<endl;
#endif
  
  unsigned int pos=1;
  string varname=pQueue.getText(0);
  
  skipWhite(pos);

  handleArrayArg(pos,tOther); //we use tOther because we don't want the args added to the constructs

  skipWhite(pos);
  if(pQueue.getToken(pos)==assignment){
    constructs.clear();
    Var v=pConf->getVars()->lookup(pQueue.getText(0));
    constructs.insert(  v.getType(),
                        "",
                        pQueue.getText(0),
                        v.getMatchType(),
                        v.getSrcType(),
                        v.getTargetType() 
                     ); 
    
    for( unsigned int i=0; i<=pos; i++ ) moveToken(out); //copy the name... literally
 
    int braces = 0;
    while((pQueue.size()>0)&&(pQueue.getToken(0)!=semi_colon)){
      updateLocals();
      constructs.expandLastValue( pQueue.getText(0) );
      ParseElem p = pQueue.pop();
      if( p.token == open_par ) braces++;
      if( p.token == close_par ) braces--;
      //if( braces != 0 ) p.token = function_argument;
      constructs.expandLastValue( p );
    }
    
    //updateLocals();
   // handleAssignment( pos, varType );    
    
    parseAssignment(out);
  }
  else if(pQueue.getToken(pos)==comparison){
    constructs.clear();
    Var v=pConf->getVars()->lookup(pQueue.getText(0));
    constructs.insert(  v.getType(),
                        "",
                        pQueue.getText(0),
                        v.getMatchType(),
                        v.getSrcType(),
                        v.getTargetType()
                     );
                      
    for( unsigned int i=0; i<=pos; i++ ) moveToken(out); //copy name etc...
    
    int braces=1;
    while((pQueue.size()>0)&&(braces>0)){
      updateLocals();
      if(pQueue.getToken(0)==open_par) braces++;
      if(pQueue.getToken(0)==close_par) braces--;
      constructs.expandLastValue( pQueue.getText(0) );
      constructs.expandLastValue( pQueue.pop() );
    }
    parseAssignment(out);
  }
  else if((pQueue.getToken(pos)==open_par)&&(varname=="printf")){
    convertPrintf(out);
  }
  else{
    //copy literally
    moveToken(out);
  }
}





/*-----------------------------------------------------------------------------
name        :convertPotential
description :adding potential variables or functions of type tOther (not to be converted but
             added to the variables list so we can identify them later on.
parameters  :ofstream& out 
return      :/
exceptions  :/
algorithm   : first check functions defined without a type (such as a constructor
              of a class )
              If it's not a function (handle function returns false)
              we remove the first word and call addVar(out,tOther) this will
              add variable(s) or a function of type tOther if they are found
-----------------------------------------------------------------------------*/
void Parser::convertPotential(ofstream& out){
  #ifdef _DEBUG_
    cerr<<"in convertPotential ";
    showQueue();
  #endif
  //check for a function without return type given
  unsigned int pos=0;
  string funcname=getVarName(pos);
  skipWhite(pos);
  if(handleFunction(pos,funcname,tOther)){
    dumpHandled(pos,out,tOther);
  }
  else{
    moveToken(out); //copy the possible type to out
    bool done=0;
    while((!done)&&(pQueue.size()>0)){
      done=addVar(out,tOther);
    }
    constructs.clear();
  }
}


/*-----------------------------------------------------------------------------
name        :checkPrintf
description :return true if there is a printf statement on the ParseQueue
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool Parser::checkPrintf(){
  unsigned int pos=0;
  if(pQueue.getText(pos)=="printf"){
    pos++;
    if(pQueue.getToken(pos)==white_space) pos++;
    if(pQueue.getToken(pos)==open_par){
      return 1;
    }
  }
  return false;
}



/*-----------------------------------------------------------------------------
name        :convertOther
description :finding printf and declarations of tOther types as well as earlier
             defined var's
parameters  :ofstream& out 
return      :/
exceptions  :/
algorithm   : check if it's a reserved word
              check if it's a var defined earlier
              find variables and functions of type tOther
              check for printf 
-----------------------------------------------------------------------------*/
void Parser::convertOther(ofstream& out){
  types varType;
  
  //check if it's a reserved word
  if(pConf->getReswords()->isReserved(pQueue.getText(0))){
    convertReserved(out);
  }
  //check if it's a variable defined earlier
  else if( ( varType = pConf->getVars()->lookup(pQueue.getText(0)).getType() ) != tError){
    //if it's not of type tOther and not a var which has to be skipped 
    //we may convert this variable  
    if( ( varType != tOther ) &&
        (!pConf->getSkipVar()->locate(pConf->getStrCurFunction(),pQueue.getText(0))) ) {
      //convert the var if necessary meaning when assignment or comparison is being done
      convertKnownVar( out );
    }
    else{
      //var may not be converted
      //check for printf statements 
      if(checkPrintf()){
        convertPrintf(out);
      }
      else{
        moveToken(out);
      }
    }
  }
  else if(checkPrintf()){
    convertPrintf(out);
  }
  else{
    //potential new var of type tOther or function
    convertPotential(out);
  }
}


/*-----------------------------------------------------------------------------
name        :convertKnownType
description :convertion of the known types : int,float, etc
parameters  :ofstream& out 
return      :/
exceptions  :/
algorithm   : repeadetly call addVar to fill up the constructs, then call
              parseVar.
              Special cases added for :
              bMain is true (copy main function literally!)
              bArrayConstruct is true we will call parseVar
-----------------------------------------------------------------------------*/
void Parser::convertKnownType(ofstream& out,types currentVarType){
  #ifdef _DEBUG_
    cerr<<"in convertKnownType at "<<pConf->getLine()<<endl;
  #endif
  bool done=0;
  while((!done)&&(pQueue.size()>0)){
    done=addVar(out,currentVarType);
  }
  if(!bMain){
    if(constructs.size()==0){
      if( currentVarType != tOther ) writeStrType(out,typeMatch.name);
    }
    else{
      if(!bArrayConstruct) parseVar(out,pQueue.poptext());
    }
  }
  else{
    //literally output whole statement line 
    //(meaning the main function and it's arguments are always copied literally)
    
    writeOriginalType(out);
    dumpq(out);
  }
}





/*-----------------------------------------------------------------------------
name        :convertStatement
description :check for a pending assignment else 
             if we find an int,float,double etc we convert it with convertKnownType
             else we call convertOther
parameters  :ofstream& out 
return      :/
exceptions  :/
algorithm   : first check if bArrayConstruct is true (used for array
              constructors ex: int a[]={1,2,3,4}; 
              
              else
              use the setVarType function to determine the possible type
              use addVar,parseVar functions for int, float,double etc....
              use the convertOther function for everything else
              
              26/3/2002 modified handling of for loops!
-----------------------------------------------------------------------------*/
void Parser::convertStatement(ofstream& out){
  
  if(bArrayConstruct){
    #ifdef _DEBUG_
      cerr<<"pending array constructor found at line "<<pConf->getLine()<<endl;
      showQueue();
    #endif
    convertKnownType(out,ArrayType);
  }
  else{
    types currentVarType=setVarType(out);
    
    if(currentVarType!=tOther){ //int float or double
      constructs.clear();
      
      //ADDED THIS 27/1/2001 for correct handling of the -nofor option!!!
      if(pConf->getSkipFor()&&pConf->getInsideForParam()){
        writeOriginalType(out);
        convertOther(out);
      }
      else{ //NORMAL CONVERSION
        convertKnownType(out,currentVarType);
      }
    }
    else{
      if (pQueue.size()>0) convertOther(out);    
    }
  }
}




/*-----------------------------------------------------------------------------
name        :showQueue
description :show queue contents on standard error stream cerr
             used for debugging purposes
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::showQueue(){
      cerr<<"queue contents: ";
      for(vector<ParseElem>::iterator p=pQueue.begin();p!=pQueue.end();p++){
        cerr<<p->text<<" | ";
      }
      cerr<<endl;
}



/*-----------------------------------------------------------------------------
name        :updateQueue
description :add last lexer token and lexer text to queue 
             if the queue contains a complete statement
             call convertStatement repeadetly until pQueue is empty or 
             does not change anymore
parameters  :/
return      :/
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
void Parser::updateQueue( int token, string text, ofstream& out){
  pQueue.add(token,text);
  
  if((pQueue.size()>0)&&(checkStatementEnd())){
    pQueue.changed=1;
    while((pQueue.changed)&&(pQueue.size()>0)){ //keep on calling convertStatement until it doesn't change anything anymore
      pQueue.changed=0;    
      convertStatement(out);
    }
    
    
    
    //just dump rest of queue literally to outfile
    #ifdef _DEBUG_
      cerr<<"dumping this queue literally:"<<endl;
      showQueue();
      cerr<<"     at line "<<pConf->getLine()<<endl;
    #endif
    dumpq(out);
  }
}


/*-----------------------------------------------------------------------------
name        :parse
description :main entry point of the parser
             we use updateQueue to convert statements and write to out
             we use updateLocals to keep track of counters such as
             the line_nr the scope depth this is done in the move and
             dumpHandled functions.
parameters  :ofstream& out
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void Parser::parse(ofstream& out,char* lexText,int lexToken){
  string strLex="";
  strLex+=lexText;
      
  pConf->updateLineNr(lexText);
  updateQueue(lexToken,lexText,out);
        
  if(bMain){
    if(lexToken==open_brace){
      bMain=0;
      setFloatPrecision(out); //the commandline precision settings (if bNoDefault == false )
      writePrecision(out);    //precision setting in xml file =pcdata in <init> tag
      writeInit(out);         //if init option was specified.
    }
  }
        
}



