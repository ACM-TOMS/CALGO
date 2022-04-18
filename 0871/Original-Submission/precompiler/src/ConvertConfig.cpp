#include "ConvertConfig.h"


#ifdef PRECOMPILER_GUI

//to make the rules available in the gui, add its name here
vector<string> ConvertConfig::DefinedOperations(){
  vector<string> names;
  names.push_back( "toStringConstructor" );
  names.push_back( "toConstructor" );
  names.push_back( "toRational" );
  names.push_back( "toString" );
  names.push_back( "toMatrix" );
  names.push_back( "typeDef" );
  names.push_back( "noChange" );
  names.push_back( "castIt" );
  //yes, we deliberately left out the callMember rule here, since it extends the
  //configuration format, it is not yet available in the current precompiler GUI!
  return names;
}

#else


/*=============================================================================
                              USER DEFINED OPERATIONS 
=============================================================================*/

string ConvertConfig::toString( const string& value){
  return "string( \""+value+"\" )";
}

string ConvertConfig::toStringConstructor(const ConvertElem& cElem,const string& value){
  return cElem.target.keyword + "( \"" + value + "\" )";
}

string ConvertConfig::toConstructor(const ConvertElem& cElem, const string& value){
  return cElem.target.keyword + "( " + value + " )";
}

string ConvertConfig::toRational(const ConvertElem& cElem, const string& value){
  return cElem.target.keyword + "( \"" + getRationalConstant( value ) + "\" )";
}


string ConvertConfig::castIt( const string& value ){
  return "toZeroInt32( "+value+ " )";
}

string ConvertConfig::callMember( const ConvertElem& cElem, const string& value ){
  return value + "." + cElem.operationArg;
}


//don't forget to include the new operation functionname in this if-else structure when specifying
//a new function!!!
string ConvertConfig::execute(const ConvertElem& cElem, const string& value){
  string retStr=value;
  
  if(cElem.hasOperation){
  
    if( cElem.operation == "toStringConstructor" ) {
      retStr=toStringConstructor( cElem, value );
    }
    else if( cElem.operation == "toConstructor" ) {
      retStr=toConstructor( cElem, value );
    }
    else if( cElem.operation == "toRational" ) {
      retStr=toRational( cElem, value );
    }
    else if( cElem.operation == "toString" ){
      retStr=toString( value );
    }
    else if( cElem.operation == "toMatrix" ){
      retStr="MpfrReal( \""+value+"\" )"; 
    }
    else if( cElem.operation == "typeDef" ){
      retStr="MpIeee( \""+value+"\" )";
    }
    else if( cElem.operation == "noChange" ){
      // do nothing
    }
    else if( cElem.operation == "castIt" ){
      retStr = castIt( value );
    }
    else if( cElem.operation == "callMember" ){
      retStr = callMember( cElem, value );
    }
    else{
      cerr<<"Warning: ConvertConfig::execute operation:"<<cElem.operation<<" of rule '"<<cElem.name<<"' is not defined, skipping"<<endl;
    }
  
  }
  else{
    retStr=cElem.target.keyword;
  }
  
  return retStr;
}



/*=============================================================================
                               TOKEN LIST
 Also look at include/defs.h and pre.lex before altering these.
=============================================================================*/


vector<ConvertElem> ConvertConfig::path(const string& srcWord, const string& target){
  vector<ConvertElem> result;
  result.clear();
  try{
    if(srcWord == "" ) return result;
  
    //check reserverd words;
    Reserved resCheck;
    if( resCheck.isReserved( srcWord ) ) return result;
  
    //check for known types
    result = fMatcher->convertLookup( srcWord, target );
    if( result.size() !=0 ) return result;
    
  }
  catch(Warning w){
    
    //catch Warning from Matcher 
    //check if it is a known variable or function
    //if not give meaningful Warning message
    
    Var v=fConf->getVars()->lookup(srcWord);
    //check if word is a known variable (for which conversion was not skipped)
    if( !fConf->getSkipVar()->locate( fConf->getStrCurFunction(), srcWord ) ){
      if( v.getTargetType() != "" ) {
        if( v.getTargetType() == target ){
          return result; // variable has same type, no conversion needed
        }
        else{ //lookup conversion rule for type of variable
          return path( v.getTargetType(), target );
        }
      }

    }
//    else if( v.getSrcType() != "" ){
//      return path( v.getSrcType(), target );
//    }
    
    //check for known function (which was not skipped)
    if( !fConf->getSkipVar()->locate( srcWord, srcWord ) ){
      Var func=fConf->getFunctions()->lookup(srcWord);
      if( func.getTargetType() == target ){
        return result; // function has same type, no conversion needed
      }
    }
    

    //give warning
    throw Warning( "Assignment or comparison to type "
                   + target + 
                   " may contain wrong conversion for '" 
                   + srcWord + "'" );       
  
  }
  
  return result;
}



/*-----------------------------------------------------------------------------
name        :path
description :give convert rules in vector which lead to target
parameters  :/
return      :string
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
vector<ConvertElem> ConvertConfig::path(const ParseElem& source, const string& target){
  
  switch(source.token){
    case sym_word : return path( source.text, target );              	   break;
    case decimal  : return fMatcher->convertLookup( "decimal", target );   break;
    case integer  : return fMatcher->convertLookup( "integer", target );   break;
    
    // case sym_eof : break; //the parser stops when receiving this
    
    case sym_other       : break; //copy literally
    case comment         : break; //copy literally 
    case directives      : break; //copy literally
    
    case quoted_string   : break; //copy literally
    
    case assignment      : break; //copy literally
    case comparison      : break; //copy literally
    case stream_op       : break; //copy literally
    case mult            : break; //copy literally
    case plus            : break; //copy literally
    case minus           : break; //copy literally
    case div             : break; //copy literally
    case double_colon    : break; //copy literally
    case pointer_op      : break; //copy literally
    case bool_op         : break; //copy literally
    case percent         : break; //copy literally
    
    case open_bracket    : break; //copy literally
    case close_bracket   : break; //copy literally
    case open_par        : break; //copy literally
    case close_par       : break; //copy literally
    case open_brace      : break; //copy literally
    case close_brace     : break; //copy literally
    case semi_colon      : break; //copy literally
    case colon           : break; //copy literally
    case white_space     : break; //copy literally
    
    case sym_and         : break; //copy literally
    case comma           : break; //copy literally
    case tilde           : break; //copy literally
    
    default:
      cerr<<"Warning: ConvertConfig::path received unidentified token from lexer."<<endl;
  }  
  vector<ConvertElem> empty;
  empty.clear();
  return empty;
}










/*=============================================================================
                         ADDITIONAL USEFULL FUNCTIONS
=============================================================================*/


/*-----------------------------------------------------------------------------
name        :strToCharp
description :convert string to char*
parameters  :string
return      :char*
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
char* ConvertConfig::strToCharp(string s){
     char* p = new char[s.length()+1];
     s.copy(p,string::npos);
     p[s.length()]=0; //terminator toevoegen
     return p;
}



/*-----------------------------------------------------------------------------
name        :getRationalConstant
description :this converts a string with a float or double ex. 3.5 to 35/10
             for the Rational constructor
parameters  :const string& strIn
return      :string
exceptions  :/
algorithm   :- convert the string to a char*
             - use algorithm from A.E. (written by Peter Kuterna) to do
               convertion
             - convert the changed char* to a string
             - return converted string
-----------------------------------------------------------------------------*/
string ConvertConfig::getRationalConstant(const string& strIn){
  string strOut="";  
  char* str=strToCharp(strIn);
    // filter exponent uit str
    char *ptr = str;
    int exp = 0, afterFrac = 0;
      
    while (*ptr != '\0') {
      if (*ptr == '.') afterFrac++;
      if (afterFrac) afterFrac++;
      if (*ptr == 'e' || *ptr == 'E' || *ptr == '^') exp = atoi(++ptr);
      else ptr++;
    }
      
    // filter significant uit str
    char *tmpStr = new char[strlen(str)+abs(exp)+afterFrac+5];
    if (tmpStr) {
      ptr = str;
      int count = 0, i=0, j=0;
      bool fract = false;
	
      while (*(ptr+i) != '\0' && *(ptr+i) != 'e' && *(ptr+i) != 'E' && *(ptr+i) != '^') {
	      if (*(ptr+i) != '.') {
	        if (fract) count++;
	        *(tmpStr+j) = *(ptr+i);
	        j++;
	      }
	      else fract = true;
	      i++;
      }
      for (i=0; i<exp; i++) {
        *(tmpStr+j) = '0';
        j++;
      }
      *(tmpStr+j) = '/';
      j++;
      *(tmpStr+j) = '1';
      j++;
      for (i=0; i<count; i++) {
        *(tmpStr+j) = '0';
        j++;
      }
      for (i=0; i>exp; i--) {
        *(tmpStr+j) = '0';	
        j++;
      }
      *(tmpStr+j) = '\0';
      strOut=string(tmpStr);
      delete [] tmpStr;
    }
  return strOut;
}

#endif



