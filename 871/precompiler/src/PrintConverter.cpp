/*=============================================================================
author        :Walter Schreppers
filename      :PrintConverter.cpp
created       :/
modified      :11/05/2001
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#include "PrintConverter.h"


/*-----------------------------------------------------------------------------
name        :init
description :used by constructors
parameters  :/
return      :/
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
void PrintConverter::init(){
  format="";
  args.clear();
  flags.clear();
  pos=0;    //position in formatstring
  argpos=0; //position in arguments
  bNumberSignFlag=0; //set to 1 if '#' is found in format string flags of printf
}



/*-----------------------------------------------------------------------------
name        :PrintConverter
description :constructor
parameters  :/
return      :PrintConverter
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
PrintConverter::PrintConverter(){
  init();
}


/*-----------------------------------------------------------------------------
name        :PrintConverter
description :constructor
parameters  :const string& format, const vector<pArgs>& args
return      :PrintConverter
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
PrintConverter::PrintConverter(const string& format,const vector<PrintArg>& args){
  init();
  this->format=format;
  for(unsigned int i=0;i<args.size();i++){
    this->args.push_back(args[i]);
  }
}





/*-----------------------------------------------------------------------------
name        :PrintConverter
description :destructor
parameters  :/
return      :PrintConverter
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
PrintConverter::~PrintConverter(){
  format="";
  args.clear();
  flags.clear();
}


/*-----------------------------------------------------------------------------
name        :defaultSettings
description :returns cout default settings for inserting after a cout which
             used format settings
parameters  :/
return      :string
exceptions  :/
algorithm   :set the following things
             cout.precision(6); -> this is cout ANSI C++ standard default precision 6
                                   for all types!
             cout.fill(' '); -> default padding is a space
             we can't or the base because it requires 2 arguments to unsetf 
             (or resetiosflags) 
-----------------------------------------------------------------------------*/
string PrintConverter::defaultSettings(){
  string retStr="cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);\n";
  if(flags.size()>0){
    retStr+="cout<<resetiosflags(";
    for(unsigned int i=0;i<flags.size()-1;i++){
      retStr+="("+flags[i]+")|";
    }
    retStr+="("+flags[flags.size()-1]+"))<<";
  }
  else{
    retStr+="cout<<";
  }
  return retStr;
}


/*-----------------------------------------------------------------------------
name        :customFlags
description :returns cout flag settings found in flags vector
             this is used by getCoutStatement
parameters  :/
return      :string
exceptions  :/
algorithm   :use setiosflags( ... contents of the vector or'ed together )
             then finally add << so that output may follow
             this is all returned as a string
-----------------------------------------------------------------------------*/
string PrintConverter::customFlags(){
  string retStr="";
  if(flags.size()>0){
    retStr+="setiosflags(";
    for(unsigned int i=0;i<flags.size()-1;i++){
      retStr+="("+flags[i]+")|";
    }
    retStr+="("+flags[flags.size()-1]+"))<<";
  }
  return retStr;
}




/*-----------------------------------------------------------------------------
name        :getCoutStatement
description :return string representing a converted printf into cout
parameters  :/
return      :string
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
string PrintConverter::getCoutStatement(){
  string retStr="cout<<";
  
  if( format.size() < 3 ){ //special case for empty format string
    retStr = "cout << \"\"";
    for( unsigned int i = 0; i<args.size();i++){
      retStr+="<<"+args[i].getArg();
    }  
    retStr += ";";
    
    return retStr;  
  }
  
  string outStr="";
  format=string(format,1,format.size()-2); //strip leading and trailing "
  string coutFormatSpec="";
  while(pos<format.size()){
    if(isFormatSpec(coutFormatSpec)){
      //copy the string and an argument after it
      retStr+="\""+outStr+"\"<<";
      outStr="";
      if(argpos<args.size()){
        if((coutFormatSpec.size()>0)||(flags.size()>0)){
          retStr+=customFlags();
          retStr+=coutFormatSpec; //stream formatting string
          retStr+=args[argpos].getArg()+
                  ";\n"+defaultSettings();
        }
        else{
          retStr+=args[argpos].getArg()+"<<";
        }
      }
      else{
        cerr<<"not enough arguments corresponding to format string of cout"<<endl;
        return retStr;
      }
      argpos++;
    }//if(isFormatSpec...
    if(pos<format.size()) outStr+=format[pos];
    pos++;
  }
  retStr+="\""+outStr+"\";"; //output the remaining text in outStr and flush
  
  return retStr;
}



/*-----------------------------------------------------------------------------
name        :operator =
description := operator overloading
parameters  :const PrintConverter& 
return      :PrintConverter&
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
PrintConverter& PrintConverter::operator=(const PrintConverter& p){
  this->args.clear();
  for(unsigned int i=0;i<p.args.size();i++){
    this->args.push_back(p.args[i]);
  }
  for(unsigned int i=0;i<p.flags.size();i++){
    this->flags.push_back(p.flags[i]);
  }
	this->pos=p.pos;
	this->argpos=p.argpos;
	this->format=p.format;
  this->bNumberSignFlag=p.bNumberSignFlag;
	return *this;
}





/*-----------------------------------------------------------------------------
name        :operator<<
description :overload << operator for use in cout ... or filestreams
             returns a cout statement that represents a converted printf
             statement.
parameters  :ostream& s, const PrintConverter& d
return      :ostream&
exceptions  :/
algorithm   :use getCoutStatement on a copy of d (since d is const)
-----------------------------------------------------------------------------*/
ostream& operator<<(ostream& s,const PrintConverter& d){
  PrintConverter pcopy=d;
  s<<pcopy.getCoutStatement();
	return s;
}


/*-----------------------------------------------------------------------------
name        :insertPercent
description :insert a string argument into args containing "%"
             so that getCoutStatement can insert is as if it were an argument
             and that padding etc can be aplied to this string.
             ex.
             printf("%05%"); has to print 0000%
             ->
             cout<<setfill('0')<<setw(5)<<"%"; also prints 0000%
             
             I've noticed my gcc doesn't apply padding to a %05%
             but according to the ANSI C/C++ standard it should!
             
parameters  :ostream& s, const PrintConverter& d
return      :ostream&
exceptions  :/
algorithm   :insert a new PrintArg in args vector 
             before position argpos
-----------------------------------------------------------------------------*/
void PrintConverter::insertPercent(){
  PrintArg a("\"%\"",pstring);

  vector<PrintArg>::iterator p=args.begin();
  unsigned int i=0;
  while((i<argpos)&&(i<args.size())){
    p++;
    i++;
  }
  args.insert(p,a);
}


/*-----------------------------------------------------------------------------
name        :conversionOperation
description :if the required conversion operation is present return true
             and set the pos to after the conversion operation
             
             if percent is given as a conversion op a "%" will be inserted
             into the arguments so that it will be put on cout (with possible
             padding etc)  
             
parameters  :/
return      :bool
exceptions  :/
algorithm   :not finished yet!!!
             the p and n conversions are not implemented -> just consume argument 
              and output it
             the d,i,u conversions do nothing -> just consume argument 
              and output it
-----------------------------------------------------------------------------*/
bool PrintConverter::conversionOperation(string& coutFormatSpec){
  unsigned int retPos=pos;
  //the required conversion operation
  if(retPos<format.size()){
    switch(format[retPos]){
      case 'd'      : break;
      case 'i'      : break;
      case 'u'      : break; 
      case 'o'      : coutFormatSpec+="oct<<"; break;
      case 'x'      : coutFormatSpec+="hex<<"; break;
      case 'X'      : coutFormatSpec+="hex<<";
                      flags.push_back("ios::uppercase"); break;
      case 'c'      : break; //arg is a char nothing to be done here
      case 's'      : break; //arg is a string nothing to be done here
      case 'p'      : break; //argument is an void*  = not common
      case 'n'      : 
        cerr<<"possible wrong conversion of printf statement with '%n' as option in the format string"<<endl;
        break; //arg is an int* the number of char's output so far are put in this arg!!! = not common
      case 'f'      : flags.push_back("ios::fixed & ios::floatfield"); break;
      case 'e'      : flags.push_back("ios::scientific & ios::floatfield"); break;
      case 'E'      : flags.push_back("ios::scientific & ios::floatfield");
                      flags.push_back("ios::uppercase");
                      break;
      case 'g'      : flags.push_back("ios::floatfield"); break;  // normally 0 & ios::floatfield but does not work in g++ 3
      case 'G'      : flags.push_back("ios::floatfield");         // normally 0 & ios::floatfield
                      flags.push_back("ios::uppercase");
                      break;
      case '%'      : insertPercent(); break; //a single percent sign must be printed
      default: return 0; //not a valid conversion operation
    }
    pos=retPos+1; //skip the conversion letter d,i,....,G
    return 1;
  }
  else{//not a valid conversion operation;
    return 0;
  }
}




/*-----------------------------------------------------------------------------
name        :getDigitString
description :scan the format string for a digit string or a '*'
             when finding such a digitstring return it and set pos to after
             last digit
             when finding a * return an argument from the args vector
             and increase argpos by 1, and pos to after '*'
parameters  :/
return      :string
exceptions  :/
algorithm   :
-----------------------------------------------------------------------------*/
string PrintConverter::getDigitString(){
  string retStr="";
  if(format[pos]=='*'){
      //eat up an argument
      retStr+=args[argpos].getArg();
      argpos++;
      pos++;
  }
  else{
    bool done=0;
    while(pos<format.size()&&!done){
      if(isdigit(format[pos])){
        retStr+=format[pos];
        pos++;
      }
      else{
        done=1;
      }
    }//while
  }//if
  return retStr;
}


/*-----------------------------------------------------------------------------
name        :getWidthStatement
description :close current cout stream and insert a
             cout.width(..); statement, solves output problems for some
             compilers where setw is ignored see also getWidth()
parameters  :/
return      :string
exceptions  :/
algorithm   :use getDigitString() to get the specified width (which might also
             consume an argument if it contains a '*'
-----------------------------------------------------------------------------*/
string PrintConverter::getWidthStatement(){
  string digits=getDigitString();
  string retStr="";
  if(digits.size()>0){
    retStr="\"\";\n"; //close previous cout;
    retStr+="cout.width("+digits+");\n"; //set the width
    retStr+="cout<<"; //open a new cout statement
  }
  return retStr;
}



/*-----------------------------------------------------------------------------
name        :getWidth
description :return setw(..) command corresponding to a printf minimum width
             empty string if no width was specified in the printf format string
parameters  :/
return      :string
exceptions  :/
algorithm   :use getDigitString() to get the specified width (which might also
             consume an argument if it contains a '*'
-----------------------------------------------------------------------------*/
string PrintConverter::getWidth(){
  string digits=getDigitString();
  string retStr="";
  if(digits.size()>0){
    retStr="setw("+digits+")<<"; //does not work properly with my g++ !
  }
  return retStr;
}



/*-----------------------------------------------------------------------------
name        :getPrecision
description :return setprecision(..) command corresponding to a printf 
             precision setting if no DigitString returns an empty string
             we set precision to 0 
parameters  :/
return      :string
exceptions  :/
algorithm   :use getDigitString() to get the specified precision (which might also
             consume an argument if it contains a '*'
-----------------------------------------------------------------------------*/
string PrintConverter::getPrecision(){
  string retStr="";
  retStr+=getDigitString();
  if(retStr.size()>0){
    retStr="setprecision("+retStr+")<<";
  }
  else{
    retStr="setprecision(0)<<";
  }
  return retStr;
}




/*-----------------------------------------------------------------------------
name        :convertFlagChars
description :get 0  or more flag chars and return the corresponding cout 
             stream settings
parameters  :
return      :string retStr
exceptions  :/
algorithm   : use case structure to set the different flags in either
              the flags vector or into the retStr if they are settings that
              require an argument;
              
              a special case is the 0 it should '0' as cout padding char
              unless the '-' was specified which means left justify and the
              zeros may not be trailing so then we skip the setfill('0')
              setting
              
              TODO:
              space ' ' is not implemented because don't know the cout equivalent!!!
-----------------------------------------------------------------------------*/
string PrintConverter::convertFlagChars(){
  string retStr="";
  bool done=0;
  bNumberSignFlag=0; 
  bool bZeroPadding=0,bLeft=0;
  
  while(pos<format.size()&&!done){
    switch(format[pos]){
      case '-': pos++; bLeft=1; flags.push_back("ios::left"); break; 
      case '+': pos++; flags.push_back("ios::showpos"); break;
      case '0': pos++; bZeroPadding=1; break;
      case ' ': pos++; break; //WHAT HERE ???
      case '.': pos++; retStr+=getPrecision();break;

      case 'l': pos++; break; //says the argument is a long nothing to be done
      case 'L': pos++; break; //same as l
      case '#': pos++; bNumberSignFlag=1; break;
      
      case '1': case'2': case '3': case '4': case '5'
      : case '6' : case '7' : case '8' : case '9'
      : case '*' : retStr+=getWidth();break;
      default : done=1;
    }
  }
  if(bZeroPadding&&!bLeft){
    retStr+="setfill('0')<<";
  }
  return retStr;
}







/*-----------------------------------------------------------------------------
name        :isFormatSpec
description :return true if a format specification is fount in src 
                change the pos to right after this specification
                if precision specification is used with .* then the args will
                also be increased. The corresponding cout format settings are
                put in the local string coutFormatSpec;
             else return false and leave pos and argpos unchanged.
             
             one exception is when % is used as type specification
             this won't eat up an arg so false is returned but we will
             change the pos so that it points to this % resulting in output
             of a single %
             
parameters  :string& coutFormatSpec
return      :bool
exceptions  :/
algorithm   :
    TODO : set the correct coutFormatStr
-----------------------------------------------------------------------------*/
bool PrintConverter::isFormatSpec(string& coutFormatSpec){
  coutFormatSpec="";
  flags.clear();
  if(format[pos]=='%'){
    pos++;
    coutFormatSpec+=convertFlagChars();
    return conversionOperation(coutFormatSpec);    
  }
  else{
    return 0;
  }
}



