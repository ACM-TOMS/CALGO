/*=============================================================================
author        :Walter Schreppers
filename      :ParseCmdLine.cpp
created       :/
modified      :25/02/2002
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#include "ParseCmdLine.h"




/*-----------------------------------------------------------------------------
name        :ParseCmdLine
description :constructor
parameters  :/
return      :ParseCmdLine
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
ParseCmdLine::ParseCmdLine(){
  init();
}



/*-----------------------------------------------------------------------------
name        :init
description :initialization
parameters  :/
return      :/
exceptions  :/
algorithm   :initialize locals
-----------------------------------------------------------------------------*/
void ParseCmdLine::init(){
  bNoFor      = 0;
  bConf       = 0;
  bInit       = 0;
  bConvert    = 0;
  bVerbose    = 0;
  bDefault    = 0;
  bNoPrintf   = 0;
  bPreParse   = 0;
  bConstants  = 0;
    
  infile      = "";
  outfile     = "";
  conffile    = "";
  convertfile = "";
  progname    = "";
  initfile    = "";
  
  fRadix        = 2;
  fPrecision    = 24;
  fExprangeLow  = -126;
  fExprangeUp   = 127;
  fRounding     = tRoundNearest;
  
  #ifdef _DEBUG_
    cerr<<"allocating the vectors for input...";
  #endif
  fMpieee.clear();
  fMpieee.push_back(DECIMAL);
  
  fRational.clear();
  fRational.push_back(RATIONAL);
  
  #ifdef _DEBUG_
    cerr<<"done"<<endl;
  #endif
}








/*-----------------------------------------------------------------------------
name        :setRoundType
description :set the rounding type
             n = tRoundNearest
             p = tRoundPlus
             m = tRoundMin
             z = tRoundZero
             return true if one of the chars n,p,m,z is given
              
parameters  :char c
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::setRoundType(char c){
  bool ok=1;
  if(c=='n'){
    fRounding=tRoundNearest;
  }
  else if(c=='p'){
    fRounding=tRoundPlus;
  }
  else if(c=='m'){
    fRounding=tRoundMin;
  }
  else if(c=='z'){
    fRounding=tRoundZero;
  }
  else{
    ok=0;
  }
  return ok;
}


/*-----------------------------------------------------------------------------
name        :addIOMpieee
description :add Mpieee output type
             p = parameter
             d = decimal
             b = binary
             y = binary representation
             h = hexidecimal representation
             r = rational
             f = flags
             return true if one of the chars p,d,b,y,h,r,f is given
parameters  :char c
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::addIOMpieee(char c){
  bool ok=1;
  switch(c){
    case('p'):  fMpieee.push_back(PARAM);      break;
    case('d'):  fMpieee.push_back(DECIMAL);    break;
    case('b'):  fMpieee.push_back(BINARY);     break;
    case('y'):  fMpieee.push_back(BINREP);     break;
    case('h'):  fMpieee.push_back(HEXREP);     break;
    case('r'):  fMpieee.push_back(RATIONAL);   break;
    case('f'):  fMpieee.push_back(FLAGS);      break;
    default:    ok=0;
  }
  return ok;
}


/*-----------------------------------------------------------------------------
name        :setIOMpieee
description :set Mpieee output type vector
             return true if a combination of chars p,d,b,y,h,r,f in 
             option string is given
parameters  :string options
return      :bool
exceptions  :/
algorithm   :walk through cp and call addIOMpieee for each char
-----------------------------------------------------------------------------*/
bool ParseCmdLine::setIOMpieee(string options){
  bool ok=0;
  if(options.size()>0) ok=1;
  unsigned int i=0;
  fMpieee.clear();
  while(ok&&(i<options.size())){
    ok=addIOMpieee(options[i]);
    i++;
  }
  return ok;
}




/*-----------------------------------------------------------------------------
name        :addIORational
description :add Rational output type
             d = decimal
             r = rational
             f = flags
             return true if one of the chars d,r,f is given
parameters  :char c
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::addIORational(char c){
  bool ok=1;
  switch(c){
    case('d'):  fRational.push_back(DECIMAL);    break;
    case('r'):  fRational.push_back(RATIONAL);   break;
    case('f'):  fRational.push_back(FLAGS);      break;
    default:    ok=0;
  }
  return ok;
}


/*-----------------------------------------------------------------------------
name        :setIORational
description :set Rational output type vector
             return true if a combination of chars p,d,b,y,h,r,f in option
             string is given
parameters  :string options
return      :bool
exceptions  :/
algorithm   :walk through cp and call addIOMpieee for each char
-----------------------------------------------------------------------------*/
bool ParseCmdLine::setIORational(string options){
  bool ok=0;
  if(options.size()>0) ok=1;
  unsigned int i=0;
  fRational.clear();
  while( ok && ( i < options.size() ) ){
    ok=addIORational(options[i]);
    i++;
  }
  return ok;
}








/*-----------------------------------------------------------------------------
name        :lpow
description :calculate a^b for longs a and b>=0
             return  a^b
parameters  :long a, long b
return      :int 
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
long ParseCmdLine::lpow(long a,long b){
  long power=1;
  for(int i=b;i>0;i--){
    power=power*a;
  }
  return power;
}


/*-----------------------------------------------------------------------------
name        :setExponentBits
description : if b>0 set the fExprangeLow = -((2^(b-1) -1) -1) 
                             fExprangeUp  = (2^(b-1) -1)
              return true if b>0
parameters  :int b
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::setExponentBits(int b){
  if(b>0){
    fExprangeLow= -( (lpow(2,b-1)-1) -1);
    fExprangeUp= lpow(2,b-1) -1;
  }
  return (b>0);
}

/*-----------------------------------------------------------------------------
name        :getUsage
description :return Usage info as a string
parameters  :/
return      :string
exceptions  :/
algorithm   :return string with info about usage of the command line
-----------------------------------------------------------------------------*/
string ParseCmdLine::getUsage() const{
  string strProg="";
  strProg+=progname;
  return "Usage: "+strProg+" [options] <inputfile> <outputfile>\n\n"+
         "Options:\n"+
         "  -nofor             don't convert variables in for loops\n"+ 
         "  -c <configfile>    specify xml file containing procedures and/or variables\n"+
         "                      which must be skipped\n"+
         "  -x <convertfile>   xml file which specifies the conversion between types\n"+
         "  -noprintf          don't convert printf statements\n"+
         "  -preparse          use preparser to make a list of all constructed \n"+
         "                      variables in <outputfile>\n"+         
         "  -constants         write position of constants 'decimal' or 'integer' to <outfile>\n"+
         "  -radix <size>      set the radix size, default size=2\n"+
         "  -precision <size>  set the precision size, default size=24\n"+
         "  -exp <low> <up>    set the exponent default low=-126, upper=127\n"+
         "  -expbits <bitsize> set the nr of bits the exponent must use\n"+
         "  -round <n|p|m|z>   set the rounding type, default = n  \n"+
         "                      n = round to nearest\n"+
         "                      p = round to plus inf\n"+
         "                      m = round to min inf\n"+
         "                      z = round to zero\n"+
         "  -outputverbose     set verbose output \n"+
         "  -outputmpieee <p|d|b|y|h|r|f>+ set the MpIeee output type, default = d \n"+
         "                      p = parameter\n"+
         "                      d = decimal\n"+
         "                      b = binary\n"+ 
         "                      y = binary representation\n"+
         "                      h = hexidecimal representation\n"+
         "                      r = rational\n"+
         "                      f = flags\n"+
         "  -outputrational <d|r|f>+ set the Rational output type, default = r  \n"+
         "                      d = decimal\n"+
         "                      r = rational\n"+
         "                      f = flags\n"+
         "  -default           insert precision settings after the main(){ \n"+
         "  -init <initfile>   insert source code contained in file after the main(){ \n\n"+
         "The options -outputmpieee and -outputrational can have more than 1\n"+
         "output type(indicated with the +)\n"+ 
         "for example: "+strProg+" in out -outputmpieee pdy\n";
}



/*-----------------------------------------------------------------------------
name        :parse
description :parse command line arguments 
parameters  :/
return      :bool
exceptions  :/
algorithm   :Parse argtable and set local vars, if commandline is
             successfully parsed true is returned.
             The if structure got somewhat longer than I had anticipated.
             But there is no easy way around it since we have to parse 
             about 30 possible options here !
-----------------------------------------------------------------------------*/
bool ParseCmdLine::parse(int argsize,char *argtable[]){
        init();
        int i = 1;
        progname=argtable[0];
        string arg;
        bool ok=1,bInf=0,bOut=0;
                
        while (i < argsize ) {
            arg="";
            arg+= argtable[i++];
            if( arg=="-nofor" ){
              bNoFor = 1;
            }
            else if( arg=="-default" ){
              bDefault = 1;
            }
            else if( arg=="-noprintf" ){
              bNoPrintf = 1; 
            }
            else if( arg=="-outputverbose" ){
              bVerbose = 1;
            }
            else if( arg=="-preparse" ){
              bPreParse = 1;
            }
            else if( arg=="-constants" ){
              bConstants = 1;
            }
            else if( arg == "-init" ){
              if( i < argsize ){
                initfile=argtable[i++];
                bInit=1;
              }
              else{
                ok=0;
              }
            }
            else if( arg == "-c" ){
                if (i < argsize){
                    conffile = argtable[i++];
                    bConf=1;
                }
                else{
                    ok=0;
                }
            }
            else if (arg=="-x"){
                if (i < argsize){
                    convertfile = argtable[i++];
                    bConvert=1;
                }
                else{
                    ok=0;
                }
            }
            else if (arg=="-radix"){
                if (i < argsize){
                    fRadix = atoi(argtable[i++]);
                }
                else{
                    ok=0;
                }
            }
            else if (arg=="-precision"){
                if (i < argsize){
                    fPrecision = atoi(argtable[i++]);
                }
                else{
                    ok=0;
                }
            }
            else if (arg=="-exp"){
                if (i < (argsize-1)){
                    fExprangeLow = atoi(argtable[i++]);
                    fExprangeUp = atoi(argtable[i++]);
                }
                else{
                    ok=0;
                }
            }
            else if (arg=="-expbits"){
                if (i < argsize){
                    int bits= atoi(argtable[i++]);
                    ok=setExponentBits(bits);
                }
                else{
                    ok=0;
                }
            }
            else if (arg=="-round"){
                if (i < argsize){
                    string cRound=argtable[i++];
                    if(cRound.size()==1) ok=setRoundType(cRound[0]);
                    else ok=0;
                }
                else{
                    ok=0;
                }
            }
            else if (arg=="-outputmpieee"){
                if (i < argsize){
                    ok=setIOMpieee(argtable[i++]);
                }
                else{
                    ok=0;
                }
            }
            else if (arg=="-outputrational"){
                if (i < argsize){
                    ok=setIORational(argtable[i++]);
                }
                else{
                    ok=0;
                }
            }
            else if(!bInf){
              bInf=1;
              infile=argtable[i-1];
            }
            else if(!bOut){
              bOut=1;
              outfile=argtable[i-1];
            }
            else{
              ok=0;
            }
        }
        
        if (!ifstream(infile)) return 0;
        return ok&&bInf&&bOut;
}




/*-----------------------------------------------------------------------------
name        :getNoFor
description :return true if -nofor option was used
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::getNoFor() const {
  return bNoFor;
}



/*-----------------------------------------------------------------------------
name        :getConfigFile
description :return true if a configfile was given 
             in command line with -c option
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::getConfigFile() const {
  return bConf;
}


/*-----------------------------------------------------------------------------
name        :getConvertFile
description :return true if a convert file was given 
             in command line with -x option
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::getConvertFile() const {
  return bConvert;
}




/*-----------------------------------------------------------------------------
name        :getInFileName
description :return input file name
parameters  :/
return      :char*
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
char* ParseCmdLine::getInFileName() const {
  return infile;
}


/*-----------------------------------------------------------------------------
name        :getOutFileName
description :return output file name
parameters  :/
return      :char*
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
char* ParseCmdLine::getOutFileName() const {
  return outfile;
}


/*-----------------------------------------------------------------------------
name        :getConfFileName
description :return config file name
parameters  :/
return      :char*
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
char* ParseCmdLine::getConfFileName() const {
  return conffile;
}



/*-----------------------------------------------------------------------------
name        :getConvertFileName
description :return convert file name
parameters  :/
return      :char*
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
char* ParseCmdLine::getConvertFileName() const {
  return convertfile;
}


/*-----------------------------------------------------------------------------
name        :getRadix
description :return radix size
parameters  :/
return      :int
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
int ParseCmdLine::getRadix() const {
  return fRadix;
}

/*-----------------------------------------------------------------------------
name        :getPrecision
description :return precision size
parameters  :/
return      :int
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
int ParseCmdLine::getPrecision() const{
  return fPrecision;
}
 
/*-----------------------------------------------------------------------------
name        :getExprangeLow
description :return exponent lower range
parameters  :/
return      :int
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
long ParseCmdLine::getExprangeLow() const{
  return fExprangeLow;
}


/*-----------------------------------------------------------------------------
name        :getExprangeUp
description :return exponent upper range
parameters  :/
return      :int
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
long ParseCmdLine::getExprangeUp() const{
  return fExprangeUp;
}
   

/*-----------------------------------------------------------------------------
name        :getRounding
description :return rounding type
parameters  :/
return      :roundType
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
roundType ParseCmdLine::getRounding() const{
  return fRounding;
}


/*-----------------------------------------------------------------------------
name        :getDefault
description : return bNoDefault
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::getDefault() const{
  return bDefault;
}

/*-----------------------------------------------------------------------------
name        :getInit
description : return bInit
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::getInit() const{
  return bInit;
}

/*-----------------------------------------------------------------------------
name        :getInitFileName
description : return initfile
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
char* ParseCmdLine::getInitFileName() const{
  return initfile;
}


/*-----------------------------------------------------------------------------
name        :getNoPrintf
description : return bNoDefault
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::getNoPrintf() const{
  return bNoPrintf;
}




/*-----------------------------------------------------------------------------
name        :getPreParse
description : return bPreParse
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::getPreParse() const{
  return bPreParse;
}


/*-----------------------------------------------------------------------------
name        :getConstants
description : return bConstants
parameters  :/
return      :bool
exceptions  :/
algorithm   :trivial
-----------------------------------------------------------------------------*/
bool ParseCmdLine::getConstants() const{
  return bConstants;
}


/*-----------------------------------------------------------------------------
name        :getIoModeStr
description :return string for setting the setIoMode in Arithmos
parameters  :/
return      :bool
exceptions  :/
algorithm   :check for verbose output and walk through the fMpieee and
             fRational vectors and catenate the settings to the return string
-----------------------------------------------------------------------------*/
string ParseCmdLine::getIoModeStr() const{
  string strIO="";
  if(bVerbose) strIO+="ARITHMOS_IO_VERBOSE|"; 
  
  unsigned int i;
  for(i=0;i<fMpieee.size();i++){
    switch(fMpieee[i]){
      case(PARAM)    : strIO+="ARITHMOS_IO_MPIEEE_PARAM|";      break;
      case(DECIMAL)  : strIO+="ARITHMOS_IO_MPIEEE_DECIMAL|";    break;
      case(BINARY)   : strIO+="ARITHMOS_IO_MPIEEE_BINARY|";     break;
      case(BINREP)   : strIO+="ARITHMOS_IO_MPIEEE_BINREP|";     break;
      case(HEXREP)   : strIO+="ARITHMOS_IO_MPIEEE_HEXREP|";     break;
      case(RATIONAL) : strIO+="ARITHMOS_IO_MPIEEE_RATIONAL|";   break;
      case(FLAGS)    : strIO+="ARITHMOS_IO_MPIEEE_FLAGS|";      break;
      default        : cerr<< "ParseCmdLine::getIoModeStr : bad iomode parameter given for an MpIeee";  break;
    }
  }
  
  for(i=0;i<fRational.size();i++){
    switch(fRational[i]){
      case(DECIMAL)  : strIO+="ARITHMOS_IO_RATIONAL_DECIMAL|";      break;
      case(RATIONAL) : strIO+="ARITHMOS_IO_RATIONAL_RATIONAL|";     break;
      case(FLAGS)    : strIO+="ARITHMOS_IO_RATIONAL_FLAGS|";        break;
      default        : cerr<< "ParseCmdLine::getIoModeStr : bad iomode parameter given for a Rational"; break;
    }
  }

  //remove the last '|'
  strIO=string(strIO,0,strIO.size()-1); 
  
  return strIO;
}

