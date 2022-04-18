/*=============================================================================
author        :Walter Schreppers
filename      :ConvertConfig.h
created       :02/04/2000
modified      :
version       :
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef CONVERTCONFIG_H
#define CONVERTCONFIG_H

#ifndef PRECOMPILER_GUI
 #include "Matcher.h"
 #include "Var.h"
 #include "Reserved.h"
 #include "ParserConf.h"
#else
 #include <string>
 #include <vector>
 using std::string;
 using std::vector;
#endif

class ConvertConfig{
  public:
    //constructor/destructor
    //======================

#ifndef PRECOMPILER_GUI
    ConvertConfig(Matcher* m, ParserConf* p){ fMatcher=m; fConf=p;}
    ~ConvertConfig(){}
	 
    //public members
    //==============
    vector<ConvertElem> path(const ParseElem&, const string&);  
    vector<ConvertElem> path( const string&, const string& );
    string execute(const ConvertElem& cElem, const string& value);    
#endif

#ifdef PRECOMPILER_GUI
    static vector<string> DefinedOperations();
#endif

#ifndef PRECOMPILER_GUI
  private:
    //locals
    //======
    
    //user defined functions
    //----------------------
    string toString           ( const string& );
    
    string toStringConstructor( const ConvertElem&, const string& );
    string toConstructor      ( const ConvertElem&, const string& );
    string toRational         ( const ConvertElem&, const string& );
    string typeDef            ( const ConvertElem&, const string& );    
    
    string castIt             ( const string& );
    string callMember         ( const ConvertElem&, const string& );
    
    //needed private members
    //----------------------
    string getRationalConstant( const string& ); //use this to convert a decimal into a string suitable for rational's
    char*  strToCharp( string s );
    

    //private members
    //===============
    Matcher* fMatcher;
    ParserConf* fConf;
#endif 

};    

#endif

