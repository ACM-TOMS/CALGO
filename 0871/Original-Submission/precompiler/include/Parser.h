/*=============================================================================
author        :Walter Schreppers
filename      :Parser.h
created       :22/05/2000
modified      :25/02/2002
version       :3
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#ifndef ROUND_TYPE
#define ROUND_TYPE
/*=============================================================================
  enumeration to set the rounding type
=============================================================================*/
  enum roundType{tRoundNearest,tRoundPlus,tRoundMin,tRoundZero };
#endif


#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "Matcher.h"

#include "ParseElem.h"
#include "ParseQueue.h"
#include "ParserConf.h"
#include "ConfigFile.h"
#include "Var.h"
#include "Variables.h"
#include "Convert.h"
#include "PrintConverter.h"

#include "defs.h"


class Parser{
	public:
	  
	  //constructor/destructor
	  //======================
	  Parser();
    Parser( ParserConf* ); 
	  virtual ~Parser();
	  
		//members
		//=======
		void init(ParserConf*);
		
		void parse(ofstream&,char*,int); //the outputfile and the yytext, token from the lexer

    void setNoDefault ( bool );
    void setInit      ( bool b, char* filename );
    void setNoPrintf  ( bool );
    void setRadix     ( int );
    void setPrecision ( int );
    void setExprange  ( long, long );
    void setRounding  ( roundType );
    void setIoModeStr ( string s );

    void writeIncludes(ofstream&);
    virtual void dumpq(ofstream&);



	protected:
	  //locals
	  //======
    Variables constructs;
    ConfigFile *skipVar;
    ParserConf *pConf;
    bool bMain, bNoDefault,
         bNoPrintf, bPrintUsed,
         bArrayConstruct, bInit;
    types ArrayType;    
    int fRadix,fPrecision;
    long fExprangeLow,fExprangeUp;
    roundType fRounding;
    string fIoModeStr;
    string strOriginalType;
    char* initFileName;
    ParseQueue pQueue;
    Matcher* varMatcher;
    MatchElem typeMatch;
	  
		//protected members
		//===============
		void writeInit( ofstream& );
    void writePrecision( ofstream& );

    string getStrType( const string& );
    string getStrType( Variables::iterator p );
		
    string getOriginalType();
    void setOriginalType( unsigned int );

    virtual void parseVar( ofstream&,string );
    virtual void parseAssignment( ofstream& );
    virtual void setFloatPrecision( ofstream& );

    virtual void moveToken ( ofstream&);  
    virtual void moveTokens( unsigned int&, ofstream& );
    virtual void copyToken ( ofstream& );  
    virtual void writeStrType( ofstream&, const string& );
    virtual void writeOriginalType( ofstream& );
    virtual void convertPrintf( ofstream& );

    
    void updateLocals();
    
    bool checkStatementEnd();
    void skipToWord( ofstream& );

    void convertForLoop( ofstream& );
    
    bool checkPrintf();
    
    void convertReturn    ( ofstream& );
    void convertReserved  ( ofstream& );    //mainly used for 'printf' and 'for'
    void convertKnownVar  ( ofstream& );    //converting earlier constructed var
    void convertPotential ( ofstream& );    //potential new var or function of type tOther 
    
    void convertKnownType ( ofstream&, types ); //used for int,float,....
    void convertOther     ( ofstream& );        //used for all other cases than int,float,...
    void convertStatement ( ofstream& );        //the main starting point of converting statements

    types   setVarType( ofstream& );
    void    skipWhite ( unsigned int& );
    string  getVarSpec( unsigned int& );
    string  getVarName( unsigned int& );
    void    handleOperatorOverloading( unsigned int&, const string& );
    bool    handleFunction( unsigned int&, const string&, types );
    string  getEvenParentheses( unsigned int& );
    void    handleArrayAssign( const types );
    void    handleAssignment( unsigned int&, const types );
    void    handleArrayArg( unsigned int&, const types );
    void    dumpHandled( unsigned int&,ofstream&,const types );
    bool    addVar( ofstream&,const types );    
    
    void updateQueue( int, string, ofstream& ); //this calls convertStatement
    void showQueue();
};

#endif

