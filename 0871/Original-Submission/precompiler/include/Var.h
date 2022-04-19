/*=============================================================================
author        :Walter Schreppers
filename      :Var.h
created       :/
modified      :27/02/2002
version       :2
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/

#ifndef VAR_TYPES
#define VAR_TYPES
  enum types { 
                tInt,tFloat,tDouble, 
                tKnown, tOther, tOpenBracket, tError,tEmpty  
             }; //tError wordt als type gezet als de variabele niet gevonden wordt met lookup
#endif



#ifndef VAR_H
#define VAR_H

#include<string>

#include "defs.h"
#include "ParseQueue.h"

using namespace std;

class Var {
  public:
    //constructor/destructor
    //======================
    Var();
    Var(const string&,const string&,types); //name,value,type
    Var(const string&,const string&,const string&,types); //name,value,spec,type
    Var(const string&,const string&,const string&,const string&,types); //name,value,spec,array arg,type
    
    Var(  const string&,
          const string&,
          const string&,
          const string&,
          types,
          const string&
       ); //name,value,spec,array arg,type,match type
    
    Var(  const string&,
          const string&,
          const string&,
          const string&,
          types,
          const string&,
          const string&,
          const string&
       ); //name,value,spec,array arg,type,match,source, target

    Var(  types vType, 
          const string& vSpec,
          const string& vName,
          const string& matchType,
          const string& src,
          const string& trg );
		
    Var(const Var&);	//copy constructor
    ~Var();
		
		
		//members
		//=======
    string getName() const;
    string getValue() const;
    string getSpec() const;
    string getArrayArg() const;
	  string getMatchType() const;
    string getTargetType() const;
    string getSrcType() const;
    
    types getType() const;

    void setName(const string&);
    void setName(char*);

    void setValue(const string&);
    void setValue(char*);
		
    void setSpec(const string&);
    void setSpec(char*);
	  	
    void setArrayArg(const string&);
    void setArrayArg(char*);
		
    void setType(types);
    void setMatchType(const string&);

		//operators
		//=========
    bool operator<(const Var& v) const {return fValue<v.fValue;}
    Var& operator=(const Var&);
    friend ostream& operator<<( ostream&, const Var& );


    //public members
    //==============
    ParseQueue valQueue; //can hold part of the parsequeue which represents it's value

	private:
    //locals
    //======
    types fType;
    string fName;
    string fValue;
    string fSpec;
    string fArrayArg;
    string fMatchType;
    string fSrcType;
    string fTargType;
		
		//private members
		//===============
    void init();
};

#endif


