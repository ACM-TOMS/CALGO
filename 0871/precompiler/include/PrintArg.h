/*=============================================================================
author        :Walter Schreppers
filename      :PrintArg.h
created       :05/05/2001
modified      :/
version       :1
copyright     :Walter Schreppers
bugreport(log):/
=============================================================================*/


#ifndef PARG_TYPE
#define PARG_TYPE
/*=============================================================
  enumeration to set argument type used in a printf statement 
=============================================================*/
  enum pargtype {pvar,pstring,pnumber,pdecnumber,unknown};
#endif




#ifndef PRINTARG_H
#define PRINTARG_H

#include <string>
using namespace std;

class PrintArg{
	public:
	  
	  //constructor/destructor
	  //======================
    PrintArg(const string&,pargtype); 
	  ~PrintArg();


		//members
		//=======
    string getArg() const;
    pargtype getType() const;

    // friends
    // =======
      

	private:
	  //locals
	  //======
	  string argument;
	  pargtype atype;

		//private members
		//===============


};

#endif



