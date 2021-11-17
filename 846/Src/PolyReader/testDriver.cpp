//testDriver.cpp
//find mixed volume, stable mixed volume and extended mixed volume of a polynomial system
//Author: Xing Li
//Date: 06/29/2000

#include "PolynomialSystemReader.h"
#include "PolynomialException.h"

#include <string>
#include <fstream>
#include <iostream>

      
int main(int argc, char** argv)
{
 	if( argc != 2 ) {
		cout <<"Usage: " << argv[0] << "  InputFile" <<endl;
	} else {
	  ifstream is( argv[1] );
	  PolynomialSystemReader ps( is );
	  ps.createSupportCoefficientFiles( "data.supp", "data.coef");
	}
	
       	return 0;
}
