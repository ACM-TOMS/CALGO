#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <stdio.h>
int main(){
	MpIeee b= MpIeee( "5" );
	int indType=  2;
	{cout<<""<<setiosflags((ios::fixed & ios::floatfield))<<b;
cout.precision(6);cout.fill(' ');cout.width(0);cout.setf(ios::dec,ios::basefield);
cout<<resetiosflags((ios::fixed & ios::floatfield))<<"";}
	indType = toZeroInt32( b ) * 2;
	return 0;
}

