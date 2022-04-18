#include <iostream>
#include <iomanip>
using namespace std;

#include "BigInt.hh"
#include "MpIeee.hh"
#include "ArithmosIO.hh"
#include <vector>



int main(){
  BigInt a= BigInt( "3" );
  MpIeee b= MpIeee( "3.2" );

  a=BigInt( MpIeee( "3.3" ) );
  b=MpIeee( "3" );
}


