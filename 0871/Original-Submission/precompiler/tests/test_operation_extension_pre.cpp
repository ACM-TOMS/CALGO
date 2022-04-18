#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"


int main(){
  int a;
  MpIeee value= MpIeee( "123.4" );
  a = value.toInt();

  return 0;
}

