#include <iostream>
#include <iomanip>
using namespace std;

#include "BigInt.hh"
#include "MpIeee.hh"
#include "ArithmosIO.hh"
#include <vector>

int main(){
  BigInt b= BigInt( "4" );
  BigInt a[5]= {BigInt( "1" ),BigInt( "2" ),BigInt( "3" ),BigInt( "4" ),BigInt( "5" )};
  a[2]=BigInt( "12" );
  if(a[3]==BigInt( "0" )) cout<<"zero"<<endl;
  BigInt b[3][3];
  b[1][2]=BigInt( "4" );
  if(b[1][2]>BigInt( "3" )) cout<<"trivial"<<endl;
  return 0;
}

