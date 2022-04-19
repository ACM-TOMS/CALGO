#include <iostream>
#include <iomanip>
using namespace std;

#include "BigInt.hh"
#include "MpIeee.hh"
#include "ArithmosIO.hh"
#include <vector>

#include <iostream>
#include <stdio.h>
using namespace std;


BigInt funca(BigInt a){
  {cout<<""<< a <<"\n";}
  return BigInt( "5" );
}


int main(){
  BigInt a= BigInt( "5" );BigInt i= BigInt( "3" );BigInt j;
  BigInt b= BigInt( MpIeee( "2.2" ) );
  BigInt u= BigInt( "5" );
  MpIeee c= (MpIeee) MpIeee( "3.2" );
  MpIeee * d= new MpIeee(MpIeee( "3.1" ));
  *d=MpIeee( "4" );
  cout<<*d<<endl;
  delete d;
  char k='b';
  a=funca(i);
  
  for(BigInt k= BigInt( "0" );k<BigInt( "10" );k++){
    {cout<<"k="<<k<<"\n";}
  }
  
  a=BigInt( "3" );
  j=BigInt( "2" );
}
