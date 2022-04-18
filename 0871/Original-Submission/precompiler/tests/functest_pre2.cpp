#include <iostream>
#include <iomanip>
using namespace std;

#include "BigInt.hh"
#include "MpIeee.hh"
#include "ArithmosIO.hh"
#include <vector>

#include <stdio.h>
BigInt f= BigInt( "1" );
BigInt pow(int  a,BigInt n){
  if(n<BigInt( "0" )) return BigInt( "0" );
  for(BigInt i= BigInt( "1" );i<n;i=i+BigInt( "1" )){
    a*=a;
  }
  return a+BigInt( "1" );
}

int main(){
  BigInt a= BigInt( "3" );BigInt b= BigInt( "1" );BigInt c= BigInt( "3" );
  if(a+b == BigInt( "3" )){
    char c;
    c=32; //c is a char here!
  }
  else{
    c=BigInt( "32" )+pow(a,b); //c is a long int
  }
  {cout<<"c="<<c<<"\n";}
  return 0;
}

