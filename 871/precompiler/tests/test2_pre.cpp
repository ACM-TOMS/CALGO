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
  BigInt a= BigInt( "3" );
  char k='b';
  for(BigInt k= BigInt( "0" );k<BigInt( "5" );k++) k=k+BigInt( "1" ); //k is an int! 
      
      //if(a>3) a=3, k=5;       // problem here: a= '3,k=5' not a='3',k='5' TODO IMPLEMENT ',' recognition!
      
  
  if(a>BigInt( "3" )) {a=BigInt( "3" );k=5;}      //this works fine k is a char here!
  k=2;  //again k is a char!

  for(BigInt k= BigInt( "0" );k<BigInt( "5" );k++){
    k=k+BigInt( "1" ); //k is an int! 
    a=k+BigInt( "1" );
  }
  
  if(a>BigInt( "3" )) {a=BigInt( "3" );k=5;}      //this works fine k is a char here!
  k=2;  //again k is a char!
  
}


