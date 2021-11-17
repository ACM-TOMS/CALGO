#include <iostream>
#include <stdio.h>
using namespace std;


unsigned int funca(int a){
  printf( "%i\n", a );
  return 5;
}


int main(){
  int a=3;
  char k='b';
  for(int k=0;k<5;k++) k=k+1; //k is an int! 
      
      //if(a>3) a=3, k=5;       // problem here: a= '3,k=5' not a='3',k='5' TODO IMPLEMENT ',' recognition!
      
  
  if(a>3) {a=3;k=5;}      //this works fine k is a char here!
  k=2;  //again k is a char!

  for(int k=0;k<5;k++){
    k=k+1; //k is an int! 
    a=k+1;
  }
  
  if(a>3) {a=3;k=5;}      //this works fine k is a char here!
  k=2;  //again k is a char!
  
}


