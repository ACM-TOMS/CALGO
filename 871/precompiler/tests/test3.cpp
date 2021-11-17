#include <iostream>
#include <stdio.h>
using namespace std;


unsigned int funca(int a){
  printf( "%i\n", a );
  return 5;
}


int main(){
  int a=5,i=3,j;
  int b=2.2;
  unsigned int u=5;
  float c=(float) 3.2;
  float* d=new float(3.1);
  *d=4;
  cout<<*d<<endl;
  delete d;
  char k='b';
  a=funca(i);
  
  for(int k=0;k<10;k++){
    printf("k=%i\n",k);
  }
  
  a=3;
  j=2;
}
